"""极端日降雨的分位、超越概率与尾部频率分析（仅依赖标准库）。

设计边界
--------
* 只用完整观测日（observed_hours == 24 且 rain_mm 非空）。部分日不补零、
  不外推为完整日雨量；其已知超阈情况只作提示，不进入频率估计。
* 尾部分布采用广义帕累托分布（GPD）的超阈（POT）模型，阈值取经验分位，
  并对多个阈值做敏感性诊断，而不是把业务讨论阈值当成拟合参数。
* GPD 参数给两套估计：L-矩（闭式、稳定，用于自助重采样）与极大似然
  （仅用于点估计的方法学对照）。两者都在合成 GPD 样本上做过恢复性验证。
* 重现期是年超越概率的倒数，是长期平均间隔，不是"每隔固定年数准时发生"。
* 资料长度不足以支持的外推（超阈样本太少、请求重现期超出可信范围）一律
  返回边界说明而不是硬填数值。
"""

import math
import random

# 允许进入尾部拟合的最少超阈样本数；低于该数只给经验结果。
MIN_EXCEEDANCES_FOR_FIT = 40
# 经验分位能直接支撑的最高重现期（约样本年数的两倍），之外只给有条件外推。
MAX_EMPIRICAL_RETURN_FACTOR = 2.0
RETURN_PERIODS = (5, 10, 20, 50, 100)
# 阈值敏感性使用的经验分位（对完整日样本）。
THRESHOLD_QUANTILES = (0.90, 0.95, 0.975)
DEFAULT_BOOTSTRAP_REPLICATES = 1000
_BOOTSTRAP_SEED = 20250909


# --------------------------------------------------------------------------- #
# 基础经验统计
# --------------------------------------------------------------------------- #
def quantile(sorted_values, q):
    """线性插值经验分位（与 numpy 'linear'/'linear' 约定一致）。"""
    if not sorted_values:
        raise ValueError("空样本无法计算分位")
    if not 0 <= q <= 1:
        raise ValueError("分位必须位于 [0,1]")
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = q * (len(sorted_values) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return sorted_values[lower]
    weight = position - lower
    return sorted_values[lower] * (1 - weight) + sorted_values[upper] * weight


def empirical_return_period(rank_desc, n_years):
    """Weibull 绘图位置：按降序第 rank_desc（1 起）的经验重现期（年）。"""
    if rank_desc < 1:
        raise ValueError("名次必须从 1 开始")
    return (n_years + 1) / rank_desc


def binomial_wilson(k, n, confidence=0.90):
    """二项比例的 Wilson 区间，返回 (low, center, high)；n=0 时给不出区间。"""
    if n <= 0:
        return None
    z = {0.90: 1.6448536269514722, 0.95: 1.959963984540054}[confidence]
    p = k / n
    denom = 1 + z * z / n
    center = (p + z * z / (2 * n)) / denom
    half = (z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / denom
    return max(0.0, center - half), p, min(1.0, center + half)


# --------------------------------------------------------------------------- #
# GPD 估计（L-矩与极大似然，均已在已知参数的合成样本上验证）
# --------------------------------------------------------------------------- #
def gpd_lmoments(exceedances):
    """用 L-矩估计 GPD 形状 xi 与尺度 alpha。

    关系（Hosking & Wallis）：tau_3 = (1-xi)/(3-xi)，
    反解 xi=(3 tau_3-1)/(1+tau_3)，alpha=lambda_2 (1-xi)(2-xi)。
    """
    y = sorted(exceedances)
    n = len(y)
    if n < 10:
        raise ValueError("L-矩 GPD 拟合至少需要 10 个超阈样本")
    b0 = sum(y) / n
    b1 = sum(i * v for i, v in enumerate(y)) / (n * n)
    b2 = sum(i * (i - 1) * v for i, v in enumerate(y)) / (n ** 3)
    lambda2 = 2 * b1 - b0
    lambda3 = 6 * b2 - 6 * b1 + b0
    if lambda2 <= 0:
        raise ValueError("L-矩尺度非正，样本不支持 GPD 拟合")
    tau3 = lambda3 / lambda2
    xi = (3 * tau3 - 1) / (1 + tau3)
    alpha = lambda2 * (1 - xi) * (2 - xi)
    if alpha <= 0:
        raise ValueError("GPD 尺度估计非正")
    return xi, alpha


def gpd_mle(exceedances):
    """GPD 极大似然：外层搜索形状 xi，内层对正确 score 做括号二分求尺度。

    对数似然 ell = -n ln alpha - (1+1/xi) sum ln(1+xi x/alpha)，
    d ell/d alpha = ((1+xi) sum x/(alpha+xi x) - n)/alpha。
    """
    y = exceedances
    n = len(y)
    if n < 10:
        raise ValueError("MLE GPD 拟合至少需要 10 个超阈样本")
    ymax = max(y)
    ybar = sum(y) / n

    def loglik(xi, alpha):
        if xi < 0 and alpha <= -xi * ymax:
            return float("-inf")
        if abs(xi) < 1e-12:
            return -n * math.log(alpha) - sum(y) / alpha
        return -n * math.log(alpha) - (1 + 1 / xi) * sum(
            math.log1p(xi * v / alpha) for v in y)

    def score_alpha(xi, alpha):
        total = 0.0
        for v in y:
            denom = alpha + xi * v
            if denom <= 1e-15:
                return float("inf")
            total += v / denom
        return ((1 + xi) * total - n) / alpha

    def solve_alpha(xi):
        if xi < 0:
            low = -xi * ymax * (1 + 1e-9) + 1e-12
        else:
            low = ybar * 1e-8
        f_low = score_alpha(xi, low)
        if not math.isfinite(f_low) or f_low <= 0:
            return None
        high = ybar * 4
        guards = 0
        while score_alpha(xi, high) >= 0 and guards < 60:
            high *= 2
            guards += 1
        if score_alpha(xi, high) >= 0:
            return None
        a, b = low, high
        for _ in range(60):
            mid = (a + b) / 2
            if score_alpha(xi, mid) > 0:
                a = mid
            else:
                b = mid
        return (a + b) / 2

    best = None
    xi = -0.3
    while xi <= 0.8 + 1e-12:
        probe = xi if abs(xi) > 1e-12 else 1e-11
        alpha = solve_alpha(probe)
        if alpha is not None:
            value = loglik(probe, alpha)
            if math.isfinite(value) and (best is None or value > best[0]):
                best = (value, probe, alpha)
        xi += 0.01
    if best is None:
        raise ValueError("MLE 未收敛，样本不支持 GPD 拟合")
    # 在最优点附近用黄金分割精修形状。
    x0 = best[1]
    low, high = max(-0.3, x0 - 0.015), min(0.8, x0 + 0.015)

    def neg_profile(shape):
        alpha = solve_alpha(shape)
        return -loglik(shape, alpha) if alpha is not None else float("inf")

    golden = (math.sqrt(5) - 1) / 2
    c = high - golden * (high - low)
    d = low + golden * (high - low)
    for _ in range(30):
        if neg_profile(c) < neg_profile(d):
            high = d
        else:
            low = c
        c = high - golden * (high - low)
        d = low + golden * (high - low)
    shape = (low + high) / 2
    return shape, solve_alpha(shape)


def gpd_return_level(threshold, xi, alpha, rate_per_year, return_period):
    """POT 重现期水平：阈值 + GPD 尾部，rate_per_year 为年均超阈次数。"""
    annual_p = 1.0 / return_period
    ratio = rate_per_year / annual_p
    if ratio <= 1:
        # 请求重现期低于平均超阈间隔，阈值本身已高于该水平。
        return threshold
    if abs(xi) < 1e-9:
        return threshold + alpha * math.log(ratio)
    return threshold + alpha / xi * (ratio ** xi - 1)


def gpd_annual_exceedance_prob(level, threshold, xi, alpha, rate_per_year):
    """给定雨量水平对应的年超越概率（POT：1-exp(-lambda*(1-F_tail))）。"""
    if level <= threshold:
        return 1.0 - math.exp(-rate_per_year)
    excess = level - threshold
    if abs(xi) < 1e-9:
        tail_survival = math.exp(-excess / alpha)
    else:
        z = 1 + xi * excess / alpha
        if z <= 0:
            # 水平超过有界尾部支撑端点，超越概率视为 0。
            return 0.0
        tail_survival = z ** (-1 / xi)
    return 1.0 - math.exp(-rate_per_year * tail_survival)


# --------------------------------------------------------------------------- #
# 自助重采样（非参数 iid，重采样超阈样本；速率按泊松过程另行处理）
# --------------------------------------------------------------------------- #
def _bootstrap_return_levels(exceedances, threshold, complete_years,
                             replicates, rng):
    n = len(exceedances)
    base_rate = n / complete_years
    levels = {period: [] for period in RETURN_PERIODS}
    for _ in range(replicates):
        resampled = [exceedances[rng.randrange(n)] for _ in range(n)]
        # 超阈次数按泊松过程重抽样，体现速率本身的抽样不确定。
        rate = rng_poisson(rng, n) / complete_years
        try:
            xi, alpha = gpd_lmoments(resampled)
        except ValueError:
            continue
        for period in RETURN_PERIODS:
            levels[period].append(
                gpd_return_level(threshold, xi, alpha, rate, period))
    intervals = {}
    for period, values in levels.items():
        if len(values) < max(20, replicates // 5):
            intervals[period] = None
            continue
        values.sort()
        intervals[period] = {
            "low": percentile(values, 0.05),
            "high": percentile(values, 0.95),
            "replicates_used": len(values),
        }
    return base_rate, intervals


def rng_poisson(rng, lam):
    """Knuth 泊松抽样，用对数空间避免小概率下溢。"""
    k = 0
    log_product = 0.0
    while log_product > -lam:
        k += 1
        log_product += math.log(rng.random())
    return k - 1


def percentile(sorted_values, q):
    return quantile(sorted_values, q)
