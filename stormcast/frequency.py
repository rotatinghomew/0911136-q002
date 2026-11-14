"""面向防汛物资储备讨论的极端降雨频率分析编排。

输入为 data-contract 约定的同代表（daily / eras / thresholds），输出可解释的
高分位、超越概率范围、阈值与年代合并的敏感性诊断，以及相邻乡镇同日超阈的
联合风险。无法由资料支撑的外推或联合估计会显式标注边界，不硬填百年值。
"""

from collections import defaultdict
from itertools import combinations
import random

from .extreme import (
    _BOOTSTRAP_SEED,
    DEFAULT_BOOTSTRAP_REPLICATES,
    MAX_EMPIRICAL_RETURN_FACTOR,
    MIN_EXCEEDANCES_FOR_FIT,
    RETURN_PERIODS,
    THRESHOLD_QUANTILES,
    _bootstrap_return_levels,
    binomial_wilson,
    empirical_return_period,
    gpd_annual_exceedance_prob,
    gpd_lmoments,
    gpd_mle,
    gpd_return_level,
    quantile,
)
from .analysis import summarize

_DAYS_PER_YEAR = 365.25


class FrequencyError(ValueError):
    pass


def _validate(tables):
    eras = {row["era_id"]: row for row in tables["eras"]}
    stations = {}
    towns = {}
    for era_id, row in eras.items():
        if row["station_id"] in stations and stations[row["station_id"]] != row["town_id"]:
            raise FrequencyError("同一站点对应了不同乡镇，资料契约不一致")
        stations[row["station_id"]] = row["town_id"]
        towns.setdefault(row["town_id"], []).append(era_id)
    threshold_map = {row["town_id"]: row["threshold_mm"] for row in tables["thresholds"]}
    return eras, stations, towns, threshold_map


def _collect_complete_days(tables, eras, stations):
    """按站点归集完整日，同时记录部分日的已知超阈提示与暴露天数。"""
    per_station = defaultdict(list)
    exposure_days = defaultdict(float)
    partial_known_exceed = defaultdict(list)
    for row in tables["daily"]:
        era = eras.get(row["era_id"])
        if era is None or era["station_id"] != row["station_id"]:
            raise FrequencyError("观测年代与站点不一致")
        if not 0 <= row["observed_hours"] <= 24:
            raise FrequencyError("日有效观测小时数不合法")
        station = row["station_id"]
        if row["observed_hours"] == 24 and row["rain_mm"] is not None:
            per_station[station].append(row)
            exposure_days[station] += 1.0
        else:
            # 部分日不进入频率估计；若已知雨量本身已超阈，保留为风险提示。
            if row["rain_mm"] is not None:
                partial_known_exceed[station].append(row)
    return per_station, exposure_days, partial_known_exceed


def _fit_one_threshold(complete_rows, threshold, complete_years, rng, replicates):
    """对单个阈值拟合并返回重现期水平、超越概率与自助区间。"""
    values = sorted(row["rain_mm"] for row in complete_rows)
    exceed_rows = [row for row in complete_rows if row["rain_mm"] > threshold]
    n_exceed = len(exceed_rows)
    result = {
        "threshold_mm": round(threshold, 3),
        "exceedances": n_exceed,
        "rate_per_year": round(n_exceed / complete_years, 4),
        "fit_status": "not_fitted",
        "return_levels": {},
        "sensitivity_note": None,
    }
    if n_exceed < MIN_EXCEEDANCES_FOR_FIT:
        result["sensitivity_note"] = (
            "超阈样本仅 %d 个，低于尾部拟合下限 %d，只给经验分位，不做分布外推"
            % (n_exceed, MIN_EXCEEDANCES_FOR_FIT))
        return result, exceed_rows

    exceedances = sorted(row["rain_mm"] - threshold for row in exceed_rows)
    try:
        xi_l, alpha_l = gpd_lmoments(exceedances)
        xi_m, alpha_m = gpd_mle(exceedances)
    except ValueError as exc:
        result["sensitivity_note"] = "GPD 拟合失败：%s；只给经验分位" % exc
        return result, exceed_rows

    rate = n_exceed / complete_years
    result["fit_status"] = "fitted"
    result["shape_lmoments"] = round(xi_l, 4)
    result["scale_lmoments"] = round(alpha_l, 4)
    result["shape_mle"] = round(xi_m, 4)
    result["scale_mle"] = round(alpha_m, 4)
    # 负形状意味着有界尾部，给出估计的支撑上端点作为可解释边界。
    if xi_l < 0:
        result["estimated_upper_endpoint_mm"] = round(threshold + alpha_l / -xi_l, 2)

    rate_boot, intervals = _bootstrap_return_levels(
        exceedances, threshold, complete_years, replicates, rng)
    empirical_limit = complete_years * MAX_EMPIRICAL_RETURN_FACTOR
    # 外推硬边界：资料无法支撑超过有效年数 5 倍的重现期，不硬填数值。
    extrapolation_ceiling = complete_years * 5

    for period in RETURN_PERIODS:
        if period > extrapolation_ceiling:
            result["return_levels"][str(period)] = {
                "return_period_years": period,
                "level_lmoments_mm": None,
                "level_mle_mm": None,
                "annual_exceedance_probability": None,
                "bootstrap_90pct_interval_mm": None,
                "support": "unsupported",
                "extrapolation_note": (
                    "%d 年重现期超出资料有效长度 %d 年的 5 倍外推上限，"
                    "现有样本无法支撑该量级估计，不提供数值，避免把模型外推当作储备依据"
                    % (period, round(complete_years))),
            }
            continue
        level_l = gpd_return_level(threshold, xi_l, alpha_l, rate, period)
        level_m = gpd_return_level(threshold, xi_m, alpha_m, rate, period)
        entry = {
            "return_period_years": period,
            "level_lmoments_mm": round(level_l, 2),
            "level_mle_mm": round(level_m, 2),
            "annual_exceedance_probability": round(
                gpd_annual_exceedance_prob(level_l, threshold, xi_l, alpha_l, rate), 5),
            "bootstrap_90pct_interval_mm": None,
            "support": "model_extrapolation" if period > empirical_limit else "within_record",
        }
        interval = intervals.get(period)
        if interval is not None:
            entry["bootstrap_90pct_interval_mm"] = {
                "low": round(interval["low"], 2),
                "high": round(interval["high"], 2),
            }
        if period > empirical_limit:
            entry["extrapolation_note"] = (
                "重现期 %d 年超过资料有效长度 %d 年的两倍，属分布外推，"
                "区间宽且对尾部形状敏感，不能当作实测保证值" % (period, round(complete_years)))
        # 有界尾部下请求水平若逼近端点，明确标注。
        if xi_l < 0 and level_l >= threshold + alpha_l / -xi_l * 0.95:
            entry["bounded_tail_note"] = "该水平接近有界尾部估计上端点，外推信息有限"
        result["return_levels"][str(period)] = entry
    return result, exceed_rows


def _town_analysis(town_id, era_ids, tables, eras, stations, threshold_mm,
                   rng, replicates, collected):
    per_station, exposure_days, partial_known = collected
    station_ids = sorted({eras[e]["station_id"] for e in era_ids})

    blocks = []
    merged_rows = []
    for station in station_ids:
        rows = per_station.get(station, [])
        station_eras = [e for e in era_ids if eras[e]["station_id"] == station]
        for era_id in station_eras:
            era_rows = [r for r in rows if r["era_id"] == era_id]
            years = len(era_rows) / _DAYS_PER_YEAR
            blocks.append({
                "era_id": era_id,
                "station_id": station,
                "instrument": eras[era_id]["instrument"],
                "complete_days": len(era_rows),
                "effective_years": round(years, 3),
                "maximum_mm": max((r["rain_mm"] for r in era_rows), default=None),
            })
        merged_rows.extend(rows)

    complete_years = sum(exposure_days[s] for s in station_ids) / _DAYS_PER_YEAR
    if not merged_rows:
        return {"town_id": town_id, "fit_status": "no_complete_data"}
    values = sorted(r["rain_mm"] for r in merged_rows)

    # 经验高分位与历史排位（不依赖任何分布假设）。
    empirical = {
        "effective_years": round(complete_years, 3),
        "complete_days": len(merged_rows),
        "quantiles": {str(q): round(quantile(values, q), 2) for q in
                      (0.90, 0.95, 0.975, 0.99)},
        "observed_maximum_mm": round(values[-1], 2),
        "top_events": [],
    }
    descending = sorted(merged_rows, key=lambda r: r["rain_mm"], reverse=True)
    for rank, row in enumerate(descending[:5], start=1):
        empirical["top_events"].append({
            "rank": rank,
            "day": row["day"],
            "station_id": row["station_id"],
            "rain_mm": round(row["rain_mm"], 2),
            "process_id": row["process_id"],
            "empirical_return_period_years": round(
                empirical_return_period(rank, complete_years), 2),
        })

    # 业务讨论阈值的经验超越频率（Wilson 区间，承认年数有限）。
    business_exceed = [r for r in merged_rows if r["rain_mm"] >= threshold_mm]
    wilson = binomial_wilson(len(business_exceed), len(merged_rows))
    annual_rate = len(business_exceed) / complete_years
    business = {
        "threshold_mm": threshold_mm,
        "exceedance_days": len(business_exceed),
        "empirical_annual_rate": round(annual_rate, 3),
        "empirical_return_period_years": round(1 / annual_rate, 2) if annual_rate else None,
        "daily_exceedance_probability_wilson90": (
            None if wilson is None else
            {"low": round(wilson[0], 6), "center": round(wilson[1], 6),
             "high": round(wilson[2], 6)}),
    }
    partial_hits = [r for s in station_ids for r in partial_known.get(s, [])
                    if r["rain_mm"] >= threshold_mm]
    if partial_hits:
        business["partial_day_known_exceedances"] = len(partial_hits)
        business["partial_day_note"] = (
            "另有 %d 个非完整观测日的已知雨量已达业务阈值，未计入频率估计，"
            "真实超越次数可能更高" % len(partial_hits))

    # 阈值敏感性：业务阈值 + 多个经验分位阈值，比较尾部形状是否稳定。
    threshold_fits = []
    candidate_thresholds = [threshold_mm]
    for q in THRESHOLD_QUANTILES:
        candidate_thresholds.append(quantile(values, q))
    # 去重并排序，避免业务阈值与经验分位过于接近时重复拟合。
    unique_thresholds = []
    for candidate in candidate_thresholds:
        if all(abs(candidate - kept) > 0.5 for kept in unique_thresholds):
            unique_thresholds.append(candidate)
    fit_shapes = []
    for threshold in sorted(unique_thresholds):
        fit, _ = _fit_one_threshold(merged_rows, threshold, complete_years,
                                    rng, replicates)
        threshold_fits.append(fit)
        if fit["fit_status"] == "fitted":
            fit_shapes.append(fit["shape_lmoments"])

    sensitivity = {"threshold_fits": threshold_fits}
    if len(fit_shapes) >= 2:
        spread = max(fit_shapes) - min(fit_shapes)
        sensitivity["shape_range_across_thresholds"] = {
            "low": round(min(fit_shapes), 4), "high": round(max(fit_shapes), 4)}
        sensitivity["threshold_sensitivity_note"] = (
            "尾部形状在不同阈值下的变动幅度 %.3f；%s"
            % (spread, "形状稳定，阈值选择对结论影响有限" if spread < 0.15
               else "形状对阈值敏感，重现期外推应谨慎，优先参考经验分位与区间"))
    else:
        sensitivity["threshold_sensitivity_note"] = (
            "可拟合的阈值不足两个，无法评估阈值敏感性")

    # 年代合并诊断：分别看各年代的超阈速率，提示平稳性假设是否成立。
    era_rates = []
    for block in blocks:
        era_rows = [r for r in merged_rows if r["era_id"] == block["era_id"]]
        exceed = sum(1 for r in era_rows if r["rain_mm"] >= threshold_mm)
        era_rates.append({
            "era_id": block["era_id"],
            "instrument": block["instrument"],
            "rate_per_year": round(exceed / block["effective_years"], 3),
        })
    rates = [e["rate_per_year"] for e in era_rates]
    era_diagnostic = {"era_rates": era_rates}
    if len(rates) >= 2:
        era_diagnostic["pooling_note"] = (
            "各年代超阈速率介于 %.2f–%.2f 次/年；合并分析假设不同仪器年代同分布，"
            "速率差异较大时该假设可疑，应分年代解读" % (min(rates), max(rates)))

    return {
        "town_id": town_id,
        "fit_status": "ok",
        "station_ids": station_ids,
        "empirical": empirical,
        "business_threshold": business,
        "era_pooling": era_diagnostic,
        "threshold_sensitivity": sensitivity,
    }


def _joint_exceedance(tables, eras, stations, threshold_map):
    """基于同一天共同超阈的经验联合概率，不假设联合分布形式。"""
    # 乡镇 -> 站点集合
    town_stations = defaultdict(set)
    for era in eras.values():
        town_stations[era["town_id"]].add(era["station_id"])
    town_ids = sorted(town_stations)
    # 以日历日为单位，重建每个乡镇当日是否有任一（完整日）站点超阈。
    day_town_exceed = defaultdict(lambda: defaultdict(bool))
    day_town_observed = defaultdict(lambda: defaultdict(bool))
    for row in tables["daily"]:
        if row["observed_hours"] != 24 or row["rain_mm"] is None:
            continue
        town = stations[row["station_id"]]
        day_town_observed[row["day"]][town] = True
        if row["rain_mm"] >= threshold_map[town]:
            day_town_exceed[row["day"]][town] = True

    pairs = []
    for town_a, town_b in combinations(town_ids, 2):
        common_days = [day for day, observed in day_town_observed.items()
                       if observed.get(town_a) and observed.get(town_b)]
        n_common = len(common_days)
        both = sum(1 for day in common_days
                   if day_town_exceed[day].get(town_a) and day_town_exceed[day].get(town_b))
        a_only_days = [day for day in common_days if day_town_exceed[day].get(town_a)]
        # 经验条件概率 P(B 超阈 | A 超阈) 与独立假设下的期望对照。
        p_a = len(a_only_days) / n_common if n_common else 0.0
        p_b = sum(1 for day in common_days if day_town_exceed[day].get(town_b)) / n_common if n_common else 0.0
        p_both = both / n_common if n_common else 0.0
        conditional = (both / len(a_only_days)) if a_only_days else None
        independence = p_a * p_b
        common_years = n_common / _DAYS_PER_YEAR
        pair = {
            "town_pair": [town_a, town_b],
            "common_complete_days": n_common,
            "common_effective_years": round(common_years, 2),
            "joint_exceedance_days": both,
            "joint_annual_rate": round(both / common_years, 3) if common_years else None,
            "joint_daily_probability": round(p_both, 6),
            "expected_under_independence": round(independence, 6),
            "dependence_ratio": round(p_both / independence, 2) if independence > 0 else None,
            "conditional_probability_b_given_a": (
                None if conditional is None else round(conditional, 4)),
        }
        if both < 20:
            pair["boundary_note"] = (
                "共同超阈仅 %d 次，联合频次的抽样不确定很大，"
                "条件概率只能作量级参考，不做联合分布外推" % both)
        elif common_years < 10:
            pair["boundary_note"] = (
                "共同完整观测仅 %.1f 年，联合概率为短期经验估计" % common_years)
        if pair["dependence_ratio"] is not None and pair["dependence_ratio"] > 1.5:
            pair["interpretation"] = (
                "同日共同超阈明显多于独立假设，相邻乡镇存在暴雨过程联动，"
                "储备调度应按可能同时受灾准备，不能把各乡镇风险简单相加")
        else:
            pair["interpretation"] = "共同超阈接近独立假设，联动证据有限"
        pairs.append(pair)

    # 三乡镇同时超阈（若适用）。
    if len(town_ids) >= 3:
        triple_days = [day for day, observed in day_town_observed.items()
                       if all(observed.get(t) for t in town_ids)]
        triple_both = sum(1 for day in triple_days
                          if all(day_town_exceed[day].get(t) for t in town_ids))
        triple_years = len(triple_days) / _DAYS_PER_YEAR
        triple = {
            "towns": town_ids,
            "common_complete_days": len(triple_days),
            "common_effective_years": round(triple_years, 2),
            "simultaneous_exceedance_days": triple_both,
            "simultaneous_annual_rate": round(
                triple_both / triple_years, 3) if triple_years else None,
            "simultaneous_daily_probability": round(
                triple_both / len(triple_days), 6) if triple_days else None,
        }
    else:
        triple = None
    return {"pairwise": pairs, "simultaneous_all": triple,
            "method_note": (
                "联合风险直接由同日共同超阈的经验频次计算，不拟合 copula 或多元分布；"
                "共同观测窗口之外的联动无法估计")}


def frequency_summary(tables, replicates=DEFAULT_BOOTSTRAP_REPLICATES, seed=_BOOTSTRAP_SEED):
    """生成完整的极端降雨频率分析（保留原经验摘要于 legacy 键）。"""
    eras, stations, towns, threshold_map = _validate(tables)
    rng = random.Random(seed)
    collected = _collect_complete_days(tables, eras, stations)

    town_results = []
    for town_id in sorted(towns):
        threshold_mm = threshold_map.get(town_id)
        if threshold_mm is None:
            raise FrequencyError("乡镇 %s 缺少业务阈值" % town_id)
        town_results.append(_town_analysis(
            town_id, towns[town_id], tables, eras, stations,
            threshold_mm, rng, replicates, collected))

    joint = _joint_exceedance(tables, eras, stations, threshold_map)

    return {
        "legacy_empirical_summary": summarize(tables),
        "town_frequency_analysis": town_results,
        "joint_exceedance": joint,
        "method": {
            "tail_distribution": "广义帕累托分布（GPD）超阈模型 POT",
            "estimators": "L-矩（主，含自助区间）与极大似然（点估计对照）",
            "bootstrap": "非参数重采样超阈样本，超阈次数按泊松过程重抽样",
            "bootstrap_coverage_caveat": (
                "自助区间刻画的是在已拟合模型附近的抽样离散度，并非经过覆盖率校准的置信带。"
                "极值理论中 GPD 形状参数在超阈样本仅百例量级时存在小样本偏差，"
                "区间可能系统性偏窄、对高重现期真值覆盖不足；因此区间应作为不确定性量级参考，"
                "配合阈值敏感性与经验分位一起解读，不能当作严格的 90% 置信保证"),
            "threshold_handling": "业务阈值仅用于讨论；尾部拟合阈值取经验分位并做敏感性诊断",
            "era_pooling": "同一乡镇不同仪器年代在平稳性假设下合并，另给分年代速率诊断",
        },
        "scope": (
            "完整日极端降雨频率分析。重现期是年超越概率的倒数，表示长期平均间隔，"
            "不代表每隔固定年数准时发生；区间反映抽样不确定而非气候变化趋势。"
            "资料长度或超阈样本不足时只给经验结果并标注边界，不硬填百年值。"),
    }
