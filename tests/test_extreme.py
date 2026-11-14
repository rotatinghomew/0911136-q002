"""验证 GPD 尾部估计、经验分位与超越概率的数值性质（纯标准库）。"""

import math
import random
import unittest

from stormcast.extreme import (
    binomial_wilson,
    empirical_return_period,
    gpd_annual_exceedance_prob,
    gpd_lmoments,
    gpd_mle,
    gpd_return_level,
    quantile,
    rng_poisson,
)


def _rgpd(rng, n, xi, alpha):
    out = []
    for _ in range(n):
        u = rng.random()
        if abs(xi) < 1e-9:
            out.append(-alpha * math.log(1 - u))
        else:
            out.append(alpha / xi * ((1 - u) ** (-xi) - 1))
    return out


class EstimatorTests(unittest.TestCase):
    def setUp(self):
        self.rng = random.Random(1234)

    def test_lmoments_recovers_known_shapes(self):
        for true_xi in (-0.2, -0.1, 0.0, 0.15, 0.4):
            sample = _rgpd(self.rng, 60000, true_xi, 10.0)
            xi, alpha = gpd_lmoments(sample)
            self.assertAlmostEqual(xi, true_xi, delta=0.04)
            self.assertAlmostEqual(alpha, 10.0, delta=0.4)

    def test_mle_recovers_known_shapes(self):
        for true_xi in (-0.2, 0.0, 0.15):
            sample = _rgpd(self.rng, 2500, true_xi, 10.0)
            xi, alpha = gpd_mle(sample)
            self.assertAlmostEqual(xi, true_xi, delta=0.06)
            self.assertAlmostEqual(alpha, 10.0, delta=0.8)

    def test_mle_and_lmoments_agree_on_exponential_limit(self):
        sample = _rgpd(self.rng, 4000, 0.0, 10.0)
        xi_l, alpha_l = gpd_lmoments(sample)
        xi_m, alpha_m = gpd_mle(sample)
        self.assertAlmostEqual(xi_l, 0.0, delta=0.08)
        self.assertAlmostEqual(xi_m, 0.0, delta=0.08)
        self.assertAlmostEqual(alpha_l, 10.0, delta=0.6)
        self.assertAlmostEqual(alpha_m, 10.0, delta=0.6)

    def test_too_few_exceedances_rejected(self):
        with self.assertRaises(ValueError):
            gpd_lmoments([1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            gpd_mle([1.0, 2.0, 3.0])

    def test_return_level_and_probability_are_consistent(self):
        threshold, xi, alpha, rate = 40.0, 0.12, 9.0, 5.0
        for period in (2, 5, 20, 50):
            level = gpd_return_level(threshold, xi, alpha, rate, period)
            annual_p = gpd_annual_exceedance_prob(
                level, threshold, xi, alpha, rate)
            # POT 模型：重现期 T 对应年超越概率 1-exp(-1/T)，而非简单 1/T。
            self.assertAlmostEqual(annual_p, 1.0 - math.exp(-1.0 / period),
                                   places=6)
            self.assertGreater(level, threshold)
        # 小重现期下 POT 概率与朴素 1/T 的差异可被明确区分，确认没用错公式。
        level2 = gpd_return_level(threshold, xi, alpha, rate, 2)
        p2 = gpd_annual_exceedance_prob(level2, threshold, xi, alpha, rate)
        self.assertAlmostEqual(p2, 1.0 - math.exp(-0.5), places=6)
        self.assertNotAlmostEqual(p2, 0.5, places=3)

    def test_exponential_tail_matches_closed_form(self):
        threshold, alpha, rate = 40.0, 8.0, 4.0
        level = gpd_return_level(threshold, 0.0, alpha, rate, 25)
        expected = threshold + alpha * math.log(rate * 25)
        self.assertAlmostEqual(level, expected, places=6)

    def test_bounded_tail_probability_is_zero_beyond_endpoint(self):
        threshold, xi, alpha, rate = 40.0, -0.2, 10.0, 4.0
        endpoint = threshold + alpha / -xi  # 90
        self.assertEqual(
            gpd_annual_exceedance_prob(endpoint + 5, threshold, xi, alpha, rate),
            0.0)

    def test_quantile_matches_linear_interpolation(self):
        values = [1.0, 2.0, 3.0, 4.0]
        self.assertEqual(quantile(sorted(values), 0.0), 1.0)
        self.assertEqual(quantile(sorted(values), 1.0), 4.0)
        self.assertAlmostEqual(quantile(sorted(values), 0.5), 2.5)

    def test_empirical_return_period_uses_weibull_position(self):
        self.assertEqual(empirical_return_period(1, 20), 21.0)
        self.assertEqual(empirical_return_period(2, 20), 10.5)

    def test_wilson_interval_covers_center_and_bounds(self):
        low, center, high = binomial_wilson(10, 100)
        self.assertLessEqual(low, center)
        self.assertLessEqual(center, high)
        self.assertAlmostEqual(center, 0.1)
        self.assertIsNone(binomial_wilson(0, 0))

    def test_poisson_sampler_matches_mean(self):
        rng = random.Random(7)
        draws = [rng_poisson(rng, 20.0) for _ in range(4000)]
        self.assertAlmostEqual(sum(draws) / len(draws), 20.0, delta=0.4)

    def test_small_sample_shape_bias_is_documented_behavior(self):
        # 现实样本量（约 150 个超阈）下 L-矩形状参数存在向下的小样本偏差。
        # 这是极值估计的已知现象，固化为测试以防止被误当回归，并支撑
        # frequency 输出中对自助区间覆盖率的免责说明。
        rng = random.Random(99)
        true_xi, true_alpha = 0.12, 9.0
        shapes = []
        for _ in range(300):
            sample = _rgpd(rng, 150, true_xi, true_alpha)
            xi, _ = gpd_lmoments(sample)
            shapes.append(xi)
        mean_shape = sum(shapes) / len(shapes)
        self.assertLess(mean_shape, true_xi)
        # 大样本下偏差应消失，证明估计器本身渐近无偏。
        large = _rgpd(random.Random(5), 60000, true_xi, true_alpha)
        xi_large, _ = gpd_lmoments(large)
        self.assertAlmostEqual(xi_large, true_xi, delta=0.04)


if __name__ == "__main__":
    unittest.main()
