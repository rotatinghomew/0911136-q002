"""验证频率分析的年代合并、边界处理、敏感性诊断与联合超越风险。"""

import copy
import json
import unittest
from pathlib import Path

from stormcast.frequency import frequency_summary


FIXTURE = json.loads(
    (Path(__file__).resolve().parents[1] / "fixtures" / "dataset.json").read_text())


def _make_tables(daily, eras=None, thresholds=None):
    eras = eras or [
        {"era_id": "e0-0", "station_id": "s0", "town_id": "t0",
         "start": "2005-01-01T00:00:00+00:00", "end": "2015-01-01T00:00:00+00:00",
         "instrument": "gauge-0"},
    ]
    thresholds = thresholds or [{"town_id": "t0", "threshold_mm": 50}]
    return {"daily": daily, "eras": eras, "thresholds": thresholds}


def _row(record_id, station, day, rain, hours=24, era="e0-0", process="p0"):
    return {"record_id": record_id, "station_id": station, "day": day,
            "rain_mm": rain, "observed_hours": hours,
            "era_id": era, "process_id": process}


class FrequencyBoundaryTests(unittest.TestCase):
    def test_short_record_marks_century_as_unsupported(self):
        # 仅 3 年完整资料：50/100 年重现期必须被判为无法支撑，而不是硬填。
        daily = []
        for i in range(3 * 365):
            day = "2010-01-01"
            daily.append(_row("r%d" % i, "s0",
                              "2010%03d" % (i + 1), (i % 7) * 1.0))
        tables = _make_tables(daily)
        summary = frequency_summary(tables, replicates=100)
        town = summary["town_frequency_analysis"][0]
        fits = town["threshold_sensitivity"]["threshold_fits"]
        # 至少有一个拟合存在，且所有 100 年条目要么不拟合要么标注 unsupported。
        statuses = {f["fit_status"] for f in fits}
        self.assertIn("not_fitted", statuses)
        for fit in fits:
            entry = fit["return_levels"].get("100")
            if entry is not None:
                self.assertEqual(entry["support"], "unsupported")
                self.assertIsNone(entry["level_lmoments_mm"])

    def test_partial_days_are_never_used_as_complete(self):
        daily = [
            _row("r0", "s0", "2010-01-01", 5.0),
            _row("r1", "s0", "2010-01-02", 200.0, hours=12),  # 部分日巨雨
            _row("r2", "s0", "2010-01-03", None, hours=0),
        ]
        tables = _make_tables(daily)
        summary = frequency_summary(tables, replicates=50)
        town = summary["town_frequency_analysis"][0]
        self.assertEqual(town["empirical"]["complete_days"], 1)
        self.assertEqual(town["empirical"]["observed_maximum_mm"], 5.0)

    def test_partial_day_known_exceedance_is_flagged_not_counted(self):
        daily = [_row("r%d" % i, "s0", "2010-%03d" % (i + 1), 3.0)
                 for i in range(400)]
        daily.append(_row("rx", "s0", "2011-02-01", 60.0, hours=12))
        tables = _make_tables(daily)
        summary = frequency_summary(tables, replicates=50)
        business = summary["town_frequency_analysis"][0]["business_threshold"]
        self.assertEqual(business["exceedance_days"], 0)
        self.assertEqual(business["partial_day_known_exceedances"], 1)

    def test_invalid_exposure_raises(self):
        tables = _make_tables([_row("r0", "s0", "2010-01-01", 1.0, hours=25)])
        with self.assertRaises(ValueError):
            frequency_summary(tables, replicates=50)

    def test_era_station_mismatch_raises(self):
        tables = _make_tables([_row("r0", "s9", "2010-01-01", 1.0)])
        with self.assertRaises(ValueError):
            frequency_summary(tables, replicates=50)

    def test_missing_town_threshold_raises(self):
        eras = [{"era_id": "e0-0", "station_id": "s0", "town_id": "t9",
                 "start": "2005-01-01T00:00:00+00:00",
                 "end": "2015-01-01T00:00:00+00:00", "instrument": "g0"}]
        tables = _make_tables([_row("r0", "s0", "2010-01-01", 1.0)], eras=eras)
        with self.assertRaises(ValueError):
            frequency_summary(tables, replicates=50)


class FrequencyContentTests(unittest.TestCase):
    def test_fixture_summary_is_deterministic_and_complete(self):
        tables = copy.deepcopy(FIXTURE["tables"])
        first = frequency_summary(tables, replicates=200)
        second = frequency_summary(copy.deepcopy(tables), replicates=200)
        self.assertEqual(json.dumps(first, sort_keys=True, ensure_ascii=False),
                         json.dumps(second, sort_keys=True, ensure_ascii=False))
        self.assertEqual(len(first["town_frequency_analysis"]), 3)
        for town in first["town_frequency_analysis"]:
            self.assertEqual(town["fit_status"], "ok")
            self.assertIn("threshold_sensitivity_note", town["threshold_sensitivity"])
            for fit in town["threshold_sensitivity"]["threshold_fits"]:
                for period, entry in fit["return_levels"].items():
                    self.assertIn(entry["support"],
                                  ("within_record", "model_extrapolation", "unsupported"))
        # 19 年记录下 100 年重现期超过 5 倍上限，必须显式 unsupported。
        for town in first["town_frequency_analysis"]:
            for fit in town["threshold_sensitivity"]["threshold_fits"]:
                self.assertEqual(fit["return_levels"]["100"]["support"],
                                 "unsupported")

    def test_scope_explicitly_warns_return_period_is_not_a_schedule(self):
        summary = frequency_summary(copy.deepcopy(FIXTURE["tables"]),
                                    replicates=100)
        self.assertIn("不代表每隔固定年数准时发生", summary["scope"])

    def test_legacy_summary_is_embedded_unchanged(self):
        from stormcast.analysis import summarize
        tables = copy.deepcopy(FIXTURE["tables"])
        summary = frequency_summary(tables, replicates=100)
        # 现行摘要的 eras 分代结果被原样内嵌，scope 保留其原始边界表述。
        self.assertEqual(summary["legacy_empirical_summary"], summarize(tables))
        self.assertIn("scope", summary["legacy_empirical_summary"])


class JointExceedanceTests(unittest.TestCase):
    def _two_town_tables(self, joint_days):
        eras = [
            {"era_id": "e0-0", "station_id": "s0", "town_id": "t0",
             "start": "2005-01-01T00:00:00+00:00",
             "end": "2015-01-01T00:00:00+00:00", "instrument": "g0"},
            {"era_id": "e1-0", "station_id": "s1", "town_id": "t1",
             "start": "2005-01-01T00:00:00+00:00",
             "end": "2015-01-01T00:00:00+00:00", "instrument": "g0"},
        ]
        thresholds = [{"town_id": "t0", "threshold_mm": 50},
                      {"town_id": "t1", "threshold_mm": 50}]
        daily = []
        n_days = 400
        import datetime
        base = datetime.date(2010, 1, 1)
        for i in range(n_days):
            day = (base + datetime.timedelta(days=i)).isoformat()
            rain0 = 60.0 if i < joint_days else 2.0
            rain1 = 60.0 if i < joint_days else 2.0
            daily.append(_row("a%d" % i, "s0", day, rain0, era="e0-0",
                              process="p%d" % i))
            daily.append(_row("b%d" % i, "s1", day, rain1, era="e1-0",
                              process="p%d" % i))
        return _make_tables(daily, eras=eras, thresholds=thresholds)

    def test_perfect_dependence_is_detected(self):
        tables = self._two_town_tables(joint_days=30)
        joint = frequency_summary(tables, replicates=50)["joint_exceedance"]
        pair = joint["pairwise"][0]
        self.assertEqual(pair["joint_exceedance_days"], 30)
        self.assertEqual(pair["conditional_probability_b_given_a"], 1.0)
        self.assertGreater(pair["dependence_ratio"], 1.5)
        self.assertIn("联动", pair["interpretation"])

    def test_joint_uses_only_common_complete_days(self):
        import datetime
        tables = self._two_town_tables(joint_days=10)
        base = datetime.date(2010, 1, 1)
        partial_days = {(base + datetime.timedelta(days=389)).isoformat(),
                        (base + datetime.timedelta(days=390)).isoformat()}
        for row in tables["daily"]:
            if row["station_id"] == "s1" and row["day"] in partial_days:
                row["observed_hours"] = 12
        joint = frequency_summary(tables, replicates=50)["joint_exceedance"]
        pair = joint["pairwise"][0]
        self.assertEqual(pair["common_complete_days"], 398)

    def test_no_copula_extrapolation_claim_in_method_note(self):
        tables = self._two_town_tables(joint_days=5)
        joint = frequency_summary(tables, replicates=50)["joint_exceedance"]
        self.assertIn("不拟合 copula", joint["method_note"])
        self.assertIn("boundary_note", joint["pairwise"][0])


class MultiEraPoolingTests(unittest.TestCase):
    def _two_era_town(self, rate_era0, rate_era1):
        import datetime
        eras = [
            {"era_id": "e0-0", "station_id": "s0", "town_id": "t0",
             "start": "2005-01-01T00:00:00+00:00",
             "end": "2010-01-01T00:00:00+00:00", "instrument": "gauge-0"},
            {"era_id": "e0-1", "station_id": "s0", "town_id": "t0",
             "start": "2010-01-01T00:00:00+00:00",
             "end": "2015-01-01T00:00:00+00:00", "instrument": "gauge-1"},
        ]
        base = datetime.date(2005, 1, 1)
        daily = []
        rid = 0
        for era_id, years, heavy_every in (
                ("e0-0", 5, rate_era0), ("e0-1", 5, rate_era1)):
            for i in range(years * 365):
                day = (base + datetime.timedelta(days=rid)).isoformat()
                rain = 60.0 if (i % heavy_every == 0) else 3.0
                daily.append(_row("r%d" % rid, "s0", day, rain, era=era_id))
                rid += 1
        return _make_tables(daily, eras=eras)

    def test_pooling_merges_eras_but_reports_separate_rates(self):
        # 两个年代超阈频率刻意不同，验证合并发生且分代速率被如实呈现。
        tables = self._two_era_town(rate_era0=40, rate_era1=15)
        town = frequency_summary(tables, replicates=100)["town_frequency_analysis"][0]
        self.assertEqual(town["station_ids"], ["s0"])
        # 合并后完整日为两个年代之和。
        self.assertEqual(town["empirical"]["complete_days"], 10 * 365)
        rates = {e["era_id"]: e["rate_per_year"]
                 for e in town["era_pooling"]["era_rates"]}
        self.assertIn("e0-0", rates)
        self.assertIn("e0-1", rates)
        self.assertNotAlmostEqual(rates["e0-0"], rates["e0-1"], places=1)
        self.assertIn("pooling_note", town["era_pooling"])
        # 两个仪器名都应出现在分代诊断里。
        instruments = {e["instrument"] for e in town["era_pooling"]["era_rates"]}
        self.assertEqual(instruments, {"gauge-0", "gauge-1"})


if __name__ == "__main__":
    unittest.main()
