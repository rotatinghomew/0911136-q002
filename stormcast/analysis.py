"""整理按观测年代划分的完整日雨量经验摘要。"""
from collections import defaultdict


def era_profile(tables):
    eras = {row["era_id"]: row for row in tables["eras"]}
    grouped = defaultdict(lambda: {"complete_days": 0, "partial_days": 0, "maximum_mm": None})
    for row in tables["daily"]:
        if row["era_id"] not in eras or eras[row["era_id"]]["station_id"] != row["station_id"]:
            raise ValueError("观测年代与站点不一致")
        if not 0 <= row["observed_hours"] <= 24:
            raise ValueError("日有效观测小时数不合法")
        profile = grouped[row["era_id"]]
        if row["observed_hours"] < 24 or row["rain_mm"] is None:
            profile["partial_days"] += 1
        else:
            profile["complete_days"] += 1
            maximum = profile["maximum_mm"]
            profile["maximum_mm"] = row["rain_mm"] if maximum is None else max(maximum, row["rain_mm"])
    return dict(sorted(grouped.items()))


def summarize(tables):
    return {"eras": era_profile(tables), "scope": "完整日经验摘要，不估计历史外重现期"}
