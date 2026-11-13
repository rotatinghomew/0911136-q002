"""验证现有领域摘要的业务边界。"""
import copy
import json
import unittest
from pathlib import Path

from stormcast.analysis import era_profile


class DomainTests(unittest.TestCase):
    def setUp(self):
        self.tables = json.loads((Path(__file__).resolve().parents[1]/'fixtures/dataset.json').read_text())['tables']

    def test_partial_day_cannot_set_complete_day_maximum(self):
        tables = copy.deepcopy(self.tables)
        tables['daily'] = [dict(tables['daily'][0],observed_hours=12,rain_mm=9999)]
        result = era_profile(tables)
        profile = next(iter(result.values()))
        self.assertEqual(0,profile['complete_days'])
        self.assertIsNone(profile['maximum_mm'])

    def test_zero_is_a_valid_complete_observation(self):
        tables = copy.deepcopy(self.tables)
        tables['daily'] = [dict(tables['daily'][0],observed_hours=24,rain_mm=0)]
        self.assertEqual(0,next(iter(era_profile(tables).values()))['maximum_mm'])

    def test_invalid_exposure_fails(self):
        tables = copy.deepcopy(self.tables)
        tables['daily'][0]['observed_hours'] = 25
        with self.assertRaises(ValueError):
            era_profile(tables)

    def test_station_cannot_borrow_another_era(self):
        tables = copy.deepcopy(self.tables)
        tables['daily'][0]['era_id'] = 'e1-0'
        with self.assertRaises(ValueError):
            era_profile(tables)
