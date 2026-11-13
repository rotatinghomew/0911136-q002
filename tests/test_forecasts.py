import unittest
from datetime import datetime, timezone

from stormcast.forecasts import Forecast, eligible


class ForecastTests(unittest.TestCase):
    def test_probability_must_be_valid_and_issued_in_time(self):
        now = datetime.now(timezone.utc)
        self.assertTrue(eligible(Forecast("p", "v", now, 60, 0.4), now))
        self.assertFalse(eligible(Forecast("p", "v", now, 60, 1.4), now))
