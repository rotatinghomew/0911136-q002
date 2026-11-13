"""预报发行与观测实况的基础表示。"""

from dataclasses import dataclass
from datetime import datetime


@dataclass(frozen=True)
class Forecast:
    process_id: str
    area_id: str
    issued_at: datetime
    lead_minutes: int
    probability: float


def eligible(forecast: Forecast, available_at: datetime) -> bool:
    return forecast.issued_at <= available_at and 0.0 <= forecast.probability <= 1.0
