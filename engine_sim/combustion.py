from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional


@dataclass
class WiebeFunction:
    """
    Wiebe function parameters to describe mass fraction burned over crank angle.
    xb(theta) = 1 - exp(-a * ((theta - theta0)/delta)^(m+1)) for theta within [theta0, theta0+delta].
    """

    start_deg: float
    duration_deg: float
    a: float = 5.0
    m: float = 2.0

    def mass_fraction_burned(self, theta_deg: float) -> float:
        if theta_deg < self.start_deg:
            return 0.0
        if theta_deg > self.start_deg + self.duration_deg:
            return 1.0
        x = (theta_deg - self.start_deg) / max(self.duration_deg, 1e-9)
        return 1.0 - math.exp(-self.a * (x ** (self.m + 1.0)))

    def dxb_dtheta(self, theta_deg: float) -> float:
        if theta_deg < self.start_deg or theta_deg > self.start_deg + self.duration_deg:
            return 0.0
        x = (theta_deg - self.start_deg) / max(self.duration_deg, 1e-9)
        # derivative wrt theta (deg)
        return (
            self.a * (self.m + 1.0) / max(self.duration_deg, 1e-9)
            * (x ** self.m)
            * math.exp(-self.a * (x ** (self.m + 1.0)))
        )