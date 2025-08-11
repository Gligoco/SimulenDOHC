from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional


@dataclass
class ValveProfile:
    """
    Defines a simple sinusoidal lift profile based on open and close angles.
    All angles are in crank degrees in [0, 720).
    - open_deg: opening angle [deg]
    - close_deg: closing angle [deg]; may be less than open_deg to indicate wrap-around
    - max_lift: maximum valve lift [m]
    - valve_diameter: valve head diameter [m]
    - num_valves: number of identical valves
    - discharge_coefficient: dimensionless flow discharge coefficient (0.5-0.9 typical)
    """

    open_deg: float
    close_deg: float
    max_lift: float
    valve_diameter: float
    num_valves: int
    discharge_coefficient: float = 0.7

    def duration(self) -> float:
        if self.close_deg >= self.open_deg:
            return self.close_deg - self.open_deg
        # wrap-around over 720 deg
        return (720.0 - self.open_deg) + self.close_deg

    def phase_fraction(self, theta_deg: float) -> Optional[float]:
        """
        Returns phase fraction f in [0,1] within the open duration if valve is open at theta_deg, else None.
        Uses wrap-around logic over 720 deg cycle.
        """
        t = theta_deg % 720.0
        if self.close_deg >= self.open_deg:
            if t < self.open_deg or t > self.close_deg:
                return None
            return (t - self.open_deg) / max(self.duration(), 1e-9)
        # wrap case
        if t < self.open_deg and t > self.close_deg:
            return None
        if t >= self.open_deg:
            return (t - self.open_deg) / max(self.duration(), 1e-9)
        # t <= close_deg
        return (t + (720.0 - self.open_deg)) / max(self.duration(), 1e-9)

    def lift(self, theta_deg: float) -> float:
        f = self.phase_fraction(theta_deg)
        if f is None:
            return 0.0
        # Smooth sinusoidal cam-like lift
        return self.max_lift * math.sin(math.pi * f)

    def effective_flow_area(self, theta_deg: float) -> float:
        """
        Returns the effective geometric flow area [m^2] of the valve curtain times discharge coefficient.
        Curtain area ~ circumference * lift per valve, times number of valves.
        Caps area by valve curtain at max lift but does not model seat throttling separately.
        """
        current_lift = self.lift(theta_deg)
        if current_lift <= 0.0:
            return 0.0
        circumference = math.pi * self.valve_diameter
        geometric_area = circumference * current_lift * self.num_valves
        return geometric_area * self.discharge_coefficient