from __future__ import annotations

import math
from dataclasses import dataclass


DEG_TO_RAD = math.pi / 180.0


@dataclass
class EngineGeometry:
    """
    Engine geometric parameters and kinematics helpers.
    - bore: Cylinder bore [m]
    - stroke: Piston stroke [m]
    - connecting_rod_length: Connecting rod length [m]
    - compression_ratio: Dimensionless compression ratio (Vmax / Vmin)
    """

    bore: float
    stroke: float
    connecting_rod_length: float
    compression_ratio: float

    def __post_init__(self) -> None:
        self.piston_area = math.pi * (self.bore ** 2) / 4.0
        self.crank_radius = self.stroke / 2.0
        self.rod_length = self.connecting_rod_length
        self.swept_volume = self.piston_area * self.stroke
        # V_clearance is the minimum volume at TDC
        self.clearance_volume = self.swept_volume / (self.compression_ratio - 1.0)

    def piston_displacement_from_tdc(self, theta_deg: float) -> float:
        """
        Returns piston displacement from TDC along cylinder axis [m] as a function of crank angle theta [deg].
        theta=0 deg corresponds to TDC between exhaust and intake (start of intake stroke).
        """
        theta = theta_deg * DEG_TO_RAD
        r = self.crank_radius
        l = self.rod_length
        # Slider-crank exact kinematics
        term = max(0.0, l**2 - (r * math.sin(theta)) ** 2)
        x = r * math.cos(theta) + math.sqrt(term)
        x_tdc = r + l
        s = x_tdc - x  # displacement from TDC
        return s

    def volume(self, theta_deg: float) -> float:
        """Instantaneous cylinder volume [m^3] at crank angle theta [deg]."""
        s = self.piston_displacement_from_tdc(theta_deg)
        return self.clearance_volume + self.piston_area * s

    def volume_and_dVdtheta(self, theta_deg: float, dtheta_deg: float = 0.1) -> tuple[float, float]:
        """
        Returns (V, dV/dtheta) where V is volume [m^3], and dV/dtheta is derivative w.r.t. crank angle [m^3/deg].
        Uses a small centered finite difference for dV/dtheta.
        """
        v = self.volume(theta_deg)
        v_plus = self.volume(theta_deg + dtheta_deg)
        v_minus = self.volume(theta_deg - dtheta_deg)
        dv_dtheta = (v_plus - v_minus) / (2.0 * dtheta_deg)
        return v, dv_dtheta

    def wall_area_estimate(self, theta_deg: float) -> float:
        """
        Very simple estimate of instantaneous in-cylinder heat transfer area [m^2]:
        piston crown area + liner area up to current piston position.
        """
        s = self.piston_displacement_from_tdc(theta_deg)
        liner_area = math.pi * self.bore * max(s, 0.0)
        return self.piston_area + liner_area