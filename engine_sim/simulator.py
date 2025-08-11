from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, List, Tuple

from .geometry import EngineGeometry
from .valve import ValveProfile
from .combustion import WiebeFunction


R_AIR = 287.0  # J/(kg.K)
GAMMA = 1.35   # specific heat ratio, approximate for hot air-fuel
CP = GAMMA * R_AIR / (GAMMA - 1.0)
CV = CP - R_AIR


@dataclass
class BoundaryConditions:
    intake_pressure_pa: float
    intake_temperature_k: float
    exhaust_pressure_pa: float
    exhaust_temperature_k: float
    wall_temperature_k: float = 450.0


@dataclass
class EngineParams:
    geometry: EngineGeometry
    rpm: float
    intake_valve: ValveProfile
    exhaust_valve: ValveProfile
    combustion: WiebeFunction
    equivalence_ratio: float = 1.0  # phi
    afr_stoich: float = 14.7        # air-fuel mass ratio stoichiometric
    heat_transfer_coeff_w_m2k: float = 150.0
    num_cycles: int = 1


def critical_pressure_ratio(gamma: float) -> float:
    return (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))


def is_choked(p_up: float, p_down: float, gamma: float) -> bool:
    return (p_down / max(p_up, 1e-9)) <= critical_pressure_ratio(gamma)


def orifice_mass_flow(area_m2: float, cd: float, p_up: float, t_up: float, p_down: float, gamma: float = GAMMA, r: float = R_AIR) -> float:
    """
    Compressible orifice mass flow [kg/s] from upstream to downstream.
    Positive direction: from upstream to downstream.
    """
    if area_m2 <= 0.0 or p_up <= 1.0:
        return 0.0

    pr = max(p_down / max(p_up, 1e-9), 1e-9)
    coef = cd * area_m2 * p_up * math.sqrt(gamma / (r * max(t_up, 1e-9)))
    if pr <= critical_pressure_ratio(gamma):
        # choked flow
        return coef * ((2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0))))
    # subsonic
    return coef * (pr ** (1.0 / gamma)) * math.sqrt((2.0 * gamma / (gamma - 1.0)) * (1.0 - pr ** ((gamma - 1.0) / gamma)))


def trapezoidal_integral(y_values: List[float], x_values: List[float]) -> float:
    total = 0.0
    for i in range(1, len(x_values)):
        dx = x_values[i] - x_values[i - 1]
        total += 0.5 * (y_values[i] + y_values[i - 1]) * dx
    return total


class EngineSimulator:
    def __init__(self, params: EngineParams, boundaries: BoundaryConditions) -> None:
        self.params = params
        self.boundaries = boundaries
        self.omega_deg_per_s = 6.0 * params.rpm  # deg/s (360 deg per rev, rpm/60 rev/s => 6*rpm deg/s)

    def run(self, dtheta_deg: float = 0.5) -> Dict[str, List[float]]:
        geometry = self.params.geometry
        total_degrees = 720.0 * self.params.num_cycles

        thetas: List[float] = []
        volumes: List[float] = []
        dv_dtheta: List[float] = []
        pressures: List[float] = []
        temperatures: List[float] = []
        masses: List[float] = []
        m_dot_cyl: List[float] = []
        heat_release_rate: List[float] = []
        q_loss_rate: List[float] = []
        intake_area: List[float] = []
        exhaust_area: List[float] = []

        # Initial state at theta=0 (TDC between exhaust and intake).
        V0, _ = geometry.volume_and_dVdtheta(0.0)
        p0 = self.boundaries.intake_pressure_pa
        T0 = self.boundaries.intake_temperature_k
        m0 = max(p0 * V0 / (R_AIR * T0), 1e-6)

        U = m0 * CV * T0
        m = m0

        # Precompute Q_total when combustion starts (based on trapped air at SOC)
        q_total = None

        theta = 0.0
        while theta <= total_degrees + 1e-9:
            V, dVdth = geometry.volume_and_dVdtheta(theta)
            thetas.append(theta)
            volumes.append(V)
            dv_dtheta.append(dVdth)

            # current p, T
            T = max(U / max(m, 1e-9) / CV, 1.0)
            p = max(m * R_AIR * T / max(V, 1e-12), 100.0)
            pressures.append(p)
            temperatures.append(T)
            masses.append(m)

            # Areas
            Ai_eff = self.params.intake_valve.effective_flow_area(theta % 720.0)
            Ae_eff = self.params.exhaust_valve.effective_flow_area(theta % 720.0)
            intake_area.append(Ai_eff)
            exhaust_area.append(Ae_eff)

            # Mass flows into cylinder (positive into cylinder)
            m_dot_in = 0.0
            m_dot_out = 0.0

            # Intake flow: from intake manifold to cylinder if p_intake > p_cyl; handle backflow
            if Ai_eff > 0.0:
                if self.boundaries.intake_pressure_pa >= p:
                    md = orifice_mass_flow(Ai_eff, self.params.intake_valve.discharge_coefficient,
                                           self.boundaries.intake_pressure_pa, self.boundaries.intake_temperature_k,
                                           p, GAMMA)
                    m_dot_in += md
                else:
                    # backflow to intake
                    md = orifice_mass_flow(Ai_eff, self.params.intake_valve.discharge_coefficient,
                                           p, T, self.boundaries.intake_pressure_pa, GAMMA)
                    m_dot_out += md

            # Exhaust flow: from cylinder to exhaust if p_cyl > p_exh; handle backflow
            if Ae_eff > 0.0:
                if p >= self.boundaries.exhaust_pressure_pa:
                    md = orifice_mass_flow(Ae_eff, self.params.exhaust_valve.discharge_coefficient,
                                           p, T, self.boundaries.exhaust_pressure_pa, GAMMA)
                    m_dot_out += md
                else:
                    md = orifice_mass_flow(Ae_eff, self.params.exhaust_valve.discharge_coefficient,
                                           self.boundaries.exhaust_pressure_pa, self.boundaries.exhaust_temperature_k,
                                           p, GAMMA)
                    m_dot_in += md

            m_dot_net = m_dot_in - m_dot_out
            m_dot_cyl.append(m_dot_net)

            # Heat release from combustion (Wiebe)
            theta_mod = theta % 720.0
            dxb_dtheta = self.params.combustion.dxb_dtheta(theta_mod)
            if q_total is None and theta_mod >= self.params.combustion.start_deg:
                # determine energy of fuel based on trapped air mass at SOC
                m_air = m
                m_fuel = (m_air / self.params.afr_stoich) * self.params.equivalence_ratio
                LHV = 44e6  # J/kg, gasoline-like
                q_total = m_fuel * LHV
            qdot_comb = 0.0 if q_total is None else (dxb_dtheta * q_total * self.omega_deg_per_s)
            heat_release_rate.append(qdot_comb)

            # Heat loss to walls (very simple model)
            Aw = geometry.wall_area_estimate(theta_mod)
            h = self.params.heat_transfer_coeff_w_m2k
            qdot_loss = h * Aw * max(T - self.boundaries.wall_temperature_k, 0.0)
            q_loss_rate.append(qdot_loss)

            # Enthalpy transport
            h_in = CP * self.boundaries.intake_temperature_k
            h_exh = CP * T  # assume outflow enthalpy at cylinder temperature

            # Work term p dV/dt
            dVdt = dVdth * self.omega_deg_per_s
            pdotV = p * dVdt

            # Integrate over dt
            dt = dtheta_deg / self.omega_deg_per_s
            dU = (qdot_comb - qdot_loss - pdotV + h_in * m_dot_in - h_exh * m_dot_out) * dt
            U = max(U + dU, 1.0)
            m = max(m + m_dot_net * dt, 1e-7)

            theta += dtheta_deg

        # Compute P-V work and IMEP
        work = trapezoidal_integral(pressures, volumes)
        displacement = self.params.geometry.swept_volume
        imep = work / displacement if displacement > 0 else 0.0

        return {
            "theta_deg": thetas,
            "volume_m3": volumes,
            "pressure_pa": pressures,
            "temperature_k": temperatures,
            "mass_kg": masses,
            "m_dot_kg_s": m_dot_cyl,
            "qdot_comb_w": heat_release_rate,
            "qdot_loss_w": q_loss_rate,
            "intake_area_m2": intake_area,
            "exhaust_area_m2": exhaust_area,
            "imep_pa": [imep],
            "work_j": [work],
        }