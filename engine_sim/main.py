from __future__ import annotations

import os
import csv

from .geometry import EngineGeometry
from .valve import ValveProfile
from .combustion import WiebeFunction
from .simulator import BoundaryConditions, EngineParams, EngineSimulator


def make_baseline_engine() -> tuple[EngineParams, BoundaryConditions]:
    geom = EngineGeometry(
        bore=0.086,              # 86 mm
        stroke=0.086,            # 86 mm (square engine)
        connecting_rod_length=0.143,  # 143 mm
        compression_ratio=10.5,
    )

    # Angles are absolute in 0..720 deg, 0 at TDC between exhaust-intake, increasing with rotation
    # Typical cam: IVO ~ 350 deg (10 deg BTDC), IVC ~ 580 deg (40 deg ABDC)
    intake = ValveProfile(
        open_deg=350.0,
        close_deg=580.0,
        max_lift=0.010,          # 10 mm
        valve_diameter=0.034,    # 34 mm
        num_valves=2,
        discharge_coefficient=0.7,
    )

    # EVO ~ 480 deg (60 deg BBDC), EVC ~ 370 deg (10 deg ATDC)
    exhaust = ValveProfile(
        open_deg=480.0,
        close_deg=370.0,  # wrap-around
        max_lift=0.009,   # 9 mm
        valve_diameter=0.029,    # 29 mm
        num_valves=2,
        discharge_coefficient=0.7,
    )

    wiebe = WiebeFunction(
        start_deg=700.0,       # 20 deg BTDC compression to TDC=720
        duration_deg=60.0,     # fast burn
        a=6.0,
        m=2.0,
    )

    params = EngineParams(
        geometry=geom,
        rpm=3000.0,
        intake_valve=intake,
        exhaust_valve=exhaust,
        combustion=wiebe,
        equivalence_ratio=1.0,
        afr_stoich=14.7,
        heat_transfer_coeff_w_m2k=150.0,
        num_cycles=2,  # run two cycles to stabilize
    )

    bc = BoundaryConditions(
        intake_pressure_pa=1.0e5,
        intake_temperature_k=300.0,
        exhaust_pressure_pa=1.1e5,
        exhaust_temperature_k=700.0,
        wall_temperature_k=420.0,
    )

    return params, bc


def run_and_save(params: EngineParams, bc: BoundaryConditions, outdir: str, label: str) -> dict:
    os.makedirs(outdir, exist_ok=True)
    sim = EngineSimulator(params, bc)
    res = sim.run(dtheta_deg=0.5)

    # Write CSV
    csv_path = os.path.join(outdir, f"{label}.csv")
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow([
            "theta_deg",
            "volume_m3",
            "pressure_pa",
            "temperature_k",
            "mass_kg",
            "m_dot_kg_s",
            "qdot_comb_w",
            "qdot_loss_w",
            "intake_area_m2",
            "exhaust_area_m2",
        ])
        for i in range(len(res["theta_deg"])):
            w.writerow([
                res["theta_deg"][i],
                res["volume_m3"][i],
                res["pressure_pa"][i],
                res["temperature_k"][i],
                res["mass_kg"][i],
                res["m_dot_kg_s"][i],
                res["qdot_comb_w"][i],
                res["qdot_loss_w"][i],
                res["intake_area_m2"][i],
                res["exhaust_area_m2"][i],
            ])

    # Write summary
    summary_path = os.path.join(outdir, f"{label}_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"IMEP [Pa]: {res['imep_pa'][0]:.0f}\n")
        f.write(f"Work per cycle [J]: {res['work_j'][0]:.2f}\n")

    return res


def main() -> None:
    outdir = "/workspace/engine_sim/outputs"
    baseline_params, bc = make_baseline_engine()

    # Baseline
    base_res = run_and_save(baseline_params, bc, outdir, label="baseline")

    # Sweeps: vary intake duration and lift
    scenarios: dict[str, dict] = {"baseline": base_res}

    for d_add in [ -20.0, 0.0, 20.0, 40.0 ]:
        p = baseline_params
        new_intake = ValveProfile(
            open_deg=p.intake_valve.open_deg,
            close_deg=(p.intake_valve.close_deg + d_add) % 720.0,
            max_lift=p.intake_valve.max_lift,
            valve_diameter=p.intake_valve.valve_diameter,
            num_valves=p.intake_valve.num_valves,
            discharge_coefficient=p.intake_valve.discharge_coefficient,
        )
        params2 = EngineParams(
            geometry=p.geometry,
            rpm=p.rpm,
            intake_valve=new_intake,
            exhaust_valve=p.exhaust_valve,
            combustion=p.combustion,
            equivalence_ratio=p.equivalence_ratio,
            afr_stoich=p.afr_stoich,
            heat_transfer_coeff_w_m2k=p.heat_transfer_coeff_w_m2k,
            num_cycles=p.num_cycles,
        )
        label = f"int_dur_{int(new_intake.duration())}deg"
        scenarios[label] = run_and_save(params2, bc, outdir, label)

    for lift_factor in [0.8, 1.0, 1.2]:
        p = baseline_params
        new_intake = ValveProfile(
            open_deg=p.intake_valve.open_deg,
            close_deg=p.intake_valve.close_deg,
            max_lift=p.intake_valve.max_lift * lift_factor,
            valve_diameter=p.intake_valve.valve_diameter,
            num_valves=p.intake_valve.num_valves,
            discharge_coefficient=p.intake_valve.discharge_coefficient,
        )
        params2 = EngineParams(
            geometry=p.geometry,
            rpm=p.rpm,
            intake_valve=new_intake,
            exhaust_valve=p.exhaust_valve,
            combustion=p.combustion,
            equivalence_ratio=p.equivalence_ratio,
            afr_stoich=p.afr_stoich,
            heat_transfer_coeff_w_m2k=p.heat_transfer_coeff_w_m2k,
            num_cycles=p.num_cycles,
        )
        label = f"int_lift_{lift_factor:.1f}x"
        scenarios[label] = run_and_save(params2, bc, outdir, label)


if __name__ == "__main__":
    main()