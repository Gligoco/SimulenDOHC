# Engine Combustion Chamber Cycle Simulation (0-D)

This project simulates a single-cylinder four-stroke engine over the full 720° cycle with configurable intake/exhaust valve timing, lift, and opening/closing angles. It uses a zero-dimensional thermodynamic model with compressible valve flows and a Wiebe heat release model.

## Quick start

1. Install dependencies:

```bash
pip install -r /workspace/engine_sim/requirements.txt
```

2. Run the simulation and generate plots:

```bash
python -m engine_sim.main
```

Plots are saved under `/workspace/engine_sim/outputs`:
- `<label>_cycle.png` for each scenario
- `overlay_pressure.png` comparing scenarios

## Notes
- This is a lumped-parameter model; it does not include full CFD or intake/exhaust wave dynamics.
- Valve timing, lift, and angles are set in `engine_sim/main.py` and can be easily modified.