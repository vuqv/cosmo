"""
Run a COSMO coarse-grained simulation from a control file (md.ini).

This is the canonical runner for the package. Use it as a CLI::

    cosmo-mdrun -f md.ini
    python -m cosmo.mdrun -f md.ini

or call :func:`mdrun` from your own script. Control-file parsing lives in
:func:`cosmo.read_simulation_config`; the build / setup / finalize steps live in
:mod:`cosmo.engine`. This runner is intentionally thin -- it just wires those
pieces together, replacing the per-example ``run_simulation.py`` copies and the
old ``cosmo.dynamics.Dynamics`` class.
"""
import argparse
import sys
import time
import warnings

import openmm as mm

try:
    from parmed.exceptions import OpenMMWarning

    warnings.filterwarnings("ignore", category=OpenMMWarning)
except Exception:  # pragma: no cover - parmed ships alongside cosmo
    pass

from cosmo import engine
from cosmo.utils.config import read_simulation_config


def mdrun():
    """Run a simulation from parameters specified in a control file.

    Usage: ``cosmo-mdrun -f md.ini`` (or ``python -m cosmo.mdrun -f md.ini``).
    """
    parser = argparse.ArgumentParser(
        prog="cosmo-mdrun",
        description="Run a COSMO coarse-grained simulation from a control file "
                    "(md.ini).")
    parser.add_argument('-input', '-f', type=str, help='simulation config file')
    # A bare `cosmo-mdrun` (no arguments) prints help, like `-h`.
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    args = parser.parse_args()

    print('-' * 70)
    print("COSMO: COarse-grained Simulation of intrinsically disordered prOteins "
          "with openMM")
    print(f"OpenMM installed version: {mm.__version__}")

    # 1. Parse the control file.
    cfg = read_simulation_config(args.input)

    # 2. Build the coarse-grained system (+ provenance dumps).
    built = engine.build_system(cfg)

    # 3. Set up the OpenMM Simulation (integrator, platform, positions/restart).
    print('Simulation started')
    start_time = time.time()
    ctx = engine.setup_simulation(cfg, built, control_file=args.input)

    # 4. Attach reporters and run.
    engine.attach_reporters(cfg, ctx.simulation, append=ctx.restart_active,
                            total_steps=cfg.md_steps)
    ctx.simulation.step(ctx.nsteps_remain)

    # 5. Save the final conformation + checkpoint and report timing.
    engine.finalize_simulation(cfg, ctx, built.topology, start_time)


if __name__ == '__main__':
    mdrun()
