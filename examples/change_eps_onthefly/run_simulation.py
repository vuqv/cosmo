#!/usr/bin/env python
"""Change the Ashbaugh-Hatch global epsilon before running.

A specialized variant of the canonical runner: it overrides the monomer-monomer
LJ energy scale (the Ashbaugh-Hatch force's global ``epsilon`` parameter) from
the ``LJ_eps`` key in ``md.ini`` before the OpenMM context is created. Everything
else -- config parsing, model build, integrator/platform, reporters, restart,
finalize -- is reused from ``cosmo.read_simulation_config`` and ``cosmo.engine``.

Usage:

    python run_simulation.py -f md.ini
"""
import argparse
import configparser
import time
import warnings

import openmm as mm

try:
    from parmed.exceptions import OpenMMWarning

    warnings.filterwarnings("ignore", category=OpenMMWarning)
except Exception:
    pass

import cosmo
from cosmo import engine


def main():
    parser = argparse.ArgumentParser(description="Run a COSMO simulation with a "
                                                 "custom Ashbaugh-Hatch epsilon.")
    parser.add_argument('-input', '-f', type=str, required=True,
                        help='simulation config file')
    args = parser.parse_args()

    print(f"OpenMM version: {mm.__version__}")

    # Standard parsing for all the usual settings...
    cfg = cosmo.read_simulation_config(args.input)

    # ...plus this example's extra key. read_simulation_config ignores unknown
    # options, so LJ_eps is read here from the same [OPTIONS] section.
    extra = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    extra.read(args.input)
    lj_eps = float(extra['OPTIONS']['LJ_eps'])
    print(f"Overriding Ashbaugh-Hatch global epsilon -> {lj_eps} kJ/mol")

    # Build the model (and write provenance files), then override epsilon before
    # the context is created. No COM remover / barostat are added for this run.
    built = engine.build_model(cfg)
    # The Ashbaugh-Hatch force exposes epsilon as global parameter index 0.
    built.ashbaugh_HatchForce.setGlobalParameterDefaultValue(0, lj_eps)

    print('Simulation started')
    start_time = time.time()
    ctx = engine.setup_simulation(cfg, built, control_file=args.input)
    engine.attach_reporters(cfg, ctx.simulation, append=ctx.restart_active,
                            total_steps=cfg.md_steps)
    ctx.simulation.step(ctx.nsteps_remain)
    engine.finalize_simulation(cfg, ctx, built.topology, start_time)


if __name__ == '__main__':
    main()
