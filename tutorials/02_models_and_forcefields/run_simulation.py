#!/usr/bin/env python
"""Run this example with COSMO's canonical runner.

All of these read the same md.ini and are equivalent:

    cosmo-mdrun -f md.ini            # console command (after `pip install -e .`)
    python -m cosmo.mdrun -f md.ini  # module form, no install needed
    python run_simulation.py -f md.ini

The boilerplate that used to live here (config parsing, model build, OpenMM
setup, reporters) now lives in ``cosmo.read_simulation_config`` + ``cosmo.engine``
and is orchestrated by ``cosmo.mdrun.mdrun``.
"""
from cosmo.mdrun import mdrun

if __name__ == "__main__":
    mdrun()
