#!/usr/bin/env python
"""Slab / NPT simulation example.

This is a plain NPT run: temperature coupling plus a Monte Carlo barostat. The
barostat is added automatically by the canonical runner when ``md.ini`` sets
``pcoupl = yes`` (with ``ref_p`` / ``frequency_p``), so no custom code is needed
here -- ``cosmo.engine`` builds the system, adds the barostat, runs and finalizes.

Equivalent ways to launch (all read the same md.ini):

    cosmo-mdrun -f md.ini
    python -m cosmo.mdrun -f md.ini
    python run_simulation.py -f md.ini
"""
from cosmo.mdrun import mdrun

if __name__ == "__main__":
    mdrun()
