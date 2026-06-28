"""COSMO simulation runner.

Exposes :func:`mdrun`, the canonical runner, available from the shell as
``cosmo-mdrun -f md.ini`` or ``python -m cosmo.mdrun -f md.ini``.
"""
from .mdrun import mdrun

__all__ = ["mdrun"]
