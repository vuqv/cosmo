"""cosmo.translation -- co-translational (protein synthesis) modeling on COSMO.

See ``PLAN.md`` in this folder. Provides the nascent-chain elongation runner
(:mod:`cosmo.translation.elongate`) and the ribosome truncation tool
(:mod:`cosmo.translation.truncate_ribosome`), available from the shell as
``cosmo-elongate`` or ``python -m cosmo.translation``.

Mirrors the sibling ``topo`` project's ``topo.translation`` package, but uses
cosmo's sequence-based IDP force field (HPS / mpipi) for the nascent chain. Build
step v1 (the nascent-only elongation loop) is implemented; build step v2 (the
rigid ribosome scenery) is planned -- see ``PLAN.md`` §6.
"""
# Re-export the high-level runner and its parameters. The CLI entry point (the
# ``elongate`` function) is intentionally NOT re-exported here, so the ``elongate``
# name keeps referring to the submodule (:mod:`cosmo.translation.elongate`).
from .elongate import (run_elongation, ElongationParams,
                       read_elongate_config, ElongateConfig)

__all__ = ["run_elongation", "ElongationParams",
           "read_elongate_config", "ElongateConfig"]
