"""``cosmo.csp`` -- the O'Brien Continuous Synthesis Protocol, ported to cosmo.

The per-codon, three-stage protein synthesis protocol of
``continuous_synthesis_v6.py`` (Yang Jiang, Dan Nissley, Ed O'Brien), expressed in
cosmo style on cosmo's **sequence-based IDP force field** (HPS / mpipi). It mirrors the
sibling ``topo`` project's ``topo.csp`` package -- same module layout, class/function
names and CLI -- but drops topo's structure-based Gō machinery (STRIDE, native
contacts, ``domain.yaml`` nscales): a length-``L`` nascent model is just
:func:`cosmo.models.buildCoarseGrainModel` on the first ``L`` residues (all cosmo
forces are sequence-local or pairwise-by-type, so the restriction is exact).

It times every residue from its mRNA codon and splits it into peptidyl-transfer /
translocation / tRNA-binding sub-stages, reusing the shared per-length MD engine
:mod:`cosmo.csp.core`.

CLI::

    cosmo-csp -f csp.ini
    python -m cosmo.csp -f csp.ini

See :mod:`cosmo.csp.protocol` (the explicit-ribosome 3-stage runner + INI),
:mod:`cosmo.csp.cylinder` (the analytic-tunnel *cylinder* runner -- a parallel model
that reuses the same kinetics; ``cosmo-cylinder``), and :mod:`cosmo.csp.kinetics`
(the timing core).
"""
from cosmo.csp.core import RunParams
from cosmo.csp.protocol import (CSPConfig, csp, read_csp_config,
                                run_continuous_synthesis)
from cosmo.csp.cylinder import (CylinderConfig, CylinderParams, cylinder,
                                read_cylinder_config, run_cylinder_synthesis)
from cosmo.csp import kinetics
from cosmo.csp import resume

__all__ = [
    "CSPConfig",
    "RunParams",
    "csp",
    "read_csp_config",
    "run_continuous_synthesis",
    # cylinder ribosome model (analytic exit tunnel)
    "CylinderConfig",
    "CylinderParams",
    "cylinder",
    "read_cylinder_config",
    "run_cylinder_synthesis",
    "kinetics",
    "resume",
]
