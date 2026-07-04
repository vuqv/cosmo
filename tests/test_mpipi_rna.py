"""
Tests for RNA support in the Mpipi model.

These pin the implementation against the paper (Joseph et al., Nat. Comput. Sci.
2021) and its SI Table 12:

* RNA residue typing: A/C/G/U map to ids 20-23, charge -0.75e, standard masses,
  and common aliases (RA/RU/RG/RC, ADE/URA/GUA/CYT) resolve to the same beads.
* Wang-Frenkel tables are extended to 24x24 with mu = 3 for every RNA-involving
  pair, nu = 1, and R_ij = 3*sigma_ij everywhere.
* The OpenMM Wang-Frenkel energy of an isolated U-U pair reproduces the analytic
  Eq. (4) to machine precision (verifies the tabulated-function wiring).
* A protein+RNA system builds and evaluates to finite energies.

Run with:  python -m pytest tests/test_mpipi_rna.py
       or:  python tests/test_mpipi_rna.py
"""
import os
import sys

import numpy as np
import openmm as mm
from openmm import unit

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import cosmo  # noqa: E402
from cosmo.parameters import model_parameters as MP  # noqa: E402

RNA_PDB = os.path.join(_HERE, "data", "rna.pdb")
COMPLEX_PDB = os.path.join(_HERE, "data", "combined.pdb")


def _build(pdb):
    return cosmo.models.buildCoarseGrainModel(
        pdb, minimize=False, model="mpipi", box_dimension=None)


def test_rna_residue_typing():
    p = MP.parameters["mpipi"]
    expected = {"A": (20, 329.2), "C": (21, 305.2), "G": (22, 345.2), "U": (23, 306.2)}
    for name, (rid, mass) in expected.items():
        assert p[name]["id"] == rid
        assert p[name]["mass"] == mass
        assert p[name]["charge"] == -0.75   # Mpipi RNA bead charge


def test_rna_aliases_resolve():
    p = MP.parameters["mpipi"]
    for alias, canon in {"RA": "A", "RU": "U", "RG": "G", "RC": "C",
                         "ADE": "A", "URA": "U", "GUA": "G", "CYT": "C"}.items():
        assert p[alias]["id"] == p[canon]["id"]
        assert p[alias]["charge"] == -0.75


def test_wf_tables_shape_and_rules():
    p = MP.parameters["mpipi"]
    eps, sig, mu, nu, rc = (p["eps_ij"], p["sigma_ij"], p["mu_ij"],
                            p["nu_ij"], p["rc_ij"])
    assert eps.shape == sig.shape == mu.shape == nu.shape == rc.shape == (24, 24)
    # symmetry and the paper's universal rules
    for M in (eps, sig, mu):
        assert np.allclose(M, M.T)
    assert np.all(nu == 1.0)                 # nu = 1 for every pair
    assert np.allclose(rc, 3.0 * sig)        # R_ij = 3 sigma_ij
    # mu = 3 for every RNA-involving pair (Table 12); protein block keeps Table 11
    assert np.all(mu[20:, :] == 3.0) and np.all(mu[:, 20:] == 3.0)
    # protein-protein special exponents preserved (V-I = 4, I-I = 11)
    order = "M G K T R A D E Y V L Q W F S H N P C I".split()
    iV, iI = order.index("V"), order.index("I")
    assert mu[iV, iI] == 4.0 and mu[iI, iI] == 11.0


def test_uu_wang_frenkel_matches_analytic(tmp_path):
    """OpenMM WF energy of an isolated U-U pair == analytic Eq. (4)."""
    pdb = tmp_path / "uu.pdb"
    pdb.write_text(
        "ATOM      1  P     U A   1       0.000   0.000   0.000  1.00  0.00           P\n"
        "ATOM      2  P     U B   1       9.000   0.000   0.000  1.00  0.00           P\n"
        "END\n"
    )
    model = _build(str(pdb))
    ctx = mm.Context(model.system, mm.VerletIntegrator(0.001),
                     mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(model.positions)
    groups = {n: i for i, n in enumerate(model.forceGroups)}
    wf_name = next(n for n in groups if "PairW" in n)  # "PairWire Energy"
    wf = ctx.getState(getEnergy=True, groups={groups[wf_name]}).getPotentialEnergy()
    wf = wf.value_in_unit(unit.kilojoules_per_mole)

    # analytic Wang-Frenkel, Eq. (4)-(5), for the U-U pair (id 23)
    p = MP.parameters["mpipi"]
    i = p["U"]["id"]
    eps, sig, mu, nu, rc = (p["eps_ij"][i, i], p["sigma_ij"][i, i], p["mu_ij"][i, i],
                            p["nu_ij"][i, i], p["rc_ij"][i, i])
    r = 0.9
    alpha = 2 * nu * (rc / sig) ** (2 * mu) * (
        (2 * nu + 1) / (2 * nu * ((rc / sig) ** (2 * mu) - 1))) ** (2 * nu + 1)
    phi = eps * alpha * ((sig / r) ** (2 * mu) - 1) * ((rc / r) ** (2 * mu) - 1) ** (2 * nu)
    assert np.isclose(wf, phi, rtol=1e-7, atol=1e-7)


def test_protein_rna_complex_builds():
    """A protein+RNA system builds and gives finite energies with RNA beads."""
    if not os.path.exists(COMPLEX_PDB):
        print(f"SKIP: {COMPLEX_PDB} not present")
        return
    model = _build(COMPLEX_PDB)
    ids = list(model.particle_type_id)
    assert any(i >= 20 for i in ids)         # at least one RNA bead present
    ctx = mm.Context(model.system, mm.VerletIntegrator(0.001),
                     mm.Platform.getPlatformByName("Reference"))
    ctx.setPositions(model.positions)
    total = ctx.getState(getEnergy=True).getPotentialEnergy().value_in_unit(
        unit.kilojoules_per_mole)
    assert np.isfinite(total)


def test_mpipi_debye_length_is_model_specific():
    assert MP.debye_length["mpipi"] == 0.795           # Joseph et al. 2021
    assert MP.debye_length.get("hps_urry", MP.DEFAULT_DEBYE_LENGTH) == 1.0


def _run_all():
    import tempfile
    import pathlib
    test_rna_residue_typing();            print("ok  test_rna_residue_typing")
    test_rna_aliases_resolve();           print("ok  test_rna_aliases_resolve")
    test_wf_tables_shape_and_rules();      print("ok  test_wf_tables_shape_and_rules")
    with tempfile.TemporaryDirectory() as d:
        test_uu_wang_frenkel_matches_analytic(pathlib.Path(d))
    print("ok  test_uu_wang_frenkel_matches_analytic")
    test_protein_rna_complex_builds();    print("ok  test_protein_rna_complex_builds")
    test_mpipi_debye_length_is_model_specific(); print("ok  test_mpipi_debye_length_is_model_specific")
    print("\nALL MPIPI-RNA TESTS PASSED")


if __name__ == "__main__":
    _run_all()
