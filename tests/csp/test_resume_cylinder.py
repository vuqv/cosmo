"""Driver-level integration tests for CYLINDER resume (cosmo.csp.cylinder).

Mirrors ``test_resume_driver.py`` but for the analytic-tunnel runner: a single MD
segment per residue, a flat ``L_<L>/`` output layout, and no ribosome/PTC geometry.
``run_length`` is stubbed, so this checks the resume *mechanism* (schedule persist/re-read,
progress log, RUNNING-drop, continuation point) deterministically without OpenMM.
``N_full`` is read from a tiny real CA-only PDB fixture (mda.Universe).

Run with ``pytest tests/csp/test_resume_cylinder.py``.
"""
import importlib

import numpy as np
import pytest

from cosmo.csp import resume as r
from cosmo.csp.cylinder import CylinderParams

# `cosmo.csp` re-exports the `cylinder` *function*, which shadows the submodule under
# `from cosmo.csp import cylinder`; import the real module object explicitly so
# monkeypatch.setattr targets the module globals run_cylinder_synthesis actually calls.
cylinder = importlib.import_module("cosmo.csp.cylinder")

N_FULL = 6


@pytest.fixture
def prot_pdb(tmp_path):
    """A tiny CA-only PDB with N_FULL ALA residues (cylinder reads N_full via mda)."""
    p = tmp_path / "prot.pdb"
    with open(p, "w") as fh:
        for i in range(N_FULL):
            fh.write(f"ATOM  {i + 1:>5d}  CA  ALA A{i + 1:>4d}    "
                     f"{i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n")
        fh.write("END\n")
    return str(p)


def _write_final_pdb(path, natoms):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        for i in range(natoms):
            fh.write(f"ATOM  {i + 1:>5d}  CA  ALA A{i + 1:>4d}    "
                     f"{i:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00  0.00           C\n")
        fh.write("END\n")


def _install_stubs(monkeypatch, crash_at_L=None):
    calls = []

    def fake_run_length(L, **kw):
        subdir = kw.get("out_subdir") or f"L_{L:03d}"
        out_root = kw["out_root"]
        if crash_at_L is not None and L == crash_at_L and subdir.startswith("L_"):
            raise RuntimeError(f"simulated crash at L={L}")
        calls.append((L, subdir))
        _write_final_pdb(out_root / subdir / "traj_final.pdb", L)
        return np.zeros((L, 3))

    monkeypatch.setattr(cylinder, "run_length", fake_run_length)
    return calls


def _params(**over):
    p = CylinderParams()
    p.uniform_codon_time = 0.005
    p.random_seed = 42
    p.max_steps_per_stage = 10
    p.device = "CPU"
    for k, v in over.items():
        setattr(p, k, v)
    return p


def _run(pdb, out_dir, params):
    cylinder.run_cylinder_synthesis(
        pdb, L0=1, L_max=5, out_root=str(out_dir), params=params)


def test_golden_run(tmp_path, monkeypatch, prot_pdb):
    _install_stubs(monkeypatch)
    _run(prot_pdb, tmp_path, _params())
    assert (tmp_path / "dwell_times.dat").is_file()
    prog = r.read_progress(tmp_path)
    assert prog.last_done_residue == 5
    for L in range(1, 6):
        assert r.cylinder_final_path(tmp_path, L).is_file()


def test_crash_then_resume_equivalence(tmp_path, monkeypatch, prot_pdb):
    # Golden reference in a separate dir (same seed/config).
    golden = tmp_path / "golden"
    _install_stubs(monkeypatch)
    _run(prot_pdb, golden, _params())
    golden_schedule = (golden / "dwell_times.dat").read_bytes()

    # Interrupted at residue 3 (L_003 left RUNNING).
    work = tmp_path / "work"
    monkeypatch.undo()
    _install_stubs(monkeypatch, crash_at_L=3)
    with pytest.raises(RuntimeError, match="simulated crash"):
        _run(prot_pdb, work, _params())
    prog = r.read_progress(work)
    assert prog.last_done_residue == 2
    assert prog.running_units() == ["L_003"]
    # Schedule persisted before the crash, byte-identical to golden.
    assert (work / "dwell_times.dat").read_bytes() == golden_schedule
    # A partial file in the RUNNING dir must be dropped on resume.
    (work / "L_003").mkdir(parents=True, exist_ok=True)
    (work / "L_003" / "junk.dcd").write_text("partial")

    # Resume (auto): drop L_003, continue to L_max.
    monkeypatch.undo()
    calls = _install_stubs(monkeypatch)
    _run(prot_pdb, work, _params())
    rerun_L = sorted({L for (L, sub) in calls if sub.startswith("L_")})
    assert rerun_L == [3, 4, 5]
    assert not (work / "L_003" / "junk.dcd").exists()
    assert r.read_progress(work).last_done_residue == 5
    # Schedule never rewritten on resume.
    assert (work / "dwell_times.dat").read_bytes() == golden_schedule


def test_presence_guard_aborts(tmp_path, monkeypatch, prot_pdb):
    _install_stubs(monkeypatch)
    _run(prot_pdb, tmp_path, _params())
    r.cylinder_final_path(tmp_path, 2).unlink()   # simulate a scratch purge of L_002
    monkeypatch.undo()
    _install_stubs(monkeypatch)
    with pytest.raises(SystemExit, match="L_002"):
        _run(prot_pdb, tmp_path, _params(resume="yes"))


def test_resume_completed_is_noop(tmp_path, monkeypatch, prot_pdb):
    _install_stubs(monkeypatch)
    _run(prot_pdb, tmp_path, _params())
    monkeypatch.undo()
    calls = _install_stubs(monkeypatch)
    _run(prot_pdb, tmp_path, _params())   # resume=auto, all DONE
    assert calls == []                    # run_length never called


def test_ini_resume_parsed(tmp_path):
    ini = tmp_path / "cylinder.ini"
    ini.write_text("[OPTIONS]\n"
                   "pdb_file = p.pdb\n"
                   "L0 = 1\n"
                   "codon_times = 0.005\n"
                   "resume = no\n")
    cfg = cylinder.read_cylinder_config(str(ini), verbose=False)
    assert cfg.params.resume == "no"

    ini.write_text("[OPTIONS]\n"
                   "pdb_file = p.pdb\n"
                   "L0 = 1\n"
                   "codon_times = 0.005\n"
                   "resume = bogus\n")
    with pytest.raises(ValueError, match="resume"):
        cylinder.read_cylinder_config(str(ini), verbose=False)
