"""Driver-level integration tests for cosmo CSP resume + consolidated layout.

Exercises the real :func:`cosmo.csp.protocol.run_continuous_synthesis` control flow --
fresh, mid-run crash, resume, presence guard, idempotent completed run -- with the heavy
MD primitives (``run_length``, ``optimal_ptc_targets``, ribosome load, anchors) stubbed.
``N_full`` is read from a tiny real CA-only PDB fixture (mda.Universe).

Run with ``pytest tests/csp/test_resume_driver.py``.
"""
import types

import numpy as np
import pytest

from cosmo.csp import protocol
from cosmo.csp import resume as r
from cosmo.csp.core import RunParams

N_FULL = 6


@pytest.fixture
def prot_pdb(tmp_path):
    """A tiny CA-only PDB with N_FULL ALA residues (mda reads N_full from it)."""
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

    monkeypatch.setattr(protocol, "load_ribosome",
                        lambda *a, **k: types.SimpleNamespace(n=10))
    monkeypatch.setattr(protocol, "read_anchor", lambda *a, **k: np.zeros(3))
    monkeypatch.setattr(protocol, "optimal_ptc_targets",
                        lambda ribo, **k: (np.array([1.0, 0.0, 0.0]),
                                           np.array([1.38, 0.0, 0.0])))

    def fake_run_length(L, **kw):
        subdir = kw["out_subdir"]
        out_root = kw["out_root"]
        outname = kw.get("outname", "traj")
        if crash_at_L is not None and L == crash_at_L and outname == "traj_s1":
            raise RuntimeError(f"simulated crash at L={L}")
        calls.append((L, subdir))
        if kw.get("persist_final", True):
            _write_final_pdb(out_root / subdir / "traj_final.pdb", L)
        return np.zeros((L, 3))

    monkeypatch.setattr(protocol, "run_length", fake_run_length)
    return calls


def _params(**over):
    p = RunParams()
    p.uniform_codon_time = 0.005
    p.random_seed = 42
    p.max_steps_per_stage = 10
    p.device = "CPU"
    for k, v in over.items():
        setattr(p, k, v)
    return p


def _run(pdb, out_dir, params):
    protocol.run_continuous_synthesis(
        pdb, "ribo.pdb", L0=1, L_max=5, out_root=str(out_dir), params=params)


def test_golden_run(tmp_path, monkeypatch, prot_pdb):
    _install_stubs(monkeypatch)
    _run(prot_pdb, tmp_path, _params())
    assert (tmp_path / "dwell_times.dat").is_file()
    prog = r.read_progress(tmp_path)
    assert prog.last_done_residue == 5
    for L in range(1, 6):
        assert r.residue_final_path(tmp_path, L).is_file()


def test_crash_then_resume_equivalence(tmp_path, monkeypatch, prot_pdb):
    golden = tmp_path / "golden"
    _install_stubs(monkeypatch)
    _run(prot_pdb, golden, _params())
    golden_schedule = (golden / "dwell_times.dat").read_bytes()

    work = tmp_path / "work"
    monkeypatch.undo()
    _install_stubs(monkeypatch, crash_at_L=3)
    with pytest.raises(RuntimeError, match="simulated crash"):
        _run(prot_pdb, work, _params())
    prog = r.read_progress(work)
    assert prog.last_done_residue == 2
    assert prog.running_units() == ["L_003"]
    assert (work / "dwell_times.dat").read_bytes() == golden_schedule
    (work / "L_003").mkdir(parents=True, exist_ok=True)
    (work / "L_003" / "junk.dcd").write_text("partial")

    monkeypatch.undo()
    calls = _install_stubs(monkeypatch)
    _run(prot_pdb, work, _params())
    rerun_L = sorted({L for (L, sub) in calls if sub.startswith("L_")})
    assert rerun_L == [3, 4, 5]
    assert not (work / "L_003" / "junk.dcd").exists()
    assert r.read_progress(work).last_done_residue == 5
    assert (work / "dwell_times.dat").read_bytes() == golden_schedule


def test_presence_guard_aborts(tmp_path, monkeypatch, prot_pdb):
    _install_stubs(monkeypatch)
    _run(prot_pdb, tmp_path, _params())
    r.residue_final_path(tmp_path, 2).unlink()
    monkeypatch.undo()
    _install_stubs(monkeypatch)
    with pytest.raises(SystemExit, match="L_002"):
        _run(prot_pdb, tmp_path, _params(resume="yes"))


def test_resume_completed_is_noop(tmp_path, monkeypatch, prot_pdb):
    _install_stubs(monkeypatch)
    _run(prot_pdb, tmp_path, _params())
    monkeypatch.undo()
    calls = _install_stubs(monkeypatch)
    _run(prot_pdb, tmp_path, _params())
    assert calls == []


def test_resume_yes_without_run_errors(tmp_path, monkeypatch, prot_pdb):
    _install_stubs(monkeypatch)
    with pytest.raises(SystemExit, match="nothing to resume"):
        _run(prot_pdb, tmp_path, _params(resume="yes"))


# --------------------------------------------------------------------------
# INI parsing of the resume key
# --------------------------------------------------------------------------
_BASE_INI = ("[OPTIONS]\n"
             "pdb_file = p.pdb\n"
             "ribosome = r.pdb\n"
             "codon_times = 0.005\n")


def _write_ini(tmp_path, extra=""):
    path = tmp_path / "csp.ini"
    path.write_text(_BASE_INI + extra)
    return str(path)


def test_ini_resume_default_auto(tmp_path):
    cfg = protocol.read_csp_config(_write_ini(tmp_path), verbose=False)
    assert cfg.params.resume == "auto"


@pytest.mark.parametrize("val", ["yes", "no", "auto", "YES"])
def test_ini_resume_valid(tmp_path, val):
    cfg = protocol.read_csp_config(_write_ini(tmp_path, f"resume = {val}\n"), verbose=False)
    assert cfg.params.resume == val.lower()


def test_ini_resume_invalid(tmp_path):
    with pytest.raises(ValueError, match="resume"):
        protocol.read_csp_config(_write_ini(tmp_path, "resume = maybe\n"), verbose=False)
