"""Deterministic unit tests for cosmo.csp.resume (the CSP resume mechanism).

Covers the parts that need no MD: schedule round-trip (write->read reproduces step
counts + PTC geometry), the cylinder single-segment schedule, the progress-log parse,
the presence guard, the RUNNING-drop, and coverage.

Run with ``pytest tests/csp/test_resume.py``.
"""
import numpy as np
import pytest

from cosmo.csp import resume as r
from cosmo.csp import kinetics
from cosmo.csp.core import RunParams


def _params():
    p = RunParams()
    p.scale_factor = 4331293.0
    p.dt_ps = 0.01
    p.time_stage_1 = 0.00034
    p.time_stage_2 = 0.004201
    p.random_seed = 12345
    p.uniform_codon_time = None
    return p


# --------------------------------------------------------------------------
# 3-stage schedule table round-trip
# --------------------------------------------------------------------------
def test_schedule_round_trip(tmp_path):
    rows = [
        r.SchedRow(L=1, codon="AUG", t_total=1.23456e-3,
                   times=(3.4e-4, 4.2e-3, 7.7e-3), steps=(10, 200, 3)),
        r.SchedRow(L=2, codon="GCU", t_total=9.87654e-3,
                   times=(3.1e-4, 5.0e-3, 1.1e-2), steps=(1, 2, 3)),
        r.SchedRow(L=3, codon="uniform", t_total=5.0e-3,
                   times=(3.4e-4, 4.2e-3, 4.6e-4), steps=(12345, 6789, 1)),
    ]
    a_target = np.array([1.123456789012345, 2.0, -3.5])
    p_target = np.array([1.503456789012345, 2.0, -3.5])
    wall_x = 1.123456789012345
    path = tmp_path / "dwell_times.dat"
    r.write_schedule(path, rows, _params(), a_target, p_target, wall_x)

    rows2, a2, p2, wall2 = r.read_schedule(path)
    assert [x.L for x in rows2] == [1, 2, 3]
    assert [x.codon for x in rows2] == ["AUG", "GCU", "uniform"]
    assert [x.steps for x in rows2] == [(10, 200, 3), (1, 2, 3), (12345, 6789, 1)]
    np.testing.assert_array_equal(a2, a_target)
    np.testing.assert_array_equal(p2, p_target)
    assert wall2 == wall_x
    for orig, got in zip(rows, rows2):
        np.testing.assert_allclose(got.times, orig.times, rtol=1e-6)


def test_schedule_round_trip_wall_none(tmp_path):
    rows = [r.SchedRow(L=5, codon="AAA", t_total=1e-3,
                       times=(1e-4, 2e-4, 3e-4), steps=(1, 1, 1))]
    path = tmp_path / "dwell_times.dat"
    r.write_schedule(path, rows, _params(), np.zeros(3), np.ones(3), None)
    _, _, _, wall2 = r.read_schedule(path)
    assert wall2 is None


def test_read_schedule_rejects_no_ptc(tmp_path):
    path = tmp_path / "legacy.dat"
    path.write_text("# no PTC header\n"
                    "   1  AUG  1.0e-3  1e-4  2e-4  3e-4  0  0  0  1  2  3\n")
    with pytest.raises(ValueError):
        r.read_schedule(path)


def test_cylinder_schedule_round_trip(tmp_path):
    rows = [
        r.CylSchedRow(L=1, codon="AUG", dwell_s=1.23456e-3, steps=137),
        r.CylSchedRow(L=2, codon="GCU", dwell_s=9.87654e-3, steps=1),
        r.CylSchedRow(L=3, codon="uniform", dwell_s=5.0e-3, steps=99999),
    ]
    path = tmp_path / "dwell_times.dat"
    r.write_cylinder_schedule(path, rows, _params())
    rows2 = r.read_cylinder_schedule(path)
    assert [x.L for x in rows2] == [1, 2, 3]
    assert [x.steps for x in rows2] == [137, 1, 99999]
    r.schedule_covers(rows2, 1, 3)
    with pytest.raises(SystemExit):
        r.schedule_covers(rows2, 1, 4)


# --------------------------------------------------------------------------
# schedule_covers
# --------------------------------------------------------------------------
def test_schedule_covers_ok():
    rows = [r.SchedRow(L, "X", 0.0, (0, 0, 0), (1, 1, 1)) for L in range(3, 8)]
    r.schedule_covers(rows, 3, 7)   # no raise


@pytest.mark.parametrize("L0,L_max", [(3, 8), (2, 7), (1, 10)])
def test_schedule_covers_mismatch(L0, L_max):
    rows = [r.SchedRow(L, "X", 0.0, (0, 0, 0), (1, 1, 1)) for L in range(3, 8)]
    with pytest.raises(SystemExit):
        r.schedule_covers(rows, L0, L_max)


# --------------------------------------------------------------------------
# progress.log
# --------------------------------------------------------------------------
def test_progress_parse_last_status(tmp_path):
    r.write_progress_header(tmp_path)
    for L in (1, 2, 3):
        r.append_progress(tmp_path, f"L_{L:03d}", "RUNNING")
        r.append_progress(tmp_path, f"L_{L:03d}", "DONE")
    r.append_progress(tmp_path, "L_004", "RUNNING")
    prog = r.read_progress(tmp_path)
    assert prog.last_done_residue == 3
    assert prog.is_done("L_002")
    assert not prog.is_done("L_004")
    assert prog.running_units() == ["L_004"]


def test_progress_last_status_wins(tmp_path):
    """A unit that goes RUNNING then DONE reads as DONE (last status wins)."""
    r.write_progress_header(tmp_path)
    r.append_progress(tmp_path, "ejection", "RUNNING")
    r.append_progress(tmp_path, "ejection", "DONE")
    prog = r.read_progress(tmp_path)
    assert prog.is_done("ejection")
    assert prog.running_units() == []


def test_progress_empty(tmp_path):
    r.write_progress_header(tmp_path)
    prog = r.read_progress(tmp_path)
    assert prog.last_done_residue == 0
    assert prog.running_units() == []


# --------------------------------------------------------------------------
# verify / drop / load
# --------------------------------------------------------------------------
def _touch_final(out_root, L):
    fp = r.residue_final_path(out_root, L)
    fp.parent.mkdir(parents=True, exist_ok=True)
    fp.write_text("ATOM\n")


def test_verify_completed_units_ok(tmp_path):
    r.write_progress_header(tmp_path)
    for L in (1, 2, 3):
        r.append_progress(tmp_path, f"L_{L:03d}", "DONE")
        _touch_final(tmp_path, L)
    prog = r.read_progress(tmp_path)
    r.verify_completed_units(tmp_path, prog, L0=1)   # no raise


def test_verify_completed_units_missing_hole(tmp_path):
    r.write_progress_header(tmp_path)
    for L in (1, 2, 3):
        r.append_progress(tmp_path, f"L_{L:03d}", "DONE")
        _touch_final(tmp_path, L)
    r.residue_final_path(tmp_path, 2).unlink()
    prog = r.read_progress(tmp_path)
    with pytest.raises(SystemExit, match="L_002"):
        r.verify_completed_units(tmp_path, prog, L0=1)


def test_drop_running_units(tmp_path):
    r.write_progress_header(tmp_path)
    for L in (1, 2):
        r.append_progress(tmp_path, f"L_{L:03d}", "DONE")
        _touch_final(tmp_path, L)
    r.append_progress(tmp_path, "L_003", "RUNNING")
    (tmp_path / "L_003").mkdir()
    (tmp_path / "L_003" / "traj_s1.dcd").write_text("partial")
    prog = r.read_progress(tmp_path)
    assert r.drop_running_units(tmp_path, prog) == ["L_003"]
    assert not (tmp_path / "L_003").exists()
    assert r.residue_final_path(tmp_path, 1).is_file()


def test_drop_running_units_none(tmp_path):
    """No RUNNING unit (clean finish) -> nothing dropped."""
    r.write_progress_header(tmp_path)
    r.append_progress(tmp_path, "L_001", "DONE")
    prog = r.read_progress(tmp_path)
    assert r.drop_running_units(tmp_path, prog) == []


def test_load_final_pdb(tmp_path):
    pdb = tmp_path / "traj_final.pdb"
    pdb.write_text(
        "ATOM      1  CA  ALA A   1      10.000  20.000  30.000  1.00  0.00           C\n"
        "ATOM      2  CA  GLY A   2      11.000  21.000  31.000  1.00  0.00           C\n"
        "END\n")
    coords = r.load_final_pdb(pdb)
    assert coords.shape == (2, 3)
    np.testing.assert_allclose(coords, [[1.0, 2.0, 3.0], [1.1, 2.1, 3.1]], atol=1e-6)


def test_load_final_pdb_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        r.load_final_pdb(tmp_path / "nope.pdb")


def test_cylinder_final_path_layout(tmp_path):
    assert r.cylinder_final_path(tmp_path, 7) == tmp_path / "L_007" / "traj_final.pdb"
    assert r.residue_final_path(tmp_path, 7) == tmp_path / "L_007" / "traj_final.pdb"


# --------------------------------------------------------------------------
# Pre-loop schedule draw == sequential (legacy in-loop) draw
# --------------------------------------------------------------------------
def test_preloop_schedule_equals_sequential_draw():
    import random
    intrinsic = [0.005] * 12
    real = list(intrinsic)
    kw = dict(time_stage_1=0.00034, time_stage_2=0.004201,
              scale_factor=4331293.0, dt_ps=0.01,
              max_steps_per_stage=50, min_steps_per_stage=1)
    L0, L_max, seed = 2, 10, 777
    rng_a = random.Random(seed)
    pre = [kinetics.stage_steps(L, intrinsic, real, rng=rng_a, **kw)
           for L in range(L0, L_max + 1)]
    rng_b = random.Random(seed)
    legacy = [kinetics.stage_steps(L, intrinsic, real, rng=rng_b, **kw)
              for L in range(L0, L_max + 1)]
    assert pre == legacy
