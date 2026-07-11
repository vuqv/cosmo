"""Layout-discovery tests for cosmo.csp.movie under the consolidated CSP layout.

Deterministic and MD-free: synthetic ``L_<L>/`` trees with marker files, checking the
per-stage segment discovery (no trajectory parsing, which happens later in the stitch).

Run with ``pytest tests/csp/test_movie_layout.py``.
"""
from cosmo.csp import movie


def _make_residue(root, L, stages=(1, 2, 3), final=True, empty=()):
    d = root / f"L_{L:03d}"
    d.mkdir()
    (d / "traj.psf").write_text("psf")
    (d / f"native_1_{L}.pdb").write_text("native")
    for s in stages:
        (d / f"traj_s{s}.dcd").write_text("" if s in empty else "dcd-bytes")
    if final:
        (d / "traj_final.pdb").write_text("ATOM\n")
    return d


def test_find_segments_consolidated(tmp_path):
    _make_residue(tmp_path, 1)
    _make_residue(tmp_path, 2)
    segs = movie.find_segments(str(tmp_path))
    assert [s[0] for s in segs] == ["L=1 s1", "L=1 s2", "L=1 s3",
                                    "L=2 s1", "L=2 s2", "L=2 s3"]
    for label, L, psf, traj in segs:
        assert psf.endswith("traj.psf")
        assert "traj_s" in traj
        assert L in (1, 2)


def test_stage3_empty_dcd_falls_back_to_final(tmp_path):
    # s1 has frames; s2 empty -> dropped; s3 empty -> falls back to the single traj_final.pdb.
    _make_residue(tmp_path, 1, stages=(1, 2, 3), final=True, empty=(2, 3))
    labels = [s[0] for s in movie.find_segments(str(tmp_path))]
    assert labels == ["L=1 s1", "L=1 s3"]


def test_zero_byte_dcds_skipped_quietly(tmp_path):
    # stages 1/2 ran fewer steps than nstout -> 0-byte DCDs; skipped, s3 via final.
    _make_residue(tmp_path, 7, stages=(1, 2, 3), final=True, empty=(1, 2))
    segs = movie.find_segments(str(tmp_path))
    assert [s[0] for s in segs] == ["L=7 s3"]


def test_flat_layout_single_segment(tmp_path):
    # Cylinder flat layout: L_<L>/traj.psf + traj.dcd -> one segment per length.
    d = tmp_path / "L_001"
    d.mkdir()
    (d / "traj.psf").write_text("psf")
    (d / "traj.dcd").write_text("dcd-bytes")
    segs = movie.find_segments(str(tmp_path))
    assert [s[0] for s in segs] == ["L=1"]


def test_residue_without_psf_skipped(tmp_path):
    d = tmp_path / "L_001"
    d.mkdir()
    (d / "traj_s1.dcd").write_text("d")  # no traj.psf -> falls to flat branch, no traj.dcd
    assert movie.find_segments(str(tmp_path)) == []
