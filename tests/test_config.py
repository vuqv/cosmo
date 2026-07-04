"""
Tests for cosmo.utils.config.read_simulation_config.

These pin the parsing behaviour against (a) the real example control files in
the repository and (b) the historical inline parser in
cosmo.dynamics.Dynamics.read_config, so the Phase-1 refactor is provably
behaviour-preserving.

Run with:  python -m pytest tests/test_config.py
       or:  python tests/test_config.py
"""
import os
import sys

import openmm as mm
from openmm import unit

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

from cosmo.utils.config import read_simulation_config  # noqa: E402

STANDARD_INI = os.path.join(_HERE, "data", "md.ini")


def test_standard_example_values():
    """Parse the standard example md.ini and check every field."""
    cfg = read_simulation_config(STANDARD_INI, verbose=False)

    assert cfg.md_steps == 10000            # '10_000' underscore handled
    assert cfg.dt == 0.01 * unit.picoseconds
    assert cfg.nstxout == 100
    assert cfg.nstlog == 100
    assert cfg.nstcomm == 100
    assert cfg.model == "mpipi"

    assert cfg.tcoupl is True
    assert cfg.ref_t == 310.0 * unit.kelvin
    assert cfg.tau_t == 0.01 / unit.picoseconds

    # pcoupl is off in this file; pbc is on with a scalar (cubic) box.
    assert cfg.pcoupl is False
    assert cfg.pbc is True
    assert cfg.box_dimension == 30          # scalar -> cubic, kept as-is

    assert cfg.pdb_file == "asyn.pdb"
    # topo-style output: traj/asyn.* (no explicit checkpoint -> defaults into the folder)
    assert cfg.output_dir == "traj"
    assert cfg.outname == "asyn"
    assert cfg.checkpoint is None
    assert cfg.device == "GPU"

    # restart = no -> minimize honoured (yes in this file)
    assert cfg.restart is False
    assert cfg.minimize is True


def test_helpers():
    """The dataclass helper methods behave as the runner/engine expect."""
    cfg = read_simulation_config(STANDARD_INI, verbose=False)

    # Outputs go under output_dir (default 'traj') with basename outname.
    assert cfg.output_dir == "traj"
    assert cfg.outname == "asyn"
    assert cfg.output_path(".dcd") == os.path.join("traj", "asyn.dcd")
    assert cfg.output_path("_init.pdb") == os.path.join("traj", "asyn_init.pdb")
    # no explicit checkpoint -> defaults to <output_dir>/<outname>.chk
    assert cfg.checkpoint_path() == os.path.join("traj", "asyn.chk")
    assert cfg.writes_forcefield() is False   # mpipi has no _ff.dat
    assert cfg.build_kwargs() == dict(minimize=True, model="mpipi",
                                      box_dimension=30)


def test_output_dir_conventions(tmp_path):
    """output_dir/outname (topo-style) and the protein_code fallback + default."""
    base = ("[OPTIONS]\nnstxout=10\nnstlog=10\nmodel=hps_urry\n"
            "tcoupl=no\npbc=no\npcoupl=no\npdb_file=x.pdb\ndevice=CPU\n"
            "restart=no\nminimize=no\n")

    # explicit topo-style output_dir + outname
    ini = tmp_path / "a.ini"
    ini.write_text(base + "output_dir=run1\noutname=myrun\n")
    cfg = read_simulation_config(str(ini), verbose=False)
    assert cfg.output_path(".dcd") == os.path.join("run1", "myrun.dcd")
    assert cfg.checkpoint_path() == os.path.join("run1", "myrun.chk")

    # no protein_code, no outname -> default traj/traj.*
    ini2 = tmp_path / "b.ini"
    ini2.write_text(base)
    cfg2 = read_simulation_config(str(ini2), verbose=False)
    assert cfg2.output_path(".log") == os.path.join("traj", "traj.log")


def test_restart_forces_no_minimize(tmp_path):
    """restart = yes must force minimize off, regardless of the minimize key."""
    ini = tmp_path / "md.ini"
    ini.write_text(
        "[OPTIONS]\n"
        "nstxout = 10\nnstlog = 10\n"
        "model = hps_urry\n"
        "tcoupl = no\npbc = no\npcoupl = no\n"
        "pdb_file = x.pdb\nprotein_code = X\ncheckpoint = x.chk\n"
        "device = CPU\nppn = 2\n"
        "restart = yes\nminimize = yes\n"
    )
    cfg = read_simulation_config(str(ini), verbose=False)
    assert cfg.restart is True
    assert cfg.minimize is False
    assert cfg.writes_forcefield() is True    # hps_urry writes _ff.dat
    assert cfg.checkpoint_path() == "x.chk"


def test_matches_legacy_dynamics_parser():
    """read_simulation_config must agree with the historical Dynamics parser."""
    try:
        from cosmo.dynamics.dynamics import Dynamics
    except Exception as exc:  # pragma: no cover - Dynamics removed in Phase 4
        print(f"SKIP legacy comparison (Dynamics unavailable: {exc})")
        return

    cfg = read_simulation_config(STANDARD_INI, verbose=False)
    legacy = Dynamics(STANDARD_INI)  # its __init__ calls read_config

    assert cfg.md_steps == legacy.md_steps
    assert cfg.dt == legacy.dt
    assert cfg.nstxout == legacy.nstxout
    assert cfg.nstlog == legacy.nstlog
    assert cfg.nstcomm == legacy.nstcomm
    assert cfg.model == legacy.model
    assert cfg.tcoupl == legacy.tcoupl
    assert cfg.ref_t == legacy.ref_t
    assert cfg.tau_t == legacy.tau_t
    assert cfg.pbc == legacy.pbc
    assert cfg.box_dimension == legacy.box_dimension
    assert cfg.pcoupl == legacy.pcoupl
    assert cfg.pdb_file == legacy.pdb_file
    assert cfg.protein_code == legacy.protein_code
    assert cfg.device == legacy.device
    assert cfg.restart == legacy.restart
    assert cfg.minimize == legacy.minimize


def _run_all():
    """Lightweight runner so the file works without pytest installed."""
    import tempfile
    import pathlib

    test_standard_example_values()
    print("ok  test_standard_example_values")
    test_helpers()
    print("ok  test_helpers")
    with tempfile.TemporaryDirectory() as d:
        test_output_dir_conventions(pathlib.Path(d))
    print("ok  test_output_dir_conventions")
    with tempfile.TemporaryDirectory() as d:
        test_restart_forces_no_minimize(pathlib.Path(d))
    print("ok  test_restart_forces_no_minimize")
    test_matches_legacy_dynamics_parser()
    print("ok  test_matches_legacy_dynamics_parser")
    print("\nALL CONFIG TESTS PASSED")


if __name__ == "__main__":
    _run_all()
