"""Resume support for the cosmo continuous-synthesis drivers (``cosmo.csp``).

A production ``cosmo-csp`` / ``cosmo-cylinder`` run is hours to days of wall time and
today survives no interruption. This module makes a run resumable with three small,
human-readable on-disk artifacts and no heavyweight simulation checkpoint (mirrors the
sibling ``topo`` project's ``topo/csp/resume.py``):

1. **The schedule** (``dwell_times.dat``). The per-residue step counts are drawn from
   the seeded generator **once, before the main loop**, and persisted. For the explicit
   protocol the file header additionally carries the deterministic-but-expensive PTC
   restraint geometry (``a_target`` / ``p_target`` / tunnel-wall plane). On resume the
   driver *reads* this table instead of re-drawing the RNG or re-running the SLSQP PTC
   solve, so the kinetic schedule and restraint geometry are pinned identical across the
   interruption. This is the immutable **plan**.

2. **The progress log** (``progress.log``). An append-only ``L_XXX RUNNING`` /
   ``L_XXX DONE`` record of how far the run got. The ``DONE`` line is the commit point.

3. **The seed conformation.** The only per-residue state that crosses a residue boundary
   besides the schedule is the previous residue's final coordinates, already written to
   disk (:func:`load_final_pdb`).

The resume unit is the **residue**. There is no serialized RNG state anywhere -- the
schedule file *is* the materialized RNG output.
"""
from __future__ import annotations

import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import openmm as mm
from openmm import unit

PROGRESS_SCHEMA = 1
PROGRESS_FILENAME = "progress.log"


# --------------------------------------------------------------------------
# Output-layout knowledge (where a completed unit's final structure lands)
# --------------------------------------------------------------------------
def residue_final_path(out_root: Path, L: int) -> Path:
    """Path of the final (stage-3) structure for residue length ``L`` (explicit CSP).

    Consolidated layout: ``<out_root>/L_<L>/traj_final.pdb`` -- the single per-residue
    final written by stage 3; the seed for residue ``L+1`` and the resume-reload target.
    """
    return Path(out_root) / f"L_{L:03d}" / "traj_final.pdb"


def cylinder_final_path(out_root: Path, L: int) -> Path:
    """Path of residue ``L``'s final structure for the cylinder runner (flat layout)."""
    return Path(out_root) / f"L_{L:03d}" / "traj_final.pdb"


def phase_final_path(out_root: Path, name: str) -> Path:
    """Path of a post-synthesis phase's final (``ejection``)."""
    return Path(out_root) / name / "traj_final.pdb"


# --------------------------------------------------------------------------
# The 3-stage schedule table (dwell_times.dat) + PTC-geometry header
# --------------------------------------------------------------------------
@dataclass
class SchedRow:
    """One residue's persisted 3-stage schedule row."""
    L: int
    codon: str
    t_total: float
    times: Tuple[float, float, float]
    steps: Tuple[int, int, int]


def _fmt_float_exact(x: float) -> str:
    """``repr`` of a Python float -- round-trips to the exact same value."""
    return repr(float(x))


def write_schedule(path: Path, rows: List[SchedRow], params,
                   a_target: np.ndarray, p_target: np.ndarray,
                   wall_x: Optional[float]) -> None:
    """Write the immutable schedule table + PTC-geometry header to ``dwell_times.dat``.

    Called **once** at a fresh start (re-read, not rewritten, on resume). Data rows are
    byte-compatible with the legacy in-loop cosmo writer; the ``#PTC`` header block is
    new and carries the restraint geometry at full float precision.
    """
    timing = "uniform" if params.uniform_codon_time is not None else "per-codon"
    with open(path, "w") as fh:
        fh.write(
            "# O'Brien continuous-synthesis per-residue dwell times (cosmo.csp)\n"
            f"#   scale_factor={params.scale_factor:g}  dt={params.dt_ps} ps  "
            f"time_stage_1={params.time_stage_1:g} s  time_stage_2={params.time_stage_2:g} s\n"
            f"#   timing={timing}  "
            f"{'ribosome_traffic=on  ' if params.ribosome_traffic else ''}"
            f"random_seed={params.random_seed}\n"
            "#   t1/t2/t3 = sampled peptidyl-transfer / translocation / tRNA-binding "
            "dwell (s); steps = clamped integration steps actually run\n")
        fh.write(f"#PTC schema {PROGRESS_SCHEMA}\n")
        fh.write("#PTC a_target " + " ".join(_fmt_float_exact(v) for v in a_target) + "\n")
        fh.write("#PTC p_target " + " ".join(_fmt_float_exact(v) for v in p_target) + "\n")
        fh.write("#PTC wall_x " + (_fmt_float_exact(wall_x) if wall_x is not None else "none") + "\n")
        fh.write(
            "# L  codon  t_invivo_total_s  t1_s  t2_s  t3_s  "
            "ns1  ns2  ns3  steps1  steps2  steps3\n")
        for r in rows:
            t1, t2, t3 = r.times
            s1, s2, s3 = r.steps
            fh.write(
                f"{r.L:4d}  {r.codon:>5s}  {r.t_total:.6e}  "
                f"{t1:.6e}  {t2:.6e}  {t3:.6e}  "
                f"{t1 * 1e9 / params.scale_factor:.6e}  "
                f"{t2 * 1e9 / params.scale_factor:.6e}  "
                f"{t3 * 1e9 / params.scale_factor:.6e}  "
                f"{s1:8d}  {s2:8d}  {s3:8d}\n")


def read_schedule(path: Path) -> Tuple[List[SchedRow], np.ndarray, np.ndarray,
                                       Optional[float]]:
    """Read a ``dwell_times.dat`` table back into rows + PTC geometry.

    Returns ``(rows, a_target, p_target, wall_x)`` (rows ascending in ``L``).
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"schedule table not found: {path}")
    a_target = p_target = None
    wall_x: Optional[float] = None
    rows: List[SchedRow] = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#PTC"):
                tok = stripped.split()
                if len(tok) >= 5 and tok[1] == "a_target":
                    a_target = np.array([float(x) for x in tok[2:5]])
                elif len(tok) >= 5 and tok[1] == "p_target":
                    p_target = np.array([float(x) for x in tok[2:5]])
                elif len(tok) >= 3 and tok[1] == "wall_x":
                    wall_x = None if tok[2] == "none" else float(tok[2])
                continue
            if stripped.startswith("#"):
                continue
            tok = stripped.split()
            if len(tok) < 12:
                raise ValueError(f"{path}: malformed schedule row: {line!r}")
            rows.append(SchedRow(
                L=int(tok[0]), codon=tok[1], t_total=float(tok[2]),
                times=(float(tok[3]), float(tok[4]), float(tok[5])),
                steps=(int(tok[9]), int(tok[10]), int(tok[11]))))
    if a_target is None or p_target is None:
        raise ValueError(f"{path}: missing '#PTC a_target'/'#PTC p_target' header lines "
                         f"(not a resume-capable schedule file).")
    rows.sort(key=lambda r: r.L)
    return rows, a_target, p_target, wall_x


def schedule_covers(rows, L0: int, L_max: int) -> None:
    """Assert the persisted schedule covers exactly ``L0..L_max`` (contiguous)."""
    have = {r.L for r in rows}
    want = set(range(L0, L_max + 1))
    if have != want:
        missing = sorted(want - have)
        extra = sorted(have - want)
        raise SystemExit(
            f"[resume] persisted schedule covers L={sorted(have)[0]}..{sorted(have)[-1]} "
            f"but this run asks for L={L0}..{L_max}"
            + (f"; missing {missing}" if missing else "")
            + (f"; unexpected {extra}" if extra else "")
            + ". Extending a run is a fresh run (the schedule is fixed at first launch).")


# --------------------------------------------------------------------------
# The cylinder schedule table (single MD segment per residue; no PTC header)
# --------------------------------------------------------------------------
@dataclass
class CylSchedRow:
    """One residue's persisted cylinder schedule row (single MD segment)."""
    L: int
    codon: str
    dwell_s: float
    steps: int


def write_cylinder_schedule(path: Path, rows: List[CylSchedRow], params) -> None:
    """Write the immutable cylinder schedule to ``dwell_times.dat`` (written once)."""
    timing = "uniform" if params.uniform_codon_time is not None else "per-codon"
    with open(path, "w") as fh:
        fh.write(
            "# cylinder continuous-synthesis per-residue dwell times (cosmo.csp.cylinder)\n"
            f"#   scale_factor={params.scale_factor:g}  dt={params.dt_ps} ps  "
            f"timing={timing}  random_seed={params.random_seed}\n"
            "#   t_dwell = sampled codon dwell (s); ns = in-silico ns; steps = integration "
            "steps actually run (single MD segment)\n"
            "# L  codon  t_dwell_s  ns  steps\n")
        for r in rows:
            ns = r.dwell_s * 1e9 / params.scale_factor
            fh.write(f"{r.L:4d}  {r.codon:>5s}  {r.dwell_s:.6e}  {ns:.6e}  {r.steps:8d}\n")


def read_cylinder_schedule(path: Path) -> List[CylSchedRow]:
    """Read a cylinder ``dwell_times.dat`` table back into rows (ascending in ``L``)."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"schedule table not found: {path}")
    rows: List[CylSchedRow] = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            tok = stripped.split()
            if len(tok) < 5:
                raise ValueError(f"{path}: malformed cylinder schedule row: {line!r}")
            rows.append(CylSchedRow(L=int(tok[0]), codon=tok[1],
                                    dwell_s=float(tok[2]), steps=int(tok[4])))
    rows.sort(key=lambda r: r.L)
    return rows


# --------------------------------------------------------------------------
# The progress log (progress.log): mutable status
# --------------------------------------------------------------------------
def progress_path(out_root: Path) -> Path:
    """Path of the progress log under ``out_root``."""
    return Path(out_root) / PROGRESS_FILENAME


def progress_exists(out_root: Path) -> bool:
    """True iff a ``progress.log`` is present under ``out_root``."""
    return progress_path(out_root).is_file()


def write_progress_header(out_root: Path) -> None:
    """Create a fresh ``progress.log`` with its schema header line (truncates any prior)."""
    with open(progress_path(out_root), "w") as fh:
        fh.write(f"# cosmo csp progress log -- schema {PROGRESS_SCHEMA}\n")


def append_progress(out_root: Path, unit: str, status: str) -> None:
    """Append one ``<unit> <status>`` line and flush (the crash-safe commit point)."""
    with open(progress_path(out_root), "a") as fh:
        fh.write(f"{unit} {status}\n")
        fh.flush()


@dataclass
class Progress:
    """Parsed ``progress.log`` -- the last status seen per unit."""
    last_status: Dict[str, str]

    def is_done(self, unit: str) -> bool:
        """True iff ``unit``'s last recorded status is ``DONE``."""
        return self.last_status.get(unit) == "DONE"

    @property
    def last_done_residue(self) -> int:
        """Highest residue length ``L`` whose unit ``L_<L>`` is ``DONE`` (0 if none)."""
        best = 0
        for unit, st in self.last_status.items():
            if st == "DONE" and unit.startswith("L_"):
                try:
                    best = max(best, int(unit[2:]))
                except ValueError:
                    continue
        return best

    def running_units(self) -> List[str]:
        """Units whose last recorded status is ``RUNNING`` (in flight at the crash)."""
        return [u for u, s in self.last_status.items() if s == "RUNNING"]


def read_progress(out_root: Path) -> Progress:
    """Parse ``progress.log`` into a :class:`Progress` (last status wins per unit)."""
    last: Dict[str, str] = {}
    with open(progress_path(out_root)) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            tok = line.split()
            if len(tok) != 2:
                continue
            last[tok[0]] = tok[1]
    return Progress(last_status=last)


# --------------------------------------------------------------------------
# Resume actions: verify, drop, reload
# --------------------------------------------------------------------------
def verify_completed_units(out_root: Path, prog: Progress, L0: int,
                           final_path_fn=residue_final_path) -> None:
    """Assert every length ``L0..last_done`` has its final structure on disk.

    A ``DONE`` residue whose directory was deleted / lost to a scratch purge would
    otherwise leave a permanent hole (resume only reloads the last ``DONE`` residue). So
    verify presence of every prior length before continuing.
    """
    for L in range(L0, prog.last_done_residue + 1):
        fp = final_path_fn(out_root, L)
        if not fp.is_file():
            raise SystemExit(
                f"[resume] L_{L:03d} is marked DONE but its final structure {fp} is "
                f"missing -- the output tree is incomplete. Refusing to resume (re-run "
                f"fresh, or restore the missing length).")


def drop_running_units(out_root: Path, prog: Progress) -> List[str]:
    """Remove the on-disk directory of every ``RUNNING`` unit; return their names."""
    dropped: List[str] = []
    for unit in prog.running_units():
        d = Path(out_root) / unit
        if d.is_dir():
            shutil.rmtree(d)
            dropped.append(unit)
    return dropped


def load_final_pdb(path: Path) -> np.ndarray:
    """Reload a written nascent final structure as an ``(N, 3)`` nm coordinate array."""
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"seed structure to resume from not found: {path}")
    pos = mm.app.PDBFile(str(path)).getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    return np.asarray(pos)


def est_walltime(total_steps: int, params) -> str:
    """Nominal cost annotation for the up-front schedule report (exact simulated time)."""
    sim_ns = total_steps * params.dt_ps * 1e-3
    return (f" (~{sim_ns:,.1f} ns simulated at dt={params.dt_ps} ps; wall-time is "
            f"nominal -- dt-halving retries and a growing chain push it higher)")
