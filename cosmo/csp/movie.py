"""Stitch the CSP per-residue/-stage trajectories into one VMD-playable movie.

The continuous-synthesis runners write one standalone trajectory per MD segment, and
each nascent length ``L`` has a **different number of beads** (``L`` CA atoms). A
single VMD molecule needs a *constant* atom count across all frames, so the per-segment
DCDs cannot just be concatenated. This tool pads every frame up to the longest length
present and parks the not-yet-synthesized beads out of view.

Two output layouts are auto-detected:

* **CSP 3-stage** (:mod:`cosmo.csp.protocol`): ``<out>/L_<L>/traj_s<1,2,3>.dcd`` + a shared ``traj.psf`` --
  the three sub-stage segments of each residue are played in order.
* **cylinder / flat** (:mod:`cosmo.csp.cylinder`): ``<out>/L_<L>/traj.*`` -- one
  segment per residue.

Any post-synthesis phase (``ejection/``) is appended after the growth sequence.
A ready-to-run ``movie.tcl`` is also written.

Mirrors the sibling ``topo`` project's ``topo/csp/movie.py``.

Usage::

    cosmo-csp-movie -o synth_out
    python -m cosmo.csp.movie -o synth_out
    vmd -e synth_out/movie.tcl
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import sys
from typing import List, Optional, Tuple

import numpy as np

# Sentinel coordinate (angstrom) for parked (not-yet-synthesized) beads. Far enough
# that the VMD selection ``x > 9000`` cleanly isolates them.
SENTINEL_A = 99999.0

# Post-synthesis phase folders, in play order (harmlessly absent if not run).
POST_PHASES = ("ejection",)


def _pick_traj(phase_dir: str, outname: str = "traj") -> Optional[str]:
    """Pick the trajectory file to read for one stage/phase, or ``None`` if absent.

    Prefer ``<outname>.dcd`` when it actually holds frames, but a coarse ``nstout``
    (relative to the steps run per stage) can leave a stage's DCD **empty** (0 bytes)
    -- every frame interval landed past the end of the short run. In that case fall
    back to ``<outname>_final.pdb``, the last-frame snapshot the runner always writes,
    so the stage still contributes its (single) final conformation instead of being
    silently dropped. Without this, a coarse-output CSP run yields a movie that skips
    most lengths.

    Mirrors the sibling ``topo`` project's ``topo/csp/movie.py._pick_traj``.

    Parameters
    ----------
    phase_dir : str
        Directory of a single stage or post-synthesis phase (e.g.
        ``<out_root>/L_<L>/stage_<s>`` or ``<out_root>/ejection``).
    outname : str, optional
        Per-stage output basename used by the runner (default ``"traj"``); the
        candidate files are ``<outname>.dcd`` and ``<outname>_final.pdb``.

    Returns
    -------
    str or None
        Path to ``<outname>.dcd`` if it exists and is non-empty, else
        ``<outname>_final.pdb`` if present, else the empty ``<outname>.dcd`` if it
        exists, else ``None`` when no trajectory file is found.
    """
    dcd = os.path.join(phase_dir, f"{outname}.dcd")
    if os.path.isfile(dcd) and os.path.getsize(dcd) > 0:
        return dcd
    final_pdb = os.path.join(phase_dir, f"{outname}_final.pdb")
    if os.path.isfile(final_pdb):
        return final_pdb
    if os.path.isfile(dcd):  # empty DCD and no final-frame fallback -- read it anyway
        return dcd
    return None


def find_segments(out_root: str,
                  outname: str = "traj") -> List[Tuple[str, int, str, str]]:
    """Return the ordered growth segments ``[(label, n_atoms, psf, dcd), ...]``.

    Auto-detects the layout per length directory ``<out_root>/L_<L>/``:

    - if it contains per-stage ``traj_s<1,2,3>.dcd`` + a shared ``traj.psf`` (the
      consolidated CSP 3-stage layout), each stage is a segment in numeric stage order
      (``L=<L> s<s>``);
    - otherwise the flat ``L_<L>/traj.*`` is a single segment (``L=<L>``; the cylinder
      runner).

    ``n_atoms`` for a length-``L`` segment is ``L`` (nascent-only output). Only segments
    with a ``.psf`` and a readable trajectory are kept. A 0-byte (headerless) stage DCD is
    skipped quietly; for stage 3 the single ``traj_final.pdb`` still contributes its
    conformation. Lengths are sorted ascending.
    """
    lengths = []
    for d in glob.glob(os.path.join(out_root, "L_*")):
        m = re.search(r"L_(\d+)$", os.path.basename(d))
        if m:
            lengths.append((int(m.group(1)), d))
    lengths.sort(key=lambda t: t[0])

    segments: List[Tuple[str, int, str, str]] = []
    for L, d in lengths:
        per_stage = sorted(
            glob.glob(os.path.join(d, "traj_s*.dcd")),
            key=lambda p: int(re.search(r"traj_s(\d+)\.dcd$", p).group(1))
            if re.search(r"traj_s(\d+)\.dcd$", p) else 0)
        psf = os.path.join(d, "traj.psf")
        if per_stage and os.path.isfile(psf):
            # Consolidated CSP layout: shared traj.psf + per-stage traj_s{1,2,3}.dcd.
            for sd in per_stage:
                sm = re.search(r"traj_s(\d+)\.dcd$", os.path.basename(sd))
                s = sm.group(1) if sm else "?"
                traj = None
                if os.path.getsize(sd) > 0:
                    traj = sd
                elif s == "3" and os.path.isfile(os.path.join(d, "traj_final.pdb")):
                    # Coarse nstout can leave a stage DCD empty; stage 3's single final
                    # still contributes (stages 1/2 have no per-stage final now, so an
                    # empty s1/s2 DCD is simply skipped).
                    traj = os.path.join(d, "traj_final.pdb")
                if traj is not None:
                    segments.append((f"L={L} s{s}", L, psf, traj))
        else:
            # Flat per-length (cylinder): L_<L>/traj.psf + traj.dcd (+ traj_final.pdb).
            psf = os.path.join(d, f"{outname}.psf")
            traj = _pick_traj(d, outname)
            if os.path.isfile(psf) and traj is not None:
                segments.append((f"L={L}", L, psf, traj))
    return segments


def find_post(out_root: str, outname: str = "traj") -> List[Tuple[str, str, str]]:
    """Return ``[(name, psf, dcd), ...]`` for present post-synthesis phases."""
    found = []
    for name in POST_PHASES:
        d = os.path.join(out_root, name)
        psf = os.path.join(d, f"{outname}.psf")
        traj = _pick_traj(d, outname)
        if os.path.isfile(psf) and traj is not None:
            found.append((name, psf, traj))
    return found


def stitch_movie(out_root: str, out_prefix: str = "movie",
                 park: str = "sentinel", outname: str = "traj",
                 ribosome_pdb: Optional[str] = None,
                 verbose: bool = True) -> Tuple[str, str, str]:
    """Stitch the per-segment DCDs into ``<out_root>/<out_prefix>.{psf,dcd}`` (+ .tcl).

    Parameters
    ----------
    out_root : str
        The synthesis run's output root (contains the ``L_<L>/`` folders).
    out_prefix : str
        Basename for the movie files (default ``movie``).
    park : {'sentinel', 'cterm'}
        Where to put not-yet-synthesized beads in each frame -- far away and hidden
        by the VMD script (``sentinel``) or stacked on the C-terminus (``cterm``).
    outname : str
        Per-segment output basename used by the runner (default ``traj``).
    ribosome_pdb : str, optional
        A CG ribosome PDB to copy next to the movie and load as static scenery.

    Returns
    -------
    (psf_path, dcd_path, tcl_path)
    """
    import MDAnalysis as mda  # heavy; import only when actually stitching

    def log(msg: str) -> None:
        if verbose:
            print(msg)

    if park not in ("sentinel", "cterm"):
        raise ValueError(f"park must be 'sentinel' or 'cterm', got {park!r}.")

    segments = find_segments(out_root, outname=outname)
    if not segments:
        raise SystemExit(
            f"no per-length/-stage trajectories found under {out_root!r} "
            f"(expected {out_root}/L_<L>/traj_s<1,2,3>.dcd + traj.psf, or a flat "
            f"{out_root}/L_<L>/{outname}.dcd + .psf).")

    # Append the post-synthesis phases (at the final length) after growth.
    for name, psf, dcd in find_post(out_root, outname=outname):
        n = len(mda.Universe(psf).atoms)
        segments.append((name, n, psf, dcd))

    out_psf = os.path.join(out_root, f"{out_prefix}.psf")
    out_dcd = os.path.join(out_root, f"{out_prefix}.dcd")
    out_tcl = os.path.join(out_root, f"{out_prefix}.tcl")

    # Movie topology = a segment with the most beads (the final length); every frame
    # is padded up to N atoms.
    N = max(n for _, n, _, _ in segments)
    psf_full = next(psf for _, n, psf, _ in segments if n == N)
    shutil.copyfile(psf_full, out_psf)
    full = mda.Universe(psf_full)
    full.load_new(np.zeros((1, N, 3), dtype=np.float32))
    log(f"Movie topology: {out_psf}  ({N} beads = final length)")
    log(f"Parking not-yet-synthesized beads: {park}")

    total_frames = 0
    skipped = 0
    with mda.Writer(out_dcd, n_atoms=N) as writer:
        for label, n, psf, dcd in segments:
            # A segment can be unreadable if its DCD is still being written (an
            # in-progress run) or was truncated by a crash -- skip it (with a warning)
            # rather than aborting the whole movie. Open + read inside the guard so a
            # premature-EOF header error or a zero-frame DCD just drops that segment.
            try:
                u = mda.Universe(psf, dcd)
                nfr = 0
                for _ in u.trajectory:
                    coords = np.empty((N, 3), dtype=np.float32)
                    coords[:n] = u.atoms.positions
                    if n < N:
                        if park == "sentinel":
                            coords[n:] = SENTINEL_A
                        else:  # stack on the current C-terminus (last real bead)
                            coords[n:] = u.atoms.positions[n - 1]
                    full.atoms.positions = coords
                    writer.write(full.atoms)
                    nfr += 1
            except Exception as exc:
                skipped += 1
                log(f"  {label:>12}: SKIPPED (unreadable/empty DCD: "
                    f"{type(exc).__name__}: {exc})")
                continue
            total_frames += nfr
            log(f"  {label:>12}: {nfr} frames")
    if skipped:
        log(f"  ({skipped} segment(s) skipped -- truncated or still being written)")

    log(f"Movie trajectory: {out_dcd}  ({total_frames} frames total)")

    ribo_name = None
    if ribosome_pdb is not None:
        ribo_name = f"{out_prefix}_ribosome.pdb"
        shutil.copyfile(ribosome_pdb, os.path.join(out_root, ribo_name))
        log(f"Static ribosome reference: {os.path.join(out_root, ribo_name)}")

    _write_tcl(out_tcl, os.path.basename(out_psf), os.path.basename(out_dcd),
               park=park, ribosome_name=ribo_name)
    log(f"VMD script: {out_tcl}")
    log("")
    log(f"View it with:  vmd -e {out_tcl}")
    return out_psf, out_dcd, out_tcl


def _write_tcl(path: str, psf_name: str, dcd_name: str, park: str,
               ribosome_name: Optional[str] = None) -> None:
    """Write a VMD script that loads the movie and grows the chain N->C."""
    if park == "sentinel":
        sel = "not (x > 9000)"
        hide_note = ("# Not-yet-synthesized beads are parked far away and hidden "
                     "each frame\n# (selection re-evaluated per frame via selupdate).")
    else:
        sel = "all"
        hide_note = ("# 'cterm' parking: future beads are stacked on the "
                     "C-terminus (no hiding).")

    ribo_block = ""
    if ribosome_name is not None:
        ribo_block = f"""
# Static ribosome scenery: a separate molecule the chain grows inside.
mol new {ribosome_name} type pdb waitfor all
mol delrep 0 top
mol representation Points 1.0
mol color ColorID 6
mol selection {{all}}
mol material Transparent
mol addrep top
mol top 0
"""

    tcl = f"""# VMD visualization of the protein synthesis movie.
# Generated by cosmo.csp.movie.
#   vmd -e {os.path.basename(path)}      (run from this folder)

mol new {psf_name} type psf waitfor all
mol addfile {dcd_name} type dcd waitfor all

mol delrep 0 top

{hide_note}
# Beads (van der Waals); colored by residue id so the growing chain is a gradient.
mol representation VDW 1.5 16.0
mol color ResID
mol selection {{{sel}}}
mol material Opaque
mol addrep top
mol selupdate 0 top on

# Backbone trace (bonds between consecutive synthesized beads).
mol representation Licorice 0.4 16.0
mol color ResID
mol selection {{{sel}}}
mol addrep top
mol selupdate 1 top on
{ribo_block}
# Fit the camera to the LAST frame (full length -> nothing parked), then rewind.
set nf [molinfo top get numframes]
if {{$nf > 0}} {{
    animate goto [expr {{$nf - 1}}]
    display resetview
    animate goto 0
}}

axes location off
display projection Orthographic
color Display Background white
animate speed 0.85

puts "Loaded $nf frames. Press Play (or run: animate forward) to watch the chain grow."
"""
    with open(path, "w") as fh:
        fh.write(tcl)


def main(argv: Optional[List[str]] = None) -> None:
    """CLI: ``cosmo-csp-movie -o <out_root>``."""
    p = argparse.ArgumentParser(
        prog="cosmo-csp-movie",
        description="Stitch the CSP per-residue/-stage trajectories "
                    "(<out_root>/L_<L>/traj_s<1,2,3>.dcd or a flat traj.dcd) -- plus any post-synthesis "
                    "phase (ejection/) -- into one "
                    "VMD-playable movie that grows the nascent chain N->C.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-o", "--out-root", required=True,
                   help="synthesis run output root (contains the L_<L>/ folders).")
    p.add_argument("--prefix", default="movie",
                   help="basename for the stitched movie files.")
    p.add_argument("--park", default="sentinel", choices=["sentinel", "cterm"],
                   help="where to put not-yet-synthesized beads each frame.")
    p.add_argument("--outname", default="traj",
                   help="per-segment output basename used by the runner.")
    p.add_argument("--ribosome", default=None,
                   help="optional CG ribosome PDB to load as static scenery.")
    p.add_argument("--tunnel", default=None, metavar="CYLINDER_INI",
                   help="cylinder.ini to draw the analytic exit tunnel (bore + PTC "
                        "cap + exit-face wall) into movie.tcl, for cylinder-synthesis "
                        "runs (cosmo-cylinder). Omit for a plain movie.")
    p.add_argument("--wall-outer", type=float, default=None,
                   help="outer radius (nm) of the drawn exit-face wall when --tunnel "
                        "is given (default: bore radius + 3 nm).")
    if argv is None and len(sys.argv) == 1:
        p.print_help()
        sys.exit(0)
    args = p.parse_args(argv)
    _psf, _dcd, tcl = stitch_movie(args.out_root, out_prefix=args.prefix, park=args.park,
                                   outname=args.outname, ribosome_pdb=args.ribosome)

    # Cylinder synthesis has no ribosome scenery; on request, append the analytic
    # tunnel (read from the run's cylinder.ini) to the generated movie.tcl.
    if args.tunnel is not None:
        from cosmo.csp.cylinder import tunnel_tcl
        with open(tcl, "a") as fh:
            fh.write("\n" + tunnel_tcl(args.tunnel, wall_outer_nm=args.wall_outer))
        print(f"Appended analytic tunnel (from {args.tunnel}) to {tcl}")


if __name__ == "__main__":
    main()
