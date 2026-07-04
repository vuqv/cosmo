"""Stitch the per-length elongation trajectories into one VMD-playable movie.

The elongation runner (:mod:`cosmo.translation.elongate`) writes one standalone
trajectory per nascent-chain length (``<out>/L_<L>/traj.dcd``), and each length
has a **different number of beads** (length ``L`` has ``L`` CA atoms). A single
VMD molecule, however, needs a *constant* atom count across all frames, so the
per-length DCDs cannot just be concatenated.

This tool builds a single fixed-width trajectory:

* atom count = the **longest** length present (``N`` = max ``L``);
* each length-``L`` frame keeps residues ``1..L`` at their simulated coordinates
  and **parks** the not-yet-synthesized residues ``L+1..N`` at a far sentinel
  coordinate (default) so they can be hidden in VMD, or stacks them on the
  C-terminus (``--park cterm``);
* frames are written in length order, so the movie is ``F`` frames of ``L=L0``,
  then ``F`` frames of ``L=L0+1``, ... (exactly the "20 frames of L=5, then 20 of
  L=6, ..." layout).

If a post-elongation phase was run (``ejection/`` or ``stallation/``; build step
v2), its frames (at the final length) are appended **after** the growth sequence,
so the movie continues into the protein's release/ejection (or stalled wiggling).

It also writes a ready-to-run VMD script (``movie.tcl``) that loads the movie and
hides the parked beads each frame, so the chain appears to grow N->C.

Mirrors the sibling ``topo`` project's ``topo/translation/make_movie.py``.

Usage::

    cosmo-elongate-movie -o synth_out
    python -m cosmo.translation.make_movie -o synth_out
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

# Sentinel coordinate (angstrom) for parked (not-yet-synthesized) beads. Far
# enough that the VMD selection ``x > 9000`` cleanly isolates them; the view is
# fit to the full-length frame (where nothing is parked), so this does not zoom
# the camera out (see the generated movie.tcl).
SENTINEL_A = 99999.0


def find_lengths(out_root: str,
                 outname: str = "traj") -> List[Tuple[int, str, str]]:
    """Return ``[(L, psf, dcd), ...]`` for each ``<out_root>/L_<L>/`` length.

    Sorted by ``L``; only lengths with both a ``.psf`` and a ``.dcd`` are kept.
    """
    items = []
    for d in glob.glob(os.path.join(out_root, "L_*")):
        m = re.search(r"L_(\d+)$", os.path.basename(d))
        if not m:
            continue
        L = int(m.group(1))
        psf = os.path.join(d, f"{outname}.psf")
        dcd = os.path.join(d, f"{outname}.dcd")
        if os.path.isfile(psf) and os.path.isfile(dcd):
            items.append((L, psf, dcd))
    items.sort(key=lambda t: t[0])
    return items


# Post-elongation phase folders, in the order they should appear after growth
# (build step v2). Harmlessly absent in v1.
POST_PHASES = ("ejection", "stallation")


def find_post(out_root: str, outname: str = "traj") -> List[Tuple[str, str, str]]:
    """Return ``[(name, psf, dcd), ...]`` for present post-elongation phases.

    Looks for ``<out_root>/ejection/`` and ``<out_root>/stallation/`` (in that
    order) -- the optional post-elongation runs written after the chain reaches its
    final length. Only phases with both a ``.psf`` and a ``.dcd`` are returned.
    """
    found = []
    for name in POST_PHASES:
        d = os.path.join(out_root, name)
        psf = os.path.join(d, f"{outname}.psf")
        dcd = os.path.join(d, f"{outname}.dcd")
        if os.path.isfile(psf) and os.path.isfile(dcd):
            found.append((name, psf, dcd))
    return found


def stitch_movie(out_root: str, out_prefix: str = "movie",
                 park: str = "sentinel", outname: str = "traj",
                 ribosome_pdb: Optional[str] = None,
                 verbose: bool = True) -> Tuple[str, str, str]:
    """Stitch per-length DCDs into ``<out_root>/<out_prefix>.{psf,dcd}`` (+ .tcl).

    Parameters
    ----------
    out_root : str
        The elongation run's output root (contains the ``L_<L>/`` folders).
    out_prefix : str
        Basename for the movie files (default ``movie``).
    park : {'sentinel', 'cterm'}
        Where to put not-yet-synthesized beads in each frame. ``sentinel`` (far
        away, hidden by the VMD script -- cleanest) or ``cterm`` (stacked on the
        C-terminus -- no VMD selection needed, but leaves a small bead cluster at
        the growing tip).
    outname : str
        Per-length output basename used by the runner (default ``traj``).
    ribosome_pdb : str, optional
        Build step v2: a CG ribosome PDB to copy next to the movie and load as
        static scenery in the generated ``movie.tcl``.

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

    items = find_lengths(out_root, outname=outname)
    if not items:
        raise SystemExit(
            f"no per-length trajectories found under {out_root!r} "
            f"(expected {out_root}/L_<L>/{outname}.dcd + .psf).")

    # Ordered playback segments: each growth length, then the post-elongation phase
    # (ejection / stallation) at full length. Each segment is (label, n_atoms,
    # psf, dcd); the post phases run at the final length, so n_atoms = max L.
    segments = [(f"L={L}", L, psf, dcd) for L, psf, dcd in items]
    for name, psf, dcd in find_post(out_root, outname=outname):
        n = len(mda.Universe(psf).atoms)
        segments.append((name, n, psf, dcd))

    out_psf = os.path.join(out_root, f"{out_prefix}.psf")
    out_dcd = os.path.join(out_root, f"{out_prefix}.dcd")
    out_tcl = os.path.join(out_root, f"{out_prefix}.tcl")

    # Movie topology = a segment with the most beads (the final length); every
    # frame is padded up to N atoms.
    N = max(n for _, n, _, _ in segments)
    psf_full = next(psf for _, n, psf, _ in segments if n == N)
    shutil.copyfile(psf_full, out_psf)
    full = mda.Universe(psf_full)
    # Give the topology-only universe a single in-memory coordinate frame so we
    # can assign per-frame positions into it before writing.
    full.load_new(np.zeros((1, N, 3), dtype=np.float32))
    log(f"Movie topology: {out_psf}  ({N} beads = final length)")
    log(f"Parking not-yet-synthesized beads: {park}")

    total_frames = 0
    with mda.Writer(out_dcd, n_atoms=N) as writer:
        for label, n, psf, dcd in segments:
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
            total_frames += nfr
            log(f"  {label:>12}: {nfr} frames")

    log(f"Movie trajectory: {out_dcd}  ({total_frames} frames total)")

    # Optional static ribosome reference (v2): copy it next to the movie so the
    # generated tcl can load it as fixed scenery the chain grows inside.
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
    # The hiding selection is only needed for the 'sentinel' parking scheme.
    if park == "sentinel":
        sel = "not (x > 9000)"
        hide_note = ("# Not-yet-synthesized beads are parked far away and hidden "
                     "each frame\n# (selection re-evaluated per frame via selupdate).")
    else:
        sel = "all"
        hide_note = ("# 'cterm' parking: future beads are stacked on the "
                     "C-terminus (no hiding).")

    # Optional: load the static ribosome (v2) as a separate molecule for context.
    ribo_block = ""
    if ribosome_name is not None:
        ribo_block = f"""
# Static ribosome scenery (v2): a separate molecule the chain grows inside.
mol new {ribosome_name} type pdb waitfor all
mol delrep 0 top
mol representation Points 1.0
mol color ColorID 6
mol selection {{all}}
mol material Transparent
mol addrep top
mol top 0
"""

    tcl = f"""# VMD visualization of the co-translational elongation movie.
# Generated by cosmo.translation.make_movie.
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
    """CLI: ``cosmo-elongate-movie -o <out_root>``."""
    p = argparse.ArgumentParser(
        prog="cosmo-elongate-movie",
        description="Stitch the per-length elongation trajectories "
                    "(<out_root>/L_<L>/traj.dcd) -- plus any post-elongation phase "
                    "(ejection/ or stallation/) -- into one VMD-playable movie that "
                    "grows the nascent chain N->C, and write a movie.tcl to view it.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("-o", "--out-root", required=True,
                   help="elongation run output root (contains the L_<L>/ folders).")
    p.add_argument("--prefix", default="movie",
                   help="basename for the stitched movie files.")
    p.add_argument("--park", default="sentinel", choices=["sentinel", "cterm"],
                   help="where to put not-yet-synthesized beads each frame.")
    p.add_argument("--outname", default="traj",
                   help="per-length output basename used by the runner.")
    p.add_argument("--ribosome", default=None,
                   help="optional CG ribosome PDB (v2) to load as static scenery "
                        "in the generated movie.tcl.")
    if argv is None and len(sys.argv) == 1:
        p.print_help()
        sys.exit(0)
    args = p.parse_args(argv)
    stitch_movie(args.out_root, out_prefix=args.prefix, park=args.park,
                 outname=args.outname, ribosome_pdb=args.ribosome)


if __name__ == "__main__":
    main()
