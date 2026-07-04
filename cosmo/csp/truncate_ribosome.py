#!/usr/bin/env python3
"""
Truncate a coarse-grained ribosome around the exit tunnel.

The ribosome is **rigid** and the non-bonded interactions are cut off at 2 nm, so
ribosome beads far from any position the nascent chain can reach exert no force on
it. Removing them is therefore exact for the nascent-chain dynamics (and removes
the bulk of the ribosome). This mirrors the sibling ``topo`` project's
``topo/translation/truncate_ribosome.py`` and the O'Brien-lab procedure (we
implement it ourselves).

It supersedes the older per-atom :func:`cosmo.utils.crop_ribosome.filter_pdb`
(``x<0`` drop, single cylinder, ``x>=x_threshold`` keep): the keep decision here
is made **per residue** (keep the whole P/R/BR or CA unit if *any* bead qualifies),
which keeps the coarse-grained beads intact and matches ``topo``.

Assumptions
-----------
* The structure is **aligned with the exit-tunnel central line on the X-axis**
  (so the tunnel line is y = z = 0; radial distance d = sqrt(y^2 + z^2)).
* The **PTC is near x = 0** and the **exit is toward +x** (true for
  ``4v9d_50S_PtR_5jte_AtR_model_cg.pdb``).

Keep rule (decided **per residue**: keep the *whole* residue/nucleotide if **any**
of its beads qualifies, so P/R/BR units stay intact):

  1. **Tunnel cylinder:** d <= ``r_cyl`` AND ``x_lo`` <= x <= ``x_exit``.
  2. **Exit half-space:** x >= ``x_exit`` (kept regardless of d).
  3. Optionally, any residue whose segID is in ``keep_segids`` (e.g. the tRNAs).

Usage
-----
    python -m cosmo.csp.truncate_ribosome \
        -i structures/4v9d_50S_PtR_5jte_AtR_model_cg.pdb \
        -o structures/4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb \
        [--r-cyl 30] [--x-lo -8] [--x-exit 58] [--keep-segids PtR,AtR]
"""
from __future__ import annotations

import argparse
from collections import OrderedDict


def _read_residues(path):
    """Group ATOM lines into residues, preserving line text and order.

    Returns list of dicts: {"key", "segid", "lines":[raw], "coords":[(x,y,z)]}.
    """
    residues = OrderedDict()
    order = []
    for line in open(path):
        if line[:6].strip() not in ("ATOM", "HETATM"):
            continue
        chain = line[21]; resseq = line[22:26]; icode = line[26]
        resname = line[17:20]; segid = line[72:76].strip()
        key = (chain, resseq, icode, resname, segid)
        if key not in residues:
            residues[key] = {"key": key, "segid": segid, "lines": [], "coords": []}
            order.append(key)
        residues[key]["lines"].append(line.rstrip("\n"))
        residues[key]["coords"].append(
            (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return [residues[k] for k in order]


def _keep_residue(res, r_cyl, x_lo, x_exit, keep_segids):
    """True if any bead satisfies the cylinder/exit rule (or segID is force-kept)."""
    if res["segid"] in keep_segids:
        return True
    for x, y, z in res["coords"]:
        d = (y * y + z * z) ** 0.5
        if (d <= r_cyl and x_lo <= x <= x_exit) or (x >= x_exit):
            return True
    return False


def truncate(input_pdb, output_pdb, r_cyl=30.0, x_lo=-8.0, x_exit=58.0,
             keep_segids=(), verbose=True):
    """Truncate ``input_pdb`` and write ``output_pdb``. Returns a stats dict."""
    keep_segids = set(keep_segids)
    residues = _read_residues(input_pdb)

    kept_lines = []
    serial = 0
    prev_chain = None
    seg_total = OrderedDict()
    seg_kept = OrderedDict()
    n_res_in = n_res_kept = n_beads_in = 0

    for res in residues:
        seg = res["segid"]
        nb = len(res["lines"])
        n_res_in += 1
        n_beads_in += nb
        seg_total[seg] = seg_total.get(seg, 0) + nb
        if not _keep_residue(res, r_cyl, x_lo, x_exit, keep_segids):
            continue
        n_res_kept += 1
        seg_kept[seg] = seg_kept.get(seg, 0) + nb
        chain = res["key"][0]
        if prev_chain is not None and chain != prev_chain:
            kept_lines.append("TER")
        prev_chain = chain
        for line in res["lines"]:
            serial += 1
            kept_lines.append("%s%5d%s" % (line[:6], serial % 100000, line[11:]))
    kept_lines.append("END")

    with open(output_pdb, "w") as fh:
        fh.write("REMARK  Truncated CG ribosome from %s\n" % input_pdb)
        fh.write("REMARK  rule: keep residue if any bead has "
                 "(d<=%.1f and %.1f<=x<=%.1f) or x>=%.1f; d=sqrt(y^2+z^2) "
                 "(tunnel = X-axis)\n" % (r_cyl, x_lo, x_exit, x_exit))
        if keep_segids:
            fh.write("REMARK  force-kept segIDs: %s\n" % ", ".join(sorted(keep_segids)))
        fh.write("\n".join(kept_lines) + "\n")

    stats = dict(beads_in=n_beads_in, beads_kept=serial,
                 residues_in=n_res_in, residues_kept=n_res_kept,
                 seg_total=seg_total, seg_kept=seg_kept)
    if verbose:
        red = 100.0 * (1 - serial / n_beads_in) if n_beads_in else 0.0
        print("Truncated %s -> %s" % (input_pdb, output_pdb))
        print("  rule: (d<=%.0f and %.0f<=x<=%.0f) or x>=%.0f   keep_segids=%s"
              % (r_cyl, x_lo, x_exit, x_exit, sorted(keep_segids) or "none"))
        print("  beads   : %d -> %d  (%.1f%% removed)" % (n_beads_in, serial, red))
        print("  residues: %d -> %d" % (n_res_in, n_res_kept))
        print("  per-segID kept/total:")
        for s in seg_total:
            print("     %-4s %6d / %-6d" % (s, seg_kept.get(s, 0), seg_total[s]))
    return stats


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="python -m cosmo.csp.truncate_ribosome",
        description="Truncate a CG ribosome around the exit tunnel "
                    "(tunnel = X-axis; PTC near x=0; exit toward +x).")
    p.add_argument("-i", "--input", required=True, help="CG ribosome PDB")
    p.add_argument("-o", "--output", required=True, help="truncated CG PDB")
    p.add_argument("--r-cyl", type=float, default=30.0,
                   help="radial cutoff from the tunnel axis (Å)")
    p.add_argument("--x-lo", type=float, default=-8.0,
                   help="low-x bound of the tunnel cylinder (PTC side, Å)")
    p.add_argument("--x-exit", type=float, default=58.0,
                   help="exit plane: beyond this x, keep all beads (Å); "
                        "58 follows O'Brien et al.")
    p.add_argument("--keep-segids", default="",
                   help="comma-separated segIDs to always keep (e.g. PtR,AtR)")
    args = p.parse_args(argv)
    keep = [s.strip() for s in args.keep_segids.split(",") if s.strip()]
    truncate(args.input, args.output, r_cyl=args.r_cyl, x_lo=args.x_lo,
             x_exit=args.x_exit, keep_segids=keep)


if __name__ == "__main__":
    main()
