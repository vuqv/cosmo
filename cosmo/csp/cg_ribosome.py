#!/usr/bin/env python3
"""
Coarse-grain a ribosome (or any protein + RNA structure) to a CG convention.

Protein is always one Cα bead per residue (residue name kept). For RNA there are
**two representations** (PLAN.md §5), selected by ``--rna-model``:

* ``topo`` (default) -- the O'Brien 3/4-bead representation:
    - ``P``  : phosphate (placed at the ``P`` atom; q = -1e),
    - ``R``  : ribose-ring centroid (C1',C2',C3',C4',O4'),
    - ``BR`` : centroid of each conjugated base ring (pyrimidine -> ``BR1``;
               purine -> ``BR1`` + ``BR2``).
* ``cosmo`` -- cosmo's own **1 bead per nucleotide** (the CA/P model): a single
    ``P`` bead at the phosphate (fallback: ribose centroid for a 5'-terminal
    nucleotide), residue name kept. Its charge/radius come from cosmo's
    ``model_parameters`` nucleic entries at load time.

Both representations interact with the nascent chain only through excluded volume +
electrostatics (PLAN.md §5), so the choice changes only the scenery's steric/charge
granularity. The PDB **segID** (cols 73-76) is preserved on every bead (``AtR`` /
``PtR`` / ``23S`` / ``5S`` / ``L*``), so molecules can be selected downstream.

Coarse-grain **before** truncation (every centroid is computed from the complete
all-atom residue). Mirrors topo's ``cg_ribosome.py`` (the ``cosmo`` rep is the
addition).

Usage
-----
    python -m cosmo.csp.cg_ribosome -i ribosome.pdb -o ribosome_cg.pdb \
        [--rna-model topo|cosmo]
"""
from __future__ import annotations

import argparse
from collections import OrderedDict

AMINO_ACIDS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}
PURINES = {"A", "G", "RA", "RG", "ADE", "GUA"}
PYRIMIDINES = {"C", "U", "RC", "RU", "CYT", "URA"}

RIBOSE_RING = ["C1'", "C2'", "C3'", "C4'", "O4'"]
BASE_RING_6 = ["N1", "C2", "N3", "C4", "C5", "C6"]   # pyrimidine ring; purine 6-ring
PURINE_RING_5 = ["C4", "C5", "N7", "C8", "N9"]       # purine 5-ring (fused)


def _centroid(coords):
    """Mean of a list of (x, y, z); None if the list is empty."""
    if not coords:
        return None
    n = len(coords)
    return (sum(c[0] for c in coords) / n,
            sum(c[1] for c in coords) / n,
            sum(c[2] for c in coords) / n)


def _parse_pdb(path):
    """Parse ATOM/HETATM records into ordered residues (first altloc kept)."""
    residues = OrderedDict()
    order = []
    with open(path) as fh:
        for line in fh:
            if line[:6].strip() not in ("ATOM", "HETATM"):
                continue
            name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21]
            resseq = line[22:26].strip()
            icode = line[26]
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except ValueError:
                continue
            segid = line[72:76].strip()
            key = (chain, resseq, icode, resname, segid)
            if key not in residues:
                residues[key] = {"chain": chain, "resseq": resseq, "icode": icode,
                                 "resname": resname, "segid": segid,
                                 "atoms": OrderedDict()}
                order.append(key)
            atoms = residues[key]["atoms"]
            if name not in atoms:          # keep first altloc only
                atoms[name] = (x, y, z)
    return [residues[k] for k in order]


def _beads_for_residue(res, rna_model, warn):
    """Return a list of (atom_name, element, (x,y,z)) CG beads for one residue."""
    rn = res["resname"]
    atoms = res["atoms"]

    # ---- protein: one Cα bead, residue name unchanged ----
    if rn in AMINO_ACIDS:
        if "CA" in atoms:
            return [("CA", "C", atoms["CA"])]
        warn(f"protein residue {rn} {res['chain']}{res['resseq']} has no CA atom; skipped")
        return []

    is_purine = rn in PURINES or ("N9" in atoms)
    is_pyrimidine = rn in PYRIMIDINES
    if not (is_purine or is_pyrimidine):
        warn(f"unrecognized residue {rn} {res['chain']}{res['resseq']} (not protein/RNA); skipped")
        return []

    # ---- RNA: cosmo 1-bead representation ----
    if rna_model == "cosmo":
        if "P" in atoms:
            return [("P", "P", atoms["P"])]
        ribose = _centroid([atoms[a] for a in RIBOSE_RING if a in atoms])
        if ribose is not None:
            return [("P", "P", ribose)]   # 5'-terminal: no phosphate -> ribose centroid
        warn(f"nucleotide {rn} {res['chain']}{res['resseq']} has no P/ribose atoms; skipped")
        return []

    # ---- RNA: topo 3/4-bead representation ----
    beads = []
    if "P" in atoms:                       # absent on a 5'-terminal nucleotide
        beads.append(("P", "P", atoms["P"]))
    ribose = [atoms[a] for a in RIBOSE_RING if a in atoms]
    c = _centroid(ribose)
    if c is None or len(ribose) < 3:
        warn(f"nucleotide {rn} {res['chain']}{res['resseq']} has too few ribose "
             f"ring atoms ({len(ribose)}/5); R bead skipped")
    else:
        beads.append(("R", "C", c))
    ring6 = _centroid([atoms[a] for a in BASE_RING_6 if a in atoms])
    if is_purine:
        ring5 = _centroid([atoms[a] for a in PURINE_RING_5 if a in atoms])
        if ring6 is not None:
            beads.append(("BR1", "C", ring6))
        if ring5 is not None:
            beads.append(("BR2", "C", ring5))
        if ring6 is None or ring5 is None:
            warn(f"purine {rn} {res['chain']}{res['resseq']} missing base-ring atoms")
    else:  # pyrimidine: single ring -> BR1
        if ring6 is not None:
            beads.append(("BR1", "C", ring6))
        else:
            warn(f"pyrimidine {rn} {res['chain']}{res['resseq']} missing base-ring atoms")
    return beads


def _pdb_atom_line(serial, name, resname, chain, resseq, icode, xyz, element, segid=""):
    """Format one standard PDB ATOM record (fixed columns; segID in 73-76)."""
    name4 = name[:4] if len(name) >= 4 else (" " + name).ljust(4)
    return ("ATOM  %5d %-4s%1s%3s %1s%4s%1s   %8.3f%8.3f%8.3f%6.2f%6.2f      %-4s%2s\n" % (
        serial % 100000, name4, " ", resname[:3], (chain or " ")[:1], resseq[:4],
        icode if icode.strip() else " ", xyz[0], xyz[1], xyz[2], 1.0, 0.0,
        segid[:4], element[:2]))


def coarse_grain(input_pdb: str, output_pdb: str, rna_model: str = "topo",
                 verbose: bool = True) -> dict:
    """Coarse-grain ``input_pdb`` and write ``output_pdb``. Returns a stats dict.

    ``rna_model`` selects the RNA representation: ``topo`` (3/4-bead P/R/BR,
    default) or ``cosmo`` (1 bead/nucleotide).
    """
    if rna_model not in ("topo", "cosmo"):
        raise ValueError(f"rna_model must be 'topo' or 'cosmo', got {rna_model!r}.")
    warnings_list = []

    def warn(msg):
        warnings_list.append(msg)

    residues = _parse_pdb(input_pdb)
    stats = dict(protein_beads=0, rna_nucleotides=0, P=0, R=0, BR=0,
                 purines=0, pyrimidines=0, residues_skipped=0)
    segids = OrderedDict()

    lines = []
    serial = 0
    prev_chain = None
    for res in residues:
        beads = _beads_for_residue(res, rna_model, warn)
        if not beads:
            stats["residues_skipped"] += 1
            continue
        if prev_chain is not None and res["chain"] != prev_chain:
            lines.append("TER\n")
        prev_chain = res["chain"]

        if res["resname"] in AMINO_ACIDS:
            stats["protein_beads"] += 1
        else:
            stats["rna_nucleotides"] += 1
            if any(b[0] == "BR2" for b in beads):
                stats["purines"] += 1
            elif rna_model == "topo":
                stats["pyrimidines"] += 1

        sid = res.get("segid", "")
        for name, element, xyz in beads:
            serial += 1
            segids[sid] = segids.get(sid, 0) + 1
            btype = name.rstrip("0123456789")   # CA/P/R/BR
            if btype in ("P", "R", "BR"):
                stats[btype] += 1
            lines.append(_pdb_atom_line(serial, name, res["resname"], res["chain"],
                                        res["resseq"], res["icode"], xyz, element,
                                        res.get("segid", "")))
    lines.append("END\n")

    with open(output_pdb, "w") as fh:
        fh.write(f"REMARK  cosmo coarse-grained model of {input_pdb}\n")
        if rna_model == "topo":
            fh.write("REMARK  protein: 1 bead/residue at CA. RNA (topo): P, R, "
                     "BR1/BR2 (3/4 beads per nucleotide).\n")
        else:
            fh.write("REMARK  protein: 1 bead/residue at CA. RNA (cosmo): 1 bead "
                     "per nucleotide at the phosphate.\n")
        fh.writelines(lines)

    if verbose:
        print(f"Coarse-grained {input_pdb} -> {output_pdb}  (rna_model={rna_model})")
        print(f"  protein Cα beads : {stats['protein_beads']}")
        print(f"  RNA nucleotides  : {stats['rna_nucleotides']} "
              f"(purines {stats['purines']}, pyrimidines {stats['pyrimidines']})")
        if rna_model == "topo":
            print(f"  RNA beads        : P={stats['P']}  R={stats['R']}  BR={stats['BR']}")
        else:
            print(f"  RNA beads        : P={stats['P']} (1 per nucleotide)")
        print(f"  total beads      : {serial}")
        print(f"  segIDs preserved : {len(segids)}")
        if stats["residues_skipped"]:
            print(f"  residues skipped : {stats['residues_skipped']}")
        if warnings_list:
            print(f"  warnings         : {len(warnings_list)} (e.g. {warnings_list[0]})")
    stats["total_beads"] = serial
    stats["segids"] = segids
    stats["warnings"] = warnings_list
    return stats


def main(argv=None):
    p = argparse.ArgumentParser(
        prog="python -m cosmo.csp.cg_ribosome",
        description="Coarse-grain a protein+RNA structure (protein: CA; RNA: "
                    "P/R/BR for rna_model=topo, or 1 bead/nucleotide for cosmo).")
    p.add_argument("-i", "--input", required=True, help="input all-atom PDB")
    p.add_argument("-o", "--output", required=True, help="output coarse-grained PDB")
    p.add_argument("--rna-model", default="topo", choices=["topo", "cosmo"],
                   help="RNA representation (default: topo = 3/4-bead P/R/BR).")
    args = p.parse_args(argv)
    coarse_grain(args.input, args.output, rna_model=args.rna_model)


if __name__ == "__main__":
    main()
