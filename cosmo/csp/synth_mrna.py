"""Synthetic mRNA preparation: build a *fastest*-, *slowest*- or *median*-codon mRNA.

A ``fastest`` / ``slowest`` / ``median`` mRNA encodes the target protein with, for every
residue, the fastest / slowest / median-dwell-time synonymous codon for that amino acid
-- a point on the "synonymous mutation" (codon-optimization) axis while the protein
sequence is held fixed (``median`` is a neutral, non-extreme reference between the two
extremes). Because the mRNA is only a *timing* input to the synthesis workflow (the
structure comes from the PDB), this is a one-shot **preparation step**: it reads the
protein's amino-acid sequence from the PDB, picks the chosen codon per residue from a
codon dwell-time table, writes ``mrna_fastest.txt`` / ``mrna_slowest.txt`` /
``mrna_median.txt`` next to the PDB, and hands that file to the normal per-codon path.

The codon -> amino-acid mapping and the per-codon dwell time both come from the
3-column dwell-time tables under ``assets/csp/codon_dwell_times/`` (``codon  time  aa``);
if a table has no amino-acid column the standard genetic code is used as a fallback.
"""
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

# The reserved mRNA keywords. A real mRNA *filename* must not be one of these.
SYNTHETIC_MRNA_MODES = ("fastest", "slowest", "median")

# Standard genetic code (RNA codon -> amino-acid 3-letter name; STOP for the three stop
# codons). Only a *fallback* -- the dwell-time tables carry the amino acid in column 3.
GENETIC_CODE: Dict[str, str] = {
    "UUU": "PHE", "UUC": "PHE", "UUA": "LEU", "UUG": "LEU",
    "CUU": "LEU", "CUC": "LEU", "CUA": "LEU", "CUG": "LEU",
    "AUU": "ILE", "AUC": "ILE", "AUA": "ILE", "AUG": "MET",
    "GUU": "VAL", "GUC": "VAL", "GUA": "VAL", "GUG": "VAL",
    "UCU": "SER", "UCC": "SER", "UCA": "SER", "UCG": "SER",
    "CCU": "PRO", "CCC": "PRO", "CCA": "PRO", "CCG": "PRO",
    "ACU": "THR", "ACC": "THR", "ACA": "THR", "ACG": "THR",
    "GCU": "ALA", "GCC": "ALA", "GCA": "ALA", "GCG": "ALA",
    "UAU": "TYR", "UAC": "TYR", "UAA": "STOP", "UAG": "STOP",
    "CAU": "HIS", "CAC": "HIS", "CAA": "GLN", "CAG": "GLN",
    "AAU": "ASN", "AAC": "ASN", "AAA": "LYS", "AAG": "LYS",
    "GAU": "ASP", "GAC": "ASP", "GAA": "GLU", "GAG": "GLU",
    "UGU": "CYS", "UGC": "CYS", "UGA": "STOP", "UGG": "TRP",
    "CGU": "ARG", "CGC": "ARG", "CGA": "ARG", "CGG": "ARG",
    "AGU": "SER", "AGC": "SER", "AGA": "ARG", "AGG": "ARG",
    "GGU": "GLY", "GGC": "GLY", "GGA": "GLY", "GGG": "GLY",
}
# Non-standard residue names normalised to their standard amino acid (protonation
# states / selenomethionine) so a real PDB still maps to a codon.
_RESNAME_ALIASES: Dict[str, str] = {
    # Histidine protonation/tautomer names across force fields all fold to HIS:
    #   CHARMM: HSD/HSE/HSP   AMBER: HID/HIE/HIP   GROMOS/GROMACS: HISD/HISE/HISH,
    #   HISA/HISB/HISH   plus HIC (methyl-His) and the bare "HIS".
    "HSD": "HIS", "HSE": "HIS", "HSP": "HIS",
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
    "HISD": "HIS", "HISE": "HIS", "HISH": "HIS",
    "HISA": "HIS", "HISB": "HIS", "HIC": "HIS",
    # Other common non-standard names:
    "MSE": "MET", "SEC": "CYS", "CYX": "CYS", "CYM": "CYS",
    "GLH": "GLU", "ASH": "ASP", "LYN": "LYS", "ARN": "ARG",
}
_AMINO_ACIDS = {aa for aa in GENETIC_CODE.values() if aa != "STOP"}


def read_codon_table(path: str) -> Dict[str, Tuple[float, str]]:
    """Read a codon dwell-time table into ``{codon: (time_seconds, amino_acid)}``.

    Format: one row per codon, ``CODON  TIME  [AMINO_ACID]`` (whitespace/tab separated;
    ``#`` comments and blank lines ignored). Codons are upper-cased RNA (``T`` -> ``U``).
    The amino acid is taken from column 3 when present, else filled from the standard
    genetic code (:data:`GENETIC_CODE`). Stops are the amino acid ``STOP``.

    Raises
    ------
    ValueError
        If a non-comment row cannot be parsed, a codon is unknown (no column-3 amino
        acid and not in the genetic code), or the table has no rows.
    """
    table: Dict[str, Tuple[float, str]] = {}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                raise ValueError(f"{path}: cannot parse codon table line: {line!r}")
            codon = parts[0].upper().replace("T", "U")
            time = float(parts[1])
            aa = parts[2].upper() if len(parts) >= 3 else GENETIC_CODE.get(codon)
            if aa is None:
                raise ValueError(
                    f"{path}: codon {codon!r} has no amino-acid column and is not in the "
                    f"standard genetic code.")
            table[codon] = (time, aa)
    if not table:
        raise ValueError(f"{path}: no codon rows found.")
    return table


def read_protein_sequence(pdb_path: str) -> List[str]:
    """Read the ordered amino-acid sequence (3-letter names) from a PDB file.

    One entry per protein residue in file order. Residues are grouped by
    ``(chain, resSeq, iCode)``; the first ``ATOM`` record of each new residue supplies
    its name. Non-standard names (HIS protonation states, ``MSE``, ...) are normalised
    to the standard amino acid.

    Raises
    ------
    ValueError
        If a residue name is not a standard amino acid, or no ``ATOM`` records exist.
    """
    seq: List[str] = []
    seen = None
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            key = (line[21], line[22:27])            # chain, resSeq + iCode
            if key == seen:
                continue
            seen = key
            resname = line[17:20].strip().upper()
            aa = _RESNAME_ALIASES.get(resname, resname)
            if aa not in _AMINO_ACIDS:
                raise ValueError(
                    f"{pdb_path}: residue {resname!r} (chain {line[21]!r} "
                    f"{line[22:27].strip()}) is not a standard amino acid; cannot build "
                    f"a fastest/slowest mRNA.")
            seq.append(aa)
    if not seq:
        raise ValueError(f"{pdb_path}: no protein ATOM records found.")
    return seq


def select_codon_per_aa(table: Dict[str, Tuple[float, str]], mode: str) -> Dict[str, str]:
    """Map each amino acid to its fastest, slowest or median-dwell-time synonymous codon.

    For every amino acid (and the ``STOP`` group), consider only its codons present in
    ``table`` and pick the smallest- (``mode='fastest'``), largest- (``'slowest'``) or
    median- (``'median'``) dwell-time codon. When the amino acid has an even number of
    synonymous codons there is no single middle, so ``median`` takes the faster
    (shorter-time) of the two central codons -- index ``(n - 1) // 2`` of the
    time-sorted list. Ties break on the codon string, so every choice is deterministic.
    """
    if mode not in SYNTHETIC_MRNA_MODES:
        raise ValueError(f"mode must be one of {SYNTHETIC_MRNA_MODES}; got {mode!r}.")
    groups: Dict[str, List[Tuple[float, str]]] = {}
    for codon, (time, aa) in table.items():
        groups.setdefault(aa, []).append((time, codon))
    chosen: Dict[str, str] = {}
    for aa, items in groups.items():
        items.sort()                                  # by (time, codon)
        if mode == "slowest":
            chosen[aa] = items[-1][1]
        elif mode == "fastest":
            chosen[aa] = items[0][1]
        else:                                         # median (lower median if even)
            chosen[aa] = items[(len(items) - 1) // 2][1]
    return chosen


def build_synthetic_codons(aa_seq: List[str], table: Dict[str, Tuple[float, str]],
                           mode: str) -> List[str]:
    """Build the codon list for a fastest/slowest/median mRNA: one chosen codon per
    residue, plus a terminating (fastest/slowest/median) stop codon.

    Raises
    ------
    ValueError
        If an amino acid in ``aa_seq`` has no codon in ``table``, or the table has no
        stop codon to terminate the mRNA.
    """
    chosen = select_codon_per_aa(table, mode)
    codons: List[str] = []
    for i, aa in enumerate(aa_seq):
        if aa not in chosen:
            raise ValueError(
                f"amino acid {aa!r} (residue #{i + 1}) has no codon in the dwell-time "
                f"table; cannot build a {mode} mRNA.")
        codons.append(chosen[aa])
    if "STOP" not in chosen:
        raise ValueError(
            "the dwell-time table has no stop codon (UAA/UAG/UGA); one is needed to "
            "terminate the synthetic mRNA.")
    codons.append(chosen["STOP"])                     # generate the stop codon too
    return codons


def write_synthetic_mrna(pdb_path: str, codon_time_table_path: str, mode: str,
                         out_dir: str = None, filename: str = None) -> str:
    """Write a fastest/slowest/median synonymous-codon mRNA for ``pdb_path`` and return its path.

    Reads the protein sequence from ``pdb_path`` and the dwell times + codon->amino-acid
    map from ``codon_time_table_path``, encodes each residue with its chosen codon
    (plus a terminating stop), and writes the raw-nucleotide mRNA to
    ``<out_dir>/mrna_<mode>.txt`` (``out_dir`` defaults to the PDB's directory).

    Parameters
    ----------
    pdb_path : str
        Protein PDB (source of the amino-acid sequence).
    codon_time_table_path : str
        Codon dwell-time table (``codon  time  amino_acid``).
    mode : {'fastest', 'slowest', 'median'}
        Pick the fastest, slowest or median-dwell-time synonymous codon per residue.
    out_dir : str, optional
        Output directory (default: the directory containing ``pdb_path``).
    filename : str, optional
        Output filename (default: ``mrna_<mode>.txt``).

    Returns
    -------
    str
        Path to the written mRNA file.
    """
    mode = mode.strip().lower()
    if mode not in SYNTHETIC_MRNA_MODES:
        raise ValueError(f"mode must be one of {SYNTHETIC_MRNA_MODES}; got {mode!r}.")
    table = read_codon_table(codon_time_table_path)
    aa_seq = read_protein_sequence(pdb_path)
    codons = build_synthetic_codons(aa_seq, table, mode)
    seq = "".join(codons)

    out_dir_p = Path(out_dir) if out_dir else Path(pdb_path).resolve().parent
    out_dir_p.mkdir(parents=True, exist_ok=True)
    out_path = out_dir_p / (filename or f"mrna_{mode}.txt")
    header = [
        f"# {mode} synonymous-codon mRNA for {Path(pdb_path).name}",
        f"# codon dwell-time table: {codon_time_table_path}",
        f"# {len(aa_seq)} residues + 1 stop codon = {len(codons)} codons; RNA (U)",
        f"# generated by cosmo.csp.synth_mrna (mrna={mode}); one codon per residue",
    ]
    body = [seq[i:i + 60] for i in range(0, len(seq), 60)]
    out_path.write_text("\n".join(header + body) + "\n")
    return str(out_path)


def main(argv=None) -> None:
    """``cosmo-make-mrna`` CLI: write a fastest/slowest/median mRNA for a protein + table."""
    ap = argparse.ArgumentParser(
        prog="cosmo-make-mrna",
        description="Write a fastest/slowest/median synonymous-codon mRNA for a protein "
                    "(preparation step for the cosmo.csp synthesis workflow).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument("--pdb", required=True, help="protein PDB (source of the AA sequence)")
    ap.add_argument("--codon-times", required=True,
                    help="codon dwell-time table (columns: codon  time  amino_acid)")
    ap.add_argument("--mode", required=True, choices=SYNTHETIC_MRNA_MODES,
                    help="pick the fastest, slowest or median synonymous codon per residue")
    ap.add_argument("--out-dir", default=None,
                    help="output directory (default: next to the PDB)")
    ap.add_argument("-o", "--outfile", default=None,
                    help="output filename (default: mrna_<mode>.txt)")
    a = ap.parse_args(argv)
    path = write_synthetic_mrna(a.pdb, a.codon_times, a.mode,
                                out_dir=a.out_dir, filename=a.outfile)
    print(f"wrote {a.mode} mRNA -> {path}")


if __name__ == "__main__":
    main()
