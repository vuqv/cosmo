#!/usr/bin/env python
"""
Align CA atoms along the X-axis and write a new PDB with only CA atoms.

Usage:
    python align_and_change_chain.py input.pdb output.pdb resid_to_set_at_origin

Behaviour:
- Reads the input PDB and keeps only CA atoms.
- Places all CA atoms along the X-axis with spacing of 3.8 Å.
- The CA with residue ID `resid_to_set_at_origin` is placed at the origin (0, 0, 0).
- All CA atoms are assigned to chain ID '8'.
"""

import sys
from typing import List, Tuple


CA_SPACING = 3.8
NEW_CHAIN_ID = "8"


def parse_ca_atoms(pdb_lines: List[str]) -> List[Tuple[int, str, str]]:
    """
    Extract CA atoms from PDB lines.

    Returns a list of (resid, line, record_name) tuples in file order.
    """
    ca_atoms = []
    for line in pdb_lines:
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        # Atom name is columns 13-16 (1-based), i.e. [12:16]
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        resid_str = line[22:26].strip()
        if not resid_str:
            continue
        try:
            resid = int(resid_str)
        except ValueError:
            continue
        record_name = line[0:6]
        ca_atoms.append((resid, line.rstrip("\n"), record_name))
    return ca_atoms


def format_ca_line(serial: int, resid: int, template_line: str, x: float, y: float, z: float) -> str:
    """
    Build a PDB ATOM line for a CA atom using fields from template_line and new coordinates.
    """
    record = template_line[0:6].strip() or "ATOM"
    atom_name = template_line[12:16]
    alt_loc = template_line[16]
    res_name = template_line[17:20]
    # Force new chain ID
    chain_id = NEW_CHAIN_ID
    i_code = template_line[26]

    # Occupancy / tempFactor (keep if present, else defaults)
    occupancy = template_line[54:60].strip() or "1.00"
    temp_factor = template_line[60:66].strip() or "0.00"

    element = template_line[76:78] if len(template_line) >= 78 else "  "
    charge = template_line[78:80] if len(template_line) >= 80 else "  "

    line = (
        f"{record:<6s}"
        f"{serial:5d} "
        f"{atom_name:^4s}"
        f"{alt_loc:1s}"
        f"{res_name:>3s} "
        f"{chain_id:1s}"
        f"{resid:4d}"
        f"{i_code:1s}"
        f"   "
        f"{x:8.3f}"
        f"{y:8.3f}"
        f"{z:8.3f}"
        f"{occupancy:>6s}"
        f"{temp_factor:>6s}"
        f"          "
        f"{element:>2s}"
        f"{charge:>2s}"
        "\n"
    )
    return line


def align_ca_atoms(input_pdb: str, output_pdb: str, resid_to_set_at_origin: int) -> None:
    # Read input PDB
    with open(input_pdb, "r") as f:
        lines = f.readlines()

    ca_atoms = parse_ca_atoms(lines)

    if not ca_atoms:
        raise ValueError("No CA atoms found in input PDB.")

    # Find index of CA with the specified residue ID
    origin_index = None
    for i, (resid, _, _) in enumerate(ca_atoms):
        if resid == resid_to_set_at_origin:
            origin_index = i
            break

    if origin_index is None:
        raise ValueError(
            f"Residue ID {resid_to_set_at_origin} not found among CA atoms."
        )

    # Generate new PDB with only CA atoms aligned along X-axis
    with open(output_pdb, "w") as out:
        serial = 1
        for i, (resid, template_line, record_name) in enumerate(ca_atoms):
            # Position relative to origin_index, spacing CA_SPACING along X
            x = (origin_index - i) * CA_SPACING
            y = 4.0
            z = 3.0
            out.write(format_ca_line(serial, resid, template_line, x, y, z))
            serial += 1

        # Write TER and END records
        out.write(f"TER   {serial:5d}      CA  {NEW_CHAIN_ID}{resid_to_set_at_origin:4d}\n")
        out.write("END\n")


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(
            "Usage: python align_and_change_chain.py input.pdb output.pdb resid_to_set_at_origin"
        )
        sys.exit(1)

    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    try:
        resid_to_set_at_origin = int(sys.argv[3])
    except ValueError:
        print("resid_to_set_at_origin must be an integer residue ID.")
        sys.exit(1)

    try:
        align_ca_atoms(input_pdb, output_pdb, resid_to_set_at_origin)
    except Exception as e:
        print(f"Error: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
