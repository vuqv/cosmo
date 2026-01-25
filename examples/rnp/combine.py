from pathlib import Path

def parse_pdb_atom_line(line: str):
    """
    Parse an ATOM/HETATM line using PDB fixed columns (with a safe fallback).
    Returns a dict with all fields needed to rewrite a clean PDB line.
    """
    rec = line[0:6].strip()
    if rec not in {"ATOM", "HETATM"}:
        return None

    # Fixed-column parse (PDB v3.x)
    try:
        atom_serial = int(line[6:11])
        atom_name   = line[12:16]
        alt_loc     = line[16:17]
        res_name    = line[17:20]
        chain_id    = line[21:22]
        res_seq     = int(line[22:26])
        ins_code    = line[26:27]
        x           = float(line[30:38])
        y           = float(line[38:46])
        z           = float(line[46:54])
        occ_str     = line[54:60].strip()
        tf_str      = line[60:66].strip()
        occupancy   = float(occ_str) if occ_str else 0.00
        tempfactor  = float(tf_str) if tf_str else 0.00
        segid       = line[72:76] if len(line) >= 76 else "    "
        element     = line[76:78] if len(line) >= 78 else "  "
        charge      = line[78:80] if len(line) >= 80 else "  "
    except Exception:
        # Fallback: whitespace split (less reliable, but helps with nonstandard formatting)
        parts = line.split()
        # Expected minimal: ATOM serial name resName chain resSeq x y z occ tf [segid] element
        # Many variants exist, so we keep this very conservative.
        atom_serial = int(parts[1])
        atom_name   = f"{parts[2]:<4}"[:4]
        alt_loc     = " "
        res_name    = f"{parts[3]:>3}"[:3]
        chain_id    = parts[4][0]
        res_seq     = int(parts[5])
        ins_code    = " "
        x, y, z     = map(float, parts[6:9])
        occupancy   = float(parts[9]) if len(parts) > 9 else 0.00
        tempfactor  = float(parts[10]) if len(parts) > 10 else 0.00
        segid       = "    "
        element     = (parts[-1] if len(parts[-1]) <= 2 else "  ").rjust(2)
        charge      = "  "

    return {
        "record": rec,
        "atom_name": atom_name,   # keep spacing
        "alt_loc": alt_loc if alt_loc else " ",
        "res_name": res_name,
        "chain_id": chain_id if chain_id else " ",
        "res_seq": res_seq,
        "ins_code": ins_code if ins_code else " ",
        "x": x, "y": y, "z": z,
        "occupancy": occupancy,
        "tempfactor": tempfactor,
        "segid": segid if segid else "    ",
        "element": element if element else "  ",
        "charge": charge if charge else "  ",
    }

def format_pdb_atom_line(atom_serial: int, atom_name: str, res_name: str, chain_id: str,
                         res_seq: int, x: float, y: float, z: float,
                         occupancy: float, tempfactor: float,
                         alt_loc: str=" ", ins_code: str=" ",
                         segid: str="    ", element: str="  ", charge: str="  ",
                         record: str="ATOM") -> str:
    """
    Write a standard PDB ATOM/HETATM line (80 columns).
    """
    # Ensure correct widths
    atom_name = f"{atom_name:<4}"[:4]
    res_name  = f"{res_name:>3}"[:3]
    chain_id  = (chain_id or " ")[:1]
    alt_loc   = (alt_loc or " ")[:1]
    ins_code  = (ins_code or " ")[:1]
    segid     = f"{segid:<4}"[:4]
    element   = f"{element:>2}"[:2]
    charge    = f"{charge:>2}"[:2]

    return (
        f"{record:<6}{atom_serial:>5} "
        f"{atom_name}{alt_loc}"
        f"{res_name} {chain_id}"
        f"{res_seq:>4}{ins_code}   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}"
        f"{occupancy:>6.2f}{tempfactor:>6.2f}          "
        f"{segid}{element}{charge}"
    )

def load_atoms(pdb_path: Path):
    atoms = []
    for line in pdb_path.read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            a = parse_pdb_atom_line(line)
            if a is not None:
                atoms.append(a)
    return atoms

def renumber_residues(atoms, new_chain_id: str, start_resseq: int = 1):
    """
    Renumber residues sequentially in the given atom list.
    A residue boundary is detected by (old_chain_id, old_res_seq, old_ins_code, res_name) changes.
    """
    new_atoms = []
    res_counter = start_resseq - 1
    prev_key = None

    for a in atoms:
        key = (a["chain_id"], a["res_seq"], a["ins_code"], a["res_name"])
        if key != prev_key:
            res_counter += 1
            prev_key = key

        b = dict(a)
        b["chain_id"] = new_chain_id
        b["res_seq"] = res_counter
        b["ins_code"] = " "
        new_atoms.append(b)

    return new_atoms

def write_combined(rna_atoms, prot_atoms, out_path: Path):
    lines = []
    atom_serial = 1

    for a in rna_atoms + prot_atoms:
        lines.append(format_pdb_atom_line(
            atom_serial=atom_serial,
            atom_name=a["atom_name"],
            res_name=a["res_name"],
            chain_id=a["chain_id"],
            res_seq=a["res_seq"],
            x=a["x"], y=a["y"], z=a["z"],
            occupancy=a["occupancy"],
            tempfactor=a["tempfactor"],
            alt_loc=a["alt_loc"],
            ins_code=a["ins_code"],
            segid=a["segid"],
            element=a["element"],
            charge=a["charge"],
            record=a["record"],
        ))
        atom_serial += 1

    lines.append("END")
    out_path.write_text("\n".join(lines) + "\n")

def main():
    rna_in = Path("rna.pdb")
    prot_in = Path("protein.pdb")
    out = Path("combined.pdb")

    rna_atoms = load_atoms(rna_in)
    prot_atoms = load_atoms(prot_in)

    rna_atoms = renumber_residues(rna_atoms, new_chain_id="A", start_resseq=1)
    prot_atoms = renumber_residues(prot_atoms, new_chain_id="B", start_resseq=1)

    write_combined(rna_atoms, prot_atoms, out)
    print(f"Wrote {out} with {len(rna_atoms) + len(prot_atoms)} atoms.")

if __name__ == "__main__":
    main()
