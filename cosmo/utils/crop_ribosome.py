#!/usr/bin/env python
"""
Crop ribosome PDB to a smaller region for translation/growing nascent chain simulations.

Filter PDB file to keep atoms based on:
1) X < 0: remove
2) 0 <= X < 60: remove if sqrt(y^2+z^2) > 40 Å
3) X >= 60: keep all

Use as executable:
  python -m cosmo.utils.crop_ribosome <input.pdb> <output.pdb> [radius] [x_threshold]

Or import:
  from cosmo.utils import filter_pdb
  filter_pdb("rib.pdb", "rib_cropped.pdb")
"""

import sys
import math


def filter_pdb(input_file, output_file, radius=40.0, x_threshold=60.0):
    """
    Filter PDB file based on X coordinate and distance from X-axis.

    Remove atoms if:
    1) X < 0
    2) 0 <= X < x_threshold AND sqrt(y^2+z^2) > radius
    3) X >= x_threshold: keep all

    Parameters
    ----------
    input_file : str
        Input PDB file path
    output_file : str
        Output PDB file path
    radius : float
        Cylinder radius in Angstroms (default: 40.0)
    x_threshold : float
        X coordinate threshold in Angstroms (default: 60.0)
    """
    kept_atoms = 0
    total_atoms = 0
    chains_with_atoms = set()  # Track chains that have kept atoms
    new_atom_index = 1  # Start reindexing from 1
    last_ter_serial = {}  # Track last atom serial for each chain in TER records

    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Handle TER records - only write if chain has kept atoms
            if line.startswith('TER'):
                # Extract chain ID from TER record (column 21, 0-indexed: 20)
                chain_id = line[21:22].strip() if len(line) > 21 else ''
                if chain_id in chains_with_atoms:
                    # Update TER serial number to last atom of this chain
                    if chain_id in last_ter_serial:
                        ter_serial = last_ter_serial[chain_id]
                        # Update serial number (columns 7-11, 0-indexed: 6-11)
                        new_line = line[:6] + f"{ter_serial:5d}" + line[11:]
                        outfile.write(new_line)
                continue

            # Write other non-ATOM lines (REMARK, CRYST1, END, etc.)
            if not line.startswith('ATOM') and not line.startswith('HETATM'):
                outfile.write(line)
                continue

            # Parse ATOM/HETATM line
            if line.startswith('ATOM') or line.startswith('HETATM'):
                total_atoms += 1

                # Extract coordinates (columns 31-54 in standard PDB format)
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except ValueError:
                    # If parsing fails, keep the line but still reindex
                    new_line = line[:6] + f"{new_atom_index:5d}" + line[11:]
                    outfile.write(new_line)
                    kept_atoms += 1
                    # Try to extract chain ID for tracking
                    chain_id = line[21:22].strip() if len(line) > 21 else ''
                    if chain_id:
                        chains_with_atoms.add(chain_id)
                        last_ter_serial[chain_id] = new_atom_index
                    new_atom_index += 1
                    continue

                # Calculate distance from X-axis (cylinder radius)
                distance_from_x_axis = math.sqrt(y**2 + z**2)

                # Extract chain ID (column 21, 0-indexed: 20)
                chain_id = line[21:22].strip() if len(line) > 21 else ''

                # Filter conditions (remove atoms if):
                # 1) X < 0: remove
                # 2) 0 <= X < 60: remove if sqrt(y^2+z^2) > 40
                # 3) X >= 60: keep all
                if x < 0:
                    # Condition 1: Remove X < 0
                    continue
                elif 0 <= x < x_threshold:
                    # Condition 2: Keep only if within cylinder
                    if distance_from_x_axis <= radius:
                        # Reindex atom serial number (columns 7-11, 0-indexed: 6-11)
                        new_line = line[:6] + f"{new_atom_index:5d}" + line[11:]
                        outfile.write(new_line)
                        kept_atoms += 1
                        chains_with_atoms.add(chain_id)
                        last_ter_serial[chain_id] = new_atom_index
                        new_atom_index += 1
                else:  # x >= x_threshold
                    # Condition 3: Keep all
                    # Reindex atom serial number (columns 7-11, 0-indexed: 6-11)
                    new_line = line[:6] + f"{new_atom_index:5d}" + line[11:]
                    outfile.write(new_line)
                    kept_atoms += 1
                    chains_with_atoms.add(chain_id)
                    last_ter_serial[chain_id] = new_atom_index
                    new_atom_index += 1

    print(f"Filtered PDB file:")
    print(f"  Input file: {input_file}")
    print(f"  Output file: {output_file}")
    print(f"  Total atoms: {total_atoms}")
    print(f"  Kept atoms: {kept_atoms} ({100*kept_atoms/total_atoms:.1f}%)")
    print(f"  Atoms reindexed: 1 to {new_atom_index - 1}")
    print(f"  Chains with kept atoms: {sorted(chains_with_atoms) if chains_with_atoms else 'None'}")
    print(f"  Filter criteria:")
    print(f"    - Remove if X < 0")
    print(f"    - Remove if 0 <= X < {x_threshold} AND distance_from_X_axis > {radius} Å")
    print(f"    - Keep all if X >= {x_threshold}")


def main():
    """CLI entry point for cropping ribosome PDB."""
    if len(sys.argv) < 3:
        print("Usage: python -m cosmo.utils.crop_ribosome <input.pdb> <output.pdb> [radius] [x_threshold]")
        print("  radius: cylinder radius in Angstroms (default: 40.0)")
        print("  x_threshold: X coordinate threshold in Angstroms (default: 60.0)")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    radius = float(sys.argv[3]) if len(sys.argv) > 3 else 40.0
    x_threshold = float(sys.argv[4]) if len(sys.argv) > 4 else 60.0

    filter_pdb(input_file, output_file, radius, x_threshold)


if __name__ == '__main__':
    main()
