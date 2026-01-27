#!/usr/bin/env python
"""
Change chain ID in a PDB file.

Usage:
    python change_chain.py <input.pdb> <output.pdb> <new_chain_id>
    python change_chain.py protein.pdb protein_new.pdb 8
"""

import sys


def change_chain_id(input_file, output_file, new_chain_id):
    """
    Change chain ID in a PDB file.
    
    Parameters
    ----------
    input_file : str
        Input PDB file path
    output_file : str
        Output PDB file path
    new_chain_id : str
        New chain ID (single character)
    """
    if len(new_chain_id) != 1:
        raise ValueError(f"Chain ID must be a single character, got: '{new_chain_id}'")
    
    changed_atoms = 0
    changed_ter = 0
    
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            # Handle ATOM and HETATM records
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # Chain ID is at column 21 (0-indexed: 20)
                if len(line) > 21:
                    old_chain = line[21]
                    # Replace chain ID
                    new_line = line[:21] + new_chain_id + line[22:]
                    outfile.write(new_line)
                    if old_chain != new_chain_id and old_chain.strip():
                        changed_atoms += 1
                else:
                    outfile.write(line)
            
            # Handle TER records
            elif line.startswith('TER'):
                # Chain ID is at column 21 (0-indexed: 20)
                if len(line) > 21:
                    old_chain = line[21]
                    # Replace chain ID
                    new_line = line[:21] + new_chain_id + line[22:]
                    outfile.write(new_line)
                    if old_chain != new_chain_id and old_chain.strip():
                        changed_ter += 1
                else:
                    outfile.write(line)
            
            # Write all other lines as-is
            else:
                outfile.write(line)
    
    print(f"Changed chain ID in PDB file:")
    print(f"  Input file: {input_file}")
    print(f"  Output file: {output_file}")
    print(f"  New chain ID: {new_chain_id}")
    print(f"  Changed ATOM/HETATM records: {changed_atoms}")
    print(f"  Changed TER records: {changed_ter}")


if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python change_chain.py <input.pdb> <output.pdb> <new_chain_id>")
        print("  new_chain_id: single character (e.g., '8', 'A', 'B')")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    new_chain_id = sys.argv[3]
    
    try:
        change_chain_id(input_file, output_file, new_chain_id)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
