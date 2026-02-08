#!/usr/bin/env python
"""
Test script: parse nascent sequence from nascent.pdb and generate mRNA sequences
(fast/slow) using codon translation time from cosmo.utils.ctf_utils.
"""

import os
import sys

# Allow imports when run from repo root or from examples/growing
_script_dir = os.path.dirname(os.path.abspath(__file__))
_repo_root = os.path.dirname(os.path.dirname(_script_dir))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)
if _script_dir not in sys.path:
    sys.path.insert(0, _script_dir)

from cosmo.utils import parse_nascent_sequence, generate_mRNA_sequence


def main():
    nascent_pdb = os.path.join(_script_dir, 'nascent.pdb')
    if not os.path.isfile(nascent_pdb):
        print(f"Error: {nascent_pdb} not found.")
        sys.exit(1)

    print("Parsing nascent sequence from", nascent_pdb)
    sequence = parse_nascent_sequence(nascent_pdb)
    print(f"Sequence: {sequence}")
    # One (res_name, res_num) per residue; take residue names in order
    aa_three = [res_name for res_name, _ in sequence]
    n_res = len(aa_three)
    print(f"Found {n_res} residues")
    print(f"Amino acid sequence (3-letter): {' '.join(aa_three[:10])}{'...' if n_res > 10 else ''}\n")

    for organism in ('yeast', 'ecoli'):
        print(f"--- {organism.upper()} ---")
        fast_mRNA, fast_time_ms = generate_mRNA_sequence(aa_three, 'fast', organism=organism)
        slow_mRNA, slow_time_ms = generate_mRNA_sequence(aa_three, 'slow', organism=organism)
        print(f"  Fast mRNA length: {len(fast_mRNA)} nt ({len(fast_mRNA)//3} codons)")
        print(f"  Fast estimated translation time: {fast_time_ms:.1f} ms")
        print(f"  Fast estimated translation time in simulation: {(10**3)*fast_time_ms/4331293:.1f} microseconds")
        print(f"  Slow mRNA length: {len(slow_mRNA)} nt ({len(slow_mRNA)//3} codons)")
        print(f"  Slow estimated translation time: {slow_time_ms:.1f} ms")
        print(f"  Slow estimated translation time in simulation: {(10**3)*slow_time_ms/4331293:.1f} microseconds")
        print(f"  Fast (first 30 nt): {fast_mRNA[:30]}...")
        print(f"  Slow (first 30 nt): {slow_mRNA[:30]}...")
        if fast_mRNA != slow_mRNA:
            print(f"  (Fast and slow mRNA differ)")
        print()

    print("Done.")


if __name__ == '__main__':
    main()
