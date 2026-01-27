#!/usr/bin/env python
"""
Align nascent protein along X-axis with:
- First CA atom at (0, 0, 0)
- Last CA atom at (-N*3.8, 0, 0) where N is the number of residues
"""

import sys
import numpy as np
from openmm.app import PDBFile


def align_nascent_protein(input_pdb, output_pdb, ca_spacing=3.8):
    """
    Align nascent protein along X-axis.
    
    Parameters
    ----------
    input_pdb : str
        Input PDB file path
    output_pdb : str
        Output PDB file path
    ca_spacing : float
        Spacing between CA atoms in Angstroms (default: 3.8)
    """
    # Load PDB file
    pdb = PDBFile(input_pdb)
    topology = pdb.topology
    positions = pdb.positions
    
    # Find all CA atoms
    ca_atoms = []
    ca_indices = []
    for atom in topology.atoms():
        if atom.name == 'CA':
            ca_atoms.append(atom)
            ca_indices.append(atom.index)
    
    if len(ca_atoms) < 2:
        raise ValueError(f"Need at least 2 CA atoms, found {len(ca_atoms)}")
    
    print(f"Found {len(ca_atoms)} CA atoms")
    
    # Get first and last CA atom positions
    first_ca_idx = ca_indices[0]
    last_ca_idx = ca_indices[-1]
    
    first_ca_pos = np.array([positions[first_ca_idx].x, 
                            positions[first_ca_idx].y, 
                            positions[first_ca_idx].z])
    last_ca_pos = np.array([positions[last_ca_idx].x, 
                           positions[last_ca_idx].y, 
                           positions[last_ca_idx].z])
    
    print(f"First CA (residue {ca_atoms[0].residue.name} {ca_atoms[0].residue.index}): {first_ca_pos}")
    print(f"Last CA (residue {ca_atoms[-1].residue.name} {ca_atoms[-1].residue.index}): {last_ca_pos}")
    
    # Target positions
    N_residues = len(ca_atoms)
    target_first = np.array([0.0, 0.0, 0.0])
    target_last = np.array([N_residues * ca_spacing, 0.0, 0.0])
    
    print(f"Target first CA: {target_first}")
    print(f"Target last CA: {target_last}")
    
    # Convert positions to numpy array
    pos_array = np.array([[p.x, p.y, p.z] for p in positions])
    
    # Step 1: Translate so first CA is at origin
    translation1 = -first_ca_pos
    pos_array += translation1
    
    # Update first and last CA positions after translation
    first_ca_pos_translated = pos_array[first_ca_idx]
    last_ca_pos_translated = pos_array[last_ca_idx]
    
    # Step 2: Calculate rotation to align along X-axis
    # Vector from first to last CA
    current_vector = last_ca_pos_translated - first_ca_pos_translated
    target_vector = target_last - target_first
    
    # Normalize vectors
    current_norm = np.linalg.norm(current_vector)
    target_norm = np.linalg.norm(target_vector)
    
    if current_norm < 1e-6:
        raise ValueError("First and last CA atoms are too close together")
    
    current_unit = current_vector / current_norm
    target_unit = target_vector / target_norm
    
    # Calculate rotation axis and angle using cross product
    cross = np.cross(current_unit, target_unit)
    cross_norm = np.linalg.norm(cross)
    
    if cross_norm < 1e-6:
        # Vectors are already aligned (or anti-aligned)
        if np.dot(current_unit, target_unit) < 0:
            # Anti-aligned, need 180 degree rotation around perpendicular axis
            # Use Y-axis as rotation axis
            rotation_axis = np.array([0.0, 1.0, 0.0])
            angle = np.pi
        else:
            # Already aligned, no rotation needed
            rotation_axis = None
            angle = 0.0
    else:
        rotation_axis = cross / cross_norm
        dot_product = np.clip(np.dot(current_unit, target_unit), -1.0, 1.0)
        angle = np.arccos(dot_product)
    
    # Apply rotation if needed
    if angle > 1e-6:
        print(f"Rotating by {np.degrees(angle):.2f} degrees around axis {rotation_axis}")
        
        # Rodrigues' rotation formula
        K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                     [rotation_axis[2], 0, -rotation_axis[0]],
                     [-rotation_axis[1], rotation_axis[0], 0]])
        
        R = (np.eye(3) + np.sin(angle) * K + 
             (1 - np.cos(angle)) * np.dot(K, K))
        
        # Rotate all positions around first CA
        pos_array_centered = pos_array - first_ca_pos_translated
        pos_array_rotated = np.dot(pos_array_centered, R.T)
        pos_array = pos_array_rotated + first_ca_pos_translated
    
    # Step 3: Scale to match target distance
    first_ca_final = pos_array[first_ca_idx]
    last_ca_final = pos_array[last_ca_idx]
    current_distance = np.linalg.norm(last_ca_final - first_ca_final)
    target_distance = np.linalg.norm(target_last - target_first)
    
    if current_distance > 1e-6:
        scale_factor = target_distance / current_distance
        print(f"Scaling by factor {scale_factor:.4f}")
        
        # Scale around first CA
        pos_array_centered = pos_array - first_ca_final
        pos_array_scaled = pos_array_centered * scale_factor
        pos_array = pos_array_scaled + first_ca_final
    
    # Step 4: Final translation to place first CA at (0,0,0)
    final_translation = target_first - pos_array[first_ca_idx]
    pos_array += final_translation
    
    # Verify final positions
    final_first = pos_array[first_ca_idx]
    final_last = pos_array[last_ca_idx]
    
    print(f"\nFinal first CA position: {final_first}")
    print(f"Final last CA position: {final_last}")
    print(f"Distance: {np.linalg.norm(final_last - final_first):.4f} Å")
    print(f"Target distance: {target_distance:.4f} Å")
    
    # Convert back to OpenMM positions
    from openmm.vec3 import Vec3
    new_positions = [Vec3(pos[0], pos[1], pos[2]) for pos in pos_array]
    
    # Write output PDB
    with open(output_pdb, 'w') as f:
        PDBFile.writeFile(topology, new_positions, f)
    
    print(f"\nAligned structure written to: {output_pdb}")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python align_nascent.py <input.pdb> <output.pdb> [ca_spacing]")
        print("  ca_spacing: spacing between CA atoms in Angstroms (default: 3.8)")
        sys.exit(1)
    
    input_pdb = sys.argv[1]
    output_pdb = sys.argv[2]
    ca_spacing = float(sys.argv[3]) if len(sys.argv) > 3 else 3.8
    
    try:
        align_nascent_protein(input_pdb, output_pdb, ca_spacing)
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
