import numpy as np
import pandas as pd
from PeptideBuilder import Geometry
import PeptideBuilder
import Bio.PDB
from mdtraj.utils.rotation import rotation_matrix_from_quaternion
import MDAnalysis as mda

fasta = 'MASNDYTQQATQSYGAYPTQPGQGYSQQSSQPYGQQSYSGYSQSTDTSGYGQSSYSSYGQSQNTGYGTQSTPQGYGSTGGYGSSQSSQSSYGQQSSYPGYGQQPAPSSTSGSYGSSSQSSSYGQPQSGSYSQQPSYGGQQQSYGQQQSYNPPQGYGQQNQYNS'

geo = Geometry.geometry(fasta[0])
geo.phi = -120
geo.psi_im1 = 150

structure = PeptideBuilder.initialize_res(geo)
for residue in fasta[1:]:
    structure = PeptideBuilder.add_residue(structure, residue, geo.phi, geo.psi_im1)


out = Bio.PDB.PDBIO()
out.set_structure(structure)
CA_coords = []
for residue in out.structure.get_residues():
    # print(residue.get_resname())
    atoms = residue.get_atoms()
    for atom in atoms:
        if atom.name == 'CA':
            CA_coords.append([atom.coord[0], atom.coord[1], atom.coord[2]])
            
CA_coords = np.array(CA_coords)

xyz = []
for atom in out.structure.get_atoms():
    xyz.append([atom.coord[0],atom.coord[1],atom.coord[2]])
xyz = np.array(xyz)
v = CA_coords[-1] - CA_coords[0]
u = np.array([0,0,1])
a = np.cross(v,u) 
a = a / np.linalg.norm(a,keepdims=True)
b = np.arccos( np.dot(v,u) / np.linalg.norm(v) )
quaternion = np.insert(np.sin(-b/2).reshape(-1,1)*a,0,np.cos(-b/2),axis=1)
newxyz = xyz - np.mean(xyz,axis=0)
newxyz = np.matmul(newxyz,rotation_matrix_from_quaternion(quaternion)) 
xyz = np.array(newxyz[0])

xyz -= np.min(xyz, axis=0)

for i, atom in enumerate(structure.get_atoms()):
    atom.coord = xyz[i]

final_structure = Bio.PDB.PDBIO()
final_structure.set_structure(structure)
out.save("example_1chain.pdb")

protein = mda.Universe('example_1chain.pdb')
n_atoms = len(list(protein.atoms))
box_len = np.max(protein.atoms.positions, axis=0)-np.min(protein.atoms.positions, axis=0)
print(box_len)
protein.dimensions = [box_len[0], box_len[1], box_len[2], 90, 90, 90]
chain_list= ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

def tile_universe(universe, n_x, n_y):
    box = universe.dimensions[:3]
    print(box)
    copied = []
    for x in range(n_x):
        for y in range(n_y):
            u_ = universe.copy()
            # move_by = box*(x, y, 1)
            move_by = np.array([10, 10, 10]) *(x, y, 0)
            print(move_by)
            u_.atoms.translate(move_by)
            u_.add_TopologyAttr('chainID', n_atoms*[chain_list[(x*n_x+y)%len(chain_list)]])
            copied.append(u_.atoms)

    new_universe = mda.Merge(*copied)
    new_box = box*(n_x, n_y, 1)
    # new_box = np.array([10, 10, 10])*(n_x, n_y, 1)
    new_universe.dimensions = list(new_box) + [90]*3
    return new_universe

tiled = tile_universe(protein, 20, 5)
all_atom = tiled.select_atoms('all')
all_atom.write('structure.pdb')