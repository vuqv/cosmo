"""
Using PyRosetta to build pdb (extended conformation) from FASTA sequence.
it is easier and more robust than ambertools since ambertools limited to about 200 residues.
"""
from pyrosetta import init, pose_from_sequence
from parmed import load_rosetta
init()
seq = 'MSDAAVDTSSEITTKDLKEKKEVVEEAENGRDAPANGNAENEENGEQEADNEVDEEEEEGGEEEEEEEEGDGEEEDGDEDEEAESATGKRAAEDDEDDDVDTKKQKTDEDD'
pose = pose_from_sequence(seq)
struct = load_rosetta(pose)
struct.save('K25.pdb') 
