# Optional: crop full ribosome PDB to a smaller region for faster simulation
python -m cosmo.utils.crop_ribosome rib_full.pdb rib.pdb

python grow_nascent.py -f md.ini --rib rib.pdb --nascent nascent.pdb --steps 50000
