import MDAnalysis as mda
from MDAnalysis.transformations.translate import center_in_box

# import nglview as nv

u = mda.Universe('asyn.psf', 'asyn_equil.dcd')

# create a transformation workflow
ag = u.select_atoms('all')
center_transformation = center_in_box(ag, point=[0, 0, 0])
# apply transform
u.trajectory.add_transformations(center_transformation)

# write new trajectories
with mda.Writer('centered.dcd', ag.n_atoms) as W:
    for ts in u.trajectory:
        W.write(ag)
