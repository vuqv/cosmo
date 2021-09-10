Hydropathy Scale are the model of protein systems based on simplifications made over classical Molecular Dynamics (MD) force fields.

The hps model is a Python library that offers flexibility to set up coarse-grained simulation of IDP using the MD framework of OpenMM toolkit.
The codebase is based on sbmOpenMM scripts.
It automates the creation of :code:`openmm.system` classes that contain the necessary force field parameters to run molecular dynamics simulations using a protein structure as the only necessary inputs.

hps is divided in three main classes:

1. :code:`geometry`
2. :code:`models`
3. :code:`system`
   
The first class, :code:`geometry`, contains methods to calculate the geometrical parameters from the input structures.
It's not useful in current need of simulation method.
The second class, :code:`models`, allows to easily set up CG models.
The third class, :code:`system`, is the main class that holds all the methods to define, modify and create CG system to be simulated with OpenMM.

The library is open-source and offers flexibility to simulate IDPs.

