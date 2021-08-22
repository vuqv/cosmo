Hydropathy Scale-Urryl (hps-urry) are representations of protein systems based on simplifications made over classical Molecular Dynamics (MD) force fields. 

The hps-urry model is a Python library that offers flexibility to set up coarse-grained simulation of IDP using the MD framework of OpenMM toolkit.
The codebase is based on sbmOpenMM.
It automates the creation of openmm.system classes that contain the necessary force field parameters to run molecular dynamics simulations using a protein structure as the only necessary inputs.

sbmOpenMM is divided in three main classes:

1. geometry
2. models
3. system
   
The first class, geometry, contains methods to calculate the geometrical parameters from the input structures. 
It's not useful in current need of simulation method.
The second class, models, allows to easily set up CG models.
The third class, system, is the main class that holds all the methods to define, modify and create CG system to be simulated with OpenMM.

The library is open-source and offers flexibility to simulate IDPs.
We test on RTX 2060 with the timestep of 30fs, we can get upto :math:`80\mu s/day`.
