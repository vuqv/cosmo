Introduction 
=========================================================

The hpsOpenMM model is a Python library that offers flexibility to set up coarse-grained simulation of IDP using the MD framework of OpenMM toolkit.
The codebase is based on sbmOpenMM scripts.
It automates the creation of :code:`openmm.system` classes that contain the necessary force field parameters to run molecular dynamics simulations using a protein structure as the only necessary inputs.

hpsOpenMM is divided in four main classes:

1. :code:`system`
2. :code:`models`
3. :code:`dynamics`
4. :code:`geometry`

------------------

:code:`system`, is the main class that holds all the methods to define,
modify and create CG system to be simulated with OpenMM.
Class inheritance from :code:`openmm.system` with some more attributes for hpsOpenMM.

------------------

:code:`models`, class contains set of models, each model contains a collection of sequence of commands
to build model, allows to easily set up CG models.

------------------

:code:`dynamics` class auto read the parameter controls, build the model and run simulation.

------------------

:code:`geometry`, contains methods to calculate the geometrical parameters from the input structures.
It's not useful in current need of simulation method.

------------------

The library is open-source and offers flexibility to simulate IDPs.

