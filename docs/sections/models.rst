The models class of sbmOpenMM contains three methods for automatic setting up predefined SBM potentials. It works by initializing a system class with the necessary force field parameters, derived from the input files, to set up one of the possible models which are detailed next:

Coarse grained, alpha-carbon (CA), model
++++++++++++++++++++++++++++++++++++++++

The coarse grained method represents the protein system as beads centered at the alpha carbons of each residue in the protein. It uses harmonic potentials to hold the covalent connectivity and geometry of the beads. Torsional geometries are modeled with a periodic torsion potential. Native contacts are represented through the use of Lennard-Jones potentials that allow to form and break non-bonded interactions, permitting complete and local unfolding of the structures.

To create a CA model, call:

sbmOpenMM.models.getCAModel(pdb_file, contacts_file)

Here, pdb_file is the path to the PDB format structure of the protein and  contacts_file is the path to the contact file containing only the CA atoms of the system. This last file should be numbered considering the CA atoms consecutively.

The force field equations are:

.. math::
	H_A = \sum_{bonds}V_{bond}+\sum_{pairwise}V_{LJ_{12-6}}+\sum_{electrostatics}V_{DH}

.. math::
        V_{bond} = \frac{k_b}{2}(r-r_0)^2


.. math::
        V_{LJ_{12-6}} = \epsilon((\frac{\sigma_{ij}}{r})^{12}-6(\frac{\sigma_{ij}}{r})^{6})

.. math::
        V_{LJ_{12}} = \epsilon_{nc}(\frac{\sigma_{nc}}{r})^{12}


Here the default values are :math:`k_b=20000\ kJ/(mol \cdot nm^2)`, :math:`k_a=40\ kJ/(mol \cdot rad^2)`, :math:`k_t=1.0\ kJ/mol`, :math:`\epsilon_{c}=1.0\ kJ/mol`, :math:`\epsilon_{nc}=1.0\ kJ/mol` and :math:`\sigma_{nc}=0.4\ nm`. The geometric parameters are set to the calculated structural values in the input structure, with :math:`r_0` the equilibrium bond distance in nanometers, :math:`\theta_0` the equilibrium angle length in radians, :math:`\phi_0` the equilibrium torsional angle in radians and :math:`\sigma_{ij}` the equilibrium contact distance in nanometers. The variable :math:`r` represents, accordingly, the current bond or (non)contact distance in nanometers, :math:`\theta` the current angle length in radians and :math:`\phi` the current torsional angle in radians.

It is possible to use a :math:`V_{LJ_{12-10-6}}` potential for the native contact interactions, defined as:

.. math::
        V_{LJ_{12-10-6}} = \epsilon_{c}(13(\frac{\sigma_{ij}}{r})^{12}-18(\frac{\sigma_{ij}}{r})^{10}+4(\frac{\sigma_{ij}}{r})^{6})

This potential gives a small energy barrier for contact formation/breaking that emulates a "desolvation effect". To use this potential as the native contact energy function, instead of the :math:`V_{LJ_{12-10}}` potential, give the option contact_force ='12-10-6' to the sbmOpenMM.models.getCAModel() method. 
 
Note that even if the units for the force constants are given in real physical units (e.g. :math:`kJ/mol`), this is just to match the variables used by OpenMM. The models are not parametrized to equate this real physical values and comparison with experiments will require further adjustment to the energy unit system. 

