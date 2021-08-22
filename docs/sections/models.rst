The models class of sbmOpenMM contains three methods for automatic setting up predefined SBM potentials. It works by initializing a system class with the necessary force field parameters, derived from the input files, to set up one of the possible models which are detailed next:

Coarse grained, alpha-carbon (CA), model
++++++++++++++++++++++++++++++++++++++++

The coarse grained method represents the protein system as beads centered at the alpha carbons of each residue in the protein. It uses harmonic potentials to hold the covalent connectivity and geometry of the beads. Torsional geometries are modeled with a periodic torsion potential. Native contacts are represented through the use of Lennard-Jones potentials that allow to form and break non-bonded interactions, permitting complete and local unfolding of the structures.

To create a CA model, call:
:code:`sbmOpenMM.models.getCAModel(pdb_file)`

Here, pdb_file is the path to the PDB format structure of the protein.

The force field equations are:

.. math::
	H_A = \sum_{bonds}V_{bond}+\sum_{i,j}\Phi_{ij}^{vdw}+\sum_{i,j}\Phi_{i,j}^{el}

The Bonded potential:
++++
.. math::
        V_{bond} = \frac{k_b}{2}(r-r_0)^2

Here the default values are :math:`k_b=20 kCal/(mol \times A^2), r_0=3.82 A^2`

The Pairwise potential:
++++++++++++++++++++++++++++++++++++++++

.. math::
        \Phi_{i,j}^{vdw}(r) = step(2^{1/6}\sigma_{ij}-r) \times \left( 4\epsilon\left[\left(\frac{\sigma_{ij}}{r}\right)^{12}- \left(\frac{\sigma_{ij}}{r}\right)^{6}\right]+(1-\mu\times\lambda_{ij}^{0}+\Delta)\times\epsilon\right)

        + \left[1-step(2^{1/6}\sigma_{ij}-r)\right]\times\left[(\mu \lambda_{ij}^{0}-\Delta)\times 4\pi \left[\left(\frac{\sigma_{ij}}{r}\right)^{12}-\left(\frac{\sigma_{ij}}{r}\right)^6\right]\right]

Since the step function behaves like: :code:`step(x) = 0 if x <= 0,and =1 otherwise`, we can separate in multiple case for short likes following:

.. math::
        \Phi_{i,j}^{vdw}(r) =  4\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12}-\left(\frac{\sigma_{ij}}{r}\right)^{6}\right]+(1-\mu	\times\lambda_{ij}^{0}+\Delta)	\times\epsilon, r\le 2^{1/6}\sigma_{ij}

        \Phi_{i,j}^{vdw}(r) = (\mu\times\lambda_{ij}^{0}-\Delta) \times \left( 4\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12}-\left(\frac{\sigma_{ij}}{r}\right)^{6}\right]\right), r > 2^{1/6}\sigma_{ij}

where, :math:`\sigma_{i,j}=\frac{\sigma_i+\sigma_j}{2}`: is the vdW radius interaction of interacting beads

:math:`\lambda_{ij}^{0}=\frac{\lambda_i+\lambda_j}{2}`: hydropathy scale interaction of residues

:math:`\mu, \Delta`: are the only free parameters in the model. In Jeetain Mittal(2021) Protein Science, he simulated for 42 IDP proteins and fit Rg vs experimental values.

In the current implementation, hydropathy scales are taken from Urry model, :math:`(\mu, \Delta) = (1, 0.08)`

The Debye-Huckle potential has following form:
++++++++++++++++++++++++++++++++++++++++
.. math::
        \Phi_{ij}^{el}(r) = \frac{q_{i}q_{j}}{4\pi\epsilon_0 D r}e^{-\kappa r}

where, :math:`q_i, q_j` are charge of residues :math:`i, j`

:math:`\epsilon_0`: Vacuum permitivity. For convenient, we precalculated the electric conversion factor
:math:`\frac{1}{4\pi\epsilon_0}= 138.935 485(9) kJ \times mol^{−1} \times nm \times e^{−2}`.

:math:`D`: dielectric constant, at 100mM monovalent salt (NaCl), it takes values of 80

:math:`\kappa`: inverse Debye length, at 100mM NaCl has values of 1nm

