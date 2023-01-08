Models 
=========================================================

The models class contains three methods for automatic setting up predefined potentials.

It works by initializing a system class with the necessary force field parameters.

Coarse grained, alpha-carbon (CA), model
++++++++++++++++++++++++++++++++++++++++

The coarse grained method represents the protein system as beads centered at the alpha carbons of each residue in the protein. 

It uses harmonic potentials to hold the covalent connectivity and geometry of the beads. 

Torsional geometries are modeled with a periodic torsion potential. 

Native contacts are represented through the use of Lennard-Jones potentials that allow to form and break non-bonded interactions, permitting complete and local unfolding of the structures.

To create a CA model, call:
:code:`hps.models.getCAModel(pdb_file, hps_scale)`

Here, pdb_file is the path to the PDB format structure of the protein.
hps_scale is hydropathy scale that are going to be used. :code:`urry` or :code:`kr`

The force field equations are:

.. math::
	H_A = \sum_{bonds}V_{bond}+\sum_{i,j}\Phi_{ij}^{vdw}+\sum_{i,j}\Phi_{i,j}^{el}

If hps_ss model is used, the Hamiltonian is:

.. math::
	H_{hps-ss} = \sum_{bonds}V_{bond}+\sum_{angle}V_{angle}+\sum_{torsion}V_{torsion}+\sum_{i,j}\Phi_{ij}^{vdw}+\sum_{i,j}\Phi_{i,j}^{el}


The Bonded potential:
++++++++++++++++++++++

.. math::
        V_{bond} = \frac{k_b}{2}(r-r_0)^2

Here the default values are:

    :math:`k_b= 8368 kJ/(mol \times nm^2),\\
    r_0=0.382 nm`

Angle Potential
+++++++++++++++
.. math::
    U_{angle}(\theta) = \frac{-1}{\gamma}
    \ln \left[ e^{ -\gamma[ k_{\alpha} (\theta-\theta_{\alpha})^2+\epsilon_{\alpha} ]} +e^{ -\gamma k_{\beta} (\theta-\theta_{\beta})^2 } \right]

Parameters:

    :math:`\gamma = 0.1 mol/kcal,\\
    \epsilon_{\alpha}=4.3 kcal/mol,\\
    \theta_{\alpha}=1.6 rad, \\
    \theta_{\beta}=2.27 rad`

Torsion Potential
++++++++++++++++++

.. math::
    U_{torsion}(\theta) = -\ln\left[ U_{torsion, \alpha}(\theta, \epsilon_d) + U_{torsion, \beta}(\theta, \epsilon_d)\right]


    U_{torsion, \alpha}(\theta, \epsilon_d)  = e^{-k_{\alpha, 1}(\theta-\theta_{\alpha,1})^2-\epsilon_d}
                                                + e^{-k_{\alpha, 2}(\theta-\theta_{\alpha,2})^4 + e_0}
                                                + e^{-k_{\alpha, 2}(\theta-\theta_{\alpha,2}+2\pi)^4 + e_0}

    U_{torsion, \beta}(\theta, \epsilon_d) = e^{-k_{\beta,1}(\theta-\theta_{\beta,1})^2+e_1+\epsilon_d}
                                           + e^{-k_{\beta,1}(\theta-\theta_{\beta,1}-2\pi)^2+e_1+\epsilon_d} \\
                                           + e^{-k_{\beta,2}(\theta-\theta_{\beta,2})^4+e_2}
                                           + e^{-k_{\beta,2}(\theta-\theta_{\beta,2}-2\pi)^4+e_2}


Parameters:

    :math:`k_{\alpha,1}=11.4 kcal/(mol \times rad^2),\\
    k_{\alpha,2}=0.15 kcal/(mol\times rad^4),\\
    \theta_{\alpha,1} = 0.9 rad,\\
    \theta_{\alpha,2}=1.02 rad,\\
    e_0 = 0.27 kcal/mol,\\
    k_{\beta,1}=1.8kcal/(mol \times rad^2),\\
    k_{\beta,2}=0.65kcal/(mol\times rad^4),\\
    \theta_{\beta,1}=-1.55 rad,\\
    \theta_{\beta,2}=-2.5 rad,\\
    e_1 = 0.14 kcal/mol,\\
    e_2 = 0.4 kcal/mol`


The Pairwise potential:
+++++++++++++++++++++++

.. math::
        \Phi_{i,j}^{vdw}(r) = step(2^{1/6}\sigma_{ij}-r) \times \left( 4\epsilon\left[\left(\frac{\sigma_{ij}}{r}\right)^{12}- \left(\frac{\sigma_{ij}}{r}\right)^{6}\right]+(1-\mu\times\lambda_{ij}^{0}+\Delta)\times\epsilon\right)

        + \left[1-step(2^{1/6}\sigma_{ij}-r)\right]\times\left[(\mu \lambda_{ij}^{0}-\Delta)\times 4\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12}-\left(\frac{\sigma_{ij}}{r}\right)^6\right]\right]

Since the step function behaves like: :code:`step(x) = 0 if x < 0,and =1 otherwise`, we can separate in multiple cases for short likes following:

.. math::
        \Phi_{i,j}^{vdw}(r) =  4\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12}-\left(\frac{\sigma_{ij}}{r}\right)^{6}\right]+(1-\mu	\times\lambda_{ij}^{0}+\Delta)	\times\epsilon, r\le 2^{1/6}\sigma_{ij}

        \Phi_{i,j}^{vdw}(r) = (\mu\times\lambda_{ij}^{0}-\Delta) \times \left( 4\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12}-\left(\frac{\sigma_{ij}}{r}\right)^{6}\right]\right), r > 2^{1/6}\sigma_{ij}

where, :math:`\sigma_{i,j}=\frac{\sigma_i+\sigma_j}{2}`: is the vdW radius interaction of interacting beads

:math:`\lambda_{ij}^{0}=\frac{\lambda_i+\lambda_j}{2}`: hydropathy scale interaction of residues

:math:`\mu, \Delta`: are the only free parameters in the model. In Jeetain Mittal(2021) Protein Science, he simulated for 42 IDP proteins and fit Rg vs experimental values.

In the current implementation, hydropathy scales are taken from Urry model, :math:`(\mu, \Delta) = (1, 0.08)`

Nonbonded exclusion rule is :code:`1-2`, for hps_kr and hps_urry which we only exclude pair of atoms in bonded.
while it is :code:`1-4` for hps-ss, which we exclude 3 bonds.

The cut-off distance for Lennard-Jone potential: :math:`2.0 nm`

The Debye-Huckle potential has following form:
++++++++++++++++++++++++++++++++++++++++++++++
.. math::
        \Phi_{ij}^{el}(r) = \frac{q_{i}q_{j}}{4\pi\epsilon_0 D r}e^{-\kappa r}

where, :math:`q_i, q_j` are charge of residues :math:`i, j`

:math:`\epsilon_0`: Vacuum permitivity. For convenient, we precalculated the electric conversion factor
:math:`\frac{1}{4\pi\epsilon_0}= 138.935 485(9) kJ \times mol^{−1} \times nm \times e^{−2}`.

:math:`D`: dielectric constant, at 100mM mono-valence salt (NaCl), it takes values of 80.
The dielectric constant here is fixed, but it can be temperature dependent as the function:
:math:`\frac{5321}{T}+233.76-0.9297T+0.1417\times 10^{-2}\times T^2 - 0.8292\times 10^{-6}\times T^3`

:math:`\kappa`: inverse Debye length, at 100mM NaCl has values of :math:`1 nm^{-1}`

The cut-off distance for Electrostatics interactions: :math:`3.5 nm`

.. autoclass:: hps.core.models

        .. automethod:: __init__
        .. automethod:: hps.core.models.buildHPSModel