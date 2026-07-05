Models
=========================================================

The :code:`models` class builds predefined coarse-grained (CG) force fields by
initialising a :code:`system` object with the right parameters. Every model is a
**one-bead-per-residue** representation: proteins are mapped to a bead at each
alpha carbon (``CA``) and nucleic acids to a bead at each phosphate (``P``).

Supported models
++++++++++++++++++++++++++++++++++++++++

.. list-table::
   :header-rows: 1
   :widths: 16 26 30 28

   * - Model
     - Short-range pairwise term
     - Components
     - Notes
   * - ``hps_urry``
     - Ashbaugh–Hatch (Urry hydropathy)
     - protein, DNA
     - Default / recommended for IDPs
   * - ``hps_kr``
     - Ashbaugh–Hatch (Kapcha–Rossy)
     - protein, RNA, phospho-protein
     - Nucleic-acid + PTM parameters
   * - ``hps_ss``
     - Ashbaugh–Hatch + bonded (angle, torsion)
     - protein
     - Adds secondary-structure bonded terms
   * - ``mpipi``
     - Wang–Frenkel
     - **protein, RNA**
     - Near-quantitative LLPS; full RNA support

To build a model, call:

.. code-block:: python

    import cosmo
    model = cosmo.models.buildCoarseGrainModel(structure_file, model='hps_urry')

where ``structure_file`` is a PDB structure and ``model`` selects one of the
force fields above. The method keeps only the ``CA``/``P`` beads, assigns
per-residue masses, charges and force-field parameters, and adds the bonded and
non-bonded forces. The sections below give the functional form and parameters of
each term.

Which model uses which term
++++++++++++++++++++++++++++++++++++++++

The models do **not** all share the same potential. Only ``hps_ss`` adds the
angle and torsion (bonded secondary-structure) terms, and only ``mpipi`` uses the
Wang–Frenkel short-range term instead of Ashbaugh–Hatch. The electrostatics term
is common to all. Each potential section below applies only to the models marked
here.

.. list-table::
   :header-rows: 1
   :widths: 34 14 14 14 14

   * - Energy term
     - ``hps_urry``
     - ``hps_kr``
     - ``hps_ss``
     - ``mpipi``
   * - Harmonic bond
     - ✓
     - ✓
     - ✓
     - ✓
   * - Gaussian angle
     - –
     - –
     - ✓
     - –
   * - Gaussian torsion
     - –
     - –
     - ✓
     - –
   * - Ashbaugh–Hatch pairwise (vdW)
     - ✓
     - ✓
     - ✓
     - –
   * - Wang–Frenkel pairwise (vdW)
     - –
     - –
     - –
     - ✓
   * - Debye–Hückel electrostatics
     - ✓
     - ✓
     - ✓
     - ✓

The corresponding Hamiltonians are therefore:

.. math::
	H_{\mathrm{hps\_urry/kr}} = \sum_{bonds}V_{bond}+\sum_{i,j}\Phi_{ij}^{vdw,AH}+\sum_{i,j}\Phi_{i,j}^{el}

.. math::
	H_{\mathrm{hps\_ss}} = \sum_{bonds}V_{bond}+\sum_{angle}V_{angle}+\sum_{torsion}V_{torsion}+\sum_{i,j}\Phi_{ij}^{vdw,AH}+\sum_{i,j}\Phi_{i,j}^{el}

.. math::
	H_{\mathrm{mpipi}} = \sum_{bonds}V_{bond}+\sum_{i,j}\Phi_{ij}^{WF}+\sum_{i,j}\Phi_{i,j}^{el}

where :math:`\Phi^{vdw,AH}` is the Ashbaugh–Hatch term and :math:`\Phi^{WF}` is the
Wang–Frenkel term documented below.


The Bonded potential:
++++++++++++++++++++++

All models connect consecutive beads with a harmonic bond:

.. math::
        V_{bond} = \frac{k_b}{2}(r-r_0)^2

The spring constant and equilibrium length are model-dependent:

.. list-table::
   :header-rows: 1
   :widths: 22 30 24 24

   * - Model
     - :math:`k_b` (kJ mol\ :sup:`-1` nm\ :sup:`-2`)
     - :math:`r_0` protein (nm)
     - :math:`r_0` nucleic (nm)
   * - ``hps_urry`` / ``hps_kr`` / ``hps_ss``
     - 8368
     - 0.382
     - 0.5
   * - ``mpipi``
     - 8030  (= 8.03 J mol\ :sup:`-1` pm\ :sup:`-2`)
     - 0.381
     - 0.5

Nucleic-acid (``P``–``P``) bonds use the nucleic equilibrium length of 0.5 nm;
protein (``CA``–``CA``) bonds use the protein value.

Rigid vs. flexible bonds (``constraints``)
------------------------------------------

By default the chain bonds are the **flexible harmonic springs** above — the
physically appropriate choice for intrinsically disordered chains, where backbone
flexibility matters. Passing ``constraints='AllBonds'`` to
:func:`~cosmo.core.models.buildCoarseGrainModel` instead makes every ``CA``/``P``
bond a **rigid distance constraint** pinned at its equilibrium length: the harmonic
bond force is not created, the fast bond-stretch vibrational mode is removed, and the
integrator can take a larger timestep. The two are mutually exclusive — a bond is
never both constrained and harmonic.

.. code-block:: python

    # flexible harmonic bonds (default)
    model = cosmo.models.buildCoarseGrainModel(structure_file, model='hps_urry')
    # rigid AllBonds constraints (larger-timestep path)
    model = cosmo.models.buildCoarseGrainModel(structure_file, model='hps_urry',
                                               constraints='AllBonds')

Constraints act **only** on the pseudo-bonds; the non-bonded potentials
(Ashbaugh–Hatch / Wang–Frenkel and the Debye–Hückel electrostatics) are unaffected.
From ``md.ini`` the option is the ``constraints`` key (``None`` / ``AllBonds``), with
``constraint_tolerance`` (default ``1e-5``) setting the integrator's relative
constraint tolerance. Unlike the sibling ``topo`` package (a Gō model that defaults
to ``AllBonds``), cosmo defaults to flexible bonds — a deliberate IDP-physics choice.

Angle Potential
+++++++++++++++

**Applies to:** ``hps_ss`` only.

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

**Applies to:** ``hps_ss`` only.

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

.. math::
    k_{\alpha,1} = 11.4 \ \mathrm{kcal}/(\mathrm{mol} \times \mathrm{rad}^2) \\
    k_{\alpha,2} = 0.15 \ \mathrm{kcal}/(\mathrm{mol} \times \mathrm{rad}^4) \\
    \theta_{\alpha,1} = 0.9 \ \mathrm{rad} \\
    \theta_{\alpha,2} = 1.02 \ \mathrm{rad} \\
    e_0 = 0.27 \ \mathrm{kcal}/\mathrm{mol} \\
    k_{\beta,1} = 1.8 \ \mathrm{kcal}/(\mathrm{mol} \times \mathrm{rad}^2) \\
    k_{\beta,2} = 0.65 \ \mathrm{kcal}/(\mathrm{mol} \times \mathrm{rad}^4) \\
    \theta_{\beta,1} = -1.55 \ \mathrm{rad} \\
    \theta_{\beta,2} = -2.5 \ \mathrm{rad} \\
    e_1 = 0.14 \ \mathrm{kcal}/\mathrm{mol} \\
    e_2 = 0.4 \ \mathrm{kcal}/\mathrm{mol}


The Pairwise potential (Ashbaugh–Hatch):
++++++++++++++++++++++++++++++++++++++++

**Applies to:** ``hps_urry``, ``hps_kr``, ``hps_ss``. The ``mpipi`` model uses the
Wang–Frenkel term instead (see *The Mpipi model* below).

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

The Mpipi model (Wang–Frenkel)
++++++++++++++++++++++++++++++++++++++++

The ``mpipi`` model (Joseph *et al.*, *Nat. Comput. Sci.* **1**, 732–743, 2021)
replaces the Ashbaugh–Hatch term with the **Wang–Frenkel** potential and
represents each amino acid or nucleotide by a single bead with a mass, charge,
and a tabulated set of pairwise parameters. Its total energy is

.. math::
        E_{\mathrm{Mpipi}} = E_{bond} + E_{elec} + E_{pair}

where :math:`E_{bond}` is the harmonic bond above and :math:`E_{elec}` is the
Debye–Hückel term below. The short-range pairwise energy between beads of types
:math:`i` and :math:`j` at separation :math:`r` is

.. math::
        \Phi_{ij}^{WF}(r) = \varepsilon_{ij}\,\alpha_{ij}
        \left[\left(\frac{\sigma_{ij}}{r}\right)^{2\mu_{ij}} - 1\right]
        \left[\left(\frac{R_{ij}}{r}\right)^{2\mu_{ij}} - 1\right]^{2\nu_{ij}}

.. math::
        \alpha_{ij} = 2\nu_{ij}\left(\frac{R_{ij}}{\sigma_{ij}}\right)^{2\mu_{ij}}
        \left[\frac{2\nu_{ij}+1}{2\nu_{ij}\left(\left(\frac{R_{ij}}{\sigma_{ij}}\right)^{2\mu_{ij}}-1\right)}\right]^{2\nu_{ij}+1}

The Wang–Frenkel potential is **finite-ranged**: it vanishes smoothly
(quadratically) at :math:`r = R_{ij}`, so no truncation/shifting is needed. The
parameters are:

* :math:`\varepsilon_{ij}, \sigma_{ij}, \mu_{ij}` — tabulated **per interacting
  pair** (not from a mixing rule), from Supplementary Tables 11 (protein–protein)
  and 12 (interactions with RNA).
* :math:`\nu_{ij} = 1` for every pair.
* :math:`R_{ij} = 3\,\sigma_{ij}` — the per-pair interaction range.
* :math:`\mu_{ij} = 2` for protein–protein pairs, **except**
  :math:`\mu_{\mathrm{V\!-\!I}} = 4` and :math:`\mu_{\mathrm{I\!-\!I}} = 11`;
  :math:`\mu_{ij} = 3` for every pair that involves an RNA bead.

Because the interaction range is per-pair, the neighbour-list cutoff is set
automatically to :math:`\max_{ij} R_{ij}` (the largest :math:`3\sigma_{ij}` in the
parameter table; ≈ 2.55 nm once RNA beads are present), so no pair is silently
truncated.

Charges follow the Mpipi convention: charged amino acids carry
:math:`q = \pm 0.75\,e`, histidine :math:`q = +0.375\,e`, and each RNA bead
:math:`q = -0.75\,e`.

RNA / nucleic-acid support
++++++++++++++++++++++++++++++++++++++++

The ``mpipi`` model supports RNA out of the box. RNA is represented by one bead
per nucleotide placed at the phosphate (``P``) atom, with the four bases
``A``, ``C``, ``G``, ``U`` carrying standard nucleotide masses and charge
:math:`-0.75\,e`:

.. list-table::
   :header-rows: 1
   :widths: 14 22 18 18

   * - Base
     - mass (g mol\ :sup:`-1`)
     - charge (:math:`e`)
     - bead id
   * - A
     - 329.2
     - −0.75
     - 20
   * - C
     - 305.2
     - −0.75
     - 21
   * - G
     - 345.2
     - −0.75
     - 22
   * - U
     - 306.2
     - −0.75
     - 23

Common residue-name aliases are accepted on input and resolve to the same beads:
``RA``/``RC``/``RG``/``RU`` and ``ADE``/``CYT``/``GUA``/``URA``. Protein–RNA and
RNA–RNA Wang–Frenkel parameters (with :math:`\mu = 3`) come from Supplementary
Table 12; the RNA bond length is the nucleic value of 0.5 nm. A worked
protein + RNA example is given in the *protein–RNA complex* tutorial.

.. note::

   The original Mpipi paper describes the RNA parameters as an *initial set*;
   the authors flag the RNA bond and angular constants in particular for future
   refinement. COSMO uses the published values.

The Debye-Huckle potential has following form:
++++++++++++++++++++++++++++++++++++++++++++++

Electrostatics is shared by all models (it acts on the per-residue charges):

.. math::
        \Phi_{ij}^{el}(r) = \frac{q_{i}q_{j}}{4\pi\epsilon_0 D r}e^{-\kappa r}

where, :math:`q_i, q_j` are charge of residues :math:`i, j`

:math:`\epsilon_0`: Vacuum permitivity. For convenient, we precalculated the electric conversion factor
:math:`\frac{1}{4\pi\epsilon_0}= 138.935 485(9) kJ \times mol^{−1} \times nm \times e^{−2}`.

:math:`D`: dielectric constant, at 100mM mono-valence salt (NaCl), it takes values of 80.
The dielectric constant here is fixed, but it can be temperature dependent as the function:
:math:`\frac{5321}{T}+233.76-0.9297T+0.1417\times 10^{-2}\times T^2 - 0.8292\times 10^{-6}\times T^3`

:math:`\kappa`: inverse Debye length. The screening length :math:`\kappa^{-1}` is
**model-dependent**:

* HPS family (``hps_urry`` / ``hps_kr`` / ``hps_ss``):
  :math:`\kappa = 1\ \mathrm{nm}^{-1}` (:math:`\kappa^{-1} = 1.0` nm, ~100 mM NaCl).
* ``mpipi``: :math:`\kappa = 1.26\ \mathrm{nm}^{-1}`
  (:math:`\kappa^{-1} = 0.795` nm, 150 mM), per Joseph *et al.* 2021.

Models not listed fall back to 1.0 nm. The value is stored in
``cosmo.parameters.model_parameters.debye_length``.

The cut-off distance for Electrostatics interactions: :math:`3.5 nm`

.. autoclass:: cosmo.core.models

        .. automethod:: __init__
        .. automethod:: cosmo.core.models.buildCoarseGrainModel