What COSMO does
===============

COSMO builds a **one-bead-per-residue, sequence-based coarse-grained model** of
intrinsically disordered proteins (and RNA) and runs Langevin dynamics on
`OpenMM <https://openmm.org/>`_. Interactions come from the **sequence** — the
hydropathy-scale HPS (Ashbaugh–Hatch) or the mpipi (Wang–Frenkel) force field plus
Debye–Hückel electrostatics — not from a folded structure, so COSMO is built for
**disordered** chains: single-chain dimensions, phase separation, and complexes.

That model powers **two complementary workflows**:

* :ref:`A. Coarse-grained simulation of IDPs <overview-simulation>` — single chains,
  multi-chain slabs / LLPS, and protein–RNA complexes.
* :ref:`B. Co-translational synthesis <overview-synthesis>` — grow the chain residue
  by residue on the ribosome under codon-resolved kinetics, and watch how a disordered
  chain extrudes and behaves *as it is made*.

Part B builds on the Part A force field, so start with A if you are new here.


.. _overview-simulation:

A. Coarse-grained simulation of IDPs
------------------------------------

Start from a sequence (a CA/P PDB) and run HPS / mpipi dynamics: single-chain radius
of gyration, temperature/pressure coupling and PBC, slab simulations of liquid–liquid
phase separation, and protein–RNA mixtures.

.. rubric:: Tutorials

* :doc:`1 · Single-chain quickstart <tutorials/01_single_chain>`
* :doc:`2 · Models & force fields <tutorials/02_models>`
* :doc:`3 · PBC, temperature & pressure <tutorials/03_pbc_temperature_pressure>`
* :doc:`4 · Restart & outputs <tutorials/04_restart_and_outputs>`
* :doc:`5 · Slab / LLPS <tutorials/05_slab_llps>`
* :doc:`6 · Protein–RNA complex <tutorials/06_protein_rna_complex>`

.. rubric:: Reference

* :doc:`Simulation control options (md.ini) <usage/simulation_control>`


.. _overview-synthesis:

B. Co-translational synthesis
-----------------------------

Grow the nascent chain **N→C, one residue at a time**, timing every residue from its
mRNA codon (the O'Brien Continuous Synthesis Protocol, ported to cosmo's IDP force
field). Two runners differ in how the ribosome exit tunnel is represented:

* **``cosmo-cylinder``** — the exit tunnel is an analytic cylindrical bore through an
  infinite wall (no explicit ribosome beads); fast, never jams, one MD segment per
  residue.
* **``cosmo-csp``** — the ribosome-based counterpart: grow the chain through an
  explicit truncated CG ribosome PDB (see :doc:`usage/ribosome_preparation`) with the
  O'Brien 12-10-6 excluded volume, in three codon-timed sub-stages per residue, then
  eject the completed chain.

**Which to use.** Reach for the **cylinder** for fast exploration of how tunnel geometry
+ codon kinetics shape a disordered chain as it extrudes, or when you have no ribosome
structure — it is the simplest starting point. Reach for the **explicit ribosome** when
the tunnel-wall charge, the real tunnel shape (constriction / vestibule), or
translocation-coupled forces matter to your question.

.. warning::

   **The two runners are comparable only in the *mean* per-residue dwell time, not in
   confinement chemistry.** The cylinder omits the ribosome's electrostatics and surface
   excluded volume and has a uniform straight bore. Do not compare extrusion/conformation
   observables (radius of gyration vs. length, contact formation) across the two runners
   without accounting for those missing terms.

Runnable proof-of-concept examples (α-synuclein) live in ``sandbox/validate/``
(``csp.ini`` and ``cylinder.ini``).

.. rubric:: Tutorials

* :doc:`7 · Continuous synthesis — analytic tunnel (cosmo-cylinder) <tutorials/07_csp_cylinder>`
* :doc:`8 · Continuous synthesis — coarse-grained E. coli ribosome (cosmo-csp) <tutorials/08_csp_cg_ribosome>`

.. rubric:: Reference

* :doc:`The ribosome structure — get or build one <usage/ribosome_preparation>`
* :doc:`Synthesis on a coarse-grained ribosome <usage/continuous_synthesis>`
* :doc:`Synthesis through an analytic tunnel <usage/cylinder_synthesis>`
* :doc:`Synthesis control options (csp.ini) <usage/synthesis_control>`
* :doc:`cosmo.csp API reference <cosmo.csp>`

.. note::

   Co-translational synthesis is the cosmo port of the sibling ``topo`` package's
   ``topo.csp``. It shares the codon kinetics, the three-stage protocol, the cylinder
   model, and the O'Brien 12-10-6 ribosome excluded volume, but the nascent chain is a
   **sequence-based IDP** (HPS / mpipi) rather than a structure-based Gō model — so
   there is no STRIDE, native-contact map, or ``domain.yaml``.
