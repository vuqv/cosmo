.. COSMO documentation master file.

Welcome to the COSMO documentation
==================================

**COSMO** (*COarse-grained Simulation of intrinsically disordered prOteins*) is a
Python library and command-line toolkit for coarse-grained molecular dynamics of
**intrinsically disordered proteins** and related biomolecules (RNA/DNA), built on
`OpenMM <https://openmm.org/>`_. From a sequence it builds a one-bead-per-residue
model under the hydropathy-scale **HPS** (Ashbaugh–Hatch) or **mpipi** (Wang–Frenkel)
force field and runs Langevin dynamics — for single-chain dimensions, liquid–liquid
phase separation, protein–RNA complexes, and **protein synthesis**.

**New here?** Read :doc:`overview` to see the two things COSMO does and jump to the
right tutorials, and :doc:`modules/introduction` for installation and package layout.

.. toctree::
   :maxdepth: 1
   :caption: Getting started

   overview
   modules/introduction

.. toctree::
   :maxdepth: 2
   :caption: Tutorials

   tutorials/index

.. toctree::
   :maxdepth: 1
   :caption: IDP simulation

   usage/simulation_control

.. toctree::
   :maxdepth: 1
   :caption: Protein synthesis

   usage/ribosome_preparation
   usage/codon_dwell_times
   usage/continuous_synthesis
   usage/cylinder_synthesis
   usage/synthesis_control

.. toctree::
   :maxdepth: 1
   :caption: Python & API reference

   modules/parameters
   modules/system
   modules/models
   cosmo.csp

.. toctree::
   :maxdepth: 2
   :caption: Full module index

   modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
