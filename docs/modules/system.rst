System
======

A class of methods and parameters for generating COSMO coarse-grained systems to be
simulated with OpenMM. It is typically constructed via
:func:`cosmo.models.buildCoarseGrainModel`, which sets bonds, (optional) angles and
torsions, Debye-Hückel electrostatics, and the short-range non-bonded force.

.. seealso::

   :doc:`models` for the force fields and the functional form of each potential,
   and :doc:`../usage/simulation_control` for how to run a simulation.

.. autoclass:: cosmo.core.system

   .. automethod:: __init__
   .. automethod:: getAtoms
   .. automethod:: getBonds
   .. automethod:: getAngles
   .. automethod:: getTorsions
   .. automethod:: setBondForceConstants
   .. automethod:: setParticlesMass
   .. automethod:: setParticlesRadii
   .. automethod:: setParticlesCharge
   .. automethod:: setParticlesHPS
   .. automethod:: setParticleTypeID
   .. automethod:: setMassPerResidueType
   .. automethod:: setRadiusPerResidueType
   .. automethod:: setChargePerResidueType
   .. automethod:: setHPSPerResidueType
   .. automethod:: setCAIDPerResidueType
   .. automethod:: addHarmonicBondForces
   .. automethod:: addGaussianAngleForces
   .. automethod:: addGaussianTorsionForces
   .. automethod:: addYukawaForces
   .. automethod:: addAshbaughHatchForces
   .. automethod:: addWangFrenkelForces
   .. automethod:: addParticles
   .. automethod:: addSystemForces
   .. automethod:: createSystemObject
   .. automethod:: checkBondDistances
   .. automethod:: checkLargeForces
   .. automethod:: coarseGrainingStructure
   .. automethod:: dumpStructure
   .. automethod:: dumpTopology
   .. automethod:: dumpForceFieldData
