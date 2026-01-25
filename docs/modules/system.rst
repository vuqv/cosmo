System 
========================================================= 

A class containing methods and parameters for generating CG systems to be simulated using the OpenMM interface.
    It offers flexibility to create default and custom CG systems and to easily modify their parameters.

.. autoclass:: cosmo.core.system

        .. automethod:: __init__
        .. automethod:: cosmo.core.system.getAtoms
        .. automethod:: cosmo.core.system.getBonds
        .. automethod:: cosmo.core.system.setBondForceConstants
        .. automethod:: cosmo.core.system.setParticlesMasses
        .. automethod:: cosmo.core.system.setParticlesRadii
        .. automethod:: cosmo.core.system.setParticlesCharge
        .. automethod:: cosmo.core.system.setParticlescosmo
        .. automethod:: cosmo.core.system.addYukawaForces
        .. automethod:: cosmo.core.system.addAshbaughHatchForces
        .. automethod:: cosmo.core.system.createSystemObject
        .. automethod:: cosmo.core.system.checkBondDistances
        .. automethod:: cosmo.core.system.checkLargeForces
        .. automethod:: cosmo.core.system.addParticles
        .. automethod:: cosmo.core.system.addSystemForces
        .. automethod:: cosmo.core.system.dumpStructure
        .. automethod:: cosmo.core.system.dumpTopology
        .. automethod:: cosmo.core.system.dumpForceFieldData
        .. automethod:: cosmo.core.system.setMassPerResidueType
        .. automethod:: cosmo.core.system.setRadiusPerResidueType
        .. automethod:: cosmo.core.system.setChargePerResidueType
        .. automethod:: cosmo.core.system.setHPSPerResidueType
        .. automethod:: cosmo.core.system._setParameters