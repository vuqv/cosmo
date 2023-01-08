System 
========================================================= 

A class containing methods and parameters for generating CG systems to be simulated using the OpenMM interface.
    It offers flexibility to create default and custom CG systems and to easily modify their parameters.

.. autoclass:: hps.core.system

        .. automethod:: __init__
        .. automethod:: hps.core.system.getAtoms
        .. automethod:: hps.core.system.getBonds
        .. automethod:: hps.core.system.setBondForceConstants
        .. automethod:: hps.core.system.setParticlesMasses
        .. automethod:: hps.core.system.setParticlesRadii
        .. automethod:: hps.core.system.setParticlesCharge
        .. automethod:: hps.core.system.setParticlesHPS
        .. automethod:: hps.core.system.addYukawaForces
        .. automethod:: hps.core.system.addAshbaughHatchForces
        .. automethod:: hps.core.system.createSystemObject
        .. automethod:: hps.core.system.checkBondDistances
        .. automethod:: hps.core.system.checkLargeForces
        .. automethod:: hps.core.system.addParticles
        .. automethod:: hps.core.system.addSystemForces
        .. automethod:: hps.core.system.dumpStructure
        .. automethod:: hps.core.system.dumpTopology
        .. automethod:: hps.core.system.dumpForceFieldData
        .. automethod:: hps.core.system.setCAMassPerResidueType
        .. automethod:: hps.core.system.setCARadiusPerResidueType
        .. automethod:: hps.core.system.setCAChargePerResidueType
        .. automethod:: hps.core.system.setCAHPSPerResidueType
        .. automethod:: hps.core.system._setParameters