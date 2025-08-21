#!/usr/bin/env python
# coding: utf-8

from typing import Any

from .system import system
from ..utils.build_nonbonded_interaction_revise import build_nonbonded_interaction

class models:
    """
    A class to hold functions for the automated generation of default cosmo models.

    Methods
    -------
    """

    @staticmethod
    def buildHPSModel(structure_file: str,
                      minimize: bool = False,
                      model: str = 'hps_urry',
                      box_dimension: Any = None):
        """
        This is a method for building a coarse-grained model for a protein system
        using the HPS (hydrophobic-polar scale) force field. The method takes as input a structure file,
        as well as optional parameters for whether to minimize the initial structure,
        the HPS scale to use (options include 'hps_urry', 'hps_ss', 'hps_kr', and 'mpipi'),
        and the dimensions of the periodic boundary conditions box.
        The method uses the hpsOpenMM system class to create an alpha-carbon only system,
        and then sets up the geometric parameters of the model (keeping only alpha carbon atoms in the topology,
        adding bonds and setting masses and charges of the alpha-carbon atoms based on residue type).
        Depending on the chosen HPS scale, it also sets the radii or atom types of the alpha-carbon atoms based on
        residue type. The method then adds bond, angle and torsional forces to the system, and uses periodic boundary
        conditions if the box_dimension parameter is provided.

        Creates an alpha-carbon only :code:`hpsOpenMM` system class object with default
        initialized parameters.

        Parameters
        ----------
        structure_file : string [required]
            Path to the input structure file.
        minimize : boolean (False)
            If True, the initial structure will undergo energy minimization.
        model : string [Optional, hps_urry]
            HPS scale. Available options are 'hps_urry', 'hps_ss', 'hps_kr', and 'mpipi'.
        box_dimension : float or array (None)
            If box_dimension is supplied, a PBC will be used.
            If a float is given, a cubic box will be used.
            If an array of (3,1) is given, a rectangular box with the given dimension will be used.
            If not specified, PBC will not be used.

        Returns
        -------
        cosmo : :code:`hpsOpenMM.system`
            Initialized hpsOpenMM.system class with default options for defining
            a coarse-grained CA force field.
        """

        # common for all model:
        print(f'Generating CA coarse-grained model for structure from file {structure_file}')
        print('')
        cosmo_model = system(structure_file, model)
        print("Checking input structure file ...")
        print("Be sure that you do not have missing residues in the initial structure. At the moment, I will not take "
              "care of that")

        # Set up geometric parameters of the model
        print('Setting up geometrical parameters ...')
        print('__________________________________________________________________')
        print('Keeping only alpha carbon atoms in topology')
        cosmo_model.getCAlphaOnly()

        print(f'There are {cosmo_model.n_chains} chain(s) in the input file.')

        # Common for all
        cosmo_model.getAtoms()
        print('Added ' + str(cosmo_model.n_atoms) + ' CA atoms')

        cosmo_model.getBonds()
        print('Added ' + str(cosmo_model.n_bonds) + ' bonds')

        print("Setting alpha-carbon masses to their average residue mass.")
        cosmo_model.setCAMassPerResidueType()

        print("Setting alpha-carbon charge to their residue charge.")
        cosmo_model.setCAChargePerResidueType()

        # difference for each model
        if model in ['hps_kr', 'hps_urry', 'hps_ss']:
            print("Setting alpha-carbon atoms radii to their statistical residue radius.")
            cosmo_model.setCARadiusPerResidueType()

            print(f"Setting hydropathy scale to their residue, Using {model} scale.")
            cosmo_model.setCAHPSPerResidueType()

        elif model in ['mpipi']:
            print(f"Setting atom type to their residue type, using {model} model.")
            cosmo_model.setCAIDPerResidueType()

        # add forces to system
        print('Adding default bond force constant:', end=' ')
        cosmo_model.setBondForceConstants()
        print('')
        print('__________________________________________________________________')

        print('Adding Forces:')
        cosmo_model.addHarmonicBondForces()
        print('Added Harmonic Bond Forces')
        print("---")

        if model == "hps_ss":
            # this model has angle bonded potential.
            # angle
            cosmo_model.getAngles()
            print(f'Added {cosmo_model.n_angles} angles ')
            cosmo_model.addGaussianAngleForces()
            print('Added Gaussian Angle Forces')
            print("---")

            # torsional
            cosmo_model.getTorsions()
            print(f'Added {cosmo_model.n_torsions} torsion angles ')
            cosmo_model.addGaussianTorsionForces()
            print('Add Gaussian Torsion Forces')
            print("---")

        if box_dimension:
            use_pbc = True
            if isinstance(box_dimension, list):
                """
                OpenMM use this to write dimension in PDB and dcd file. Require one-argument, so zip box dimension into 
                one variable.
                Rectangular box, given parameter is array of three number
                """
                cosmo_model.topology.setPeriodicBoxVectors(
                    ((box_dimension[0], 0, 0), (0, box_dimension[1], 0), (0, 0, box_dimension[2])))
            else:
                # cubic box, given parameter is single float
                cosmo_model.topology.setPeriodicBoxVectors(
                    ((box_dimension, 0, 0), (0, box_dimension, 0), (0, 0, box_dimension)))

            unit_cell = cosmo_model.topology.getPeriodicBoxVectors()
            # use this to write coordinate in PBC box. requires 3 numbers, unzip to 3
            cosmo_model.system.setDefaultPeriodicBoxVectors(*unit_cell)

        else:
            use_pbc = False

        cosmo_model.addYukawaForces(use_pbc)
        print('Added Yukawa Force')
        print("---")

        if model in ['hps_kr', 'hps_urry', 'hps_ss']:
            cosmo_model.addAshbaughHatchForces(use_pbc)
            print('Added PairWise Force')
            print("---")
        elif model in ['mpipi']:
            cosmo_model.addWangFrenkelForces(use_pbc)
            print('Added Wang-Frenkel Force')
            print("---")
        print('')
        print('__________________________________________________________________')

        # Generate the system object and add previously generated forces

        print('Creating System Object:')
        # print('______________________')
        cosmo_model.createSystemObject(minimize=minimize, check_bond_distances=True)
        print('cosmo system Object created')
        print('')

        return cosmo_model

    def buildEDModel(structure_file: str,
                     minimize: bool = False,
                     model: str = 'ed',
                     domain_def: str = None,
                     stride_output_file: str = None,
                     box_dimension: Any = None):
        """
        This is a method for building a coarse-grained model for a protein system
        using the ED force field. The method takes as input a structure file,
        as well as optional parameters for whether to minimize the initial structure,
        the ED scale to use (options include 'ed'),
        and the dimensions of the periodic boundary conditions box.
        """

        # common for all model:
        print(f'Generating CA coarse-grained model for structure from file {structure_file}')
        print('')
        cosmo_model = system(structure_file, model)
        print("Checking input structure file ...")
        print("Be sure that you do not have missing residues in the initial structure. At the moment, I will not take "
              "care of that")

        # Set up geometric parameters of the model
        print('Setting up geometrical parameters ...')
        print('__________________________________________________________________')
        print('Keeping only alpha carbon atoms in topology')
        cosmo_model.getCAlphaOnly()

        print(f'There are {cosmo_model.n_chains} chain(s) in the input file.')

        # set particle's properties
        # Common for all
        cosmo_model.getAtoms()
        print('Added ' + str(cosmo_model.n_atoms) + ' CA atoms')

        cosmo_model.getBonds()
        print('Added ' + str(cosmo_model.n_bonds) + ' bonds')

        print("Setting alpha-carbon masses to their average residue mass.")
        cosmo_model.setCAMassPerResidueType()

        print("Setting alpha-carbon charge to their residue charge.")
        cosmo_model.setCAChargePerResidueType()

        
        
        # set particle interactions
        # add forces to system
        print('Adding default bond force constant:', end=' ')
        cosmo_model.setBondForceConstants()
        print('')
        print('__________________________________________________________________')

        # all models have bonded interactions
        print('Adding Forces:')
        cosmo_model.addHarmonicBondForces()
        print('Added Harmonic Bond Forces')
        print("---")


        # this model has angle bonded potential.
        # angle
        cosmo_model.getAngles()
        print(f'Added {cosmo_model.n_angles} angles ')
        cosmo_model.addGaussianAngleForces()
        print('Added Gaussian Angle Forces')
        print("---")

        #add Periodic Torsion angle for Ed model
        cosmo_model.getTorsions()
        cosmo_model.addPeriodicTorsionForce()



        if box_dimension:
            use_pbc = True
            if isinstance(box_dimension, list):
                """
                OpenMM use this to write dimension in PDB and dcd file. Require one-argument, so zip box dimension into 
                one variable.
                Rectangular box, given parameter is array of three number
                """
                cosmo_model.topology.setPeriodicBoxVectors(
                    ((box_dimension[0], 0, 0), (0, box_dimension[1], 0), (0, 0, box_dimension[2])))
            else:
                # cubic box, given parameter is single float
                cosmo_model.topology.setPeriodicBoxVectors(
                    ((box_dimension, 0, 0), (0, box_dimension, 0), (0, 0, box_dimension)))

            unit_cell = cosmo_model.topology.getPeriodicBoxVectors()
            # use this to write coordinate in PBC box. requires 3 numbers, unzip to 3
            cosmo_model.system.setDefaultPeriodicBoxVectors(*unit_cell)

        else:
            use_pbc = False

        cosmo_model.addYukawaForces(use_pbc)
        print('Added Yukawa Force')
        print("---")

        # non-bonded interaction
        print("Building non-bonded interactions for ED model...")
        try:
            distance_matrix, energy_matrix = build_nonbonded_interaction(
                structure_file,
                domain_def,
                stride_output_file
            )
            print(f"Built non-bonded interaction matrices: {distance_matrix.shape}, {energy_matrix.shape}")
            
            # Store the matrices in the cosmo_model object for later use
            cosmo_model.distance_matrix = distance_matrix
            cosmo_model.energy_matrix = energy_matrix
            
            # You can also add them to the system if needed
            cosmo_model.addCustomNonBondedForce(distance_matrix, energy_matrix, use_pbc)
            
        except Exception as e:
            print(f"Warning: Could not build non-bonded interactions: {e}")
            print("Continuing with default non-bonded interactions...")
            distance_matrix = None
            energy_matrix = None




        print('')
        print('__________________________________________________________________')

        # Generate the system object and add previously generated forces

        print('Creating System Object:')
        # print('______________________')
        cosmo_model.createSystemObject(minimize=minimize, check_bond_distances=True)
        print('cosmo system Object created')
        print('')

        return cosmo_model
