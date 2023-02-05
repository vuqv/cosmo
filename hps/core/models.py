#!/usr/bin/env python
# coding: utf-8

from typing import Any

from .system import system


class models:
    """
    A class to hold functions for the automated generation of default hps models.

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
        hps : :code:`hpsOpenMM.system`
            Initialized hpsOpenMM.system class with default options for defining
            a coarse-grained CA force field.
        """

        # common for all model:
        print(f'Generating CA-hps model for structure from file {structure_file}')
        print('')
        hps = system(structure_file, model)
        print("Checking input structure file ...")
        print("Be sure that you do not have missing residues in the initial structure. At the moment, I will not take "
              "care of that")

        # Set up geometric parameters of the model
        print('Setting up geometrical parameters ...')
        print('__________________________________________________________________')
        print('Keeping only alpha carbon atoms in topology')
        hps.getCAlphaOnly()

        print(f'There are {hps.n_chains} chain(s) in the input file.')

        # Common for all
        hps.getAtoms()
        print('Added ' + str(hps.n_atoms) + ' CA atoms')

        hps.getBonds()
        print('Added ' + str(hps.n_bonds) + ' bonds')

        print("Setting alpha-carbon masses to their average residue mass.")
        hps.setCAMassPerResidueType()

        print("Setting alpha-carbon charge to their residue charge.")
        hps.setCAChargePerResidueType()

        # difference for each model
        if model in ['hps_kr', 'hps_urry', 'hps_ss']:
            print("Setting alpha-carbon atoms radii to their statistical residue radius.")
            hps.setCARadiusPerResidueType()

            print(f"Setting hydropathy scale to their residue, Using {model} scale.")
            hps.setCAHPSPerResidueType()

        elif model in ['mpipi']:
            print(f"Setting atom type to their residue type, using {model} model.")
            hps.setCAIDPerResidueType()

        # add forces to system
        print('Adding default bond force constant:', end=' ')
        hps.setBondForceConstants()
        print('')
        print('__________________________________________________________________')

        print('Adding Forces:')
        hps.addHarmonicBondForces()
        print('Added Harmonic Bond Forces')
        print("---")

        if model == "hps_ss":
            # this model has angle bonded potential.
            # angle
            hps.getAngles()
            print(f'Added {hps.n_angles} angles ')
            hps.addGaussianAngleForces()
            print('Added Gaussian Angle Forces')
            print("---")

            # torsional
            hps.getTorsions()
            print(f'Added {hps.n_torsions} torsion angles ')
            hps.addGaussianTorsionForces()
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
                hps.topology.setPeriodicBoxVectors(
                    ((box_dimension[0], 0, 0), (0, box_dimension[1], 0), (0, 0, box_dimension[2])))
            else:
                # cubic box, given parameter is single float
                hps.topology.setPeriodicBoxVectors(
                    ((box_dimension, 0, 0), (0, box_dimension, 0), (0, 0, box_dimension)))

            unit_cell = hps.topology.getPeriodicBoxVectors()
            # use this to write coordinate in PBC box. requires 3 numbers, unzip to 3
            hps.system.setDefaultPeriodicBoxVectors(*unit_cell)

        else:
            use_pbc = False

        hps.addYukawaForces(use_pbc)
        print('Added Yukawa Force')
        print("---")

        if model in ['hps_kr', 'hps_urry', 'hps_ss']:
            hps.addAshbaughHatchForces(use_pbc)
            print('Added PairWise Force')
            print("---")
        elif model in ['mpipi']:
            hps.addWangFrenkelForces(use_pbc)
            print('Added Wang-Frenkel Force')
            print("---")
        print('')
        print('__________________________________________________________________')

        # Generate the system object and add previously generated forces

        print('Creating System Object:')
        # print('______________________')
        hps.createSystemObject(minimize=minimize, check_bond_distances=True)
        print('OpenMM system Object created')
        print('')

        return hps
