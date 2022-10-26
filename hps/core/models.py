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
                      hps_scale: str = 'hps_urry',
                      box_dimension: Any = None):
        """
        Creates an alpha-carbon only :code:`hpsOpenMM` system class object with default
        initialized parameters.

        Initializes a coarse-grained, carbon alpha (CA), hpsOpenMM system class
        from a structure and a contact file defining the native contacts for the
        coarse grained model.

        The system creation steps are:

        1) Add the geometrical parameters for the model.
        2) Add the default force field parameters for the model.
        3) Create the default force objects.
        4) Create the OpenMM system class.

        The method can be used to generate an initialized hpsOpenMM system class, that only
        contains the geometrical parameters, by passing the option default_parameters as False.


        Parameters
        ----------
        structure_file : string [requires]
            Path to the input structure file.
        minimize : boolean (False)
            If True the initial structure will undergo the energy minimization.
        hps_scale : string ('hps_urry')
            HPS scale. There are three options correspond to two scale:
                * 'hps_urry': using Urry scale (default).
                * 'hps_ss': hps_urry with angle and torsion potential.
                * 'hps_kr': using Kapcha-Rossy scale.
        box_dimension : float or array (None)
            If box_dimension is supplied, then will use PBC.
            if float is given, then use cubic box
            if an array of (3,1) is given, then use rectangular box with the given dimension
            if not specify: do not use PBC

        Returns
        -------
        hps : :code:`hpsOpenMM.system`
            Initialized hpsOpenMM.system class with default options for defining
            a coarse-grained CA force field.
        """

        print('Generating CA hps for structure file ' + structure_file)
        print('')
        # Set up geometric parameters of the model
        print('Setting up geometrical parameters:')
        print('_________________________________')
        hps = system(structure_file, hps_scale)
        print('Keeping only alpha carbon atoms in topology')
        hps.getCAlphaOnly()
        hps.getAtoms()
        print('Added ' + str(hps.n_atoms) + ' CA atoms')

        print("Setting alpha-carbon masses to their average residue mass.")
        hps.setCAMassPerResidueType()

        print("Setting alpha-carbon atoms radii to their statistical residue radius.")
        hps.setCARadiusPerResidueType()

        print("Setting alpha-carbon charge to their residue charge.")
        hps.setCAChargePerResidueType()

        print(f"Setting hydropathy scale to their residue, Using {hps_scale} scale.")
        hps.setCAHPSPerResidueType()

        hps.getBonds()
        print('Added ' + str(hps.n_bonds) + ' bonds')

        if hps_scale == "hps_ss":
            hps.getAngles()
            print(f'Added {hps.n_angles} angles ')

            hps.getTorsions()
            print(f'Added {hps.n_torsions} torsion angles ')

        print('Adding default bond force constant...')
        # measured in unit of kj/mol/nm^2 (k_bond is set to 20kcal/mol/A^2)
        hps.setBondForceConstants(8368.0)
        print('')
        print('_________________________________')

        print('Adding Forces:')
        hps.addHarmonicBondForces()
        print('Added Harmonic Bond Forces')

        if hps_scale == "hps_ss":
            hps.addGaussianAngleForces()
            print('Added Gaussian Angle Forces')

            hps.addGaussianTorsionForces()
            print('Add Gaussian Torsion Forces')

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

        hps.addAshbaughHatchForces(use_pbc)
        print('Added PairWise Force')
        print('')
        print('_________________________________')

        # Generate the system object and add previously generated forces

        print('Creating System Object:')
        # print('______________________')
        hps.createSystemObject(minimize=minimize, check_bond_distances=True)
        print('OpenMM system Object created')
        print('')

        return hps
