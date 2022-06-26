#!/usr/bin/env python
# coding: utf-8

from .system import system


class models:
    """
    A class to hold functions for the automated generation of default hps models.

    Methods
    -------

    buildHPSModel(structure_file, kwarg**)
        Creates an alpha-carbon only :code:`hpsOpenMM` system class object with default
        initialized parameters.

    """

    def buildHPSModel(structure_file,
                      minimize=False,
                      hps_scale='hps_kr',
                      box_dimension=None,
                      forcefield_file=None):
        """
        Initialises a coarse-grained, carbon alpha (CA), hpsOpenMM system class
        from a structure and a contact file defining the native contacts for the
        coarse grained model.

        The system creation steps are:

        1) Add the geometrical parameters for the model.
        2) Add the default force field parameters for the model.
        3) Create the default force objects.
        4) Create the OpenMM system class.

        The method can be used to generate an initialised hpsOpenMM system class, that only
        contains the geometrical parameters, by passing the option default_parameters as False.

        Finally, a forcefield file can be given in order to read the forcefield parameters from it.

        Parameters
        ----------
        structure_file : string [requires]
            Path to the input structure file.
        minimize : boolean (False)
            If True the initial structure will undergo the energy minimization.
        hps_scale : string ('kr')
            HPS scale. There are two options correspond to two scale:
            'hps_urry': using Urry scale
            'hps_kr': using Kapcha-Rossy scale (default).
            In the future will add more scale like, Tesei scale, HPS-T which hps is temperature dependent.
        box_dimension : float or array (None)
            If box_dimension is supplied, then will use PBC.
            if float is given, then use cubic box
            if an array of (3,1) is given, then use rectangular box with the given dimension
            if not specify: do not use PBC
        forcefield_file : string (None)
            Path to the input forcefield file.

        Returns
        -------
        hps : :code:`hpsOpenMM.system`
            Initialised hpsOpenMM.system class with default options for defining
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
        if forcefield_file is None:

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

        elif forcefield_file is not None:
            print(
                'Forcefield file given. Bonds, angles, torsions and native contacts definitions will be read from it!')

        # print('')
        print('Adding default bond force constant...')
        # measured in unit of kj/mol/nm^2= 20kcal/mol/A^2
        hps.setBondForceConstants(8368.0)
        print('')
        print('_________________________________')

        print('Adding Forces:')
        hps.addHarmonicBondForces()
        print('Added Harmonic Bond Forces')

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
