#!/usr/bin/env python
# coding: utf-8

from .system import system


class models:
    """
    A class to hold functions for the automated generation of default hps models.

    Methods
    -------

    getCAModel(structure_file, kwarg**)
        Creates an alpha-carbon only hpsOpenMM system class object with default
        initialized parameters.

    """

    def getCAModel(structure_file,
                   create_system=True,
                   minimize=False,
                   residue_masses=True,
                   residue_radii=True,
                   residue_charge=True,
                   residue_hps=True,
                   hps_scale='kr',
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

        The method can be used to generate an initialized hpsOpenMM system class, that only
        contains the geometrical parameters, by passing the option default_parameters as False.
        This is useful to store the geometrical values of bonds, angles, dihedrals, etc. in
        order to add custom parameters and forces.

        The method can also be created without the initialisation of the forces classes by
        setting default_forces to False. This allows to load the default forcefield parameters
        and to modified them before creating the OpenMM force objects.

        The method can also be stopped before creating the OpenMM system class using create_system=False.

        Finally, a forcefield file can be given in order to read the forcefield parameters from it. You can give
        None to the contact file argument if you provide a forcefield file. Contacts definitions will be extracted from this
        file even if you give a path to a contact file.

        Parameters
        ----------
        structure_file : string [requires]
            Path to the input structure file.
        create_system : boolean (True)
            If True the function will call the createSystemObject() method
            to create an OpenMM system object. If modifications to the default
            forcefield are necessary this option should be given False.
        minimize : boolean (False)
            If True the initial structure will undergo the energy minimization.
        residue_masses : boolean (True)
            Set each alpha carbon atom mass to its average amino acid residue mass.
        residue_radii : boolean (True)
            Set each alpha carbon atom radius to its statistical amino acid residue radius.
        residue_charge : boolean (True)
            set charge for atoms in the system
        residue_hps : boolean (True)
            set HPS scale for atoms in the system
        hps_scale : string ('kr')
            HPS scale. There are two options correspond to two scale:
            'urry': using Urry scale
            'kr': using Kapcha-Rossy scale (default).
        box_dimension : float or array (None)
            If box_dimension is supplied, then will use PBC.
            if float is given, then use cubic box
            if an array of (3,1) is given, then use rectangular box with the given dimension
            if not specify: do not use PBC
        forcefield_file : string (None)
            Path to the input forcefield file.

        Returns
        -------
        hps : hpsOpenMM.system
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
        if forcefield_file is None:
            if residue_masses:
                print("Setting alpha-carbon masses to their average residue mass.")
                hps.setCAMassPerResidueType()

            if residue_radii:
                print("Setting alpha-carbon atoms radii to their statistical residue radius.")
                hps.setCARadiusPerResidueType()
            else:
                print('Setting default vdw radii (0.4 nm) for all atoms in the system')
                hps.setParticlesRadii(0.4)

            if residue_charge:
                print("Setting alpha-carbon charge to their residue charge.")
                hps.setCAChargePerResidueType()

            if residue_hps:
                print(f"Setting hydropathy scale to their residue, Using {hps_scale} scale.")
                hps.setCAHPSPerResidueType()

            hps.getBonds()
            print('Added ' + str(hps.n_bonds) + ' bonds')

        elif forcefield_file is not None:
            print(
                'Forcefield file given. Bonds, angles, torsions and native contacts definitions will be read from it!')

        # print('')
        print('Adding default bond force constant...')
        hps.setBondParameters(8368.0)
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
                hps.topology.setPeriodicBoxVectors(((box_dimension[0], 0, 0), (0, box_dimension[1], 0), (0, 0, box_dimension[2])))
            else:
                # cubic box, given parameter is single float
                hps.topology.setPeriodicBoxVectors(((box_dimension, 0, 0), (0, box_dimension, 0), (0, 0, box_dimension)))

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
