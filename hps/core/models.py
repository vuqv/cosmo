#!/usr/bin/env python
# coding: utf-8

from simtk.openmm import *
from simtk.openmm.app import *

from .system import system


class models:
    """
    A class to hold functions for the automated generation of default SBM models.

    Methods
    -------

    getCAModel(structure_file, kwarg**)
        Creates an alpha-carbon only sbmOpenMM system class object with default
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
                   forcefield_file=None):
        """
        Initialises a coarse-grained, carbon alpha (CA), sbmOpenMM system class
        from a structure and a contact file defining the native contacts for the
        coarse grained model.

        The system creation steps are:

        1) Add the geometrical parameters for the model.
        2) Add the default force field parameters for the model.
        3) Create the default force objects.
        4) Create the OpenMM system class.

        The method can be used to generate an initialized sbmOpenMM system class, that only
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
        structure_file : string
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
        forcefield_file : string
            Path to the input forcefield file.

        Returns
        -------
        sbm : sbmOpenMM.system
            Initialized sbmOpenMM.system class with default options for defining
            a coarse-grained CA force field.
        """

        print('Generating CA SBM for structure file ' + structure_file)
        print('')
        # Set up geometric parameters of the model
        print('Setting up geometrical parameters:')
        print('_________________________________')
        sbm = system(structure_file)
        print('Keeping only alpha carbon atoms in topology')
        sbm.getCAlphaOnly()
        sbm.getAtoms()
        print('Added ' + str(sbm.n_atoms) + ' CA atoms')
        if forcefield_file is None:
            if residue_masses:
                print("Setting alpha-carbon masses to their average residue mass.")
                sbm.setCAMassPerResidueType()

            if residue_radii:
                print("Setting alpha-carbon atoms radii to their statistical residue radius.")
                sbm.setCARadiusPerResidueType()
            else:
                print('Setting default vdw radii (0.4 nm) for all atoms in the system')
                sbm.setParticlesRadii(0.4)

            if residue_charge:
                print("Setting alpha-carbon charge to their residue charge.")
                sbm.setCAChargePerResidueType(hps_scale)

            if residue_hps:
                print("Setting hydropathy scale to their residue.")
                if hps_scale == 'urry':
                    print("Using Urry scale.")
                    sbm.setCAHPSUrryPerResidueType()
                if hps_scale == 'kr':
                    print("Using Kapcha-Rossy scale.")
                    sbm.setCAHPSKRPerResidueType()

            sbm.getBonds()
            print('Added ' + str(sbm.n_bonds) + ' bonds')

        elif forcefield_file is not None:
            print(
                'Forcefield file given. Bonds, angles, torsions and native contacts definitions will be read from it!')

        # print('')
        print('Adding default bond force constant...')
        sbm.setBondParameters(8368.0)
        print('')
        print('_________________________________')

        print('Adding Forces:')
        sbm.addHarmonicBondForces()
        print('Added Harmonic Bond Forces')

        sbm.addYukawaForces()
        print('Added Yukawa Force')

        sbm.addAshbaughHatchForces()
        print('Added PairWise Force')
        print('')
        print('_________________________________')

        # Generate the system object and add previously generated forces

        print('Creating System Object:')
        # print('______________________')
        sbm.createSystemObject(minimize=minimize, check_bond_distances=True)
        print('OpenMM system Object created')
        print('')

        return sbm
