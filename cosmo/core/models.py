#!/usr/bin/env python
# coding: utf-8

from typing import Any

from .system import system
import openmm.unit as unit

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
                      box_dimension: Any = None,
                      PTC_atom_index: int = None):
        """
        This is a method for building a coarse-grained model for protein and
        nucleic-acid systems using the HPS (hydrophobic-polar scale) force field.
        The method takes as input a structure file, as well as optional parameters
        for whether to minimize the initial structure, the HPS scale to use
        (options include 'hps_urry', 'hps_ss', 'hps_kr', and 'mpipi'), and the
        dimensions of the periodic boundary conditions box. The method uses the
        COSMO system class to keep CA atoms for proteins and P atoms for RNA/DNA,
        then sets up the geometric parameters of the model (adding bonds and
        setting masses and charges of CA/P atoms based on residue type). Depending
        on the chosen HPS scale, it also sets the radii, hydropathy scale, or atom
        types of CA/P atoms based on residue type. The method then adds bond,
        angle, and torsional forces to the system, and uses periodic boundary
        conditions if the box_dimension parameter is provided.

        Creates a CA/P coarse-grained :code:`COSMO` system class object with
        default initialized parameters.

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
        cosmo : :code:`COSMO.system`
            Initialized COSMO.system class with default options for defining
            a coarse-grained CA/P force field.
        """

        # common for all model:
        print(f'Generating CA/P coarse-grained model for structure from file {structure_file}')
        print('')
        cosmo_model = system(structure_file, model)
        print("Checking input structure file ...")
        print("Be sure that you do not have missing residues in the initial structure. At the moment, I will not take "
              "care of that")

        # Set up geometric parameters of the model
        print('Setting up geometrical parameters ...')
        print('-' * 70)
        print('Keeping only alpha carbon and phosphate P atoms in topology')
        cosmo_model.coarseGrainingStructure()

        print(f'There are {cosmo_model.n_chains} chain(s) in the input file.')

        # Common for all
        cosmo_model.getAtoms()
        ca_atoms = sum(1 for atom in cosmo_model.atoms if atom.name == 'CA')
        p_atoms = sum(1 for atom in cosmo_model.atoms if atom.name == 'P')
        print(f'Added {cosmo_model.n_atoms} atoms ({ca_atoms} CA, {p_atoms} P)')

        # Get atom index of ribosome and nascent chain
        ribosome_atom_indices = [atom.index for atom in cosmo_model.atoms if atom.residue.chain.id != '8']
        nascent_atom_indices = [atom.index for atom in cosmo_model.atoms if atom.residue.chain.id == '8']
        untranslated_atom_indices = [i for i in nascent_atom_indices if i > PTC_atom_index]
        translated_atom_indices = [i for i in nascent_atom_indices if i <= PTC_atom_index]
        print(f"There are {len(untranslated_atom_indices)} untranslated atoms and {len(translated_atom_indices)} translated atoms")
        print(f"Untranslated atom indices: {untranslated_atom_indices}")
        print(f"Translated atom indices: {translated_atom_indices}")
        # Get bonds, optionally excluding bonds of specified chains. In the synthesis pipeline, we want to exclude bonds of the ribosome chain.
        nascent_chain_id = '8'
        except_chains = [chain.id for chain in cosmo_model.chains if chain.id != nascent_chain_id]
        cosmo_model.getBonds(except_chains=except_chains)
        print('Added ' + str(cosmo_model.n_bonds) + ' bonds')

        print("Setting CA/P masses to their average residue mass.")
        cosmo_model.setMassPerResidueType()


        print("Setting CA/P charge to their residue charge.")
        cosmo_model.setChargePerResidueType()

        # difference for each model
        if model in ['hps_kr', 'synthesis_kr', 'hps_urry', 'hps_ss']:
            print("Setting CA/P atom radii to their statistical residue radius.")
            cosmo_model.setRadiusPerResidueType()

            print(f"Setting CA/P hydropathy scale to their residue, Using {model} scale.")
            cosmo_model.setHPSPerResidueType()

        elif model in ['mpipi']:
            print(f"Setting CA/P atom type to their residue type, using {model} model.")
            cosmo_model.setCAIDPerResidueType()

        # add forces to system
        print('Adding default bond force constant:', end=' ')
        cosmo_model.setBondForceConstants()
        print('')
        print('-' * 70)

        print('Adding Forces:')
        cosmo_model.addHarmonicBondForces()
        print('Added Harmonic Bond Forces')
        print('-' * 30)

        if model == "hps_ss":
            # this model has angle bonded potential.
            # angle
            cosmo_model.getAngles()
            print(f'Added {cosmo_model.n_angles} angles ')
            cosmo_model.addGaussianAngleForces()
            print('Added Gaussian Angle Forces')
            print('-' * 30)

            # torsional
            cosmo_model.getTorsions()
            print(f'Added {cosmo_model.n_torsions} torsion angles ')
            cosmo_model.addGaussianTorsionForces()
            print('Add Gaussian Torsion Forces')
            print('-' * 30)

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
        #TODO: 
        """
        For Nonbonded forces, we need to exclude all nonbonded interactions within ribosome,
        And for nascent chain, exclude all nonbonded interactions that involves untranslated regions.
        """
        exclusion_list = []
        # exclude ribosome-ribosome interactions
        for i in ribosome_atom_indices:
            for j in ribosome_atom_indices:
                if i < j:
                    exclusion_list.append((i, j))
        # exclude ribosome-untranslated interactions
        for i in ribosome_atom_indices:
            for j in untranslated_atom_indices:
                if i < j:
                    exclusion_list.append((i, j))
        # exclude untranslated-untranslated interactions
        for i in untranslated_atom_indices:
            for j in untranslated_atom_indices:
                if i < j:
                    exclusion_list.append((i, j))
        # exclude translated-untranslated interactions
        for i in translated_atom_indices:
            for j in untranslated_atom_indices:
                if i < j:
                    exclusion_list.append((i, j))
        cosmo_model.addYukawaForces(use_pbc, exclusion_list)
        print('Added Yukawa Force')
        print('-' * 30)

        if model in ['hps_kr', 'synthesis_kr', 'hps_urry', 'hps_ss']:
            cosmo_model.addAshbaughHatchForces(use_pbc, exclusion_list)
            print('Added PairWise Force')
            print('-' * 30)
        elif model in ['mpipi']:
            cosmo_model.addWangFrenkelForces(use_pbc)
            print('Added Wang-Frenkel Force')
            print('-' * 30)
        print('')
        print('-' * 70)

        # Generate the system object and add previously generated forces

        print('Creating System Object:')
        # print('______________________')
        cosmo_model.createSystemObject(minimize=minimize, check_bond_distances=True)

        # all atoms of ribosome has mass 0 since they are frozen
        for atom_index in ribosome_atom_indices:
            # print(f"Setting mass of atom {atom_index} to 0.0")
            cosmo_model.system.setParticleMass(atom_index, 0.0 * unit.dalton)
        for atom_index in untranslated_atom_indices:
            cosmo_model.system.setParticleMass(atom_index, 0.0 * unit.dalton)

        print('cosmo system Object created')
        print('')

        return cosmo_model
