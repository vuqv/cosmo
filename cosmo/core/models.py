#!/usr/bin/env python
# coding: utf-8

from typing import Any

from .system import system
import openmm.unit as unit

class models:
    """
    A class to hold functions for the automated generation of default COSMO models.

    Methods
    -------
    """

    @staticmethod
    def buildCoarseGrainModel(structure_file: str,
                              minimize: bool = False,
                              model: str = 'hps_urry',
                              box_dimension: Any = None,
                              frozen_indices: list = None,
                              except_chains: list = None,
                              nb_exclusions: list = None,
                              check_forces: bool = True):
        """
        Build a CA/P coarse-grained model for protein and nucleic-acid systems
        using an HPS (hydropathy-scale) force field.

        Creates an alpha-carbon (protein) / phosphate (nucleic-acid) system with
        harmonic bonds, Yukawa electrostatics, and a short-range non-bonded
        potential whose form depends on the chosen model: Ashbaugh-Hatch (LJ 12-6)
        for the HPS-scale models, Wang-Frenkel for ``mpipi``. The ``hps_ss`` model
        additionally adds Gaussian angle and torsion bonded terms.

        Parameters
        ----------
        structure_file : str
            Path to the input structure file.
        minimize : bool, optional (default: False)
            If True, run energy minimization on the initial structure.
        model : str, optional (default: 'hps_urry')
            HPS scale. Available options: 'hps_urry', 'hps_ss', 'hps_kr', 'mpipi'.
        box_dimension : float or array, optional
            If set, use PBC (cubic if float, rectangular if [x, y, z]). If not
            specified, PBC is not used.
        frozen_indices : list, optional
            List of atom indices to freeze (mass set to zero).
        except_chains : list, optional
            List of chain IDs to exclude from bonding.
        nb_exclusions : list, optional
            List of pairs of atom indices to exclude from nonbonded forces.
        check_forces : bool, optional (default: True)
            If True, run the build-time large-force / initial-energy check on the
            input structure. Set False when restarting from a checkpoint, where
            the input-structure energy is irrelevant (the loaded state, not the
            PDB geometry, is what gets simulated).

        Returns
        -------
        cosmo_model : cosmo.core.system
            Initialized CA/P coarse-grained system ready for simulation.
        """

        print('')
        print('=' * 66)
        print('[ System build ]')
        print('=' * 66)
        print(f'Building CA/P coarse-grained model (model={model}) from {structure_file}')

        cosmo_model = system(structure_file, model)

        # Build CA/P topology: keep only alpha-carbon (protein) and phosphate
        # (nucleic-acid) beads, then collect atoms and bonds.
        cosmo_model.coarseGrainingStructure()
        cosmo_model.getAtoms()
        cosmo_model.getBonds(except_chains=except_chains)

        ca_atoms = sum(1 for atom in cosmo_model.atoms if atom.name == 'CA')
        p_atoms = sum(1 for atom in cosmo_model.atoms if atom.name == 'P')
        print(f'  chains={cosmo_model.n_chains}  atoms={cosmo_model.n_atoms} '
              f'({ca_atoms} CA, {p_atoms} P)  bonds={cosmo_model.n_bonds}')

        # Per-residue particle properties (mass, charge, and -- depending on the
        # model family -- excluded-volume radius + hydropathy, or the mpipi atom
        # type).
        cosmo_model.setMassPerResidueType()
        cosmo_model.setChargePerResidueType()
        if model in ('hps_kr', 'hps_urry', 'hps_ss'):
            cosmo_model.setRadiusPerResidueType()
            cosmo_model.setHPSPerResidueType()
            print(f'  per-residue: mass, charge, radius, hydropathy ({model} scale)')
        elif model in ('mpipi',):
            cosmo_model.setCAIDPerResidueType()
            print(f'  per-residue: mass, charge, atom type ({model})')

        # Set particle interactions and add forces to the system. Force-group
        # insertion order below is significant (each force gets a group index in
        # this order) and must be kept in sync with the reporter / benchmark.
        cosmo_model.setBondForceConstants()
        cosmo_model.addHarmonicBondForces()

        if model == 'hps_ss':
            # hps_ss adds bonded angle + torsion terms on top of hps_urry.
            cosmo_model.getAngles()
            cosmo_model.addGaussianAngleForces()
            cosmo_model.getTorsions()
            cosmo_model.addGaussianTorsionForces()
            print(f'  bonded extras: {cosmo_model.n_angles} angles, '
                  f'{cosmo_model.n_torsions} torsions (Gaussian)')

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

        cosmo_model.addYukawaForces(use_pbc, nb_exclusions=nb_exclusions)

        # Short-range non-bonded potential: Ashbaugh-Hatch (LJ 12-6) for the
        # HPS-scale models, Wang-Frenkel for mpipi.
        if model in ('hps_kr', 'hps_urry', 'hps_ss'):
            cosmo_model.addAshbaughHatchForces(use_pbc, nb_exclusions=nb_exclusions)
            nb_name = 'Ashbaugh-Hatch (PairWise)'
        elif model in ('mpipi',):
            cosmo_model.addWangFrenkelForces(use_pbc, nb_exclusions=nb_exclusions)
            nb_name = 'Wang-Frenkel'

        extra = ', Gaussian angle+torsion' if model == 'hps_ss' else ''
        print(f'  forces: harmonic bond, Yukawa, {nb_name}{extra}  '
              f'(pbc={use_pbc})')

        # Generate the system object and add the previously generated forces. The
        # bond-distance check always runs (it validates the built geometry); the
        # large-force / initial-energy check is skipped when check_forces is False
        # (e.g. restarting from a checkpoint, where the loaded state -- not the
        # input PDB geometry -- is what gets simulated).
        cosmo_model.createSystemObject(minimize=minimize, check_bond_distances=True,
                                       check_large_forces=check_forces)

        if frozen_indices:
            print(f'  freezing {len(frozen_indices)} atoms from moving')
            for idx in frozen_indices:
                cosmo_model.system.setParticleMass(idx, 0.0 * unit.dalton)

        print('System build complete')
        print('')

        return cosmo_model

    # Deprecated alias. Retained only so the frozen ``examples/growing/`` scripts
    # keep working until they are rewritten (see TODO). New code should call
    # ``buildCoarseGrainModel`` (the name used by the sibling ``topo`` project).
    @staticmethod
    def buildHPSModel(*args, **kwargs):
        """Deprecated alias for :meth:`buildCoarseGrainModel`."""
        return models.buildCoarseGrainModel(*args, **kwargs)
