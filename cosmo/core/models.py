#!/usr/bin/env python
# coding: utf-8

from typing import Any

from .system import system

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
                              constraints: Any = None,
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
        constraints : str or None, optional (default: None)
            Bond treatment. ``None`` (or ``'None'``) keeps **flexible harmonic bonds**
            -- the physically appropriate default for intrinsically disordered chains,
            where backbone flexibility matters. ``'AllBonds'`` makes every bond a **rigid
            distance constraint** (pinned at its equilibrium length; no harmonic bond
            force is created), which removes the fast bond-stretch mode and lets the
            integrator take a larger timestep. Constraints act only on the CA/P
            pseudo-bonds; the non-bonded potentials (Ashbaugh-Hatch / Wang-Frenkel and,
            with a ribosome, the 12-10-6 excluded volume) are unaffected.
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
        cosmo_model.getBonds()

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

        # Resolve the bond-constraint mode. Accepted: None / 'None' / 'none' (flexible
        # harmonic bonds -- the default for IDPs) or 'AllBonds' (rigid distance
        # constraints). Rigid and flexible are mutually exclusive, so a bond is never
        # both constrained and harmonic. Unlike topo (Gō model, default 'AllBonds'),
        # cosmo defaults to flexible bonds -- a genuine IDP-physics choice, not a mirror
        # gap -- but supports 'AllBonds' for users who want the larger-timestep path.
        if constraints is None or str(constraints).strip().lower() == 'none':
            use_constraints = False
        elif str(constraints).strip().lower() == 'allbonds':
            use_constraints = True
        else:
            raise ValueError(
                f"Invalid constraints option: {constraints!r}. Expected 'AllBonds' or None.")
        cosmo_model.use_bond_constraints = use_constraints

        # Set particle interactions and add forces to the system. Force-group
        # insertion order below is significant (each force gets a group index in
        # this order) and must be kept in sync with the reporter / benchmark.
        cosmo_model.setBondForceConstants()
        # Only add the harmonic bond force when bonds are flexible. With rigid bonds
        # (constraints='AllBonds') the distance is pinned by a constraint added in
        # createSystemObject (after particles exist), so a harmonic term would be
        # redundant. Force-group order is preserved: the flexible path is unchanged, and
        # the rigid path simply omits the (absent) 'Harmonic Bond Energy' group.
        if not use_constraints:
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

        cosmo_model.addYukawaForces(use_pbc)

        # Short-range non-bonded potential: Ashbaugh-Hatch (LJ 12-6) for the
        # HPS-scale models, Wang-Frenkel for mpipi.
        if model in ('hps_kr', 'hps_urry', 'hps_ss'):
            cosmo_model.addAshbaughHatchForces(use_pbc)
            nb_name = 'Ashbaugh-Hatch (PairWise)'
        elif model in ('mpipi',):
            cosmo_model.addWangFrenkelForces(use_pbc)
            nb_name = 'Wang-Frenkel'

        extra = ', Gaussian angle+torsion' if model == 'hps_ss' else ''
        bond_term = 'rigid bond constraints (AllBonds)' if use_constraints else 'harmonic bond'
        print(f'  forces: {bond_term}, Yukawa, {nb_name}{extra}  '
              f'(pbc={use_pbc})')

        # Generate the system object and add the previously generated forces. The
        # bond-distance check always runs (it validates the built geometry); the
        # large-force / initial-energy check is skipped when check_forces is False
        # (e.g. restarting from a checkpoint, where the loaded state -- not the
        # input PDB geometry -- is what gets simulated).
        cosmo_model.createSystemObject(minimize=minimize, check_bond_distances=True,
                                       check_large_forces=check_forces)

        print('System build complete')
        print('')

        return cosmo_model
