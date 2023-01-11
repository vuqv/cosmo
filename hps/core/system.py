#!/usr/bin/env python
# coding: utf-8

from collections import OrderedDict

import numpy as np
import openmm
import parmed as pmd

from ..parameters import model_parameters


# from openmm import *
# from openmm.app import *


class system:
    """
    A class containing methods and parameters for generating CG systems to be simulated using the OpenMM interface.
    It offers flexibility to create default and custom CG systems and to easily modify their parameters.

    Parameters
    ----------
        structure_path : string [requires]
            Name of the input PDB or CIF file
        model: 'hps_kr','hps_urry' [optional, default='hps_urry']
            Hydropathy scale. Currently, there are two models are supported.

    Attributes
    ----------

    structure : :code:`openmm.app.pdbfile.PDBFile or openmm.app.pdbxfile.PDBxFile`
        Object that holds the information of OpenMM PDB or CIF parsing methods.
    topology : :code:`openmm.app.topology.Topology`
        OpenMM topology of the model.
    positions : :code:`unit.quantity.Quantity`
        Atomic positions of the model.
    particles_mass : :code:`float or list`
        Mass of each particle. If float then uniform masses are given to all
        particles. If list per-particle masses are assigned.
    particles_charge : :code:`list`
        Charge of each particle.
    rf_sigma : :code:`float`
        Sigma parameter used in the pairwise force object.
        This is vdw Radius of beads
    atoms : :code:`list`
        A list of the current atoms in the model. The items are :code:`openmm.app.topology.atoms`
        initialised classes.
    n_atoms : :code:`int`
        Total numer of atoms in the model.
    bonds : :code:`collections.OrderedDict`
        A dict that uses bonds (2-tuple of :code:`openmm.app.topology.bonds` objects)
        present in the model as keys and their forcefield properties as values.
    bonds_indexes : :code:`list`
        A list containing the zero-based indexes of the atoms defining the bonds in the model.
    n_bonds : :code:`int`
        Total number of bonds in the model.
    bonded_exclusions_index : :code:`int`
        Exclusion rule for nonbonded force. =1 for hps_kr and hps_urry, =3 for hps_ss
    harmonicBondForce : :code:`openmm.HarmonicBondForce`
        Stores the OpenMM :code:`HarmonicBondForce` initialised-class. Implements
        a harmonic bond potential between pairs of particles, that depends
        quadratically on their distance.
    n_angles : :code:`int`
        Total number of angles in the model.
    gaussianAngleForce : :code:`openmm.CustomAngleForce`
        Stores the OpenMM :code:`CustomAngleForce` initialised-class. Implements
        a Gaussian angle bond potential between pairs of three particles.
    n_torsions : :code:`int`
        Total number of torsion angles in the model.
    gaussianTorsionForce : :code:`openmm.CustomTorsionForce`
        Stores the OpenMM :code:`CustomTorsionForce` initialised-class. Implements
        a Gaussian torsion angle bond potential between pairs of four particles.
    yukawaForce : :code:`openmm.CustomNonbondedForce`
        Stores the OpenMM :code:`CustomNonbondedForce` initialized-class.
        Implements the Debye-Huckle potential.
    ashbaugh_HatchForce : :code:`openmm.CustomNonbondedForce`
        Stores the OpenMM :code:`CustomNonbondedForce` initialized-class. Implements the pairwise short-range
        potential.
    forceGroups : :code:`collections.OrderedDict`
        A dict that uses force names as keys and their corresponding force
        as values.
    system : :code:`openmm.System`
        Stores the OpenMM System initialised class. It stores all the forcefield
        information for the hps model.


    Methods
    -------


    loadForcefieldFromFile()
        Loads forcefield parameters from a force field file written with
        the :code:`dumpForceFieldData()` method.

    """

    def __init__(self, structure_path, model):
        """
        Initialises the hps OpenMM system class.

        Parameters
        ----------
        structure_path : string [requires]
            Name of the input PDB or CIF file
        model: 'hps_kr','hps_urry', or 'hps_ss' [optional, default='hps_urry']
            Hydropathy scale. Currently, there are three models are supported.
        Returns
        -------
        None
        """

        # Define structure object attributes

        self.structure_path = structure_path
        # Recognize format of input structure file
        if structure_path.endswith('.pdb'):
            self.structure = openmm.app.PDBFile(structure_path)
        elif structure_path.endswith('.cif'):
            self.structure = openmm.app.pdbxfile.PDBxFile(structure_path)
        else:
            raise ValueError(
                'Structure file extension not recognized. It must end with .pdb or .cif accordingly.')
        self.topology = self.structure.topology
        self.positions = self.structure.positions

        self.model = model
        # particle properties
        self.particles_mass = None
        self.rf_sigma = None  # particle vdw radius
        self.particles_charge = None
        self.particles_hps = None  # hydropathy scale of particles
        self.particle_type_id = None

        # Define geometric attributes
        self.atoms = []
        self.n_atoms = None

        self.chains = []
        self.n_chains = None

        self.bonds = OrderedDict()
        self.bonds_indexes = []
        self.n_bonds = None
        self.bond_length = model_parameters.parameters[model]["bond_length"]
        self.bondedTo = None
        self.harmonicBondForce = None

        self.angles = OrderedDict()
        self.angles_indexes = []
        self.n_angles = None
        self.gaussianAngleForce = None

        self.torsions = OrderedDict()
        self.torsions_indexes = []
        self.n_torsions = None
        self.gaussianTorsionForce = None

        # Exclusion rule for nonbonded forces
        self.bonded_exclusions_index = model_parameters.parameters[model]["bonded_exclusions_index"]

        # Parameters for PairWise potential Force
        self.ashbaugh_HatchForce = None

        # instance for Wang-Frenkel potential, alternative to Ashbaugh_Hatch
        self.wang_Frenkel_Force = None

        # Define parameters for DH potential Force
        self.yukawaForce = None

        self.forceGroups = OrderedDict()

        # Initialise an OpenMM system class instance
        self.system = openmm.System()

    def getCAlphaOnly(self) -> None:
        """
        Filter in only alpha carbon atoms from the input structure and updates
        the topology object to add new bonds between them. Used specially for
        creating alpha-carbon (CA) coarse-grained models.

        Keeps in the :code:`hps system` only the alpha carbon atoms from the :code:`OpenMM topology`.

        Parameters
        ----------


        Returns
        -------
        None
        """

        # save all non C-alpha atoms
        atoms_to_remove = []
        # oldIndex = []
        for a in self.topology.atoms():
            if a.name != 'CA':
                atoms_to_remove.append(a)

        # Remove all non C-alpha atoms
        modeller_topology = openmm.app.modeller.Modeller(self.topology, self.positions)
        modeller_topology.delete(atoms_to_remove)
        self.topology = modeller_topology.getTopology()
        self.positions = modeller_topology.getPositions()

        chains = []
        for chain in self.topology.chains():
            chains.append(chain)
            # self.n_chains += 1

        self.n_chains = 0
        for chain in chains:
            self.chains.append(chain)
            self.n_chains += 1

        """
            Update system atoms, then add bonds between C-alpha atoms of the same chain.
            openMM load input PDB file and separate chain if encounter TER instruction. It doesn't matter if two chains
            have the same chain name.
            This may have a vulnerable is that if two CA atom on the same chain not consecutive, like residue
            8 and 10 are listed consecutively on input file but residue 9 was missing. The code will add bond between
            residue 8 and 10 which cause large bond. The better solution is check the input file carefully and the logic
            condition in code has nothing to do.
        """
        atoms = list(self.topology.atoms())
        for i in np.arange(1, len(atoms)):
            if atoms[i].residue.chain == atoms[i - 1].residue.chain:
                self.topology.addBond(atoms[i - 1], atoms[i])

    def getAtoms(self):
        """
        Reads atoms from topology, adds them to the main class and sorts them
        into a dictionary to store their forcefield properties.

        After getCAlphaOnly, C-alpha atoms are stored on :code:`self.topology only`.
        We need to add them to atoms attribute and system also.
        Adds :code:`atoms` in the :code:`OpenMM topology` instance to the :code:`hpsOpenMM system` class.

        Parameters
        ----------

        Returns
        -------
        None
        """

        # Get Atoms From Topology
        atoms = []
        for atom in self.topology.atoms():
            atoms.append(atom)

        # Sort atoms by index
        atoms = sorted(atoms, key=lambda x: x.index)

        # Add atoms to hps object
        self.n_atoms = 0
        for atom in atoms:
            self.atoms.append(atom)
            self.n_atoms += 1

    def getBonds(self, except_chains=None):
        """
        Reads bonds from topology, adds them to the main class and sorts them
        into a dictionary to store their forcefield properties.

        Adds :code:`bonds` in the :code:`OpenMM topology` instance to the :code:`hpsOpenMM system` class.

        Parameters
        ----------
        except_chains: String [optional]

        Returns
        -------
        None
        """

        if isinstance(except_chains, str):
            except_chains = list(except_chains)

        # Get Bonds From Topology
        bonds = []
        for bond in self.topology.bonds():
            if except_chains is not None:
                if bond[0].residue.chain.id not in except_chains:
                    if bond[1].residue.chain.id not in except_chains:
                        if bond[0].index > bond[1].index:
                            bonds.append((bond[1], bond[0]))
                        else:
                            bonds.append((bond[0], bond[1]))
            else:
                if bond[0].index > bond[1].index:
                    bonds.append((bond[1], bond[0]))
                else:
                    bonds.append((bond[0], bond[1]))

        # Sort bonds by index of first atom
        bonds = sorted(bonds, key=lambda x: x[0].index)

        # Add bonds to hps object
        self.n_bonds = 0
        for bond in bonds:
            bond_length = self.bond_length * openmm.unit.nanometer
            self.bonds[bond] = (bond_length, None)
            self.n_bonds += 1

            # Store bond indexes
            self.bonds_indexes.append((bond[0].index,
                                       bond[1].index))

        # Record which atoms are bonded to each other
        self.bondedTo = {}
        for atom in self.topology.atoms():
            self.bondedTo[atom] = []
        for bond in self.topology.bonds():
            self.bondedTo[bond[0]].append(bond[1])
            self.bondedTo[bond[1]].append(bond[0])

    def getAngles(self):
        """
        Reads bonds from topology, adds them to the main class and sorts them
        into a dictionary to store their forcefield properties.

        Angles are built from bond (bondedTo instance).
        This function modify :code:`self.angles`, :code:`self.angles_indexes` and :code:`self.n_angles`
        """
        unique_angles = set()
        for bond in self.bonds:
            for atom in self.bondedTo[bond[0]]:
                if atom != bond[1]:
                    if atom.index < bond[1].index:
                        unique_angles.add((atom, bond[0], bond[1]))
                    else:
                        unique_angles.add((bond[1], bond[0], atom))
            for atom in self.bondedTo[bond[1]]:
                if atom != bond[0]:
                    if atom.index > bond[0].index:
                        unique_angles.add((bond[0], bond[1], atom))
                    else:
                        unique_angles.add((atom, bond[1], bond[0]))

        # sort angles by index of first atom
        unique_angles = sorted(list(unique_angles), key=lambda x: x[0].index)

        # add angles to hps object
        self.n_angles = 0
        for angle in unique_angles:
            self.angles[angle] = None
            self.n_angles += 1
            self.angles_indexes.append((angle[0].index, angle[1].index, angle[2].index))

    def getTorsions(self):
        """
        Reads bonds from topology, adds them to the main class and sorts them
        into a dictionary to store their forcefield properties.

        Torsion Angles are built from angles and bondedTo instance.
        This function modify :code:`self.torsions`, :code:`self.torsions_indexes` and :code:`self.n_torsions`
        """
        unique_torsions = set()
        for angle in self.angles:
            for atom in self.bondedTo[angle[0]]:
                if atom not in angle:
                    if atom.index < angle[2].index:
                        unique_torsions.add((atom, angle[0], angle[1], angle[2]))
                    else:
                        unique_torsions.add((angle[2], angle[1], angle[0], atom))

            for atom in self.bondedTo[angle[2]]:
                if atom not in angle:
                    """ 
                        it is stupid for checking atom.index as following line but it happens in all-atom model, when
                        chain of atoms my be branched. just keep it here.
                    """
                    if atom.index > angle[0].index:
                        unique_torsions.add((angle[0], angle[1], angle[2], atom))
                    else:
                        unique_torsions.add((atom, angle[2], angle[1], angle[0]))

        # sort dihedral angles by the first atom index
        unique_torsions = sorted(list(unique_torsions), key=lambda x: x[0].index)

        # load parameter eps_di from parameter file
        params = model_parameters.parameters[self.model]
        # add dihedral angle to hps object
        self.n_torsions = 0
        for torsion in unique_torsions:
            """
            when there are two residues bonded to first atom of torsion mean that this torsion angle is not the first
            angle. One residues will be in list of torsion atoms and the other is out of the list.
            Atom not in the list of torsion atom is preceding of residue (i).
            
            Same as last atom of torsion angle, atom not in torsion atom list is succeeding atom, (i+4)
            
            When there is only one atom bonded to first or last atom in torsion angle, means this is the first or 
            the last torsion angle.
            """
            if len(self.bondedTo[torsion[0]]) > 1:
                for atom in self.bondedTo[torsion[0]]:
                    if atom not in torsion:
                        eps_di_preceding = params[atom.residue.name]['eps_di']
                        weight_preceding = 1
            else:
                eps_di_preceding, weight_preceding = 0, 0

            if len(self.bondedTo[torsion[3]]) > 1:
                for atom in self.bondedTo[torsion[3]]:
                    if atom not in torsion:
                        eps_di_succeeding = params[atom.residue.name]['eps_di']
                        weight_succeeding = 1
            else:
                eps_di_succeeding, weight_succeeding = 0, 0

            # using weighting rule 1-1001-1: (i-1), i, (i+3), (i+4) are 1 and (i+1), (i+2) are 0
            eps_d = (eps_di_preceding + params[torsion[0].residue.name]['eps_di'] + params[torsion[3].residue.name][
                'eps_di'] + eps_di_succeeding) / (weight_succeeding + weight_preceding + 2)

            self.torsions[torsion] = (eps_d, None)  # add torsion angle parameters
            self.n_torsions += 1
            self.torsions_indexes.append((torsion[0].index, torsion[1].index, torsion[2].index, torsion[3].index))

    """ Functions for setting force specific parameters """

    def setBondForceConstants(self, bond_force_constant):
        """
        Change the forcefield parameters for bonded terms.

        Set the harmonic bond constant force parameters. The input can be
        a float, to set the same parameter for all force interactions, or
        a list, to define a unique parameter for each force interaction.

        Parameters
        ----------
        bond_force_constant : float or list
            Parameter(s) to set up for the harmonic bond forces.

        Returns
        -------
        None
        """

        system._setParameters(self.bonds, bond_force_constant)

    def setParticlesMass(self, particles_mass):
        """
        Change the mass parameter for each atom in the system.

        Set the masses of the particles in the system. The input can be a
        float, to set the same mass for all particles, or a list, to define
        a unique mass for each particle.

        Parameters
        ----------
        particles_mass : float or list
            Mass(es) values to add for the particles in the hpsOpenMM system class.

        Returns
        -------
        None
        """

        self.particles_mass = particles_mass

    def setParticlesRadii(self, particles_radii):
        """
        Change the excluded volume radius parameter for each atom in the system.

        Set the radii of the particles in the system. The input can be a
        float, to set the same radius for all particles, or a list, to define
        a unique radius for each particle.

        Parameters
        ----------
        particles_radii : float or list
            Radii values to add for the particles in the hpsOpenMM system class.

        Returns
        -------
        None
        """

        self.rf_sigma = particles_radii

    def setParticlesCharge(self, particles_charge):
        """
        Set the charge of the particles in the system. The input can be a
        float, to set the same charge for all particles, or a list, to define
        a unique charge for each particle.

        Parameters
        ----------
        particles_charge : float or list
            Charge values to add for the particles in the hpsOpenMM system class.

        Returns
        -------
        None
        """

        self.particles_charge = particles_charge

    def setParticlesHPS(self, particles_hps):
        """
        Set the hydropathy scale of the particles in the system. The input can be a
        float, to set the same hydropathy for all particles, or a list, to define
        a unique hydropathy for each particle.

        Parameters
        ----------
        particles_hps : float or list
            HPS scale values to add for the particles in the hpsOpenMM system class.

        Returns
        -------
        None
        """

        self.particles_hps = particles_hps

    def setParticleTypeID(self, particle_id):
        self.particle_type_id = particle_id

    """ Functions for creating force objects with defined parameters """

    def addHarmonicBondForces(self) -> None:
        """
        Creates a harmonic bonded force term for each bond in the main
        class using their defined forcefield parameters.

        Creates an :code:`openmm.HarmonicBondForce()` object with the bonds and
        parameters set up in the "bonds" dictionary attribute. The force object
        is stored at the :code:`harmonicBondForce` attribute.

        openMM uses harmonic bond that has energy term of form:

        .. math::
            E= \\frac{1}{2}k(r-r_0)^2

        The force parameters must be contained in self.bonds as follows:

        self.bonds is a dictionary:
            - The keys are 2-tuples for two atom items in :code:`self.topology.atoms` attribute.
            - The values are a 2-tuple of parameters in the following order:
                - first  -> bond0 (quantity)
                - second -> k (float) (measured in unit of kj/mol/nm^2)


        Returns
        -------
        None
        """

        self.harmonicBondForce = openmm.HarmonicBondForce()
        for bond in self.bonds:
            self.harmonicBondForce.addBond(bond[0].index,
                                           bond[1].index,
                                           self.bonds[bond][0],
                                           self.bonds[bond][1])

    def addGaussianAngleForces(self) -> None:
        """
        Add Gaussian functional form of angle.
        Note that in openMM log is neutral logarithm.

        Angle potential take Gaussian functional form in hps-ss model.

        .. math::
            U_{angle}(\\theta) = \\frac{-1}{\gamma}
                \\ln{[e^{-\gamma[k_\\alpha( \\theta-\\theta_\\alpha)^2+\\epsilon_\\alpha]}
                +e^{-\\gamma k_\\beta(\\theta-\\theta_\\beta)^2}]}

        Angle potential is taken from reference:
        """

        gamma = 0.0239 / openmm.unit.kilojoule_per_mole  # 0.1 mol/kcal
        eps_alpha = 17.9912 * openmm.unit.kilojoule_per_mole  # 4.3 kcal/mol
        theta_alpha = 1.6 * openmm.unit.radian
        theta_beta = 2.27 * openmm.unit.radian
        k_alpha = 445.1776 * openmm.unit.kilojoule_per_mole / openmm.unit.radian ** 2
        k_beta = 110.0392 * openmm.unit.kilojoule_per_mole / openmm.unit.radian ** 2

        energy_function = '(-1 / gamma) * log(exp(-gamma * (k_alpha * (theta - theta_alpha) ^ 2 + eps_alpha)) ' \
                          '+ exp(-gamma * k_beta * (theta - theta_beta) ^ 2))'
        self.gaussianAngleForce = openmm.CustomAngleForce(energy_function)
        self.gaussianAngleForce.addGlobalParameter('gamma', gamma)
        self.gaussianAngleForce.addGlobalParameter('eps_alpha', eps_alpha)
        self.gaussianAngleForce.addGlobalParameter('theta_alpha', theta_alpha)
        self.gaussianAngleForce.addGlobalParameter('theta_beta', theta_beta)
        self.gaussianAngleForce.addGlobalParameter('k_alpha', k_alpha)
        self.gaussianAngleForce.addGlobalParameter('k_beta', k_beta)

        for angle in self.angles:
            self.gaussianAngleForce.addAngle(angle[0].index, angle[1].index, angle[2].index)

    def addGaussianTorsionForces(self) -> None:
        """
        Torsion potential in hps-ss model takes the form:

        .. math::
            U_{torsion}(\\theta) = -\\ln\\left[ U_{torsion, \\alpha}(\\theta, \\epsilon_d) + U_{torsion, \\beta}(\\theta, \\epsilon_d)\\right]

        where,


        .. math::

            U_{torsion, \\alpha}(\\theta, \\epsilon_d)  &= e^{-k_{\\alpha, 1}(\\theta-\\theta_{\\alpha,1})^2-\\epsilon_d}
                                                        + e^{-k_{\\alpha, 2}(\\theta-\\theta_{\\alpha,2})^4 + e_0}
                                                        + e^{-k_{\\alpha, 2}(\\theta-\\theta_{\\alpha,2}+2\\pi)^4 + e_0} \\


            U_{torsion, \\beta}(\\theta, \\epsilon_d) &= e^{-k_{\\beta,1}(\\theta-\\theta_{\\beta,1})^2+e_1+\\epsilon_d}
                                                    + e^{-k_{\\beta,1}(\\theta-\\theta_{\\beta,1}-2\\pi)^2+e_1+\\epsilon_d} \\

                                                    &+ e^{-k_{\\beta,2}(\\theta-\\theta_{\\beta,2})^4+e_2}
                                                    + e^{-k_{\\beta,2}(\\theta-\\theta_{\\beta,2}-2\\pi)^4+e_2}

        """

        k_alpha1 = 47.6976 * openmm.unit.kilojoule_per_mole / openmm.unit.radian ** 2  # 11.4 kcal/mol/rad^2
        k_alpha2 = 0.6276 * openmm.unit.kilojoule_per_mole / openmm.unit.radian ** 4  # 0.15 kcal/mol/rad^4
        theta_alpha1 = 0.9 * openmm.unit.radian
        theta_alpha2 = 1.02 * openmm.unit.radian
        e_0 = 1.12968 * openmm.unit.kilojoule_per_mole  # 0.27 kcal/mol

        k_beta1 = 7.5312 * openmm.unit.kilojoule_per_mole / openmm.unit.radian ** 2  # 1.8 kcal/mol/rad^2
        k_beta2 = 2.7196 * openmm.unit.kilojoule_per_mole / openmm.unit.radian ** 4  # 0.65 kcal/mol/rad^4
        theta_beta1 = -1.55 * openmm.unit.radian
        theta_beta2 = -2.5 * openmm.unit.radian
        e_1 = 0.58576 * openmm.unit.kilojoule_per_mole  # 0.14 kcal/mol
        e_2 = 1.6736 * openmm.unit.kilojoule_per_mole  # 0.4 kcal/mol
        # pi = np.pi

        energy_function_alpha = 'exp(-k_alpha1*(theta-theta_alpha1)^2 - esp_d)'
        energy_function_alpha += '+exp(-k_alpha2*(theta-theta_alpha2)^4 + e_0)'
        energy_function_alpha += '+exp(-k_alpha2*(theta-theta_alpha2+2*pi)^4 + e_0)'

        energy_function_beta = '+exp(-k_beta1*(theta-theta_beta1)^2 + e_1 + esp_d)'
        energy_function_beta += '+exp(-k_beta1*(theta-theta_beta1-2*pi)^2 + e_1 + esp_d)'
        energy_function_beta += '+exp(-k_beta2*(theta-theta_beta2)^4 + e_2)'
        energy_function_beta += '+exp(-k_beta2*(theta-theta_beta2-2*pi)^4 + e_2)'

        energy_function = '-log(' + energy_function_alpha + energy_function_beta + ')'
        self.gaussianTorsionForce = openmm.CustomTorsionForce(energy_function)
        self.gaussianTorsionForce.addGlobalParameter("k_alpha1", k_alpha1)
        self.gaussianTorsionForce.addGlobalParameter("theta_alpha1", theta_alpha1)
        self.gaussianTorsionForce.addGlobalParameter("k_alpha2", k_alpha2)
        self.gaussianTorsionForce.addGlobalParameter("theta_alpha2", theta_alpha2)
        self.gaussianTorsionForce.addGlobalParameter("e_0", e_0)

        self.gaussianTorsionForce.addGlobalParameter("k_beta1", k_beta1)
        self.gaussianTorsionForce.addGlobalParameter("theta_beta1", theta_beta1)
        self.gaussianTorsionForce.addGlobalParameter("k_beta2", k_beta2)
        self.gaussianTorsionForce.addGlobalParameter("theta_beta2", theta_beta2)
        self.gaussianTorsionForce.addGlobalParameter("e_1", e_1)
        self.gaussianTorsionForce.addGlobalParameter("e_2", e_2)
        self.gaussianTorsionForce.addGlobalParameter("pi", np.pi)

        self.gaussianTorsionForce.addPerTorsionParameter("esp_d")

        for torsion in self.torsions:
            self.gaussianTorsionForce.addTorsion(torsion[0].index, torsion[1].index, torsion[2].index, torsion[3].index,
                                                 (self.torsions[torsion][0],))

    def addYukawaForces(self, use_pbc: bool) -> None:
        """
        Creates a nonbonded force term for electrostatic interaction DH potential.

        Creates an :code:`openmm.CustomNonbondedForce()` object with the parameters
        sigma and epsilon given to this method. The custom non-bonded force
        is initialized with the formula:

        .. math::
            energy = f \\times \\frac{q_1q_2}{\epsilon_r \\times r}\\times e^{(-r/lD)}


        where :math:`f=\\frac{1}{4\\pi\\epsilon_0}=138.935458` is the factor for short to convert dimensionless
        in calculation to :math:`kj.nm/(mol\\times e^2)` unit.

        :math:`\\epsilon_r=80`: Dielectric constant of water at 100mM mono-valent ion

        The force object is stored at the :code:`yukawaForce` attribute.

        Parameters
        ----------
        use_pbc: (bool) whether use PBC, cutoff periodic boundary condition

        Returns
        -------
        None
        """

        lD = 1.0 * openmm.unit.nanometer
        electric_factor = 138.935458 * openmm.unit.kilojoule_per_mole * openmm.unit.nanometer / openmm.unit.elementary_charge ** 2
        yukawa_cutoff = 3.5 * openmm.unit.nanometer
        epsilon_r = 80.0

        energy_function = 'factor*charge1*charge2/epsilon_r/r*exp(-r/lD)'
        self.yukawaForce = openmm.CustomNonbondedForce(energy_function)
        self.yukawaForce.addGlobalParameter('factor', electric_factor)
        self.yukawaForce.addGlobalParameter('epsilon_r', epsilon_r)
        self.yukawaForce.addGlobalParameter('lD', lD)
        self.yukawaForce.addPerParticleParameter('charge')
        if use_pbc:
            self.yukawaForce.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        else:
            self.yukawaForce.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)

        self.yukawaForce.setCutoffDistance(yukawa_cutoff)

        if isinstance(self.particles_charge, float):
            for i in np.arange(len(self.atoms)):
                self.yukawaForce.addParticle((self.particles_charge,))

        # in the case each atom has different sigma para.
        elif isinstance(self.particles_charge, list):
            assert self.n_atoms == len(self.particles_charge)
            for i, atom in enumerate(self.atoms):
                self.yukawaForce.addParticle((self.particles_charge[i],))

        # set exclusions rule
        bonded_exclusions = [(b[0].index, b[1].index) for b in list(self.topology.bonds())]
        self.yukawaForce.createExclusionsFromBonds(bonded_exclusions, self.bonded_exclusions_index)

    def addAshbaughHatchForces(self, use_pbc: bool) -> None:
        """
        Creates a nonbonded force term for pairwise interaction (customize LJ 12-6 potential).

        Creates an :code:`openmm.CustomNonbondedForce()` object with the parameters
        sigma and epsilon given to this method. The custom non-bonded force
        is initialized with the formula: (note: hps here is :math:`\lambda_{ij}^{0}` in the paper)

        Unlike :code:`BondForce` class, where we specify index for atoms pair to add bond, it means
        that number of bondForces may differ from number of particle.
        :code:`NonBondedForce` is added to all particles, hence we don't need to pass the :code:`atom index`.

        .. math::
            \\Phi_{i,j}^{vdw}(r) = step(2^{1/6}\\sigma_{ij}-r) \\times
            \\left( 4\\epsilon\\left[\\left(\\frac{\\sigma_{ij}}{r}\\right)^{12}-
            \\left(\\frac{\\sigma_{ij}}{r}\\right)^{6}\\right]+(1-\\lambda_{ij})\\epsilon\\right)

            + \\left[1-step(2^{1/6}\\sigma_{ij}-r)\\right]\\times\\left[(\\lambda_{ij})\\times 4\\epsilon
            \\left[\\left(\\frac{\\sigma_{ij}}{r}\\right)^{12}-\\left(\\frac{\\sigma_{ij}}{r}\\right)^6\\right]\\right]



        Here, :math:`\\sigma= \\frac{(\\sigma_1+\\sigma_2)}{2}; \\lambda_{ij}^{0}=\\frac{(\\lambda_i+\\lambda_j)}{2};
        \\epsilon = 0.8368 kj/mol`

        The force object is stored at the :code:`ashbaugh_HatchForce` attribute.

        epsilon : float
            Value of the epsilon constant in the energy function.
        sigma : float or list
            Value of the sigma constant (in nm) in the energy function. If float the
            same sigma value is used for every particle. If list a unique
            parameter is given for each particle.
        cutoff : float
            The cutoff distance (in nm) being used for the non-bonded interactions.

        Parameters
        ----------

        use_pbc : bool. Whether use PBC, cutoff periodic boundary condition

        Returns
        -------
        None
        """
        epsilon = 0.8368 * openmm.unit.kilojoule_per_mole
        ashbaugh_Hatch_cutoff = 2.0 * openmm.unit.nanometer

        energy_function = 'step(2^(1/6)*sigma - r) *'
        energy_function += '(4*epsilon* ((sigma/r)^12-(sigma/r)^6) + (1-hps)*epsilon )'
        energy_function += '+(1-step(2^(1/6)*sigma-r)) * (hps*4*epsilon*((sigma/r)^12-(sigma/r)^6));'
        energy_function += 'sigma=0.5*(sigma1+sigma2);'
        energy_function += 'hps=0.5*(hps1+hps2)'
        self.ashbaugh_HatchForce = openmm.CustomNonbondedForce(energy_function)
        self.ashbaugh_HatchForce.addGlobalParameter('epsilon', epsilon)
        self.ashbaugh_HatchForce.addPerParticleParameter('sigma')
        self.ashbaugh_HatchForce.addPerParticleParameter('hps')
        #
        if use_pbc:
            self.ashbaugh_HatchForce.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        else:
            self.ashbaugh_HatchForce.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)

        self.ashbaugh_HatchForce.setCutoffDistance(ashbaugh_Hatch_cutoff)

        if isinstance(self.rf_sigma, float):
            for i in np.arange(len(self.atoms)):
                self.ashbaugh_HatchForce.addParticle((self.rf_sigma, self.particles_hps[i],))

        # in the case each atom has different sigma para.
        elif isinstance(self.rf_sigma, list):
            assert self.n_atoms == len(self.rf_sigma) and self.n_atoms == len(self.particles_hps)
            for i, atom in enumerate(self.atoms):
                self.ashbaugh_HatchForce.addParticle((self.rf_sigma[i], self.particles_hps[i],))

        # set exclusions rule
        bonded_exclusions = [(b[0].index, b[1].index) for b in list(self.topology.bonds())]
        self.ashbaugh_HatchForce.createExclusionsFromBonds(bonded_exclusions, self.bonded_exclusions_index)

    def add_Wang_Frenkel_Forces(self, use_pbc: bool):
        """
        MPIPI model. using TabulatedFunction for pair interaction.
        More information about TabulatedFUnction can be found here:
        http://docs.openmm.org/7.2.0/api-c++/generated/OpenMM.Discrete2DFunction.html
        """
        wang_frenkel_cutoff = 2.5 * openmm.unit.nanometer

        table_eps = np.array([[0.165536, 0.284583, 0.136758, 0.147114, 0.717619,
                               0.186280, 0.306152, 0.324582, 0.959705, 0.094437,
                               0.105776, 0.502105, 1.233991, 0.902083, 0.211635,
                               0.851406, 0.488298, 0.247358, 0.227764, 0.083592],
                              [0.284583, 0.403630, 0.255806, 0.266161, 0.836666,
                               0.305327, 0.425199, 0.443630, 1.078752, 0.213484,
                               0.224823, 0.621152, 1.353039, 1.021130, 0.330682,
                               0.970454, 0.607345, 0.366405, 0.346812, 0.202639],
                              [0.136758, 0.255806, 0.079986, 0.118336, 0.510189,
                               0.157502, 0.002063, 0.002201, 0.429701, 0.065660,
                               0.076998, 0.473328, 0.438466, 0.482315, 0.182858,
                               0.365556, 0.459520, 0.218581, 0.198987, 0.054815],
                              [0.147114, 0.266161, 0.118336, 0.128691, 0.699197,
                               0.167858, 0.287729, 0.306160, 0.941283, 0.076015,
                               0.087354, 0.483683, 1.215565, 0.883661, 0.193213,
                               0.832984, 0.469876, 0.228936, 0.209342, 0.065170],
                              [0.717619, 0.836666, 0.510189, 0.699197, 0.376209,
                               0.738363, 0.006368, 0.006502, 2.543119, 0.646520,
                               0.657859, 1.054188, 2.893893, 2.280037, 0.763718,
                               0.519808, 1.040381, 0.799441, 0.779847, 0.635675],
                              [0.186280, 0.305327, 0.157502, 0.167858, 0.738363,
                               0.207024, 0.326896, 0.345326, 0.980449, 0.115181,
                               0.126520, 0.522849, 1.254736, 0.922827, 0.232379,
                               0.872151, 0.509046, 0.268102, 0.248509, 0.104336],
                              [0.306152, 0.425199, 0.002063, 0.287729, 0.006368,
                               0.326896, 0.330938, 0.344590, 1.100321, 0.235049,
                               0.246387, 0.642721, 1.374603, 1.042699, 0.352247,
                               0.007355, 0.628914, 0.387974, 0.368380, 0.224208],
                              [0.324582, 0.443630, 0.002201, 0.306160, 0.006502,
                               0.345326, 0.344590, 0.358242, 1.118751, 0.253479,
                               0.264818, 0.661151, 1.393034, 1.061129, 0.370677,
                               0.007494, 0.647344, 0.406404, 0.386811, 0.242639],
                              [0.959705, 1.078752, 0.429701, 0.941283, 2.543119,
                               0.980449, 1.100321, 1.118751, 1.753874, 0.888606,
                               0.899945, 1.296274, 2.028156, 1.696252, 1.005804,
                               1.645576, 1.282467, 1.041527, 1.021934, 0.877761],
                              [0.094437, 0.213484, 0.065660, 0.076015, 0.646520,
                               0.115181, 0.235049, 0.253479, 0.888606, 0.023338,
                               0.034677, 0.431006, 1.162888, 0.830984, 0.140536,
                               0.780308, 0.417199, 0.176259, 0.156666, 0.012493],
                              [0.105776, 0.224823, 0.076998, 0.087354, 0.657859,
                               0.126520, 0.246387, 0.264818, 0.899945, 0.034677,
                               0.046016, 0.442345, 1.174227, 0.842323, 0.151875,
                               0.791646, 0.428538, 0.187598, 0.168004, 0.023832],
                              [0.502105, 0.621152, 0.473328, 0.483683, 1.054188,
                               0.522849, 0.642721, 0.661151, 1.296274, 0.431006,
                               0.442345, 0.838674, 1.570556, 1.238652, 0.548204,
                               1.187976, 0.824867, 0.583927, 0.564334, 0.420161],
                              [1.233991, 1.353039, 0.438466, 1.215565, 2.893893,
                               1.254736, 1.374603, 1.393034, 2.028156, 1.162888,
                               1.174227, 1.570556, 2.302443, 1.970538, 1.280086,
                               1.919858, 1.556753, 1.315814, 1.296220, 1.152048],
                              [0.902083, 1.021130, 0.482315, 0.883661, 2.280037,
                               0.922827, 1.042699, 1.061129, 1.696252, 0.830984,
                               0.842323, 1.238652, 1.970538, 1.638630, 0.948182,
                               1.587954, 1.224845, 0.983905, 0.964312, 0.820139],
                              [0.211635, 0.330682, 0.182858, 0.193213, 0.763718,
                               0.232379, 0.352247, 0.370677, 1.005804, 0.140536,
                               0.151875, 0.548204, 1.280086, 0.948182, 0.257734,
                               0.897506, 0.534397, 0.293457, 0.273864, 0.129691],
                              [0.851406, 0.970454, 0.365556, 0.832984, 0.519808,
                               0.872151, 0.007355, 0.007494, 1.645576, 0.780308,
                               0.791646, 1.187976, 1.919858, 1.587954, 0.897506,
                               0.113872, 1.174168, 0.933229, 0.913635, 0.769463],
                              [0.488298, 0.607345, 0.459520, 0.469876, 1.040381,
                               0.509046, 0.628914, 0.647344, 1.282467, 0.417199,
                               0.428538, 0.824867, 1.556753, 1.224845, 0.534397,
                               1.174168, 0.811064, 0.570124, 0.550531, 0.406358],
                              [0.247358, 0.366405, 0.218581, 0.228936, 0.799441,
                               0.268102, 0.387974, 0.406404, 1.041527, 0.176259,
                               0.187598, 0.583927, 1.315814, 0.983905, 0.293457,
                               0.933229, 0.570124, 0.329185, 0.309591, 0.165419],
                              [0.227764, 0.346812, 0.198987, 0.209342, 0.779847,
                               0.248509, 0.368380, 0.386811, 1.021934, 0.156666,
                               0.168004, 0.564334, 1.296220, 0.964312, 0.273864,
                               0.913635, 0.550531, 0.309591, 0.289997, 0.145825],
                              [0.083592, 0.202639, 0.054815, 0.065170, 0.635675,
                               0.104336, 0.224208, 0.242639, 0.877761, 0.012493,
                               0.023832, 0.420161, 1.152048, 0.820139, 0.129691,
                               0.769463, 0.406358, 0.165419, 0.145825, 0.001653]])
        table_eps_ravel = table_eps.ravel().tolist()

        table_sigma = np.array([[0.646795, 0.557618, 0.656778, 0.617823, 0.664396,
                                 0.586850, 0.614146, 0.631818, 0.659004, 0.632057,
                                 0.648344, 0.636519, 0.675573, 0.653821, 0.593894,
                                 0.639251, 0.618802, 0.613446, 0.609429, 0.649639],
                                [0.557618, 0.469511, 0.567134, 0.528442, 0.576639,
                                 0.498013, 0.525925, 0.543637, 0.571317, 0.541239,
                                 0.557893, 0.548630, 0.587924, 0.566122, 0.505252,
                                 0.551541, 0.530902, 0.524999, 0.520872, 0.558033],
                                [0.656778, 0.567134, 0.667134, 0.627940, 0.673819,
                                 0.596669, 0.623699, 0.641358, 0.668408, 0.643354,
                                 0.659256, 0.645969, 0.684969, 0.663228, 0.603618,
                                 0.648660, 0.628255, 0.623106, 0.619115, 0.661800],
                                [0.617823, 0.528442, 0.627940, 0.588906, 0.635202,
                                 0.557795, 0.584987, 0.602654, 0.629806, 0.603697,
                                 0.619809, 0.607329, 0.646378, 0.624625, 0.564801,
                                 0.610053, 0.589610, 0.584320, 0.580321, 0.621462],
                                [0.664396, 0.576639, 0.673819, 0.635202, 0.683905,
                                 0.604909, 0.632986, 0.650709, 0.678624, 0.647780,
                                 0.664461, 0.655828, 0.695255, 0.673424, 0.612191,
                                 0.658837, 0.638103, 0.631985, 0.627828, 0.664491],
                                [0.586850, 0.498013, 0.596669, 0.557795, 0.604909,
                                 0.527007, 0.554510, 0.572198, 0.599537, 0.571467,
                                 0.587918, 0.576988, 0.616122, 0.594352, 0.534119,
                                 0.579777, 0.559268, 0.553738, 0.549687, 0.588726],
                                [0.614146, 0.525925, 0.623699, 0.584987, 0.632986,
                                 0.554510, 0.582352, 0.600058, 0.627647, 0.597874,
                                 0.614509, 0.605004, 0.644244, 0.622456, 0.561730,
                                 0.607876, 0.587280, 0.581455, 0.577340, 0.614699],
                                [0.631818, 0.543637, 0.641358, 0.602654, 0.650709,
                                 0.572198, 0.600058, 0.617767, 0.645370, 0.615504,
                                 0.632146, 0.622721, 0.661969, 0.640180, 0.579422,
                                 0.625600, 0.604999, 0.599152, 0.595036, 0.632333],
                                [0.659004, 0.571317, 0.668408, 0.629806, 0.678624,
                                 0.599537, 0.627647, 0.645370, 0.673363, 0.642355,
                                 0.659039, 0.650525, 0.690005, 0.668159, 0.606828,
                                 0.653569, 0.632798, 0.626630, 0.622465, 0.659055],
                                [0.632057, 0.541239, 0.643354, 0.603697, 0.647780,
                                 0.571467, 0.597874, 0.615504, 0.642355, 0.626600,
                                 0.639093, 0.619973, 0.658912, 0.637178, 0.578166,
                                 0.622610, 0.602261, 0.597433, 0.593567, 0.671801],
                                [0.648344, 0.557893, 0.659256, 0.619809, 0.664461,
                                 0.587918, 0.614509, 0.632146, 0.659039, 0.639093,
                                 0.653407, 0.636649, 0.675594, 0.653860, 0.594697,
                                 0.639294, 0.618936, 0.614021, 0.610125, 0.660772],
                                [0.636519, 0.548630, 0.645969, 0.607329, 0.655828,
                                 0.576988, 0.605004, 0.622721, 0.650525, 0.619973,
                                 0.636649, 0.627785, 0.667141, 0.645328, 0.584255,
                                 0.630743, 0.610059, 0.604035, 0.599885, 0.636715],
                                [0.675573, 0.587924, 0.684969, 0.646378, 0.695255,
                                 0.616122, 0.644244, 0.661969, 0.690005, 0.658912,
                                 0.675594, 0.667141, 0.706655, 0.684798, 0.623414,
                                 0.670209, 0.649416, 0.643217, 0.639053, 0.675603],
                                [0.653821, 0.566122, 0.663228, 0.624625, 0.673424,
                                 0.594352, 0.622456, 0.640180, 0.668159, 0.637178,
                                 0.653860, 0.645328, 0.684798, 0.662955, 0.601640,
                                 0.648367, 0.627603, 0.621441, 0.617279, 0.653879],
                                [0.593894, 0.505252, 0.603618, 0.564801, 0.612191,
                                 0.534119, 0.561730, 0.579422, 0.606828, 0.578166,
                                 0.594697, 0.584255, 0.623414, 0.601640, 0.541267,
                                 0.587063, 0.566533, 0.560932, 0.556848, 0.595253],
                                [0.639251, 0.551541, 0.648660, 0.610053, 0.658837,
                                 0.579777, 0.607876, 0.625600, 0.653569, 0.622610,
                                 0.639294, 0.630743, 0.670209, 0.648367, 0.587063,
                                 0.633778, 0.613018, 0.606865, 0.602702, 0.639315],
                                [0.618802, 0.530902, 0.628255, 0.589610, 0.638103,
                                 0.559268, 0.587280, 0.604999, 0.632798, 0.602261,
                                 0.618936, 0.610059, 0.649416, 0.627603, 0.566533,
                                 0.613018, 0.592335, 0.586308, 0.582164, 0.618999],
                                [0.613446, 0.524999, 0.623106, 0.584320, 0.631985,
                                 0.553738, 0.581455, 0.599152, 0.626630, 0.597433,
                                 0.614021, 0.604035, 0.643217, 0.621441, 0.560932,
                                 0.606865, 0.586308, 0.580616, 0.576515, 0.614386],
                                [0.609429, 0.520872, 0.619115, 0.580321, 0.627828,
                                 0.549687, 0.577340, 0.595036, 0.622465, 0.593567,
                                 0.610125, 0.599885, 0.639053, 0.617279, 0.556848,
                                 0.602702, 0.582164, 0.576515, 0.572436, 0.610587],
                                [0.649639, 0.558033, 0.661800, 0.621462, 0.664491,
                                 0.588726, 0.614699, 0.632333, 0.659055, 0.671801,
                                 0.660772, 0.636715, 0.675603, 0.653879, 0.595253,
                                 0.639315, 0.618999, 0.614386, 0.610587, 0.692168]])
        table_sigma_ravel = table_sigma.ravel().tolist()

        table_nu = np.array([[1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.],
                             [1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.,
                              1., 1., 1., 1.]])
        table_nu_ravel = table_nu.ravel().tolist()

        table_mu = np.array([[2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 4.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 2.],
                            [2., 2., 2., 2., 2., 2., 2., 2., 2., 4., 2., 2., 2.,
                             2., 2., 2., 2., 2., 2., 11.]])
        table_mu_ravel = table_mu.ravel().tolist()

        table_rc = np.array([[1.940385, 1.672854, 1.970335, 1.853468, 1.993189,
                              1.760551, 1.842438, 1.895453, 1.977011, 1.896171,
                              1.945033, 1.909556, 2.026719, 1.961463, 1.781682,
                              1.917752, 1.856406, 1.840339, 1.828288, 1.948918],
                             [1.672854, 1.408533, 1.701401, 1.585327, 1.729916,
                              1.494040, 1.577776, 1.630912, 1.713950, 1.623718,
                              1.673679, 1.645890, 1.763771, 1.698366, 1.515757,
                              1.654623, 1.592705, 1.574997, 1.562617, 1.674099],
                             [1.970335, 1.701401, 2.001402, 1.883820, 2.021458,
                              1.790006, 1.871096, 1.924073, 2.005225, 1.930062,
                              1.977767, 1.937908, 2.054906, 1.989684, 1.810854,
                              1.945981, 1.884764, 1.869318, 1.857346, 1.985399],
                             [1.853468, 1.585327, 1.883820, 1.766718, 1.905605,
                              1.673384, 1.754960, 1.807962, 1.889419, 1.811090,
                              1.859427, 1.821986, 1.939133, 1.873875, 1.694402,
                              1.830159, 1.768830, 1.752959, 1.740963, 1.864385],
                             [1.993189, 1.729916, 2.021458, 1.905605, 2.051715,
                              1.814728, 1.898958, 1.952126, 2.035872, 1.943340,
                              1.993382, 1.967484, 2.085765, 2.020273, 1.836572,
                              1.976511, 1.914308, 1.895955, 1.883483, 1.993474],
                             [1.760551, 1.494040, 1.790006, 1.673384, 1.814728,
                              1.581022, 1.663531, 1.716595, 1.798612, 1.714402,
                              1.763755, 1.730964, 1.848366, 1.783057, 1.602356,
                              1.739330, 1.677803, 1.661215, 1.649060, 1.766179],
                             [1.842438, 1.577776, 1.871096, 1.754960, 1.898958,
                              1.663531, 1.747057, 1.800175, 1.882942, 1.793621,
                              1.843527, 1.815012, 1.932732, 1.867367, 1.685190,
                              1.823628, 1.761839, 1.744364, 1.732021, 1.844096],
                             [1.895453, 1.630912, 1.924073, 1.807962, 1.952126,
                              1.716595, 1.800175, 1.853301, 1.936109, 1.846512,
                              1.896437, 1.868164, 1.985908, 1.920540, 1.738266,
                              1.876800, 1.814997, 1.797457, 1.785109, 1.896998],
                             [1.977011, 1.713950, 2.005225, 1.889419, 2.035872,
                              1.798612, 1.882942, 1.936109, 2.020090, 1.927065,
                              1.977116, 1.951575, 2.070015, 2.004477, 1.820484,
                              1.960707, 1.898395, 1.879889, 1.867395, 1.977165],
                             [1.896171, 1.623718, 1.930062, 1.811090, 1.943340,
                              1.714402, 1.793621, 1.846512, 1.927065, 1.879800,
                              1.917280, 1.859919, 1.976735, 1.911535, 1.734499,
                              1.867831, 1.806782, 1.792299, 1.780701, 2.015403],
                             [1.945033, 1.673679, 1.977767, 1.859427, 1.993382,
                              1.763755, 1.843527, 1.896437, 1.977116, 1.917280,
                              1.960222, 1.909947, 2.026781, 1.961581, 1.784091,
                              1.917883, 1.856809, 1.842064, 1.830374, 1.982317],
                             [1.909556, 1.645890, 1.937908, 1.821986, 1.967484,
                              1.730964, 1.815012, 1.868164, 1.951575, 1.859919,
                              1.909947, 1.883355, 2.001422, 1.935983, 1.752766,
                              1.892230, 1.830178, 1.812104, 1.799656, 1.910146],
                             [2.026719, 1.763771, 2.054906, 1.939133, 2.085765,
                              1.848366, 1.932732, 1.985908, 2.070015, 1.976735,
                              2.026781, 2.001422, 2.119964, 2.054395, 1.870242,
                              2.010627, 1.948249, 1.929652, 1.917160, 2.026810],
                             [1.961463, 1.698366, 1.989684, 1.873875, 2.020273,
                              1.783057, 1.867367, 1.920540, 2.004477, 1.911535,
                              1.961581, 1.935983, 2.054395, 1.988865, 1.804920,
                              1.945102, 1.882808, 1.864322, 1.851836, 1.961636],
                             [1.781682, 1.515757, 1.810854, 1.694402, 1.836572,
                              1.602356, 1.685190, 1.738266, 1.820484, 1.734499,
                              1.784091, 1.752766, 1.870242, 1.804920, 1.623801,
                              1.761190, 1.699598, 1.682795, 1.670544, 1.785760],
                             [1.917752, 1.654623, 1.945981, 1.830159, 1.976511,
                              1.739330, 1.823628, 1.876800, 1.960707, 1.867831,
                              1.917883, 1.892230, 2.010627, 1.945102, 1.761190,
                              1.901335, 1.839054, 1.820595, 1.808106, 1.917944],
                             [1.856406, 1.592705, 1.884764, 1.768830, 1.914308,
                              1.677803, 1.761839, 1.814997, 1.898395, 1.806782,
                              1.856809, 1.830178, 1.948249, 1.882808, 1.699598,
                              1.839054, 1.777005, 1.758925, 1.746493, 1.856996],
                             [1.840339, 1.574997, 1.869318, 1.752959, 1.895955,
                              1.661215, 1.744364, 1.797457, 1.879889, 1.792299,
                              1.842064, 1.812104, 1.929652, 1.864322, 1.682795,
                              1.820595, 1.758925, 1.741847, 1.729545, 1.843158],
                             [1.828288, 1.562617, 1.857346, 1.740963, 1.883483,
                              1.649060, 1.732021, 1.785109, 1.867395, 1.780701,
                              1.830374, 1.799656, 1.917160, 1.851836, 1.670544,
                              1.808106, 1.746493, 1.729545, 1.717309, 1.831760],
                             [1.948918, 1.674099, 1.985399, 1.864385, 1.993474,
                              1.766179, 1.844096, 1.896998, 1.977165, 2.015403,
                              1.982317, 1.910146, 2.026810, 1.961636, 1.785760,
                              1.917944, 1.856996, 1.843158, 1.831760, 2.076504]])
        table_rc_ravel = table_rc.ravel().tolist()

        # number of atom types in model. currently with protein, there are 20.
        n_atom_types = table_sigma.shape[0]

        # eps, sigma, nu, mu, rc: load from tabular table
        energy_function = 'eps * 2*nu*(rc/sigma)^(2*mu) * ((2*nu+1)/(2*nu*((rc/sigma)^(2*mu)-1)) )^(2*nu+1)'
        energy_function += '* ((sigma/r)^(2*mu)-1 ) * ((rc/r)^(2*mu)-1)^(2*nu);'
        energy_function += 'eps = eps_table(id1, id2); sigma = sigma_table(id1, id2);'
        energy_function += 'nu = nu_table(id1, id2); mu = mu_table(id1, id2);'
        energy_function += 'rc=rc_table(id1, id2)'

        self.wang_Frenkel_Force = openmm.CustomNonbondedForce(energy_function)
        self.wang_Frenkel_Force.addTabulatedFunction('eps_table', openmm.Discrete2DFunction(n_atom_types, n_atom_types,
                                                                                            table_eps_ravel))
        self.wang_Frenkel_Force.addTabulatedFunction('sigma_table',
                                                     openmm.Discrete2DFunction(n_atom_types, n_atom_types,
                                                                               table_sigma_ravel))
        self.wang_Frenkel_Force.addTabulatedFunction('nu_table', openmm.Discrete2DFunction(n_atom_types, n_atom_types,
                                                                                           table_nu_ravel))
        self.wang_Frenkel_Force.addTabulatedFunction('mu_table', openmm.Discrete2DFunction(n_atom_types, n_atom_types,
                                                                                           table_mu_ravel))
        self.wang_Frenkel_Force.addTabulatedFunction('rc_table', openmm.Discrete2DFunction(n_atom_types, n_atom_types,
                                                                                           table_rc_ravel))
        self.wang_Frenkel_Force.addPerParticleParameter('id')

        for i, atom in enumerate(self.atoms):
            self.wang_Frenkel_Force.addParticle((self.particle_type_id[i],))

        if use_pbc:
            self.wang_Frenkel_Force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        else:
            self.wang_Frenkel_Force.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)

        self.wang_Frenkel_Force.setCutoffDistance(wang_frenkel_cutoff)



        # set exclusion rule
        bonded_exclusions = [(b[0].index, b[1].index) for b in list(self.topology.bonds())]
        self.wang_Frenkel_Force.createExclusionsFromBonds(bonded_exclusions, self.bonded_exclusions_index)

    """ Functions for creating OpenMM system object """

    def createSystemObject(self, check_bond_distances: bool = True, minimize: bool = False,
                           check_large_forces: bool = True, force_threshold: float = 10.0,
                           bond_threshold: float = 0.5) -> None:
        """
        Creates OpenMM system object adding particles, masses and forces.
        It also groups the added forces into Force-Groups for the hpsReporter
        class.

        Creates an :code:`openmm.System()` object using the force field parameters
        given to the 'system' class. It adds particles, forces and
        creates a force group for each force object. Optionally the method
        can check for large bond distances (default) and minimize the atomic
        positions if large forces are found in any atom (default False).

        Parameters
        ----------
        minimize : boolean (False)
            Whether to minimize the system if large forces are found.
        check_bond_distances : boolean (True)
            Whether to check for large bond distances.
        check_large_forces : boolean (False)
            Whether to print force summary of force groups
        force_threshold : float (10.0)
            Threshold to check for large forces.
        bond_threshold : float (0.5)
            Threshold to check for large bond distances.

        Returns
        -------
        None
        """

        if check_bond_distances:
            # Check for large bond_distances
            self.checkBondDistances(threshold=bond_threshold)

        # Add particles to system
        self.addParticles()

        # Add created forces into the system
        self.addSystemForces()

        # Create force group for each added force
        for i, name in enumerate(self.forceGroups):
            self.forceGroups[name].setForceGroup(i)

        if minimize:
            check_large_forces = True

        if check_large_forces:
            # Check for high forces in atoms and minimize the system if necessary
            self.checkLargeForces(minimize=minimize, threshold=force_threshold)

    def checkBondDistances(self, threshold: float = 0.5) -> None:
        """
        Searches for large bond distances for the atom pairs defined in
        the 'bonds' attribute. It raises an error when large bonds are found.

        Parameters
        ----------
        threshold : (float, default=0.5 nm)
            Threshold to check for large bond distances.

        Returns
        -------
        None
        """
        print('Checking large bonds ...')
        if isinstance(threshold, float):
            threshold = threshold * openmm.unit.nanometer

        for b in self.bonds:
            if self.bonds[b][0] >= threshold:
                print('Problem with distance between atoms: ' + b[0].name + ' and ' + b[1].name)
                r1 = b[0].residue.name + '_' + str(b[0].residue.id)
                if b[0].residue == b[1].residue:
                    print('of residue: ' + r1)
                else:
                    r2 = b[1].residue.name + '_' + str(b[1].residue.id)
                    print('of residues: ' + r1 + ' and ' + r2 + ', respectively.')
                raise ValueError('The bond distance between them ' + str(self.bonds[b][0]) +
                                 'nm is larger than ' + str(threshold) + ' nm. Please check your input structure.')
            # else:
        print(f'All bonds seem to be OK (less than threshold: {threshold})')
        print('')

    def checkLargeForces(self, minimize: bool = False, threshold: float = 10) -> None:
        """
        Prints the hps system energies of the input configuration of the
        system. It optionally checks for large forces acting upon all
        particles in the hps system and iteratively minimizes the system
        configuration until no forces larger than a threshold are found.

        Parameters
        ----------
        threshold : (float, default=10)
            Threshold to check for large forces.
        minimize : (bool, default= False)
            Whether to iteratively minimize the system until all forces are lower or equal to
            the threshold value.

        Returns
        -------
        None
        """

        # minimized = False
        print('__________________________________________________________________')
        print('Potential Energy from initial structure (input structure):')

        # Define test simulation to extract forces
        integrator = openmm.LangevinIntegrator(1 * openmm.unit.kelvin, 1 / openmm.unit.picosecond,
                                               0.0005 * openmm.unit.picoseconds)
        sim = openmm.app.Simulation(self.topology, self.system, integrator)
        sim.context.setPositions(self.positions)
        state = sim.context.getState(getForces=True, getEnergy=True)

        # Print initial state of the system
        print(f'The Potential Energy of the system is : {state.getPotentialEnergy()}')
        for i, n in enumerate(self.forceGroups):
            energy = sim.context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(
                openmm.unit.kilojoules_per_mole)
            print('The ' + n.replace('Force', 'Energy') + ' is: ' + str(energy) + ' kJ/mol')
        print('')

        if minimize:
            print('__________________________________________________________________')
            print('Perform energy minimization ...')
            # Find if there is an acting force larger than threshold
            # minimize the system until forces have converged
            forces = [np.linalg.norm([f[0]._value, f[1]._value, f[2]._value]) for f in state.getForces()]
            print(state.getForces())
            prev_force = None
            tolerance = 10

            while np.max(forces) > threshold:

                # Write atom with the largest force if not reported before
                if np.max(forces) != prev_force:
                    atom = self.atoms[np.argmax(forces)]
                    residue = atom.residue
                    print(f'Large force {np.max(forces):.3f} kJ/(mol nm) found in:')
                    print(f'Atom: {atom.index} {atom.name}')
                    print(f'Residue: {residue.name} {residue.index}')
                    print(forces[np.argmax(forces)])
                    print(f'Minimising system with energy tolerance of {tolerance:.1f} kJ/mol')
                    print('_______________________')
                    print('')

                sim.minimizeEnergy(tolerance=tolerance * openmm.unit.kilojoule / openmm.unit.mole)
                # minimized = True
                state = sim.context.getState(getForces=True)
                prev_force = np.max(forces)
                forces = [np.linalg.norm([f.x, f.y, f.z]) for f in state.getForces()]
                if tolerance > 1:
                    tolerance -= 1
                elif tolerance > 0.1:
                    tolerance -= 0.1
                elif tolerance == 0.1:
                    raise ValueError('The system could not be minimized at the requested convergence\n' +
                                     'Try to increase the force threshold value to achieve convergence.')

            print(f'All forces are less than {threshold:.2f} kJ/mol/nm')
            print('______________________')
            state = sim.context.getState(getPositions=True, getEnergy=True)
            print('Potential Energy After minimisation:')
            print(f'The Potential Energy of the system (after minimized) is : {state.getPotentialEnergy()}')
            for i, n in enumerate(self.forceGroups):
                energy = sim.context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(
                    openmm.unit.kilojoules_per_mole)
                print('The ' + n.replace('Force', 'Energy') + ' is: ' + str(energy) + ' kJ/mol')

            print('')
            self.positions = state.getPositions()
            print('Saving minimized positions')
            print('__________________________________________________________________')
            print('')

    def addParticles(self) -> None:
        """
        Add particles to the system OpenMM class instance.

        Add a particle to the system for each atom in it. The mass
        of each particle is set up with the values in the :code:`particles_mass`
        attribute.

        """

        # Set same mass for each atom
        if isinstance(self.particles_mass, float):
            for i in np.arange(len(self.atoms)):
                self.system.addParticle(self.particles_mass)

        # Set unique masses for each atom
        if isinstance(self.particles_mass, list):
            assert len(self.particles_mass) == len(self.atoms)
            for i in np.arange(len(self.particles_mass)):
                self.system.addParticle(self.particles_mass[i])

    def addSystemForces(self) -> None:
        """
        Add forces to the system OpenMM class instance. It also save
        names for the added forces to include them in the reporter class.

        Adds generated forces to the system, also adding
        a force group to the :code:`forceGroups` attribute dictionary.

        """

        if self.harmonicBondForce is not None:
            print("harmonicBondForce is not None and then Added it...")
            self.system.addForce(self.harmonicBondForce)
            self.forceGroups['Harmonic Bond Energy'] = self.harmonicBondForce

        if self.gaussianAngleForce is not None:
            print("gaussianAngleForce is not None and then Added it...")
            self.system.addForce(self.gaussianAngleForce)
            self.forceGroups['Gaussian Angle Energy'] = self.gaussianAngleForce

        if self.gaussianTorsionForce is not None:
            print("gaussianTorsionForce is not None and then Added it...")
            self.system.addForce(self.gaussianTorsionForce)
            self.forceGroups['Gaussian Torsion Energy'] = self.gaussianTorsionForce

        if self.yukawaForce is not None:
            print("yukawaForce is not None and then Added it...")
            self.system.addForce(self.yukawaForce)
            self.forceGroups['Yukawa Energy'] = self.yukawaForce

        if self.ashbaugh_HatchForce is not None:
            print("ashbaugh_HatchForce is not None and then Added it...")
            self.system.addForce(self.ashbaugh_HatchForce)
            self.forceGroups['PairWise Energy'] = self.ashbaugh_HatchForce

        if self.wang_Frenkel_Force is not None:
            print("Wang-Frenkel is not None and then Added it...")
            self.system.addForce(self.wang_Frenkel_Force)
            self.forceGroups['PairWire Energy'] = self.wang_Frenkel_Force

    def dumpStructure(self, output_file: str) -> None:
        """
        Writes a structure file of the system in its current state.

        Writes a PDB file containing the currently defined CG system atoms and its positions.

        Parameters
        ----------
        output_file : string
            name of the PDB output file.

        Returns
        -------
        None
        """

        self.structure.writeFile(self.topology, self.positions, file=open(output_file, 'w'))

    def dumpTopology(self, output_file: str) -> None:
        """
        Writes a topology file of the system in PSF format, this is used for visualization and post-analysis.

        Writes a file containing the current topology in the
        hpsOpenMM system. This file contains topology of system, used in visualization and analysis.

        Here, we used :code:`parmed` to load :code:`openMM topology` and :code:`openMM system` to create
        :code:`Structure` object in :code:`parmed`.
        Because parmed doesn't automatically recognize :code:`charge`, :code:`mass` of atoms by their name.
        We need to set :code:`charge`, :code:`mass` back to residues properties.

        Parameters
        ----------
        output_file : string [requires]
            name of the output PSF file.

        Returns
        -------
        None
        """

        top = pmd.openmm.load_topology(self.topology, self.system)
        for i, a in enumerate(top.atoms):
            a.mass, a.charge = self.particles_mass[i], self.particles_charge[i]

        """
        Parmed know exactly what chain is but it doesn't use to write psf. Instead, it writes segid (1-letter)
        to distinguish between segment, and PSF file does not has chain identifier.
        What we will do here is copy chain ID to segID so that Parmed can write a meaningful psf.
        """
        for r in top.residues:
            r.segid = r.chain
        top.save(f'{output_file}', overwrite=True)

    def dumpForceFieldData(self, output_file: str) -> None:
        """
        Writes to a file the parameters of the forcefield.

        Writes a file containing the current forcefield parameters in the
        CG system.

        Parameters
        ----------
        output_file : string [requires]
            name of the output file.

        Returns
        -------
        None
        """

        with open(output_file, 'w') as ff:

            ff.write('#### CG Force Field Parameters ####\n')
            ff.write('\n')
            if self.atoms != OrderedDict():
                ff.write('[atoms]\n')
                ff.write(
                    '# %2s %3s %9s %9s %9s \t %14s\n' % ('atom', 'mass', 'exc_radius', 'charge', 'hps', 'atom_name'))

                for i, atom in enumerate(self.atoms):

                    if isinstance(self.particles_mass, list):
                        mass = self.particles_mass[i]
                    elif isinstance(self.particles_mass, float):
                        mass = self.particles_mass
                    if isinstance(self.rf_sigma, list):
                        sigma = self.rf_sigma[i]
                    elif isinstance(self.rf_sigma, float):
                        sigma = self.rf_sigma
                    if isinstance(self.particles_charge, list):
                        charge = self.particles_charge[i]
                    elif isinstance(self.particles_charge, float):
                        charge = self.particles_charge
                    if isinstance(self.particles_hps, list):
                        hps = self.particles_hps[i]
                    elif isinstance(self.particles_hps, float):
                        hps = self.particles_hps
                    res = atom.residue

                    ff.write('%4s %5s %9.3f %9.3f %9.3f\t# %12s\n' % (atom.index + 1,
                                                                      mass,
                                                                      sigma,
                                                                      charge,
                                                                      hps,
                                                                      atom.name + '_' + res.name + '_' + str(
                                                                          res.index + 1)))
            if self.bonds != OrderedDict():
                ff.write('\n')
                ff.write('[bonds]\n')
                ff.write('# %2s %3s %-6s %-s \t\t self.bonds[bond][1] %12s %12s\n' % (
                    'at1', 'at2', 'b0', 'k', 'at1_name', 'at2_name'))
                for bond in self.bonds:
                    atom1 = bond[0]
                    atom2 = bond[1]
                    res1 = atom1.residue
                    res2 = atom2.residue
                    ff.write('%4s %4s %.4f %s \t# %12s %12s\n' % (atom1.index + 1,
                                                                  atom2.index + 1,
                                                                  self.bonds[bond][0]._value,
                                                                  self.bonds[bond][1],
                                                                  atom1.name + '_' + res1.name + '_' + str(
                                                                      res1.index + 1),
                                                                  atom2.name + '_' + res2.name + '_' + str(
                                                                      res2.index + 1)))

    def setCAMassPerResidueType(self):
        """
        Sets alpha carbon atoms to their average residue mass. Used specially for
        modifying alpha-carbon (CA) coarse-grained models.

        Sets the masses of the alpha carbon atoms to the average mass
        of its amino acid residue.

        Parameters
        ----------

        Returns
        -------
        None
        """
        # Load mass parameters from parameters package
        params = model_parameters.parameters[self.model]
        masses = []
        for r in self.topology.residues():
            if r.name in params:
                masses.append(params[r.name]['mass'])
            else:
                raise ValueError('Residue ' + r.name + ' not found in masses dictionary.')

        self.setParticlesMass(masses)

    def setCARadiusPerResidueType(self):
        """
        Sets alpha carbon atoms to their average residue mass. Used specially for
        modifying alpha-carbon (CA) coarse-grained models.

        Sets the excluded volume radii of the alpha carbon atoms
        to characteristic radii of their corresponding amino acid
        residue.

        Parameters
        ----------

        Returns
        -------
        None
        """

        # Load radii from parameters package
        params = model_parameters.parameters[self.model]

        radii = []

        for r in self.topology.residues():
            if r.name in params:
                radii.append(params[r.name]['radii'])
            else:
                raise ValueError('Residue ' + r.name + ' not found in radii dictionary.')

        self.setParticlesRadii(radii)

    def setCAChargePerResidueType(self):
        """
        Sets the charge of the alpha carbon atoms
        to characteristic charge of their corresponding amino acid
        residue.

        Parameters
        ----------

        Returns
        -------
        None
        """

        # Load charge from parameters package
        params = model_parameters.parameters[self.model]
        charge = []

        for r in self.topology.residues():
            if r.name in params:
                charge.append(params[r.name]['charge'])
            else:
                raise ValueError('Residue ' + r.name + ' not found in charge dictionary.')

        self.setParticlesCharge(charge)

    def setCAHPSPerResidueType(self):
        """
        Sets alpha carbon atoms to their residue hydropathy scale. Used specially for
        modifying alpha-carbon (CA) coarse-grained models.

        Sets the HPS model of the alpha carbon atoms using corresponding scale.

        Parameters
        ----------

        Returns
        -------
        None
        """

        # Load hydropathy scale from parameters package
        params = model_parameters.parameters[self.model]

        hps = []

        for r in self.topology.residues():
            if r.name in params:
                hps.append(params[r.name]['hps'])
            else:
                raise ValueError('Residue ' + r.name + ' not found in hps dictionary.')

        self.setParticlesHPS(hps)

    def setCAIDPerResidueType(self):
        params = model_parameters.parameters['mpipi']
        atom_type = []

        for r in self.topology.residues():
            if r.name in params:
                atom_type.append(params[r.name]['id'])
            else:
                raise ValueError(f'Residue {r.name} not found in model parameter')

        self.setParticleTypeID(atom_type)

    # User-hidden functions #
    @staticmethod
    def _setParameters(term, parameters):
        """
        General function to set up or change force field parameters.
        protected method, can be called only inside class system.

        Parameters
        ----------
        term : dict
            Dictionary object containing the set of degrees of freedom
            (DOF) to set up attributes to (e.g. :code:`bonds` attribute)

        parameters : integer or float or list
            Value(s) for the specific forcefield parameters. If integer
            or float, sets up the same value for all the DOF in terms.
            If a list is given, sets a unique parameter for each DOF.

        Returns
        -------
        None
        """

        if isinstance(parameters, int):
            parameters = float(parameters)

        # Set constant parameter for each item in FF term
        if isinstance(parameters, float):
            for item in term:
                """
                For example, use in set bond force constant:
                in the :code:`getBond()`, we set :code:`bond_length` to the first item of tuple, 
                so lack of bond force constant.
                In :code:`setBondForceConstants(8368.0)`: parameter passed is 8368.0, 
                :code:`setBondForceConstants(param)` calls to this function. Term is self.bonds, item looks like: 
                :code:`item: (<Atom 0 (CA) of chain 0 residue 0 (MET)>, <Atom 1 (CA) of chain 0 residue 1 (ASP)>)`
                :code:`term[item]` (before calling this function): :code:`(Quantity(value=0.382, unit=nanometer), None)`
                we take all values except the last (bond force constant which is currently :code:`None` = no parameters:
                :code:`(Quantity(value=0.382, unit=nanometer),)`: a tuple
                after calling this function:
                :code:`term[item]: (Quantity(value=0.382, unit=nanometer), 8368.0)`
                :code:`[tuple+tuple = tuple]`
                """
                term[item] = term[item][:-1] + (parameters,)

        # Set unique parameter for each item in FF term
        if isinstance(parameters, list):
            assert len(parameters) == len(list(term.keys()))
            for i, item in enumerate(term):
                term[item] = term[item][:-1] + (parameters[i],)


"""
End of class system. Happy simulating ...
"""
