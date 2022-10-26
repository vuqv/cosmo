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
        hps_scale: 'hps_kr','hps_urry' [optional, default='hps_kr']
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
    atoms : :code:`list`
        A list of the current atoms in the model. The items are :code:`simtk.openmm.app.topology.Atom`
        initialised classes.
    n_atoms : :code:`int`
        Total numer of atoms in the model.
    bonds : :code:`collections.OrderedDict`
        A dict that uses bonds (2-tuple of :code:`simtk.openmm.app.topology.Atom` objects)
        present in the model as keys and their forcefield properties as values.
    bonds_indexes : :code:`list`
        A list containing the zero-based indexes of the atoms defining the bonds in the model.
    n_bonds : :code:`int`
        Total number of bonds in the model.
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
    rf_sigma : :code:`float`
        Sigma parameter used in the pairwise force object.
        This is vdw Radius of beads

    Methods
    -------


    loadForcefieldFromFile()
        Loads forcefield parameters from a force field file written with
        the :code:`dumpForceFieldData()` method.

    """

    def __init__(self, structure_path, hps_scale):
        """
        Initialises the hps OpenMM system class.

        Parameters
        ----------
        structure_path : string [requires]
            Name of the input PDB or CIF file
        hps_scale: 'hps_kr','hps_urry' [optional, default='hps_kr']
            Hydropathy scale. Currently, there are two models are supported.
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
        self.particles_mass = None
        # Define geometric attributes
        self.atoms = []
        self.n_atoms = None

        self.bonds = OrderedDict()
        self.bonds_indexes = []
        self.n_bonds = None

        self.angles = OrderedDict()
        self.angles_indexes = []
        self.n_angles = None
        self.gaussianAngleForce = None

        self.torsions = OrderedDict()
        self.torsions_indexes = []
        self.n_torsions = None
        self.gaussianTorsionForce = None

        self.hps_scale = hps_scale
        self.bond_length = model_parameters.parameters[hps_scale]["bond_length"]
        self.bondedTo = None
        self.bonded_exclusions_index = model_parameters.parameters[hps_scale]["bonded_exclusions_index"]

        # Define force attributes
        self.harmonicBondForce = None

        self.rf_sigma = None

        # PairWise potential
        self.ashbaugh_HatchForce = None
        self.epsilon = 0.8368 * openmm.unit.kilojoule_per_mole
        self.cutoff_Ashbaugh_Hatch = 2.0 * openmm.unit.nanometer
        self.particles_hps = None

        # Define parameter for DH potential
        self.yukawaForce = None
        self.particles_charge = None
        self.lD = 1.0 * openmm.unit.nanometer
        self.electric_factor = 138.935458 * openmm.unit.kilojoule_per_mole * openmm.unit.nanometer / openmm.unit.elementary_charge ** 2
        self.yukawa_cutoff = 3.5 * openmm.unit.nanometer
        self.epsilon_r = 80.0

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
        for i in range(1, len(atoms)):
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
        params = model_parameters.parameters[self.hps_scale]
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

    def setParticlesMasses(self, particles_mass):
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

        energy_function = 'factor*charge1*charge2/epsilon_r/r*exp(-r/lD)'
        self.yukawaForce = openmm.CustomNonbondedForce(energy_function)
        self.yukawaForce.addGlobalParameter('factor', self.electric_factor)
        self.yukawaForce.addGlobalParameter('epsilon_r', self.epsilon_r)
        self.yukawaForce.addGlobalParameter('lD', self.lD)
        self.yukawaForce.addPerParticleParameter('charge')
        if use_pbc:
            self.yukawaForce.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        else:
            self.yukawaForce.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)

        self.yukawaForce.setCutoffDistance(self.yukawa_cutoff)

        if isinstance(self.particles_charge, float):
            for i in range(len(self.atoms)):
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

        energy_function = 'step(2^(1/6)*sigma - r) *'
        energy_function += '(4*epsilon* ((sigma/r)^12-(sigma/r)^6) + (1-hps)*epsilon )'
        energy_function += '+(1-step(2^(1/6)*sigma-r)) * (hps*4*epsilon*((sigma/r)^12-(sigma/r)^6));'
        energy_function += 'sigma=0.5*(sigma1+sigma2);'
        energy_function += 'hps=0.5*(hps1+hps2)'
        self.ashbaugh_HatchForce = openmm.CustomNonbondedForce(energy_function)
        self.ashbaugh_HatchForce.addGlobalParameter('epsilon', self.epsilon)
        self.ashbaugh_HatchForce.addPerParticleParameter('sigma')
        self.ashbaugh_HatchForce.addPerParticleParameter('hps')
        #
        if use_pbc:
            self.ashbaugh_HatchForce.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        else:
            self.ashbaugh_HatchForce.setNonbondedMethod(openmm.NonbondedForce.CutoffNonPeriodic)

        self.ashbaugh_HatchForce.setCutoffDistance(self.cutoff_Ashbaugh_Hatch)

        if isinstance(self.rf_sigma, float):
            for i in range(len(self.atoms)):
                self.ashbaugh_HatchForce.addParticle((self.rf_sigma,))

        # in the case each atom has different sigma para.
        elif isinstance(self.rf_sigma, list):
            assert self.n_atoms == len(self.rf_sigma) and self.n_atoms == len(self.particles_hps)
            for i, atom in enumerate(self.atoms):
                self.ashbaugh_HatchForce.addParticle((self.rf_sigma[i], self.particles_hps[i],))

        # set exclusions rule
        bonded_exclusions = [(b[0].index, b[1].index) for b in list(self.topology.bonds())]
        self.ashbaugh_HatchForce.createExclusionsFromBonds(bonded_exclusions, self.bonded_exclusions_index)

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
        print('_________________________________')
        print('Energy from initial structure (input structure):')

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
            print('The ' + n.replace('Force', 'Energy') + ' is: ' + str(energy) + ' kj/mol')
        print('')

        if minimize:
            print('_________________________________')
            print('Perform energy minimization ...')
            # Find if there is an acting force larger than threshold
            # minimize the system until forces have converged
            forces = [np.linalg.norm([f[0]._value, f[1]._value, f[2]._value]) for f in state.getForces()]
            prev_force = None
            tolerance = 10

            while np.max(forces) > threshold:

                # Write atom with the largest force if not reported before
                if np.max(forces) != prev_force:
                    atom = self.atoms[np.argmax(forces)]
                    residue = atom.residue
                    print('Large force %.3f kj/(mol nm) found in:' % np.max(forces))
                    print(f'Atom: {atom.index} {atom.name}')
                    print(f'Residue: {residue.name} {residue.index}')
                    print('Minimising system with energy tolerance of %.1f kj/mol' % tolerance)
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

            print(f'All forces are less than {threshold:.2f} kj/mol/nm')
            print('______________________')
            state = sim.context.getState(getPositions=True, getEnergy=True)
            print('After minimisation:')
            print(f'The Potential Energy of the system (after minimized) is : {state.getPotentialEnergy()}')
            for i, n in enumerate(self.forceGroups):
                energy = sim.context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(
                    openmm.unit.kilojoules_per_mole)
                print('The ' + n.replace('Force', 'Energy') + ' is: ' + str(energy) + ' kj/mol')

            print('')
            self.positions = state.getPositions()
            print('Saving minimized positions')
            print('______________________')
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
            for i in range(len(self.atoms)):
                self.system.addParticle(self.particles_mass)

        # Set unique masses for each atom
        if isinstance(self.particles_mass, list):
            assert len(self.particles_mass) == len(self.atoms)
            for i in range(len(self.particles_mass)):
                self.system.addParticle(self.particles_mass[i])

    def addSystemForces(self) -> None:
        """
        Add forces to the system OpenMM class instance. It also save
        names for the added forces to include them in the reporter class.

        Adds generated forces to the system, also adding
        a force group to the :code:`forceGroups` attribute dictionary.

        """

        if self.harmonicBondForce is not None:
            self.system.addForce(self.harmonicBondForce)
            self.forceGroups['Harmonic Bond Energy'] = self.harmonicBondForce

        if self.gaussianAngleForce is not None:
            self.system.addForce(self.gaussianAngleForce)
            self.forceGroups['Gaussian Angle Energy'] = self.gaussianAngleForce

        if self.gaussianTorsionForce is not None:
            self.system.addForce(self.gaussianTorsionForce)
            self.forceGroups['Gaussian Torsion Energy'] = self.gaussianTorsionForce

        if self.yukawaForce is not None:
            self.system.addForce(self.yukawaForce)
            self.forceGroups['Yukawa Energy'] = self.yukawaForce

        if self.ashbaugh_HatchForce is not None:
            self.system.addForce(self.ashbaugh_HatchForce)
            self.forceGroups['PairWise Energy'] = self.ashbaugh_HatchForce

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
        params = model_parameters.parameters[self.hps_scale]
        masses = []
        for r in self.topology.residues():
            if r.name in params:
                masses.append(params[r.name]['mass'])
            else:
                raise ValueError('Residue ' + r.name + ' not found in masses dictionary.')

        self.setParticlesMasses(masses)

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
        params = model_parameters.parameters[self.hps_scale]

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
        params = model_parameters.parameters[self.hps_scale]
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
        params = model_parameters.parameters[self.hps_scale]

        hps = []

        for r in self.topology.residues():
            if r.name in params:
                hps.append(params[r.name]['hps'])
            else:
                raise ValueError('Residue ' + r.name + ' not found in hps dictionary.')

        self.setParticlesHPS(hps)

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
