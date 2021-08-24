#!/usr/bin/env python
# coding: utf-8

import re
from collections import OrderedDict

import numpy as np
import parmed as pmd
from simtk import unit
from simtk.openmm import *
from simtk.openmm.app import *

from ..parameters import ca_parameters


class system:
    """
    A class containing methods and parameters for generating CG systems to be simulated using the OpenMM interface.
    It offers flexibility to create default and custom CG systems and to easily modify their parameters.

    Attributes
    ----------
    structure_path : :code:`string`
        Path to the pdb or cif input file.
    structure : :code:`openmm.app.pdbfile.PDBFile or openmm.app.pdbxfile.PDBxFile`
        Object that holds the information of OpenMM PDB or CIF parsing methods.
    topology : :code:`openmm.app.topology.Topology`
        OpenMM topology of the model.
    positions : :code:`unit.quantity.Quantity`
        Atomic positions of the model.
    particles_mass : :code:`float or list`
        Mass of each particle. If float then uniform masses are given to all
        particles. If list per-particle masses are assigned.
    model_type : :code:`string`
        String representing the model type: All-atom (AA), alpha-carbon (CA)
        or multi-basin variants (AA-MB, CA-MB).
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
        an harmonic bond potential between pairs of particles, that depends
        quadratically on their distance.
    yukawaForce : :code:`openmm.CustomNonbondedForce`
        Stores the OpenMM :code:`CustomNonbondedForce` initialized-class.
        Implements the Debye-Huckle potential.
    pairWiseForce : :code:`openmm.CustomNonbondedForce`
        Stores the OpenMM :code:`CustomNonbondedForce` initialized-class. Implements the pairwise short-range
        potential.
    forceGroups : :code:`collections.OrderedDict`
        A dict that uses force names as keys and their corresponding force
        as values.
    system : :code:`openmm.System`
        Stores the OpenMM System initialised class. It stores all the forcefield
        information for the SBM model.
    rf_sigma : :code:`float`
        Sigma parameter used in the pairwise force object.
        This is vdw Radius of beads
    exclusion_NB : Matrix of size NxN
        This is a pairwise interaction matrix.
        initialize by :code:`np.ones((N, N))`, N is the number of beads in the system.
        By default, all atoms interact via pairwise. but this is not true, we should
        exclude atom pairs that covalently bonded by when adding Harmonic bonds, we set matrix element of :code:`(i,j)=
        (j,i)=0`, where i and j are atom indices in bond.

    Methods
    -------
    removeHydrogens()
        Remove hydrogens from the input structure by using a regexpression pattern.
        Used specially for creating all atom (AA) models.
    getCAlphaOnly()
        Filter in only alpha carbon atoms from the input structure and updates
        the topology object to add new bonds between them. Used specially for
        creating alpha-carbon (CA) coarse-grained models.
    getAtoms()
        Reads atoms from topology, adds them to the main class and sorts them
        into a dictionary to store their forcefield properties.
    getBonds()
        Reads bonds from topology, adds them to the main class and sorts them
        into a dictionary to store their forcefield properties.
    setBondParameters()
        Change the forcefield parameters for bonded terms.
    setParticlesMasses()
        Change the mass parameter for each atom in the system.
    setParticlesRadii()
        Change the excluded volume radius parameter for each atom in the system.
    addHarmonicBondForces()
        Creates an harmonic bonded force term for each bond in the main
        class using their defined forcefield parameters.
    addYukawaForces()
        Creates a nonbonded force term for electrostatic interaction DH potential.
    addPairWiseForces()
        Creates a nonbonded force term for pairwise interaction (customize LJ 12-6 potential).
    createSystemObject()
        Creates OpenMM system object adding particles, masses and forces.
        It also groups the added forces into Force-Groups for the sbmReporter
        class.
    addParticles()
        Add particles to the system OpenMM class instance.
    addSystemForces()
        Add forces to the system OpenMM class instance. It also save
        names for the added forces to include them in the reporter class.
    dumpStructure()
        Writes a structure file of the system in its current state.
    dumpTopology()
        Writes a topology file of the system in PSF format, this is used for visualization and post-analysis.
    dumpForceFieldData()
        Writes to a file the parameters of the forcefield.
    loadForcefieldFromFile()
        Loads forcefield parameters from a force field file written with
        the :code:`dumpForceFieldData()` method.
    setCAMassPerResidueType()
        Sets alpha carbon atoms to their average residue mass. Used specially for
        modifying alpha-carbon (CA) coarse-grained models.
    setCARadiusPerResidueType()
        Sets alpha carbon atoms to their average residue mass. Used specially for
        modifying alpha-carbon (CA) coarse-grained models.
    setCAHPSPerResidueType()
        Sets alpha carbon atoms to their residue Urry Hydropathy scale. Used specially for
        modifying alpha-carbon (CA) coarse-grained models.
    """

    def __init__(self, structure_path, particles_mass=1.0):
        """
        Initialises the SBM OpenMM system class.

        Parameters
        ----------
        structure_path : string
            Name of the input PDB or CIF file
        particles_mass : float or list
            mass of all the particles in the system.

        Returns
        -------
        None
        """

        # Define structure object attributes
        self.structure_path = structure_path
        # Recognize format of input structure file
        if structure_path.endswith('.pdb'):
            self.structure = PDBFile(structure_path)
        elif structure_path.endswith('.cif'):
            self.structure = pdbxfile.PDBxFile(structure_path)
        else:
            raise ValueError(
                'Structure file extension not recognized. It must end with .pdb or .cif accordingly.')
        self.topology = self.structure.topology
        self.positions = self.structure.positions
        self.particles_mass = particles_mass
        self.model_type = None

        # Define geometric attributes
        self.atoms = []
        self.n_atoms = None
        self.bonds = OrderedDict()
        self.bonds_indexes = []
        self.n_bonds = None
        self.bond_length = 0.382

        # Define force attributes
        self.harmonicBondForce = None
        self.rf_sigma = None
        # self.rf_epsilon = None
        self.rf_cutoff = None
        self.exclusion_NB = None
        self.particles_hps = None

        # self.exclusions = []
        # PairWise potential
        self.pairWiseForce = None
        self.muy = 1
        self.delta = 0.08
        self.epsilon = 0.8368

        # Define parameter for DH potential
        self.yukawaForce = None
        self.particles_charge = None
        self.lD = 1.0 * unit.nanometer
        self.electric_factor = 138.935458
        self.yukawa_cutoff = 4.0 * self.lD
        self.epsilon_r = 80.0

        self.forceGroups = OrderedDict()

        # Initialise an OpenMM system class instance
        self.system = openmm.System()

    def removeHydrogens(self, except_chains=None):
        """
        Removes all hydrogen atoms in the topology from the SBM system.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if isinstance(except_chains, str):
            except_chains = list(except_chains)

        # save all hydrogen atoms
        atomsToRemove = []
        _hydrogen = re.compile("[123 ]*H.*")
        for a in self.topology.atoms():
            if except_chains is not None:
                if a.residue.chain.id not in except_chains:
                    if _hydrogen.match(a.name):
                        atomsToRemove.append(a)
            else:
                if _hydrogen.match(a.name):
                    atomsToRemove.append(a)

        # Remove all hydrogen atoms
        modeller_topology = modeller.Modeller(self.topology, self.positions)
        modeller_topology.delete(atomsToRemove)
        self.topology = modeller_topology.getTopology()
        self.positions = modeller_topology.getPositions()

    def getCAlphaOnly(self):
        """
        Keeps in the SBM system only the alpha carbon atoms from the OpenMM topology.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # save all non C-alpha atoms
        atomsToRemove = []
        oldIndex = []
        for a in self.topology.atoms():
            if a.name != 'CA':
                atomsToRemove.append(a)

        # Remove all non C-alpha atoms
        modeller_topology = modeller.Modeller(self.topology, self.positions)
        modeller_topology.delete(atomsToRemove)
        self.topology = modeller_topology.getTopology()
        self.positions = modeller_topology.getPositions()

        # Update system atoms
        atoms = list(self.topology.atoms())

        # Add bonds between C-alpha atoms of the same chain
        for i in range(1, len(atoms)):
            current_chain = atoms[i].residue.chain
            if i == 1:
                previous_chain = current_chain

            # Add bond only when atoms belong to the same chain
            if current_chain == previous_chain:
                self.topology.addBond(atoms[i - 1], atoms[i])
            previous_chain = current_chain

        self.model_type = 'CA'

    def getAtoms(self):
        """
        Adds atoms in the OpenMM topology instance to the sbmOpenMM system class.

        Parameters
        ----------
        None

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

        # Add atoms to sbm object
        self.n_atoms = 0
        for atom in atoms:
            self.atoms.append(atom)
            self.n_atoms += 1

        """
        Prepare for pair-wire excluded interaction. If two atoms covalently bonded then is not appears in pair-wire
        nonbonded interactions (still present in electrostatic interactions)
        
        """
        self.exclusion_NB = np.ones((self.n_atoms, self.n_atoms))

    def getBonds(self, except_chains=None):
        """
        Adds bonds in the OpenMM topology instance to the sbmOpenMM system class.

        Parameters
        ----------
        None

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

        # Add bonds to sbm object
        self.n_bonds = 0
        for bond in bonds:
            p1 = self.positions[bond[0].index]
            p2 = self.positions[bond[1].index]
            # bond_length = geometry.bond(p1, p2)
            bond_length = self.bond_length * unit.nanometer
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

    """ Functions for setting force specific parameters """

    def setBondParameters(self, bond_parameters):
        """
        Set the harmonic bond constant force parameters. The input can be
        a float, to set the same parameter for all force interactions, or
        a list, to define a unique parameter for each force interaction.

        Parameters
        ----------
        bond_parameters : float or list
            Parameter(s) to set up for the harmonic bond forces.

        Returns
        -------
        None
        """

        system._setParameters(self.bonds, bond_parameters)

    def setParticlesMasses(self, particles_mass):
        """
        Set the masses of the particles in the system. The input can be a
        float, to set the same mass for all particles, or a list, to define
        a unique mass for each particle.

        Parameters
        ----------
        particles_mass : float or list
            Mass(es) values to add for the particles in the sbmOpenMM system class.

        Returns
        -------
        None
        """

        self.particles_mass = particles_mass

    def setParticlesRadii(self, particles_radii):
        """
        Set the radii of the particles in the system. The input can be a
        float, to set the same radius for all particles, or a list, to define
        a unique radius for each particle.

        Parameters
        ----------
        particles_radii : float or list
            Radii values to add for the particles in the sbmOpenMM system class.

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
            Charge values to add for the particles in the sbmOpenMM system class.

        Returns
        -------
        None
        """

        self.particles_charge = particles_charge

    def setParticlesHPS(self, particles_hps):
        """
        Set the charge of the particles in the system. The input can be a
        float, to set the same charge for all particles, or a list, to define
        a unique charge for each particle.

        Parameters
        ----------
        particles_hps : float or list
            HPS scale values to add for the particles in the sbmOpenMM system class.

        Returns
        -------
        None
        """

        self.particles_hps = particles_hps

    """ Functions for creating force objects with defined parameters """

    def addHarmonicBondForces(self):
        """
        Creates an :code:`openmm.HarmonicBondForce()` object with the bonds and
        parameters set up in the "bonds" dictionary attribute. The force object
        is stored at the :code:`harmonicBondForce` attribute.

        The force parameters must be contained in self.bonds as follows:
        self.bonds is a dictionary:

            - The keys are 2-tuples for two atom items in :code:`self.topology.atoms` attribute.
            - The values are a 2-tuple of parameters in the following order:

                - first  -> bond0 (quantity)
                - second -> k (float)

        Parameters
        ----------
        None

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
            # two atoms in bond will have nonbonded interaction matrix is 0 - do not interactions via nonbonded
            self.exclusion_NB[bond[0].index][bond[1].index] = 0.0
            self.exclusion_NB[bond[1].index][bond[0].index] = 0.0

    def addYukawaForces(self):
        """
        Creates an :code:`openmm.CustomNonbondedForce()` object with the parameters
        sigma and epsilon given to this method. The custom non-bonded force
        is initialized with the formula:

        .. math::
            energy = f \\times \\frac{q_1q_2}{\epsilon_r \\times r}\\times e^{(-r/lD)}


        where :math:`f=\\frac{1}{4\\pi\\epsilon_0}=138.935458` is the factor for short to convert dimensionless
        in calculation to :math:`kj.nm/(mol\\times e^2)` unit.

        The force object is stored at the :code:`yukawaForce` attribute.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        energy_function = 'factor*charge1*charge2/epsilon_r/r*exp(-r/lD)'
        self.yukawaForce = CustomNonbondedForce(energy_function)
        self.yukawaForce.addGlobalParameter('factor', self.electric_factor)
        self.yukawaForce.addGlobalParameter('epsilon_r', self.epsilon_r)
        self.yukawaForce.addGlobalParameter('lD', self.lD)
        self.yukawaForce.addPerParticleParameter('charge')

        if isinstance(self.particles_charge, float):
            for atom in self.atoms:
                self.yukawaForce.addParticle((self.particles_charge,))

        # in the case each atoms have different sigma para.
        elif isinstance(self.particles_charge, list):
            assert self.n_atoms == len(self.particles_charge)
            for i, atom in enumerate(self.atoms):
                self.yukawaForce.addParticle((self.particles_charge[i],))

        # self.yukawaForce.setCutoffDistance(self.yukawa_cutoff)

    def addPairWiseForces(self):
        """
        Creates an :code:`openmm.CustomNonbondedForce()` object with the parameters
        sigma and epsilon given to this method. The custom non-bonded force
        is initialized with the formula: (note: hps here is :math:`\lambda_{ij}^{0}` in the paper)

        Unlike :code:`BondForce` class, where we specify index for atoms pair to add bond, it means
        that number of bondForces may differ from number of particle.
        :code:`NonBondedForce` is add to all particles, hence we don't need to pass the :code:`atom index`.

        .. math::
            \\Phi_{i,j}^{vdw}(r) = step(2^{1/6}\\sigma_{ij}-r) \\times
            \\left( 4\\epsilon\\left[\\left(\\frac{\\sigma_{ij}}{r}\\right)^{12}-
            \\left(\\frac{\\sigma_{ij}}{r}\\right)^{6}\\right]+(1-\\mu\\times\\lambda_{ij}^{0}+\\Delta)\\epsilon\\right)

            + \\left[1-step(2^{1/6}\\sigma_{ij}-r)\\right]\\times\\left[(\\mu \\lambda_{ij}^{0}-\\Delta)\\times 4\\pi
            \\left[\\left(\\frac{\\sigma_{ij}}{r}\\right)^{12}-\\left(\\frac{\\sigma_{ij}}{r}\\right)^6\\right]\\right]



        Here, :math:`\\sigma= \\frac{(\\sigma_1+\\sigma_2)}{2}; \\lambda_{ij}^{0}=\\frac{(\\lambda_i+\\lambda_j)}{2};
        \\mu= 1;\\Delta= 0.08; \\epsilon = 0.8368 kj/mol`

        The
        The force object is stored at the :code:`pairWiseForce` attribute.

        Parameters
        ----------
        epsilon : float
            Value of the epsilon constant in the energy function.
        sigma : float or list
            Value of the sigma constant (in nm) in the energy function. If float the
            same sigma value is used for every particle. If list a unique
            parameter is given for each particle.
        cutoff : float
            The cutoff distance (in nm) being used for the nonbonded interactions.
        Returns
        -------
        None
        """
        nbMatrix_LST = self.exclusion_NB.ravel().tolist()

        energy_function = 'nb*(step(2^(1/6)*sigma - r) *'
        energy_function += '(4*epsilon* ((sigma/r)^12-(sigma/r)^6) + (1-muy*hps+delta)*epsilon )'
        energy_function += '+(1-step(2^(1/6)-r)) * ((muy*hps-delta)*4*epsilon*((sigma/r)^12-(sigma/r)^6)));'
        energy_function += 'nb=nb_matrix(idx1, idx2);'
        energy_function += 'sigma=0.5*(sigma1+sigma2);'
        energy_function += 'hps=0.5*(hps1+hps2)'
        self.pairWiseForce = CustomNonbondedForce(energy_function)
        self.pairWiseForce.addTabulatedFunction('nb_matrix',
                                                Discrete2DFunction(self.n_atoms, self.n_atoms, nbMatrix_LST))
        self.pairWiseForce.addGlobalParameter('muy', self.muy)
        self.pairWiseForce.addGlobalParameter('delta', self.delta)
        self.pairWiseForce.addGlobalParameter('epsilon', self.epsilon)
        self.pairWiseForce.addPerParticleParameter('idx')
        self.pairWiseForce.addPerParticleParameter('sigma')
        self.pairWiseForce.addPerParticleParameter('hps')

        if isinstance(self.rf_sigma, float):
            for atom in self.atoms:
                self.pairWiseForce.addParticle((self.rf_sigma,))

        # in the case each atoms have different sigma para.
        elif isinstance(self.rf_sigma, list):
            assert self.n_atoms == len(self.rf_sigma) and self.n_atoms == len(self.particles_hps)
            for i, atom in enumerate(self.atoms):
                self.pairWiseForce.addParticle((atom.index, self.rf_sigma[i], self.particles_hps[i],))
                """Or can be add as follow"""
                # self.pairWiseForce.addParticle([atom.index, self.rf_sigma[i], self.particles_hps[i]])
        """set cut off for nonbonded interaction"""
        # self.pairWiseForce.setCutoffDistance(self.rf_cutoff)

    """ Functions for creating OpenMM system object """

    def createSystemObject(self, check_bond_distances=True, minimize=False, check_large_forces=True,
                           force_threshold=10.0, bond_threshold=0.39):
        """
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
        bond_threshold : float (0.39)
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

    def checkBondDistances(self, threshold=0.39):
        """
        Searches for large bond distances for the atom pairs defined in
        the 'bonds' attribute. It raises an error when large bonds are found.

        Parameters
        ----------
        threshold : float
            Threshold to check for large bond distances.

        Returns
        -------
        None
        """
        print('checking large bonds ...')
        if isinstance(threshold, float):
            threshold = threshold * unit.nanometer

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

    def checkLargeForces(self, threshold=1, minimize=False):
        """
        Prints the SBM system energies of the input configuration of the
        system. It optionally checks for large forces acting upon all
        particles in the SBM system and iteratively minimizes the system
        configuration until no forces larger than a threshold are found.

        Parameters
        ----------
        threshold : float
            Threshold to check for large forces.
        minimize : float
            Whether to iteratively minimize the system until all forces are lower or equal to
            the threshold value.

        Returns
        -------
        None
        """

        # minimized = False
        print('_________________________________')
        print('Energy minimization:')

        # Define test simulation to extract forces
        integrator = LangevinIntegrator(1 * unit.kelvin, 1 / unit.picosecond, 0.0005 * unit.picoseconds)
        simulation = Simulation(self.topology, self.system, integrator)
        simulation.context.setPositions(self.positions)
        state = simulation.context.getState(getForces=True, getEnergy=True)

        # Print initial state of the system
        print('The Potential Energy of the system is : %s' % state.getPotentialEnergy())
        for i, n in enumerate(self.forceGroups):
            energy = simulation.context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(
                unit.kilojoules_per_mole)
            print('The ' + n.replace('Force', 'Energy') + ' is: ' + str(energy) + ' kj/mol')
        print('')

        if minimize:
            # Find if there is an acting force larger than threshold
            # minimize the system until forces have converged
            forces = [np.linalg.norm([f[0]._value, f[1]._value, f[2]._value]) for f in state.getForces()]
            prev_force = None
            tolerance = 10

            while np.max(forces) > threshold:

                # Write atom with largest force if not reported before
                if np.max(forces) != prev_force:
                    atom = self.atoms[np.argmax(forces)]
                    residue = atom.residue
                    print('Large force %.3f kj/(mol nm) found in:' % np.max(forces))
                    print(f'Atom: {atom.index} {atom.name}')
                    print(f'Residue: {residue.name} {residue.index}')
                    print('Minimising system with energy tolerance of %.1f kj/mol' % tolerance)
                    print('')

                simulation.minimizeEnergy(tolerance=tolerance * unit.kilojoule / unit.mole)
                # minimized = True
                state = simulation.context.getState(getForces=True)
                prev_force = np.max(forces)
                forces = [np.linalg.norm([f.x, f.y, f.z]) for f in state.getForces()]
                if tolerance > 1:
                    tolerance -= 1
                elif tolerance > 0.1:
                    tolerance -= 0.1
                elif tolerance == 0.1:
                    raise ValueError('The system could no be minimized at the requested convergence\n' +
                                     'Try to increase the force threshold value to achieve convergence.')

            print('______________________')
            state = simulation.context.getState(getPositions=True, getEnergy=True)
            print('After minimisation:')
            print('The Potential Energy of the system is : %s' % state.getPotentialEnergy())
            for i, n in enumerate(self.forceGroups):
                energy = simulation.context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(
                    unit.kilojoules_per_mole)
                print('The ' + n.replace('Force', 'Energy') + ' is: ' + str(energy) + ' kj/mol')
            print('All forces are less than %.2f kj/mol/nm' % threshold)
            print('')
            print('______________________')
            print('Saving minimized positions')
            print('')
            self.positions = state.getPositions()

    def addParticles(self):
        """
        Add a particle to the system for each atom in it. The mass
        of each particle is set up with the values in the :code:`particles_mass`
        attribute.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # Set same mass for each atom
        if isinstance(self.particles_mass, float):
            for atom in self.atoms:
                self.system.addParticle(self.particles_mass)

        # Set unique masses for each atom
        if isinstance(self.particles_mass, list):
            assert len(self.particles_mass) == len(self.atoms)
            for i in range(len(self.particles_mass)):
                self.system.addParticle(self.particles_mass[i])

    def addSystemForces(self):
        """
        Adds generated forces to the system, also adding
        a force group to the :code:`forceGroups` attribute dictionary.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        if self.harmonicBondForce is not None:
            self.system.addForce(self.harmonicBondForce)
            self.forceGroups['Harmonic Bond Energy'] = self.harmonicBondForce

        if self.yukawaForce is not None:
            self.system.addForce(self.yukawaForce)
            self.forceGroups['Yukawa Energy'] = self.yukawaForce

        if self.pairWiseForce is not None:
            self.system.addForce(self.pairWiseForce)
            self.forceGroups['PairWise Energy'] = self.pairWiseForce

    def dumpStructure(self, output_file):
        """
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

    def dumpTopology(self, output_file):
        """
        Writes a file containing the current topology in the
        sbmOpenMM system. This file contains topology of system, used in visualization and analysis.

        Here, we used :code:`parmed` to load openMM topology, openMM system to create Structure object in parmed.
        Because parmed automatically recognizes charge, mass of atoms by their name.
        We need to set charge, mass back to residues properties.

        Parameters
        ----------
        output_file : string
            name of the output PSF file.

        Returns
        -------
        None
        """

        top = pmd.openmm.load_topology(self.topology, self.system)
        for i, a in enumerate(top.atoms):
            a.mass, a.charge = self.particles_mass[i], self.particles_charge[i]
        top.save(f'{output_file}', overwrite=True)

    def dumpForceFieldData(self, output_file):
        """
        Writes a file containing the current forcefield parameters in the
        CG system.

        Parameters
        ----------
        output_file : string
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

    def loadForcefieldFromFile(self, forcefield_file):
        """
        Loads force field parameters from a force field file written by the
        :code:`dumpForceFieldData()` method into the sbmOpenMM system.

        NOTE: I still keep this function since I think I will use it later.
        I have not customized this function yet so do not use :code:`loadForcefieldFromFile()` function in simulation.
        Just keep here for prototype and customize later...

        Parameters
        ----------
        forcefield_file : string
            path to the force field file.

        Returns
        -------
        None
        """

        # Open forcefield file
        with open(forcefield_file, 'r') as ff:
            print('Reading Forcefield parameters from file ' + forcefield_file + ':')
            print('________________________________________' + '_' * len(forcefield_file))

            # Initialize ff-file section booleans to 0
            atoms = False
            bonds = False

            # Iterate for all lines in forcefield file
            for i, line in enumerate(ff):

                # Ignore comment and empty lines.
                if not line.startswith('#') and line.split() != []:

                    # Turn off all sections when a new is being reading.
                    if line.startswith('['):
                        atoms = False
                        bonds = False

                    else:

                        # Remove comments at the end of line
                        line = line[:line.find('#')]
                        ls = line.split()

                        # Reading [atoms] section
                        if atoms:
                            if not isinstance(self.particles_mass, list):
                                self.particles_mass = [1.0 for m in range(self.n_atoms)]
                            if not isinstance(self.rf_sigma, list):
                                self.rf_sigma = [0 for s in range(self.n_atoms)]
                            if len(ls) > 3:
                                raise ValueError('More than three parameters given in [atoms] section at line ' + str(
                                    i) + ':\n' + line)
                            if len(ls) < 3:
                                raise ValueError('Less than three parameters given in [atoms] section at line ' + str(
                                    i) + ':\n' + line)

                            # Check if ff atom index is the same as openmm atom index
                            assert int(ls[0]) - 1 == self.atoms[int(ls[0]) - 1].index

                            # Save mass and sigma values into list
                            self.particles_mass[int(ls[0]) - 1] = float(ls[1])
                            self.rf_sigma[int(ls[0]) - 1] = float(ls[2])

                        # Reading [bonds] section
                        if bonds:
                            if len(ls) > 4:
                                raise ValueError('More than four parameters given in [bonds] section at line ' + str(
                                    i) + ':\n' + line)
                            if len(ls) < 4:
                                raise ValueError('Less than four parameters given in [bonds] section at line ' + str(
                                    i) + ':\n' + line)
                            at1 = self.atoms[int(ls[0]) - 1]
                            at2 = self.atoms[int(ls[1]) - 1]
                            bond_length = float(ls[2]) * unit.nanometer
                            k = float(ls[3])
                            self.bonds[(at1, at2)] = (bond_length, k)

                    # Select which section is being reading by changing its boolean to 1
                    if line.startswith('[atoms]'):
                        print('Reading atom parameters')
                        self.bonds = OrderedDict()
                        atoms = True
                    elif line.startswith('[bonds]'):
                        print('Reading bond parameters')
                        self.bonds = OrderedDict()
                        bonds = True

        print('Done reading Forcefield paramters')
        print('')

    def setCAMassPerResidueType(self):
        """
        Sets the masses of the alpha carbon atoms to the average mass
        of its amino acid residue.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        # Load mass parameters from parameters package
        aa_masses = ca_parameters.aa_masses

        masses = []

        for r in self.topology.residues():
            if r.name in aa_masses:
                masses.append(aa_masses[r.name])
            else:
                raise ValueError('Residue ' + r.name + ' not found in masses dictionary.')

        self.setParticlesMasses(masses)

    def setCARadiusPerResidueType(self):
        """
        Sets the excluded volume radii of the alpha carbon atoms
        to characteristic radii of their corresponding amino acid
        residue.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # Load radii from parameters package
        aa_radii = ca_parameters.aa_radii

        radii = []

        for r in self.topology.residues():
            if r.name in aa_radii:
                radii.append(aa_radii[r.name])
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
        None

        Returns
        -------
        None
        """

        # Load charge from parameters package
        aa_charge = ca_parameters.aa_charge

        charge = []

        for r in self.topology.residues():
            if r.name in aa_charge:
                charge.append(aa_charge[r.name])
            else:
                raise ValueError('Residue ' + r.name + ' not found in radii dictionary.')

        self.setParticlesCharge(charge)

    def setCAHPSPerResidueType(self):
        """
        Sets the HPS of the alpha carbon atoms.
        The current implementation is using Urry scale.

        Parameters
        ----------
        None

        Returns
        -------
        None
        """

        # Load hydropathy scale from parameters package
        aa_hps = ca_parameters.aa_hps

        hps = []

        for r in self.topology.residues():
            if r.name in aa_hps:
                hps.append(aa_hps[r.name])
            else:
                raise ValueError('Residue ' + r.name + ' not found in radii dictionary.')

        self.setParticlesHPS(hps)

    ## User-hidden functions ##

    def _setParameters(term, parameters):
        """
        General function to set up or change force field parameters.

        Parameters
        ----------
        term : dict
            Dictionary object containing the set of degrees of freedom
            (DOF) to set up attributes to (e.g. :code:`bonds` attribute)

        parameters : integer or float or list
            Value(s) for the specific forcefield parameters. If integer
            or float, sets up the same value for all the DOF in terms. If
            list, sets a unique parameter for each DOF.

        Returns
        -------
        None
        """

        if isinstance(parameters, int):
            parameters = float(parameters)

        # Set constant parameter for each item in FF term
        if isinstance(parameters, float):
            for item in term:
                term[item] = term[item][:-1] + (parameters,)

        # Set unique parameter for each item in FF term
        if isinstance(parameters, list):
            assert len(parameters) == len(list(term.keys()))
            for i, item in enumerate(term):
                term[item] = term[item][:-1] + (parameters[i],)
