import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import pandas as pd
import re
from collections import defaultdict
import yaml
from typing import Dict, List, Tuple, Set, Optional
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Constants
KCAL_TO_KJ = 4.184
DEFAULT_CUTOFF = 4.5
LOCAL_SEPARATION = 2
SIGMA_SCALE_FACTOR = 2**(1/6)
DISTANCE_TO_NM = 10.0

# Energy parameters (in kcal/mol, converted to kJ/mol)
ENERGY_PARAMS = {
    'hydrogen_bond': 0.75 * KCAL_TO_KJ,
    'backbone_sidechain': 0.37 * KCAL_TO_KJ,
    'non_native': 0.000132 * KCAL_TO_KJ,
    'yang_shift': 0.6
}

def get_residue_mapping(universe: mda.Universe) -> Tuple[Dict[int, int], Dict[int, str], int]:
    """
    Get residue mapping from universe.
    
    Args:
        universe: MDAnalysis universe object
        
    Returns:
        Tuple of (resid_to_index, index_to_resname, n_residues)
    """
    residues = universe.select_atoms("protein").residues
    n_residues = len(residues)
    resid_to_index = {res.resid: idx for idx, res in enumerate(residues)}
    index_to_resname = {idx: res.resname for idx, res in enumerate(residues)}
    return resid_to_index, index_to_resname, n_residues

def parse_hydrogen_bonds(stride_output_file: str) -> List[Tuple]:
    """
    Parse hydrogen bonds from STRIDE output file.
    
    Args:
        stride_output_file: Path to STRIDE output file
        
    Returns:
        List of hydrogen bond pairs
    """
    try:
        with open(stride_output_file, "r") as f:
            lines = f.readlines()
    except FileNotFoundError:
        logger.error(f"STRIDE output file not found: {stride_output_file}")
        raise
    
    pattern = re.compile(
        r"(?:DNR|ACC)\s+(\w+)\s+-\s+(\d+)\s+\d+\s+->\s+(\w+)\s+-\s+(\d+)\s+\d+"
    )
    
    hb_pairs = []
    seen_pairs = set()
    
    for line in lines:
        if line.startswith(("DNR", "ACC")):
            match = pattern.search(line)
            if match:
                res1, res1_pdb, res2, res2_pdb = match.groups()
                
                if line.startswith("DNR"):
                    donor = (res1, int(res1_pdb))
                    acceptor = (res2, int(res2_pdb))
                else:  # ACC
                    acceptor = (res1, int(res1_pdb))
                    donor = (res2, int(res2_pdb))
                
                # Avoid duplicates
                pair = tuple(sorted([donor, acceptor]))
                if pair not in seen_pairs:
                    seen_pairs.add(pair)
                    hb_pairs.append(pair)
    
    return hb_pairs

def build_hb_contact_matrix(hb_pairs: List[Tuple], n_residues: int) -> np.ndarray:
    """
    Build hydrogen bond contact matrix.
    
    Args:
        hb_pairs: List of hydrogen bond pairs
        n_residues: Number of residues
        
    Returns:
        Hydrogen bond contact matrix
    """
    hb_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)
    
    # Count hydrogen bonds between residue pairs
    pair_counts = defaultdict(int)
    for donor, acceptor in hb_pairs:
        pair = tuple(sorted([donor, acceptor]))
        pair_counts[pair] += 1
    
    # Fill contact matrix
    for (res1, res2), count in pair_counts.items():
        i = res1[1] - 1  # 0-based index
        j = res2[1] - 1  # 0-based index
        hb_contact_matrix[i, j] = count
        hb_contact_matrix[j, i] = count
    
    return hb_contact_matrix

def get_hb_contact_matrix(stride_output_file: str, n_residues: int) -> np.ndarray:
    """
    Get hydrogen bond contact matrix from STRIDE output.
    
    Args:
        stride_output_file: Path to STRIDE output file
        n_residues: Number of residues
        
    Returns:
        Hydrogen bond contact matrix
    """
    hb_pairs = parse_hydrogen_bonds(stride_output_file)
    return build_hb_contact_matrix(hb_pairs, n_residues)

def get_bs_contact_matrix(u: mda.Universe, cutoff: float = DEFAULT_CUTOFF) -> np.ndarray:
    """
    Get Backbone-Sidechain Contact Matrix.
    
    Args:
        u: MDAnalysis universe
        cutoff: Distance cutoff for contacts
        
    Returns:
        Backbone-sidechain contact matrix
    """
    resid_to_index, _, n_residues = get_residue_mapping(u)
    
    backbone = u.select_atoms('protein and backbone and not name H*')
    sidechain = u.select_atoms("protein and not backbone and not name H*")
    
    dists_bs = distance_array(backbone.positions, sidechain.positions)
    
    # Build directional residue-residue contact list
    contacts_bs = set()
    for i, atom1 in enumerate(backbone):
        for j, atom2 in enumerate(sidechain):
            if (dists_bs[i, j] <= cutoff and 
                abs(atom1.resid - atom2.resid) > LOCAL_SEPARATION):
                pair = (atom1.resid, atom2.resid)  # preserve direction: bb → sc
                contacts_bs.add(pair)
    
    # Build asymmetric matrix first
    bs_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)
    for contact in sorted(contacts_bs, key=lambda x: x[0]):
        bs_contact_matrix[resid_to_index[contact[0]], 
                        resid_to_index[contact[1]]] = 1
    
    # Build symmetric count matrix
    bs_symmetric_count = np.zeros((n_residues, n_residues), dtype=int)
    for i in range(n_residues):
        for j in range(n_residues):
            if abs(i - j) > LOCAL_SEPARATION:
                count = bs_contact_matrix[i, j] + bs_contact_matrix[j, i]
                bs_symmetric_count[i, j] = count
                bs_symmetric_count[j, i] = count
    
    return bs_symmetric_count

def get_ss_contact_matrix(u: mda.Universe, cutoff: float = DEFAULT_CUTOFF) -> np.ndarray:
    """
    Get sidechain-sidechain contact matrix.
    
    Args:
        u: MDAnalysis universe
        cutoff: Distance cutoff for contacts
        
    Returns:
        Sidechain-sidechain contact matrix
    """
    resid_to_index, _, n_residues = get_residue_mapping(u)
    sidechain = u.select_atoms("protein and not backbone and not name H*")
    
    dists_ss = distance_array(sidechain.positions, sidechain.positions)
    ss_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)
    
    # Build residue-residue contact list
    contacts_ss = set()
    for i, atom1 in enumerate(sidechain):
        for j, atom2 in enumerate(sidechain):
            if (dists_ss[i, j] <= cutoff and 
                abs(atom1.resid - atom2.resid) > LOCAL_SEPARATION):
                pair = tuple(sorted([atom1.resid, atom2.resid]))
                contacts_ss.add(pair)
    
    # Fill contact matrix
    for contact in sorted(contacts_ss, key=lambda x: x[0]):
        ss_contact_matrix[resid_to_index[contact[0]], 
                        resid_to_index[contact[1]]] = 1
        ss_contact_matrix[resid_to_index[contact[1]], 
                        resid_to_index[contact[0]]] = 1
    
    return ss_contact_matrix

def load_bt_potential(bt_file: str = 'bt_potential.csv') -> pd.DataFrame:
    """
    Load BT potential from CSV file.
    
    Args:
        bt_file: Path to BT potential CSV file
        
    Returns:
        BT potential matrix in kJ/mol
    """
    try:
        df = pd.read_csv(bt_file, index_col=0)
        return KCAL_TO_KJ * np.abs(df - ENERGY_PARAMS['yang_shift'])
    except FileNotFoundError:
        logger.error(f"BT potential file not found: {bt_file}")
        raise

def get_ss_interaction_energy(u: mda.Universe, bt_file: str = 'bt_potential.csv') -> np.ndarray:
    """
    Calculate sidechain-sidechain interaction energy matrix.
    
    Args:
        u: MDAnalysis universe
        bt_file: Path to BT potential CSV file
        
    Returns:
        Sidechain-sidechain interaction energy matrix
    """
    eps_ss = load_bt_potential(bt_file)
    _, index_to_resname, n_residues = get_residue_mapping(u)
    
    sc_interaction_energy = np.zeros((n_residues, n_residues))
    for i in range(n_residues):
        for j in range(n_residues):
            sc_interaction_energy[i, j] = eps_ss.loc[index_to_resname[i], 
                                                   index_to_resname[j]]
    
    return sc_interaction_energy

def parse_residue_list(residue_items: List) -> List[int]:
    """
    Parse residue list from various formats.
    
    Args:
        residue_items: List of residue items (int or str)
        
    Returns:
        List of residue numbers
    """
    residues = []
    for item in residue_items:
        if isinstance(item, int):
            residues.append(item)
        elif isinstance(item, str):
            if '-' in item:
                start, end = map(int, item.split('-'))
                residues.extend(range(start, end + 1))
            else:
                residues.append(int(item))
    return residues

def read_yaml_config(filepath: str) -> Tuple[Dict, Dict, Dict]:
    """
    Read and parse YAML configuration file.
    
    Args:
        filepath: Path to YAML configuration file
        
    Returns:
        Tuple of (domain_to_residues, intra_strengths, inter_strengths)
    """
    try:
        with open(filepath, 'r') as f:
            config = yaml.safe_load(f)
    except FileNotFoundError:
        logger.error(f"Domain configuration file not found: {filepath}")
        raise
    
    intra = config['intra_domains']
    inter = config['inter_domains']
    n_residues = int(config['n_residues'])
    
    domain_to_residues = {}
    intra_strengths = {}
    all_residues = set()
    
    # Parse intra-domain configurations
    for domain, values in intra.items():
        raw_residues = values['residues']
        residues = parse_residue_list(raw_residues)
        domain_to_residues[domain] = residues
        intra_strengths[domain] = values['strength']
        all_residues.update(residues)
    
    # Handle unassigned residues
    full_residues = set(range(1, n_residues + 1))
    unassigned_residues = sorted(full_residues - all_residues)
    if unassigned_residues:
        domain_to_residues['X'] = unassigned_residues
        intra_strengths['X'] = 1.0
    
    # Parse inter-domain configurations
    inter_strengths = {}
    for pair_str, strength in inter.items():
        d1, d2 = pair_str.strip().split('-')
        inter_strengths[(d1, d2)] = strength
        inter_strengths[(d2, d1)] = strength  # ensure symmetry
    
    # Add inter-domain interactions for domain X
    if 'X' in domain_to_residues:
        for other in domain_to_residues:
            if other != 'X':
                inter_strengths[('X', other)] = 1.0
                inter_strengths[(other, 'X')] = 1.0
    
    return domain_to_residues, intra_strengths, inter_strengths

def get_scaling_ss_matrix(domain_def: str) -> np.ndarray:
    """
    Build scaling matrix for sidechain-sidechain interactions.
    
    Args:
        domain_def: Path to domain definition YAML file
        
    Returns:
        Scaling matrix
    """
    domain_to_residues, intra_strengths, inter_strengths = read_yaml_config(domain_def)
    
    # Build residue to domain mapping
    residue_to_domain = {}
    residue_list = []
    for domain, residues in domain_to_residues.items():
        for res in residues:
            residue_to_domain[res] = domain
            residue_list.append(res)
    
    residue_list = sorted(set(residue_list))
    res_to_idx = {res: i for i, res in enumerate(residue_list)}
    n = len(residue_list)
    matrix = np.zeros((n, n))
    
    # Fill scaling matrix
    for i_res in residue_list:
        i_idx = res_to_idx[i_res]
        dom_i = residue_to_domain[i_res]
        for j_res in residue_list:
            j_idx = res_to_idx[j_res]
            dom_j = residue_to_domain[j_res]
            
            if dom_i == dom_j:
                matrix[i_idx, j_idx] = intra_strengths[dom_i]
            else:
                key = (dom_i, dom_j)
                if key in inter_strengths:
                    matrix[i_idx, j_idx] = inter_strengths[key]
                elif (dom_j, dom_i) in inter_strengths:
                    matrix[i_idx, j_idx] = inter_strengths[(dom_j, dom_i)]
                else:
                    matrix[i_idx, j_idx] = 0.0
    
    return matrix

def calculate_sigma_values(binary_contact_matrix: np.ndarray, ca_distances: np.ndarray, 
                          n_residues: int) -> List[float]:
    """
    Calculate sigma values for repulsive interactions.
    
    Args:
        binary_contact_matrix: Binary contact matrix
        ca_distances: CA-CA distance matrix
        n_residues: Number of residues
        
    Returns:
        List of sigma values
    """
    sigma = []
    for i in range(n_residues):
        not_in_contact_with_i = [
            j for j in range(n_residues) 
            if abs(i - j) > LOCAL_SEPARATION and binary_contact_matrix[i, j] == 0
        ]
        if not_in_contact_with_i:
            distance_to_i = ca_distances[i, not_in_contact_with_i]
            sigma.append(SIGMA_SCALE_FACTOR * np.min(distance_to_i))
        else:
            sigma.append(0.0)  # fallback value
    return sigma

def build_nonbonded_interaction(pdb_file: str, domain_def: str, stride_output_file: str) -> Tuple[np.ndarray, np.ndarray]:
    """
    Build non-bonded interaction matrices.
    
    Args:
        pdb_file: Path to PDB file
        domain_def: Path to domain definition YAML file
        stride_output_file: Path to STRIDE output file
        
    Returns:
        Tuple of (distance_matrix, energy_matrix)
    """
    logger.info("Loading protein structure...")
    u = mda.Universe(pdb_file)
    resid_to_index, index_to_resname, n_residues = get_residue_mapping(u)
    
    logger.info("Building hydrogen bond contact matrix...")
    hb_contact_matrix = get_hb_contact_matrix(stride_output_file, n_residues)
    hb_interaction_energy = ENERGY_PARAMS['hydrogen_bond'] * hb_contact_matrix
    
    logger.info("Building backbone-sidechain contact matrix...")
    bs_contact_matrix = get_bs_contact_matrix(u, cutoff=DEFAULT_CUTOFF)
    bs_interaction_energy = bs_contact_matrix * ENERGY_PARAMS['backbone_sidechain']
    
    logger.info("Building sidechain-sidechain contact matrix...")
    scaling_matrix = get_scaling_ss_matrix(domain_def)
    ss_contact_matrix = get_ss_contact_matrix(u, cutoff=DEFAULT_CUTOFF)
    ss_interaction_energy = get_ss_interaction_energy(u)
    
    # Element-wise multiplication
    scaled_ss_interaction_energy = scaling_matrix * ss_contact_matrix * ss_interaction_energy
    
    # Total interaction energy for native contacts
    eps_ij = hb_interaction_energy + bs_interaction_energy + scaled_ss_interaction_energy
    
    # Build binary contact matrix
    contact_matrix = hb_contact_matrix + bs_contact_matrix + ss_contact_matrix
    binary_contact_matrix = (contact_matrix > 0).astype(int)
    
    # Calculate distance matrix
    ca_atoms = u.select_atoms('protein and name CA')
    ca_distances = distance_array(ca_atoms, ca_atoms)
    distance_matrix = np.zeros_like(ca_distances)
    
    # Set distances for native contacts
    contact_mask = binary_contact_matrix == 1
    distance_matrix[contact_mask] = ca_distances[contact_mask]
    
    # Calculate sigma values for non-native contacts
    sigma = calculate_sigma_values(binary_contact_matrix, ca_distances, n_residues)
    
    # Set distances and energies for non-native contacts
    for i in range(n_residues):
        for j in range(n_residues):
            if binary_contact_matrix[i, j] == 0:
                distance_matrix[i, j] = 0.5 * (sigma[i] + sigma[j])
                eps_ij[i, j] = ENERGY_PARAMS['non_native']
    
    # Convert to nm for OpenMM compatibility
    distance_matrix /= DISTANCE_TO_NM
    
    logger.info("Non-bonded interaction matrices built successfully")
    return distance_matrix, eps_ij

# Main execution
if __name__ == "__main__":
    try:
        R_ij, eps_ij = build_nonbonded_interaction('2ww4.pdb', 'domain.yaml', 'stride.dat')
        logger.info(f"Distance matrix shape: {R_ij.shape}")
        logger.info(f"Energy matrix shape: {eps_ij.shape}")
    except Exception as e:
        logger.error(f"Error building non-bonded interactions: {e}")
        raise