import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import pandas as pd
import re
from collections import defaultdict
import yaml

def get_hb_contact_matrix(stride_output_file, n_residues):
    with open(stride_output_file, "r") as f:
        lines = f.readlines()

    # Pattern to extract donor and acceptor (using PDB residue numbers)
    pattern = re.compile(
        r"(?:DNR|ACC)\s+(\w+)\s+-\s+(\d+)\s+\d+\s+->\s+(\w+)\s+-\s+(\d+)\s+\d+"
    )

    # Count H-bonds between donor and acceptor PDB residues
    hb_counts = defaultdict(int)

    seen_pairs = set()
    hb_counts = defaultdict(int)

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

                # Avoid duplicates: donor–acceptor and acceptor–donor are considered the same
                pair = tuple([donor, acceptor])

                # Count only once
                if pair not in seen_pairs:
                    seen_pairs.add(pair)
                    hb_counts[pair] += 1


    # Optional: deduplicate bidirectional pairs by treating (A, B) same as (B, A)
    dedup_counts = defaultdict(int)
    for (donor, acceptor), count in hb_counts.items():
        pair = tuple(sorted([donor, acceptor]))  # sort to avoid direction
        dedup_counts[pair] += count

    # Convert results to DataFrame
    hb_df = pd.DataFrame([
        {
            "Residue_1": f"{res1[0]}-{res1[1]}",
            "Residue_2": f"{res2[0]}-{res2[1]}",
            "Num_Hydrogen_Bonds": count
        }
        for (res1, res2), count in dedup_counts.items()
    ])


    hb_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)
    for (res1, res2), count in dedup_counts.items():
        i = res1[1] - 1 # 0-based index
        j = res2[1] - 1 # 0-based index
        hb_contact_matrix[i,j] = count
        hb_contact_matrix[j,i] = count

    return hb_contact_matrix


def get_bs_contact_matrix(u, cutoff=4.5):
    """
    Get Backbone-Sidechain Contact Matrix
    """
    # Build mapping from resid to 0-based index
    residues = u.select_atoms("protein").residues
    resid_to_index = {res.resid: idx for idx, res in enumerate(residues)}
    # Initialize contact matrix
    n_residues = len(residues)

    bs_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)
    
    backbone = u.select_atoms('protein and backbone and not name H*')
    sidechain = u.select_atoms("protein and not backbone and not name H*")

    dists_bs = distance_array(backbone.positions, sidechain.positions)
    
    # Build directional residue–residue contact list
    contacts_bs = set()
    for i, atom1 in enumerate(backbone):
        for j, atom2 in enumerate(sidechain):
            if dists_bs[i, j] <= cutoff and abs(atom1.resid - atom2.resid) > 2:
                pair = (atom1.resid, atom2.resid)  # preserve direction: bb → sc
                contacts_bs.add(pair)

    # Sort by backbone residue index (first value in tuple)
    sorted_contacts_bs = sorted(contacts_bs, key=lambda x: x[0])

    # assign bs_contact_matrix
    for contact in sorted_contacts_bs:
        bs_contact_matrix[resid_to_index[contact[0]], resid_to_index[contact[1]]] = 1    

    """
    count number of bb-sc interactions between residues i and j.
    if bb(i) <-> sc(j) but sc(i) <-/-> bb(j) then i and j are consider to form 1 bb-sc contact
    if bb(i) <-> sc(j) and sc(i) <-> bb(j) then i and j is considered to form 2 bb-sc contacts
    """
    # Now build symmetric count:
    # For each residue pair (i, j), count how many directional bb-sc contacts exist
    bs_symmetric_count = np.zeros((n_residues, n_residues), dtype=int)

    for i in range(n_residues):
        for j in range(n_residues):
            if abs(i - j) > 2:  # skip local pairs if needed
                # Count how many directional bb–sc links between i and j
                count = bs_contact_matrix[i, j] + bs_contact_matrix[j, i]
                bs_symmetric_count[i, j] = count
                bs_symmetric_count[j, i] = count

    return bs_symmetric_count

def get_ss_contact_matrix(u, cutoff=4.5):
    """
    Function to get sidechain-sidechain contact Matrix
    """
    # Build mapping from resid to 0-based index
    residues = u.select_atoms("protein").residues
    resid_to_index = {res.resid: idx for idx, res in enumerate(residues)}
    

    # Initialize contact matrix
    n_residues = len(residues)
    sidechain = u.select_atoms("protein and not backbone and not name H*")
    
    dists_ss = distance_array(sidechain.positions, sidechain.positions)

    ss_contact_matrix = np.zeros((n_residues, n_residues), dtype=int)

    # Build residue–residue contact list
    contacts_ss = set()
    for i, atom1 in enumerate(sidechain):
        for j, atom2 in enumerate(sidechain):
            if dists_ss[i, j] <= cutoff and abs(atom1.resid - atom2.resid) > 2:
                pair = tuple(sorted([ atom1.resid, atom2.resid]))
                contacts_ss.add(pair)
                
    sorted_contact_ss = sorted(contacts_ss, key=lambda x:x[0])

    # assign ss_contact_matrix
    for contact in sorted_contact_ss:
        ss_contact_matrix[resid_to_index[contact[0]], resid_to_index[contact[1]]] = 1
        ss_contact_matrix[resid_to_index[contact[1]], resid_to_index[contact[0]]] = 1

    return ss_contact_matrix

def get_ss_interaction_energy(u):
    # construct SC-SC interaction strength
    # SC-SC was set ti BT potential and Yang shift by -0.6 (line 668-670)
    # he commented is in kT unit however, it looks like in kcal/mol unit.
    # read BT-potential     
    df = pd.read_csv('bt_potential.csv', index_col=0)
    eps_ss = 4.184*np.abs(df-0.6) #shift by 0.6, then convert to kj/mol
    
    residues = u.select_atoms("protein").residues
    n_residues = len(residues)
    index_to_resname = {idx: res.resname for idx, res in enumerate(residues)}
    
    sc_interaction_energy = np.zeros((n_residues, n_residues))
    for i in range(n_residues):
        for j in range(n_residues):
            sc_interaction_energy[i, j] = eps_ss.loc[index_to_resname[i], index_to_resname[j]]
    
    return sc_interaction_energy

def parse_residue_list(residue_items):
    """
    Just read the raw residue index, not shifted yet.
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

def read_yaml_config(filepath):
    with open(filepath, 'r') as f:
        config = yaml.safe_load(f)

    intra = config['intra_domains']
    inter = config['inter_domains']
    n_residues = int(config['n_residues'])

    domain_to_residues = {}
    intra_strengths = {}
    #all-residues present in domain definition
    all_residues = set()

    for domain, values in intra.items():
        raw_residues = values['residues']
        residues = parse_residue_list(raw_residues)
        domain_to_residues[domain] = residues
        intra_strengths[domain] = values['strength']
        
        all_residues.update(residues)

    # Infer full range of residues (assume minimum is 1)    
    full_residues = set(range(1, n_residues+1))
    unassigned_residues = sorted(full_residues - all_residues)
    # assign unassigned residues to domain 'X' (assume X is not presented in the domain def yet.)
    if unassigned_residues:
        domain_to_residues['X'] = unassigned_residues
        intra_strengths['X'] = 1.0

    inter_strengths = {}
    for pair_str, strength in inter.items():
        d1, d2 = pair_str.strip().split('-')
        inter_strengths[(d1, d2)] = strength
        inter_strengths[(d2, d1)] = strength  # ensure symmetry
        
    # add inter-domain interactions for domain X
    if 'X' in domain_to_residues:
        for other in domain_to_residues:
            if other != 'X':
                inter_strengths[('X', other)] = 1.0
                inter_strengths[(other, 'X')] = 1.0

    return domain_to_residues, intra_strengths, inter_strengths

def get_scaling_ss_matrix(domain_def):
    """
    This function works on python matrix so need to shift the residue indices by 1.
    """
    domain_to_residues, intra_strengths, inter_strengths = read_yaml_config(domain_def)
    residue_to_domain = {}
    residue_list = []

    for domain, residues in domain_to_residues.items():
        for res in residues:
            residue_to_domain[res] = domain
            residue_list.append(res)

    # shifted by 1 since python is 0-based index, residues is 1-based index
    residue_list = sorted(set(residue_list))
    res_to_idx = {res: i for i, res in enumerate(residue_list)}
    n = len(residue_list)
    matrix = np.zeros((n, n))

    for i_res in residue_list:
        i_idx = res_to_idx[i_res]
        dom_i = residue_to_domain[i_res]
        for j_res in residue_list:
            j_idx = res_to_idx[j_res]
            dom_j = residue_to_domain[j_res]

            if dom_i == dom_j:
                # intra-domain interaction
                matrix[i_idx, j_idx] = intra_strengths[dom_i]
            else:
                # inter-domain interaction
                # the else-if statement test for symetry, if (dom_i, dom_j) are not in inter_strength, it also test
                # for (dom_j, dom_i)
                key = (dom_i, dom_j)
                if key in inter_strengths:
                    matrix[i_idx, j_idx] = inter_strengths[key]
                elif (dom_j, dom_i) in inter_strengths:
                    matrix[i_idx, j_idx] = inter_strengths[(dom_j, dom_i)]
                else:
                    matrix[i_idx, j_idx] = 0.0

    return matrix #, res_to_idx

def build_nonbonded_interaction(pdb_file, domain_def, stride_output_file):
    
    u = mda.Universe('2ww4.pdb')

    # Build mapping from resid to 0-based index
    residues = u.select_atoms("protein").residues
    resid_to_index = {res.resid: idx for idx, res in enumerate(residues)}
    # Initialize contact matrix
    n_residues = len(residues)
    
    eps_hb = 0.75*4.184 # kj/mol, convert from kcal/mol
    eps_bs = 0.37*4.184 # kj/mol, convert from kcal/mol
    eps_nn = 0.000132 * 4.184 # non-native interaction energy, kj/mol
    
    hb_contact_matrix = get_hb_contact_matrix(stride_output_file, n_residues)
    hb_interaction_energy = eps_hb*hb_contact_matrix
    
    bs_contact_matrix = get_bs_contact_matrix(u, cutoff=4.5)
    bs_interaction_energy = bs_contact_matrix*eps_bs
    
    scaling_matrix = get_scaling_ss_matrix(domain_def)
    ss_contact_matrix = get_ss_contact_matrix(u, cutoff=4.5)
    ss_interaction_energy = get_ss_interaction_energy(u)
    
    # element-wise multiple
    scaled_ss_intertaction_energy = scaling_matrix*ss_contact_matrix*ss_interaction_energy
    
    # total interaction energy for native contact, non-native contact is still 0, set at last
    eps_ij = hb_interaction_energy + bs_interaction_energy + scaled_ss_intertaction_energy
    
    contact_matrix = hb_contact_matrix + bs_contact_matrix + ss_contact_matrix 
    binary_contact_matrix = (contact_matrix > 0).astype(int)

    # set the R_ij for native contacts
    ca_atoms = u.select_atoms('protein and name CA')
    ca_distances = distance_array(ca_atoms, ca_atoms)

    # Mask where contact exists
    contact_mask = binary_contact_matrix == 1

    # Initialize output
    distance_matrix = np.zeros_like(ca_distances)

    # Copy distances where contacts exist
    distance_matrix[contact_mask] = ca_distances[contact_mask]
    
    """
    The repulsive radius of each residue is the minimum distance of residues not in contact with the current residue, then scaled by 2**(1/6)
    """
    sigma = []
    for i in range(n_residues):
        not_in_contact_with_i = [j for j in range(n_residues) if abs(i - j) > 2 and binary_contact_matrix[i, j] == 0]
        distance_to_i = ca_distances[i, not_in_contact_with_i]
        sigma.append(2**(1/6)*np.min(distance_to_i))

    # set distance and interaction energy for non-native contact    
    for i in range(n_residues):
        for j in range(n_residues):
            if binary_contact_matrix[i, j] == 0:
                distance_matrix[i, j] = 0.5 * (sigma[i] + sigma[j])
                eps_ij[i,j] = eps_nn
    distance_matrix /= 10 #convert to nm for consistent with openMM
    return distance_matrix, eps_ij

R_ij, eps_ij = build_nonbonded_interaction('2ww4.pdb', 'domain.yaml', 'stride.dat')