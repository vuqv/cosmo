import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np
import pandas as pd
import re
from collections import defaultdict

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


def get_bs_contact_matrix_symmetric(u, cutoff=4.5):
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

def get_bs_contact_matrix_asymmetric(u, cutoff=4.5):
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

    return bs_contact_matrix

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


def build_nonbonded_interaction(pdb_file, stride_output_file):
    
    u = mda.Universe('2ww4.pdb')

    # Build mapping from resid to 0-based index
    residues = u.select_atoms("protein").residues
    resid_to_index = {res.resid: idx for idx, res in enumerate(residues)}
    # Initialize contact matrix
    n_residues = len(residues)
    
    hb_contact_matrix = get_hb_contact_matrix(stride_output_file, n_residues)
    
    bs_contact_matrix = get_bs_contact_matrix_symmetric(u, cutoff=4.5)
    

    ss_contact_matrix = get_ss_contact_matrix(u, cutoff=4.5)
    
    
    contact_matrix = hb_contact_matrix + bs_contact_matrix + ss_contact_matrix 
    binary_contact_matrix = (contact_matrix > 0).astype(int)



R_ij, eps_ij = build_nonbonded_interaction('2ww4.pdb', 'domain.yaml', 'stride.dat')