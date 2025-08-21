import numpy as np
import yaml
import re

def parse_residue_list(residue_items):
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

    domain_to_residues = {}
    intra_strengths = {}

    for domain, values in intra.items():
        raw_residues = values['residues']
        residues = parse_residue_list(raw_residues)
        domain_to_residues[domain] = residues
        intra_strengths[domain] = values['strength']

    inter_strengths = {}
    for pair_str, strength in inter.items():
        d1, d2 = pair_str.strip().split('-')
        inter_strengths[(d1, d2)] = strength

    return domain_to_residues, intra_strengths, inter_strengths

def generate_interaction_matrix(domain_to_residues, intra_strengths, inter_strengths):
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

    return matrix, res_to_idx


domain_to_residues, intra_strengths, inter_strengths = read_yaml_config("domain.yaml")
matrix, res_idx = generate_interaction_matrix(domain_to_residues, intra_strengths, inter_strengths)

print(matrix[res_idx[1], res_idx[22]])   # Within domain A
print(matrix[res_idx[1], res_idx[60]])   # Between A and B
