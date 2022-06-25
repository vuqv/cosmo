"""
Dictionary contains parameters for hps model.
First level is the model name
"""

bond_length = dict(kr=0.38, urry=0.382, nucleic=0.5)
"""
Bond length for various models. nucleic acids is only worked on KR scale
"""

parameters = {
    "kr": {
        "ALA": { #H
            "mass": 71.08,
            "radii": 0.504,
            "charge": 0.0,
            "hps": 0.730
        },
        "ARG": {
            "mass": 156.20,
            "radii": 0.656,
            "charge": 1.0,
            "hps": 0.000
        },
        "ASN": {
            "mass": 114.10,
            "radii": 0.568,
            "charge": 0.0,
            "hps": 0.432
        },
        "ASP": {
            "mass": 115.10,
            "radii": 0.558,
            "charge": -1.0,
            "hps": 0.378
        },
        "CYS": {
            "mass": 103.10,
            "radii": 0.548,
            "charge": 0.0,
            "hps": 0.595
        },
        "GLU": {
            "mass": 129.10,
            "radii": 0.592,
            "charge": -1.0,
            "hps": 0.459
        },
        "GLN": {
            "mass": 128.10,
            "radii": 0.602,
            "charge": 0.0,
            "hps": 0.514
        },
        "GLY": {
            "mass": 57.05,
            "radii": 0.450,
            "charge": 0.0,
            "hps": 0.649
        },
        "HIS": { #A
            "mass": 137.10,
            "radii": 0.608,
            "charge": 0.5,
            "hps": 0.514
        },
        "ILE": { #H
            "mass": 113.20,
            "radii": 0.618,
            "charge": 0.0,
            "hps": 0.973
        },
        "LEU": { #H
            "mass": 113.20,
            "radii": 0.618,
            "charge": 0.0,
            "hps": 0.973
        },
        "LYS": {
            "mass": 128.20,
            "radii": 0.636,
            "charge": 1.0,
            "hps": 0.514
        },
        "ALY": {
            "mass": 170.20,
            "radii": 0.681,
            "charge": 0.0,
            "hps": 0.7567
        },  # N-epsilon acetyllysine
        "MET": {#H
            "mass": 131.20,
            "radii": 0.618,
            "charge": 0.0,
            "hps": 0.838
        },
        "PHE": { #A
            "mass": 147.20,
            "radii": 0.636,
            "charge": 0.0,
            "hps": 1.000
        },
        "PRO": {
            "mass": 97.12,
            "radii": 0.556,
            "charge": 0.0,
            "hps": 1.000
        },
        "SER": {
            "mass": 87.08,
            "radii": 0.518,
            "charge": 0.0,
            "hps": 0.595
        },
        "SEP": {
            "mass": 165.03,
            "radii": 0.636,
            "charge": -2.0,
            "hps": 0.162
        },  # phosposerine: phosphorylation of serine
        "THR": {
            "mass": 101.10,
            "radii": 0.562,
            "charge": 0.0,
            "hps": 0.676
        },
        "TPO": {
            "mass": 179.05,
            "radii": 0.662,
            "charge": -2.0,
            "hps": 0.0081
        },  # Phosphothreonine:  phosphorylation of Threonine
        "TRP": { #A
            "mass": 186.20,
            "radii": 0.678,
            "charge": 0.0,
            "hps": 0.946
        },
        "TYR": { #A
            "mass": 163.20,
            "radii": 0.646,
            "charge": 0.0,
            "hps": 0.865
        },
        "PTR": {
            "mass": 241.15,
            "radii": 0.738,
            "charge": -2.0,
            "hps": 0.189
        },  # phosphotyrosine
        "VAL": {#H
            "mass": 99.07,
            "radii": 0.586,
            "charge": 0.0,
            "hps": 0.892
        },
        # Parameters for RNA (KR scale from OPLS-AA forcefield)
        "A": {
            "mass": 329.20,
            "radii": 0.844,
            "charge": -1.0,
            "hps": -0.054
        },
        "C": {
            "mass": 305.2,
            "radii": 0.822,
            "charge": -1.0,
            "hps": -0.027
        },
        "G": {
            "mass": 345.2,
            "radii": 0.851,
            "charge": -1.0,
            "hps": -0.189,
        },
        "U": {
            "mass": 306.2,
            "radii": 0.817,
            "charge": -1.0,
            "hps": -0.027
        }
    },
    "urry": {
        "ALA": {
            "mass": 71.08,
            "radii": 0.504,
            "charge": 0.0,
            "hps": 0.522942
        },
        "ARG": {
            "mass": 156.20,
            "radii": 0.656,
            "charge": 1.0,
            "hps": 0.478824
        },
        "ASN": {
            "mass": 114.10,
            "radii": 0.568,
            "charge": 0.0,
            "hps": 0.508236
        },
        "ASP": {
            "mass": 115.10,
            "radii": 0.558,
            "charge": -1.0,
            "hps": 0.214119
        },
        "CYS": {
            "mass": 103.10,
            "radii": 0.548,
            "charge": 0.0,
            "hps": 0.56706
        },
        "GLU": {
            "mass": 129.10,
            "radii": 0.592,
            "charge": -1.0,
            "hps": -0.08
        },
        "GLN": {
            "mass": 128.10,
            "radii": 0.602,
            "charge": 0.0,
            "hps": 0.478824
        },
        "GLY": {
            "mass": 57.05,
            "radii": 0.450,
            "charge": 0.0,
            "hps": 0.49353
        },
        "HIS": {
            "mass": 137.10,
            "radii": 0.608,
            "charge": 0.0,
            "hps": 0.684707
        },
        "ILE": {
            "mass": 113.20,
            "radii": 0.618,
            "charge": 0.0,
            "hps": 0.625883
        },
        "LEU": {
            "mass": 113.20,
            "radii": 0.618,
            "charge": 0.0,
            "hps": 0.640589
        },
        "LYS": {
            "mass": 128.20,
            "radii": 0.636,
            "charge": 1.0,
            "hps": 0.302354
        },
        "MET": {
            "mass": 131.20,
            "radii": 0.618,
            "charge": 0.0,
            "hps": 0.596471
        },
        "PHE": {
            "mass": 147.20,
            "radii": 0.636,
            "charge": 0.0,
            "hps": 0.74353
        },
        "PRO": {
            "mass": 97.12,
            "radii": 0.556,
            "charge": 0.0,
            "hps": 0.678824
        },
        "SER": {
            "mass": 87.08,
            "radii": 0.518,
            "charge": 0.0,
            "hps": 0.508236
        },
        "THR": {
            "mass": 101.10,
            "radii": 0.562,
            "charge": 0.0,
            "hps": 0.508236
        },
        "TRP": {
            "mass": 186.20,
            "radii": 0.678,
            "charge": 0.0,
            "hps": 0.92
        },
        "TYR": {
            "mass": 163.20,
            "radii": 0.646,
            "charge": 0.0,
            "hps": 0.817059
        },
        "VAL": {
            "mass": 99.07,
            "radii": 0.586,
            "charge": 0.0,
            "hps": 0.584707
        }
    }
}
"""
* HPS-Kr scale was taken from:
Dignon, G. L., Zheng, W., Kim, Y. C., Best, R. B., ; Mittal, J. (2018). 
Sequence determinants of protein phase behavior from a coarse-grained model. 
PLoS Computational Biology, 1–23. 
https://doi.org/10.1101/238170

* Parameter for Nucleic acids (KR scale):
Regy, R. M., Dignon, G. L., Zheng, W., Kim, Y. C., Mittal, J. (2020). 
Sequence dependent phase separation of protein-polynucleotide mixtures elucidated using molecular simulations. 
Nucleic Acids Research, 48(22), 12593–12603. 
https://doi.org/10.1093/nar/gkaa1099

* Phosphorylation version of some residues for KR scale are taken from:
Perdikari, T. M., Jovic, N., Dignon, G. L., Kim, Y. C., Fawzi, N. L.,  Mittal, J. (2021). 
A predictive coarse-grained model for position-specific effects of post-translational modifications. 
Biophysical Journal, 120(7), 1187–1197. 
https://doi.org/10.1016/j.bpj.2021.01.034

-----------------------------------------------
Note on hps (lambda) in urry scale:

# These parameters were shifted by 0.08 from original parameters directly.

in the original paper: 
Regy, R. M., Thompson, J., Kim, Y. C., ; Mittal, J. (2021). 
Improved coarse-grained model for studying sequence dependent phase separation of disordered proteins.
Protein Science, 30(7), 1371–1379. https://doi.org/10.1002/pro.4094

lambda_ij = muy*lambda0_ij -delta
where muy=1, delta= 0.08 is the optimal set for the set of 42 proteins they studied.
lambda_ij = lambda0_ij - 0.08 = 0.5*(lambda_i+lambda_j) - 0.08 = 0.5(lambda_i -0.08 + lambda_j-0.08)

In both version, KR and Urry, we can tune directly lambda parameter in Urry by 0.08 so we can use only one equation for
two model (choose parameter when passing hps_scale parameter)
"""