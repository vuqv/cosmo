"""
This module contains bead's parameters to use in simulation.
They are mass, radius, hydropathy scale, charge
"""

aa_masses = {'ALA': 71.08, 'ARG': 156.20, 'ASN': 114.10,
             'ASP': 115.10, 'CYS': 103.10, 'GLU': 129.10,
             'GLN': 128.10, 'GLY': 57.05, 'HIS': 137.10,
             'ILE': 113.20, 'LEU': 113.20, 'LYS': 128.20,
             'MET': 131.20, 'PHE': 147.20, 'PRO': 97.12,
             'SER': 87.08, 'THR': 101.10, 'TRP': 186.20,
             'TYR': 163.20, 'VAL': 99.07}
"""
Mass of beads. Unit in amu
"""

aa_radii = {'ALA': 0.504, 'ARG': 0.656, 'ASN': 0.568,
            'ASP': 0.558, 'CYS': 0.548, 'GLU': 0.592,
            'GLN': 0.602, 'GLY': 0.450, 'HIS': 0.608,
            'ILE': 0.618, 'LEU': 0.618, 'LYS': 0.636,
            'MET': 0.618, 'PHE': 0.636, 'PRO': 0.556,
            'SER': 0.518, 'THR': 0.562, 'TRP': 0.667,
            'TYR': 0.646, 'VAL': 0.586}
"""
vdw radius of beads. Unit in Angstrom
"""

aa_hps = {'ALA': 0.602942, 'ARG': 0.558824, 'ASN': 0.588236,
          'ASP': 0.294119, 'CYS': 0.64706, 'GLU': 0.,
          'GLN': 0.558824, 'GLY': 0.57353, 'HIS': 0.764707,
          'ILE': 0.705883, 'LEU': 0.720589, 'LYS': 0.382354,
          'MET': 0.676471, 'PHE': 0.82353, 'PRO': 0.758824,
          'SER': 0.588236, 'THR': 0.588236, 'TRP': 1.0,
          'TYR': 0.897059, 'VAL': 0.664707}
"""
Hydropathy scale (Urry scale)
Regy, R. M., Thompson, J., Kim, Y. C., & Mittal, J. (2021). 
Improved coarse-grained model for studying sequence dependent phase separation of disordered proteins. 
Protein Science, 30(7), 1371–1379. https://doi.org/10.1002/pro.4094
"""

aa_charge = {'ALA': 0.0, 'ARG': 1.0, 'ASN': 0.0,
             'ASP': -1.0, 'CYS': 0.0, 'GLU': -1.0,
             'GLN': 0.0, 'GLY': 0.0, 'HIS': 0.0,
             'ILE': 0.0, 'LEU': 0.0, 'LYS': 1.0,
             'MET': 0.0, 'PHE': 0.0, 'PRO': 0.0,
             'SER': 0.0, 'THR': 0.0, 'TRP': 0.0,
             'TYR': 0.0, 'VAL': 0.0}
"""
Charge of beads, assigns to it alpha-carbon atoms.
Unit in elementary charge
"""