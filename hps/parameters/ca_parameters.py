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
Mass of the Coarse-grained beads are set to residue masses.
Unit in amu
"""

aa_radii = {'ALA': 0.504, 'ARG': 0.656, 'ASN': 0.568,
            'ASP': 0.558, 'CYS': 0.548, 'GLU': 0.592,
            'GLN': 0.602, 'GLY': 0.450, 'HIS': 0.608,
            'ILE': 0.618, 'LEU': 0.618, 'LYS': 0.636,
            'MET': 0.618, 'PHE': 0.636, 'PRO': 0.556,
            'SER': 0.518, 'THR': 0.562, 'TRP': 0.678,
            'TYR': 0.646, 'VAL': 0.586}
"""
vdw radius of the coarse-grained beads. 
Unit in nanometer
"""
aa_hps_urry = {'ALA': 0.522942, 'ARG': 0.478824, 'ASN': 0.508236,
               'ASP': 0.214119, 'CYS': 0.56706, 'GLU': -0.08,
               'GLN': 0.478824, 'GLY': 0.49353, 'HIS': 0.684707,
               'ILE': 0.625883, 'LEU': 0.640589, 'LYS': 0.302354,
               'MET': 0.596471, 'PHE': 0.74353, 'PRO': 0.678824,
               'SER': 0.508236, 'THR': 0.508236, 'TRP': 0.92,
               'TYR': 0.817059, 'VAL': 0.584707}

aa_hps_kr = {'ALA': 0.730, 'ARG': 0.000, 'ASN': 0.432,
             'ASP': 0.378, 'CYS': 0.595, 'GLU': 0.459,
             'GLN': 0.514, 'GLY': 0.649, 'HIS': 0.514,
             'ILE': 0.973, 'LEU': 0.973, 'LYS': 0.514,
             'MET': 0.838, 'PHE': 1.000, 'PRO': 1.000,
             'SER': 0.595, 'THR': 0.676, 'TRP': 0.946,
             'TYR': 0.865, 'VAL': 0.892}
# aa_hps = {'ALA': 0.602942, 'ARG': 0.558824, 'ASN': 0.588236,
#           'ASP': 0.294119, 'CYS': 0.64706, 'GLU': 0.,
#           'GLN': 0.558824, 'GLY': 0.57353, 'HIS': 0.764707,
#           'ILE': 0.705883, 'LEU': 0.720589, 'LYS': 0.382354,
#           'MET': 0.676471, 'PHE': 0.82353, 'PRO': 0.758824,
#           'SER': 0.588236, 'THR': 0.588236, 'TRP': 1.0,
#           'TYR': 0.897059, 'VAL': 0.664707}
"""
Hydropathy scale (Urry scale)
Here we use muy=1 and delta = 0.08. Directly scale the hps.
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
Charge of the coarse-grain beads, assigns to it alpha-carbon atoms.
Charges are determined at neutral pH (pH=7).
Setting to the values from, where ARG, LYS =1 and ASP, GLU = -1,
Other residues have charge of 0.
Regy, R. M., Thompson, J., Kim, Y. C., & Mittal, J. (2021). 
Improved coarse-grained model for studying sequence dependent phase separation of disordered proteins. 
Protein Science, 30(7), 1371–1379. https://doi.org/10.1002/pro.4094.
Unit in elementary charge
"""
