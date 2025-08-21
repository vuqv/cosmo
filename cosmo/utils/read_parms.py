import numpy as np
import json
import os

"""
It loads dihedral angle parameters from a JSON file named karanicolas_dihe_parm.json located in a subdirectory called data, 
relative to the current Python file's location.

__file__: The path to the current Python script.
os.path.dirname(__file__): Gets the directory containing the script.
'data', 'karanicolas_dihe_parm.json': Specifies the subdirectory and JSON filename.
os.path.join(...): Constructs the full path to the JSON file.
"""
# Load dihedral parameters from JSON file
def load_dihedral_parameters():
    json_path = os.path.join(os.path.dirname(__file__), 'data', 'karanicolas_dihe_parm.json')
    with open(json_path, 'r') as f:
        return json.load(f)