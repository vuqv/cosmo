"""
COSMO: COarse-grained Simulation of intrinsically disordered prOteins with openMM

COSMO is a Python library to run CG simulations using OpenMM toolkit.
The library offers flexibility for creating CG model that can be customised to implemented different potential model.

Considering an input structure, the library automatizes the creation of forces to specify it.

The library offers methods to tailor forcefield parameter for each force term.

hpsOpenMM is divided in three main classes:

    1. geometry

    2. models

    3. system

The first class, geometry, contains methods to calculate the geometrical parameters from the input structures.

The library is open-source and offers flexibility.
"""
__all__ = ['system', 'models', 'dynamics', 'build_structure']

from .core import geometry
from .core import models
from .core import system
from .dynamics import dynamics
from .parameters import model_parameters
from .reporter import cosmoReporter
from .utils import build_structure
