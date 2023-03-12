"""
core package of the hpsOpenMM package that contains the main cosmo classes.

The cosmo.core package contains the three hpsOpenMM main classes:

    1. geometry

    2. models

    3. system

The first class, geometry, contains methods to calculate the geometrical parameters from the input structures.

The second class, models, allows to easily set up predefined CG models.
The final class, system, is the main class that holds all the methods to define, modify and create CG to be simulated
with OpenMM.

"""
__all__ = ['system', 'models']
from .geometry import geometry
from .models import models
from .system import system
