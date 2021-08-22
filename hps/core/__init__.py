"""
core package of the HPS-Urry package that contains the main HPS-Urry classes.

The hps_urry.core package contains the three sbmOpenMM main classes:

    1. geometry

    2. models

    3. system

The first class, geometry, contains methods to calculate the geometrical parameters from the input structures.

The second class, models, allows to easily set up predefined CG models.
The third class, system, is the main class that holds all the methods to define, modify and create CG to be simulated with OpenMM.

"""

from .geometry import geometry
from .models import models
from .system import system
