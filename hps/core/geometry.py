#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from simtk import unit


class geometry:
    """
    A class to hold methods for calculating geometrical values 
    given sets of atom coordinates.

    Parameters
    ----------
    None

    Methods
    -------
    position2Array:
        convert position to array
    bond:
        get bond between two point
    """

    def position2Array(position, output_unit):
        """Converts an OpenMM position object quantity into a numpy array.

        Parameters
        ----------
        position : simtk.unit.quantity.Quantity
            Array containing quantity objects [e.g. (x,y,z) array returned
            from positions].
        output_unit : openmm.unit
            Unit in which to return the items of the array.

        Returns
        -------
        numpy.ndarray
            A numpy array containing the quantity values converted to floats.
        """

        return np.array([c.value_in_unit(output_unit) for c in position])

    def bond(coord1, coord2):
        """Calculate the distance length between two (x,y,z) quantity coordinates.

        Parameters
        ----------
        coord1 : simtk.unit.quantity.Quantity array
            Vector for the first coordinate.
        coord2 : simtk.unit.quantity.Quantity array
            Vector for the second coordinate.

        Returns
        -------
        simtk.unit.quantity.Quantity
            Quantity (value and unit) of the distance length in nanometers.
        """

        coord1 = geometry.position2Array(coord1, unit.nanometer)

        coord2 = geometry.position2Array(coord2, unit.nanometer)

        bond_length = np.linalg.norm(coord2 - coord1)

        return bond_length * unit.nanometer
