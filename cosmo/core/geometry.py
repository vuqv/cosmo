#!/usr/bin/env python
# coding: utf-8

import numpy as np
import openmm.unit as unit


class geometry:
    """
    A class to hold methods for calculating geometrical values 
    given sets of atom coordinates.

    Methods
    -------

    """

    @staticmethod
    def position2Array(position: unit.quantity.Quantity, output_unit):
        """Converts an OpenMM position object quantity into a numpy array.

        Parameters
        ----------
        position : openmm.unit.quantity.Quantity
            Array containing quantity objects [e.g. (x,y,z) array returned
            from positions].
        output_unit : openmm.unit.nanometer
            Unit in which to return the items of the array.

        Returns
        -------
        numpy.ndarray
            A numpy array containing the quantity values in unit of nm, converted to float.
        """

        return np.array([c.value_in_unit(output_unit) for c in position])

    @staticmethod
    def bond(coord1: unit.quantity.Quantity, coord2: unit.quantity.Quantity):
        """Calculate the distance length between two (x,y,z) quantity coordinates.

        Parameters
        ----------
        coord1 : openmm.unit.quantity.Quantity array
            Vector for the first coordinate.
        coord2 : openmm.unit.quantity.Quantity array
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
