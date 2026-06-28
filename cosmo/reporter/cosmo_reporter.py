#!/usr/bin/env python
# coding: utf-8

import re
import time

from openmm.app.statedatareporter import StateDataReporter
from openmm import unit

from cosmo.core import system


class cosmoReporter(StateDataReporter):
    """
    A :class:`~openmm.app.StateDataReporter` that writes a clean, fixed-width log.

    Compared with OpenMM's reporter this adds:

    * ``precision`` -- floating-point columns are written with a fixed number of
      decimals instead of full ``repr`` precision.
    * ``width`` -- every column is padded to a fixed width so the columns line up
      (each column uses ``max(len(header), width)``, so long headers still fit).
    * ``sbmObject`` -- optionally append one energy column per force group of a
      COSMO system object.

    The header line is written as a ``#`` comment, and columns are separated by
    the ``separator`` string (use two spaces for a readable, machine-parsable log;
    see :func:`readOpenMMReporterFile`).

    For more customization reporter, check this:
    http://docs.openmm.org/latest/userguide/application/04_advanced_sim_examples.html#extracting-and-reporting-forces-and-other-data

    Methods
    -------
    """

    def __init__(self, file, reportInterval, sbmObject=None, precision=4, width=14, **kwargs):
        """
        Initialises the COSMO reporter.

        Parameters
        ----------
        reportInterval : int
            The interval (in time steps) at which to write frames
        sbmObject : cosmo.core.system, optional
            If given, append one energy column per force group.
        precision : int or None, optional (default: 4)
            Decimal places for floating-point columns. ``None`` keeps OpenMM's
            full ``repr`` precision.
        width : int or None, optional (default: 14)
            Minimum total width of each column (right-justified) for aligned
            output. ``None`` disables fixed-width formatting (plain separator
            join, OpenMM-style quoted header).
        **kwargs : openMM StateDataReporter arguments

        Returns
        -------
        initialized StateDataReporter class.

        """
        super(cosmoReporter, self).__init__(file, reportInterval, **kwargs)
        self._sbmObject = sbmObject
        self._precision = precision
        self._width = width
        self._col_widths = None

    def _constructHeaders(self):
        """
        Build headers for the StateDataReporter class.

        Build headers for the StateDataReporter class. It builds the headers
        for the force groups contained in the sbmOpenMM system instance.

        Parameters
        ----------
        None

        Returns
        -------
        headers : list
            List with strings representing the headers to be written to the report file.
        """

        headers = super()._constructHeaders()
        if isinstance(self._sbmObject, system):
            for i, n in enumerate(self._sbmObject.forceGroups):
                headers.append(n + ' (kJ/mol)')

        return headers

    def _constructReportValues(self, simulation, state):
        """
        Calculates the energies for the force groups in the COSMO system instance.

        Parameters
        ----------
        None

        Returns
        -------
        values : list
            List with floats representing the values to be written to the report file.
        """

        values = super()._constructReportValues(simulation, state)

        if isinstance(self._sbmObject, system):
            for i, n in enumerate(self._sbmObject.forceGroups):
                values.append(
                    simulation.context.getState(getEnergy=True, groups={i}).getPotentialEnergy().value_in_unit(
                        unit.kilojoules_per_mole))

        return values

    def _value_to_str(self, v):
        """Format one value: floats to fixed precision, everything else via ``str``."""
        if self._precision is not None and isinstance(v, float):
            return f'{v:.{self._precision}f}'
        return str(v)

    def _format_header(self, headers):
        if self._width is None:
            # OpenMM-style quoted, separator-joined comment line.
            return '#"%s"' % ('"' + self._separator + '"').join(headers)
        line = self._separator.join(h.rjust(w) for h, w in zip(headers, self._col_widths))
        # Mark the header as a comment without shifting the columns: replace the
        # first leading space with '#', or prepend '# ' if there is none.
        return ('#' + line[1:]) if line[:1] == ' ' else ('# ' + line)

    def _format_row(self, fields):
        if self._width is None:
            return self._separator.join(fields)
        return self._separator.join(f.rjust(w) for f, w in zip(fields, self._col_widths))

    def report(self, simulation, state):
        """Generate a report (fixed-width when ``width`` is set)."""
        if not self._hasInitialized:
            self._initializeConstants(simulation)
            headers = self._constructHeaders()
            if self._width is not None:
                # Column width = max(header length, requested width) so headers and
                # data align even when a header is longer than `width`.
                self._col_widths = [max(len(h), self._width) for h in headers]
            if not self._append:
                print(self._format_header(headers), file=self._out)
            try:
                self._out.flush()
            except AttributeError:
                pass
            self._initialClockTime = time.time()
            self._initialSimulationTime = state.getTime()
            self._initialSteps = simulation.currentStep
            self._hasInitialized = True

        self._checkForErrors(simulation, state)
        values = [self._value_to_str(v) for v in self._constructReportValues(simulation, state)]
        print(self._format_row(values), file=self._out)
        try:
            self._out.flush()
        except AttributeError:
            pass


def readOpenMMReporterFile(reporter_file):
    """
    Read a fixed-width COSMO/OpenMM log into a ``{column_name: [values...]}`` dict.

    The first line is the ``#``-commented header. Columns are separated by runs of
    two or more spaces (so multi-word headers like ``Potential Energy (kJ/mole)``
    stay intact), and numeric data rows are split on whitespace.

    Parameters
    ----------
    reporter_file : str
        Path to the reporter output file.
    """
    with open(reporter_file, 'r') as ef:
        lines = [ln.rstrip('\n') for ln in ef if ln.strip()]

    header = lines[0].lstrip('#').strip()
    names = [c.strip().strip('"') for c in re.split(r'\s{2,}', header)]
    data = {name: [] for name in names}
    for line in lines[1:]:
        if line.lstrip().startswith('#'):
            continue
        fields = line.split()
        for name, value in zip(names, fields):
            try:
                data[name].append(float(value))
            except ValueError:
                data[name].append(value)   # non-numeric column (e.g. Time Remaining)
    return data
