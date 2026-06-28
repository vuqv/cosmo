Running simulations
=========================================================

A simulation is driven by a control file (``md.ini``). The canonical runner
parses it, builds the coarse-grained model, runs the dynamics and writes the
trajectory / log / checkpoint / final structure. Launch it as a console command
or a module:

.. code-block:: bash

    cosmo-mdrun -f md.ini
    python -m cosmo.mdrun -f md.ini

The runner is intentionally thin: control-file parsing lives in
:func:`cosmo.read_simulation_config` (returning a :class:`cosmo.SimulationConfig`)
and the build / setup / run / finalize steps live in :mod:`cosmo.engine`, so a
custom workflow can reuse the same pieces.

Runner
------------------------------------------------------------

.. automodule:: cosmo.mdrun.mdrun
    :members:

Configuration
------------------------------------------------------------

.. autoclass:: cosmo.utils.config.SimulationConfig
    :members:

.. autofunction:: cosmo.utils.config.read_simulation_config

Engine
------------------------------------------------------------

.. automodule:: cosmo.engine
    :members:
