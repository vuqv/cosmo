Simulation control options
==========================

Simulation parameters are read from an ``.ini`` file (e.g. ``md.ini``) by
:func:`cosmo.read_simulation_config`, which returns a
:class:`~cosmo.utils.config.SimulationConfig`. The section title ``[OPTIONS]`` is
required.

* Comments: inline or on their own line, starting with ``;`` or ``#``.
* Keyword and value are separated by ``=`` or ``:``.
* Every option below has a default **except** the ones marked *required*; you
  only need to set the options you want to change.

Running a simulation
--------------------

The runner lives in the package as :mod:`cosmo.mdrun`. Once COSMO is installed
(``pip install -e .`` from the repo root) any of these are equivalent:

.. code-block:: bash

    cosmo-mdrun -f md.ini             # installed console command
    python -m cosmo.mdrun -f md.ini   # module form
    python run_simulation.py -f md.ini  # tutorial shim (calls cosmo.mdrun)

Example ``md.ini``:

.. code-block::

        [OPTIONS]
        md_steps = 500_000   ; number of steps (underscores allowed)
        dt = 0.01            ; time step in ps
        nstxout = 1000       ; steps between trajectory (DCD) frames
        nstlog = 1000        ; steps between log writes
        nstchk = 1000        ; steps between checkpoint writes (defaults to nstxout)
        nstcomm = 100        ; center-of-mass motion removal; off by default, omit for multi-chain runs
        log_precision = 4    ; decimal places for float columns in the .log
        log_width = 14       ; min column width for aligned, fixed-width .log
        ; force-field model: hps_urry, hps_kr, hps_ss, or mpipi
        model = hps_urry
        ; bond treatment: None (flexible harmonic bonds, default) or AllBonds (rigid)
        ; constraints = AllBonds
        ; constraint_tolerance = 1e-5   ; only used with constraints = AllBonds

        ; temperature coupling
        tcoupl = yes
        ref_t = 310          ; Kelvin
        tau_t = 0.01         ; ps^-1

        ; pressure coupling (requires pbc = yes)
        pcoupl = no
        ref_p = 1
        frequency_p = 25

        ; periodic boundary condition
        pbc = yes
        box_dimension = 30   ; cubic 30 nm; or [30, 30, 60] for a box

        ; input
        pdb_file = asyn.pdb
        ; init_position = traj/asyn_final.pdb  ; optional: start from these coords
        ; output  (all files -> <output_dir>/<outname>.*, default traj/traj.*)
        output_dir = traj
        outname = asyn
        ; hardware
        device = GPU
        ppn = 4
        ; restart
        restart = no
        minimize = yes


Parameter summary
+++++++++++++++++

"Required = yes" means the run cannot proceed without it. Options with a default
may be omitted. ``—`` in the *Default* column means there is no default (the
option is either required, or only meaningful in a specific mode noted in the
description).

.. list-table::
   :header-rows: 1
   :widths: 20 14 14 14 38

   * - Option
     - Type
     - Required
     - Default
     - Description
   * - ``md_steps``
     - int
     - no
     - ``1000``
     - Total number of integration steps. Underscores are allowed (``500_000``).
   * - ``dt``
     - float [ps]
     - no
     - ``0.01``
     - Integration time step.
   * - ``nstxout``
     - int
     - no
     - ``10``
     - Steps between writing the trajectory (DCD) frames. Also the default checkpoint write frequency (see ``nstchk``).
   * - ``nstlog``
     - int
     - no
     - ``10``
     - Steps between writing the energy/temperature log.
   * - ``nstchk``
     - int
     - no
     - ``nstxout``
     - Steps between checkpoint (``.chk``) writes. Independent of ``nstxout``; when omitted it falls back to ``nstxout``.
   * - ``nstcomm``
     - int
     - no
     - ``none`` (off)
     - Steps between center-of-mass motion removals; off by default (set ``none`` or omit to disable). COM removal suits a single chain but couples the drift of independent chains, so leave it off (or use a large value) for multi-chain / slab runs.
   * - ``log_precision``
     - int
     - no
     - ``4``
     - Decimal places for floating-point columns (energies, time, temperature, ...) in the ``.log``. Set to ``none`` for OpenMM's full ``repr`` precision.
   * - ``log_width``
     - int
     - no
     - ``14``
     - Minimum width (characters) of each ``.log`` column, right-justified, so columns line up. Each column uses ``max(header_length, log_width)``. Set to ``none`` to disable fixed-width formatting.
   * - ``model``
     - str
     - no
     - ``hps_urry``
     - Force-field model: ``hps_urry``, ``hps_kr``, ``hps_ss``, or ``mpipi``. See :doc:`../tutorials/02_models`.
   * - ``tcoupl``
     - bool
     - no
     - ``yes``
     - Langevin thermostat on/off. (NVE is not used.)
   * - ``ref_t``
     - float [K]
     - if ``tcoupl = yes``
     - ``—``
     - Reference temperature. Read only when ``tcoupl = yes``.
   * - ``tau_t``
     - float [ps⁻¹]
     - if ``tcoupl = yes``
     - ``—``
     - Friction coefficient coupling the system to the heat bath. Read only when ``tcoupl = yes``.
   * - ``pcoupl``
     - bool
     - no
     - ``no``
     - Monte Carlo barostat on/off. Requires ``pbc = yes``.
   * - ``ref_p``
     - float [bar]
     - if ``pcoupl = yes``
     - ``1``
     - Reference pressure. Read only when ``pcoupl = yes``.
   * - ``frequency_p``
     - int [steps]
     - if ``pcoupl = yes``
     - ``25``
     - Barostat move attempt frequency. Read only when ``pcoupl = yes``.
   * - ``pbc``
     - bool
     - no
     - ``no``
     - Periodic boundary conditions on/off.
   * - ``box_dimension``
     - float or [x,y,z] [nm]
     - if ``pbc = yes``
     - ``—``
     - Box size: a scalar ``L`` gives a cubic ``L×L×L`` box; a list ``[x, y, z]`` a rectangular box.
   * - ``pdb_file``
     - str
     - **yes**
     - ``—``
     - Input structure (``.pdb`` / ``.cif``) used to build the model (topology, force field) and the initial coordinates.
   * - ``init_position``
     - str
     - no
     - ``—``
     - Optional starting coordinates for a fresh run (e.g. a previous run's ``_final.pdb``). When set, these coordinates are used as-is (no origin shift) instead of ``pdb_file``'s; the atom count must match the built system. Ignored on restart (coordinates come from the checkpoint).
   * - ``output_dir``
     - str
     - no
     - ``traj``
     - Folder for all generated files; created if missing. One run = one self-contained folder.
   * - ``outname``
     - str
     - no
     - ``traj``
     - Basename for generated files: ``<output_dir>/<outname>.dcd``, ``.log``, ``.psf``, ``.chk``, ``_init.pdb``, ``_final.pdb`` (and ``_ff.dat`` for HPS models).
   * - ``protein_code``
     - str
     - no
     - ``—``
     - Legacy output prefix. When set and ``outname`` is not given, it is used as ``outname``. Prefer ``outname`` in new control files.
   * - ``checkpoint``
     - str
     - no
     - ``<output_dir>/<outname>.chk``
     - Explicit checkpoint path override. Normally leave unset so it lands in the run folder.
   * - ``device``
     - str
     - no
     - ``CPU``
     - Compute platform: ``CPU`` or ``GPU`` (CUDA).
   * - ``ppn``
     - int
     - no
     - ``1``
     - Number of CPU threads. Read only when ``device = CPU``.
   * - ``restart``
     - bool
     - no
     - ``no``
     - Restart from ``checkpoint`` instead of the PDB coordinates. Forces ``minimize = no``.
   * - ``minimize``
     - bool
     - no
     - ``yes``
     - Energy-minimize the input structure before dynamics. Forced ``no`` when ``restart = yes``.
   * - ``constraints``
     - str
     - no
     - ``None``
     - Bond treatment: ``None`` (flexible harmonic bonds, the default) or ``AllBonds`` (rigid distance constraints, so a larger ``dt`` can be used). See the note below.
   * - ``constraint_tolerance``
     - float
     - no
     - ``1e-5``
     - Integrator relative constraint tolerance. Only meaningful when ``constraints = AllBonds``; harmless otherwise.

.. note::

   Boolean options accept ``yes``/``no``, ``true``/``false``, ``1``/``0``
   (parsed by :func:`cosmo.utils.config.strtobool`).


Notes on individual options
+++++++++++++++++++++++++++

Output layout (``output_dir`` / ``outname``)
    Every generated file is written to ``<output_dir>/<outname><suffix>``, so a run
    is one self-contained folder (default ``traj/``): ``traj.dcd`` (trajectory),
    ``traj.log`` (state log), ``traj.psf`` (topology), ``traj.chk`` (checkpoint),
    ``traj_init.pdb`` (the built/minimized starting structure), ``traj_final.pdb``
    (last conformation), ``traj_ff.dat`` (per-residue force-field dump, HPS
    models only), and ``traj_runinfo.log`` (run provenance: software versions,
    hardware, GPU, timing). ``output_dir`` is created automatically if missing. To keep
    several runs side by side, point each at its own folder (e.g.
    ``output_dir = runs/FUS_T300``) or change ``outname``. ``protein_code`` is the
    legacy prefix: when set and ``outname`` is unset it becomes the basename.

Output frequency and log formatting (``nstxout`` / ``nstlog`` / ``nstchk`` / ``nstcomm`` / ``log_precision`` / ``log_width``)
    ``nstxout`` controls trajectory (DCD) frames; ``nstchk`` controls checkpoint
    (``.chk``) writes independently and defaults to ``nstxout`` when omitted.
    ``nstlog`` controls how often a row is appended to the ``.log``.
    ``nstcomm`` is the center-of-mass motion removal interval — **off by default**
    (omit it, or set ``none``). Enable it for a single chain; leave it off for
    multi-chain / slab runs, where COM removal would couple the drift of
    independent chains. The log reporter is
    :class:`cosmo.cosmoReporter`, which writes a fixed-width, aligned ``.log``:
    ``log_precision`` sets the decimals for float columns (default ``4``;
    ``none`` for full precision) and ``log_width`` sets the minimum column width
    (default ``14``; ``none`` to disable padding). Columns are separated by two
    spaces and the header is a ``#`` comment, so the log stays both human-readable
    and machine-parsable (see :func:`cosmo.reporter.readOpenMMReporterFile`).

Force-field model (``model``)
    Selects the coarse-grained potential: ``hps_urry`` (default, Urry hydropathy
    scale), ``hps_kr`` (Kapcha–Rossky scale; broad residue coverage including RNA
    and phosphorylated residues, but less accurate), ``hps_ss`` (``hps_urry`` plus
    bonded angle and torsion terms), or ``mpipi`` (Wang–Frenkel non-bonded
    potential). All models use flexible harmonic bonds by default (see the
    ``constraints`` option to switch to rigid constraints). See
    :doc:`../tutorials/02_models`.

Bond treatment (``constraints`` / ``constraint_tolerance``)
    By default (``constraints`` unset or ``None``) the CA/P chain bonds are
    **flexible harmonic springs** — the physically appropriate choice for disordered
    chains. Setting ``constraints = AllBonds`` instead makes every bond a **rigid
    distance constraint** pinned at its equilibrium length: the harmonic bond force is
    not created, the fast bond-stretch mode is removed, and a larger ``dt`` can be
    used. Constraints act **only** on the pseudo-bonds — the non-bonded potentials
    (Ashbaugh–Hatch / Wang–Frenkel and Debye–Hückel electrostatics) are untouched.
    ``constraint_tolerance`` (default ``1e-5``) sets the integrator's relative
    constraint tolerance and is read only in the ``AllBonds`` case. Unlike the sibling
    ``topo`` package (a Gō model that defaults to ``AllBonds``), cosmo defaults to
    flexible bonds — a deliberate IDP-physics choice.

Temperature / pressure coupling
    ``ref_t`` and ``tau_t`` are only consumed when ``tcoupl = yes`` (and are then
    required); ``ref_p`` and ``frequency_p`` only when ``pcoupl = yes``. Pressure
    coupling additionally requires ``pbc = yes`` (asserted at parse time).

Periodic boundary conditions
    Turning ``pbc`` on affects the non-bonded forces and how coordinates are
    written to the PDB/DCD (handled internally). ``box_dimension`` must be given
    when ``pbc = yes``: a scalar ``L`` gives a cubic ``L×L×L`` box and a list
    ``[x, y, z]`` a rectangular box (e.g. an elongated box for slab/LLPS runs).

Hardware (``device`` / ``ppn``)
    ``device = GPU`` runs on CUDA (mixed precision, device 0). ``device = CPU``
    uses ``ppn`` threads; ``ppn`` is ignored on GPU.

``restart`` and ``minimize``
    ``restart = yes`` loads positions **and** velocities from ``checkpoint`` and
    continues; reporters append to the existing log/trajectory, and ``minimize``
    is forced off. With ``restart = no`` you may choose ``minimize``; on a fresh
    start velocities are drawn from the Boltzmann distribution at ``ref_t``. The
    checkpoint path defaults to ``<output_dir>/<outname>.chk`` and can be
    overridden with ``checkpoint`` (e.g. to restart from a differently named run).

Starting coordinates (``init_position``)
    On a fresh run (``restart = no``) the initial coordinates come from
    ``pdb_file`` by default. Set ``init_position`` to start instead from another
    structure — typically a previous run's ``<outname>_final.pdb`` — to chain runs
    together without restarting from a checkpoint. These coordinates are used
    as-is (no origin shift) and must have the same atom count as the built system.
    Ignored on a restart, where the checkpoint supplies positions and velocities.
