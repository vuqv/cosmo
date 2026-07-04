"""
Shared build / setup / run / finalize machinery for the COSMO runner.

This is the reusable core that the CLI runner (:mod:`cosmo.mdrun`) orchestrates.
Splitting it out of the old monolithic ``cosmo.dynamics.Dynamics.run`` (and the
copy-pasted ``run_simulation.py`` scripts) means a custom driver can reuse the
exact same steps::

    cfg   = cosmo.read_simulation_config("md.ini")
    built = cosmo.engine.build_system(cfg)
    ctx   = cosmo.engine.setup_simulation(cfg, built)
    cosmo.engine.attach_reporters(cfg, ctx.simulation, append=ctx.restart_active,
                                  total_steps=cfg.md_steps)
    ctx.simulation.step(ctx.nsteps_remain)
    cosmo.engine.finalize_simulation(cfg, ctx, built.topology, start_time)

The design mirrors the sibling ``topo`` project's ``engine.py``.
"""
import time
from dataclasses import dataclass

import numpy as np
import openmm as mm
from openmm import unit

from .core import models
from .reporter import cosmoReporter
from .utils import runinfo


@dataclass
class RunContext:
    """Everything the runner needs after the OpenMM ``Simulation`` is set up.

    Attributes
    ----------
    simulation : openmm.app.Simulation
        The configured simulation (positions/velocities set, or checkpoint loaded).
    restart_active : bool
        True if this run resumed from a checkpoint.
    done_steps : int
        Steps already completed (0 for a fresh run; the checkpoint's step count
        on a restart).
    nsteps_remain : int
        Steps still to run to reach ``cfg.md_steps``.
    checkpoint : str
        Resolved checkpoint path for this run.
    runinfo_path : str
        Path of the ``_runinfo.log`` provenance file for this run.
    """
    simulation: object
    restart_active: bool
    done_steps: int
    nsteps_remain: int
    checkpoint: object = None
    runinfo_path: object = None


def build_model(cfg):
    """Build the coarse-grained model and write its provenance files.

    Builds the model via :func:`cosmo.models.buildCoarseGrainModel` and dumps the
    force-field data (HPS models only), the initial coordinates and the
    topology. It does **not** add the MD extras (COM remover / barostat) -- that
    is :func:`add_md_forces` -- so a specialized runner can modify the system
    (e.g. tweak a force's global parameter, add a restraint) in between.

    Parameters
    ----------
    cfg : cosmo.SimulationConfig
        Parsed control file.

    Returns
    -------
    cosmo.core.system.system
        The built model (carries ``.system``, ``.topology`` and ``.positions``).
    """
    # All outputs go to cfg.output_dir (default traj/); create it up front so the
    # provenance dumps below and every later writer land in the run folder.
    cfg.prepare_output_dir()

    hps_model = models.buildCoarseGrainModel(cfg.pdb_file, **cfg.build_kwargs())

    # Provenance dumps (same set and naming as the historical runner).
    if cfg.writes_forcefield():
        # dumpForceFieldData only supports the sigma/epsilon-per-residue format
        # used by the HPS-scale models.
        hps_model.dumpForceFieldData(cfg.output_path('_ff.dat'))
    hps_model.dumpStructure(cfg.output_path('_init.pdb'))
    hps_model.dumpTopology(cfg.output_path('.psf'))

    return hps_model


def add_md_forces(cfg, built):
    """Add the standard MD extras to a built system, in place.

    Adds a center-of-mass motion remover (when ``cfg.nstcomm`` is set) and, when
    pressure coupling is requested, a Monte Carlo barostat. Both are added to
    ``built.system`` before any OpenMM ``Context`` is created.

    Parameters
    ----------
    cfg : cosmo.SimulationConfig
        Parsed control file.
    built : cosmo.core.system.system
        Output of :func:`build_model`.
    """
    # Remove center-of-mass motion (energy-neutral; does not affect benchmarks).
    # Skipped when nstcomm is 0/None so a runner can opt out.
    if cfg.nstcomm:
        built.system.addForce(mm.CMMotionRemover(cfg.nstcomm))

    # Pressure coupling: add a barostat. read_simulation_config already asserts
    # pbc is on when pcoupl is on.
    if cfg.pcoupl:
        barostat = mm.MonteCarloBarostat(cfg.ref_p, cfg.ref_t, cfg.frequency_p)
        built.system.addForce(barostat)


def build_system(cfg):
    """Build the model, write provenance files and add the MD extras.

    Convenience composition of :func:`build_model` + :func:`add_md_forces`, used
    by the canonical runner. Returns the ``cosmo.system`` model ready for
    :func:`setup_simulation`.
    """
    built = build_model(cfg)
    add_md_forces(cfg, built)
    return built


def setup_simulation(cfg, built, control_file=None, shift_positions=True):
    """Create the OpenMM ``Simulation`` and prime it for stepping.

    On a fresh run the input coordinates are (optionally) shifted into the
    positive octant (min corner at the origin), positions are set and velocities
    are drawn from ``ref_t``. On a restart the checkpoint is loaded and the
    remaining step count is computed.

    Parameters
    ----------
    cfg : cosmo.SimulationConfig
        Parsed control file.
    built : cosmo.core.system.system
        Output of :func:`build_model` / :func:`build_system`.
    control_file : str, optional
        Path to the control file, recorded in the run-metadata header.
    shift_positions : bool, optional (default: True)
        Shift the input coordinates so the minimum corner is at the origin.
        Set False when the absolute coordinate frame is meaningful (e.g. a
        co-translational run that restrains atoms relative to a fixed point).
        Ignored on restart (coordinates come from the checkpoint) and when
        ``init_position`` is given (those coordinates are used as-is).

    Returns
    -------
    RunContext
    """
    platform, properties = cfg.make_platform()
    integrator = mm.LangevinIntegrator(cfg.ref_t, cfg.tau_t, cfg.dt)
    simulation = mm.app.Simulation(built.topology, built.system, integrator,
                                   platform, properties)

    checkpoint = cfg.checkpoint_path()

    if cfg.restart:
        simulation.loadCheckpoint(checkpoint)
        done_steps = simulation.context.getState().getStepCount()
        print(f"Restarting from step: {done_steps}")
        nsteps_remain = cfg.md_steps - done_steps
        restart_active = True
        coord_source = f"checkpoint ({checkpoint})"
        vel_source = f"checkpoint ({checkpoint})"
    else:
        if cfg.init_position:
            # Explicit starting coordinates (e.g. a previous run's _final.pdb).
            # Used as-is -- no shift, since the absolute frame is deliberate.
            init_pos = mm.app.PDBFile(cfg.init_position).getPositions()
            if len(init_pos) != built.system.getNumParticles():
                raise SystemExit(
                    f"init_position '{cfg.init_position}' has {len(init_pos)} atoms but "
                    f"the system has {built.system.getNumParticles()}; they must match.")
            built.positions = init_pos
            coord_source = f"init_position ({cfg.init_position})"
        else:
            if shift_positions:
                # Shift coordinates so the minimum corner sits at the origin.
                xyz = np.array(built.positions / unit.nanometer)
                xyz[:, 0] -= np.amin(xyz[:, 0])
                xyz[:, 1] -= np.amin(xyz[:, 1])
                xyz[:, 2] -= np.amin(xyz[:, 2])
                built.positions = xyz * unit.nanometer
            coord_source = f"pdb_file ({cfg.pdb_file})"
        simulation.context.setPositions(built.positions)
        simulation.context.setVelocitiesToTemperature(cfg.ref_t)
        done_steps = 0
        nsteps_remain = cfg.md_steps
        restart_active = False
        vel_source = f"Boltzmann distribution at {cfg.ref_t}"

    # Record run provenance (package versions, hardware, GPU, timing) to a
    # side-channel <outname>_runinfo.log -- does not affect the simulation.
    runinfo_path = cfg.output_path('_runinfo.log')
    runinfo.write_run_start(
        runinfo_path,
        control_file=control_file,
        checkpoint_file=checkpoint,
        restart=restart_active,
        steps_planned=nsteps_remain,
        simulation=simulation,
        use_gpu=(cfg.device == 'GPU'),
        ppn=cfg.ppn,
        coord_source=coord_source,
        vel_source=vel_source,
    )
    print(f"Writing run metadata to {runinfo_path}")

    return RunContext(simulation=simulation, restart_active=restart_active,
                      done_steps=done_steps, nsteps_remain=nsteps_remain,
                      checkpoint=checkpoint, runinfo_path=runinfo_path)


def attach_reporters(cfg, simulation, append=False, total_steps=None):
    """Attach checkpoint, trajectory and log reporters to ``simulation``.

    Replaces any existing reporters. Files are written to
    ``<output_dir>/<outname>.{chk,dcd,log}``.

    Parameters
    ----------
    cfg : cosmo.SimulationConfig
        Parsed control file.
    simulation : openmm.app.Simulation
        The simulation to attach reporters to.
    append : bool, optional
        Append to existing trajectory/log files (used on restart).
    total_steps : int, optional
        Total step count advertised to the StateDataReporter (for the
        remaining-time/speed estimate). Defaults to ``cfg.md_steps``.
    """
    if total_steps is None:
        total_steps = cfg.md_steps

    # Checkpoint frequency (nstchk) is independent of the trajectory frequency
    # (nstxout); it falls back to nstxout when not set in the config.
    simulation.reporters = []
    simulation.reporters.append(
        mm.app.CheckpointReporter(cfg.checkpoint_path(), cfg.checkpoint_interval()))
    simulation.reporters.append(
        mm.app.DCDReporter(cfg.output_path('.dcd'), cfg.nstxout,
                           enforcePeriodicBox=bool(cfg.pbc), append=append))
    # cosmoReporter writes a clean, fixed-width log: each float column uses
    # log_precision decimals and every column is padded to log_width characters so
    # the columns line up. Columns are separated by two spaces (aligned and still
    # machine-parsable via cosmo.reporter.readOpenMMReporterFile).
    simulation.reporters.append(
        cosmoReporter(
            cfg.output_path('.log'), cfg.nstlog,
            precision=cfg.log_precision, width=cfg.log_width,
            step=True, time=True,
            potentialEnergy=True, kineticEnergy=True, totalEnergy=True,
            temperature=True, remainingTime=True, speed=True,
            totalSteps=total_steps, separator='  ', append=append))


def finalize_simulation(cfg, ctx, topology, start_time, write_pdb=None):
    """Write the final conformation + checkpoint and report wall-clock time.

    Parameters
    ----------
    cfg : cosmo.SimulationConfig
        Parsed control file.
    ctx : RunContext
        The run context from :func:`setup_simulation`.
    topology : openmm.app.Topology
        Topology to pair with the final coordinates (``built.topology``).
    start_time : float
        ``time.time()`` captured when stepping began.
    write_pdb : callable, optional
        Custom final-frame writer ``write_pdb(topology, positions, path)``. When
        None (default) the standard ``openmm.app.PDBFile.writeFile`` is used.
        Pass e.g. ``cosmo.utils.write_pdb_with_chain_ids`` to preserve chain IDs.
    """
    simulation = ctx.simulation
    last_frame = simulation.context.getState(
        getPositions=True, enforcePeriodicBox=bool(cfg.pbc)).getPositions()
    final_path = cfg.output_path('_final.pdb')
    if write_pdb is None:
        with open(final_path, 'w') as fh:
            mm.app.PDBFile.writeFile(topology, last_frame, fh)
    else:
        write_pdb(topology, last_frame, final_path)
    simulation.saveCheckpoint(cfg.checkpoint_path())

    # Close the run-metadata side channel with timing + final-state info.
    if ctx.runinfo_path is not None:
        runinfo.write_run_end(ctx.runinfo_path, simulation=simulation,
                              start_epoch=start_time, final_structure=final_path)
    print(f"Finished in {(time.time() - start_time):.2f} seconds")
