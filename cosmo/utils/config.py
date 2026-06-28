"""
Read a simulation control file (``md.ini``) into a structured configuration.

The control file is an INI file with a single ``[OPTIONS]`` section. Parsing it
used to be a ~100-line block copy-pasted into every ``run_simulation.py`` (and
re-implemented in :class:`cosmo.dynamics.Dynamics`). This module centralizes that
logic so a runner can simply do::

    import cosmo

    cfg = cosmo.read_simulation_config("md.ini")
    model = cosmo.models.buildCoarseGrainModel(cfg.pdb_file, **cfg.build_kwargs())
    integrator = openmm.LangevinIntegrator(cfg.ref_t, cfg.tau_t, cfg.dt)

All parsed values are returned on a :class:`SimulationConfig` dataclass, with
OpenMM units already applied where appropriate (``dt``, ``ref_t``, ``tau_t``,
``ref_p``).

The design mirrors the sibling ``topo`` project's ``utils/config.py``: a parsing
function that returns a dataclass carrying both the settings and a few small
helpers (``build_kwargs``, ``output_path``, ``checkpoint_path``,
``make_platform``) that the runner and engine consume.
"""
import configparser
from dataclasses import dataclass, field
from json import loads
from pathlib import Path
from typing import Any, List, Optional

import openmm as mm
from openmm import unit


def strtobool(value):
    """Convert a truthy/falsy string to 0/1.

    Local replacement for ``distutils.util.strtobool`` (distutils was removed in
    Python 3.12). Accepts the same vocabulary: true/yes/on/1 -> 1,
    false/no/off/0 -> 0 (case-insensitive); anything else raises ValueError.
    """
    v = str(value).strip().lower()
    if v in ("y", "yes", "t", "true", "on", "1"):
        return 1
    if v in ("n", "no", "f", "false", "off", "0"):
        return 0
    raise ValueError(f"invalid truth value {value!r}")


@dataclass
class SimulationConfig:
    """Parsed contents of a COSMO simulation control file (``md.ini``).

    Defaults match the historical behaviour of
    :class:`cosmo.dynamics.Dynamics`. Options absent from the file fall back to
    these defaults.
    """

    # run control
    md_steps: int = 1000
    # default_factory: a Quantity is mutable as far as dataclasses are concerned.
    dt: Any = field(default_factory=lambda: 0.01 * unit.picoseconds)
    nstxout: int = 10
    nstlog: int = 10
    nstcomm: int = 100
    model: str = 'hps_urry'

    # log formatting (cosmoReporter fixed-width columns). precision=None keeps
    # OpenMM's full repr precision; width=None disables fixed-width alignment.
    log_precision: Optional[int] = 4
    log_width: Optional[int] = 14

    # temperature coupling
    tcoupl: bool = True
    ref_t: Any = 300.0          # kelvin when tcoupl is on
    tau_t: Any = None           # ps^-1 when tcoupl is on

    # pressure coupling (requires pbc)
    pcoupl: bool = False
    ref_p: Any = 1.0            # bar when pcoupl is on
    frequency_p: int = 25

    # periodic boundary conditions
    pbc: bool = False
    # Either a scalar L (cubic box) or an [x, y, z] list, matching what
    # cosmo.models.buildCoarseGrainModel accepts; None when pbc is off.
    box_dimension: Optional[Any] = None

    # input
    pdb_file: Optional[str] = None

    # output: all generated files go to <output_dir>/<outname><suffix>, so a run
    # is one self-contained folder (default traj/traj.dcd, traj.log, traj.psf,
    # traj.chk, traj_init.pdb, traj_final.pdb, ...). This mirrors the sibling
    # `topo` project. `protein_code` is the legacy COSMO output prefix; when set
    # (and `outname` is not) it is used as `outname` for backward compatibility.
    output_dir: str = 'traj'
    outname: str = 'traj'
    protein_code: Optional[str] = None
    checkpoint: Optional[str] = None   # explicit override; defaults to <output_dir>/<outname>.chk

    # hardware
    device: str = 'CPU'
    ppn: int = 1

    # restart / minimize
    restart: bool = False
    minimize: bool = True

    # the path this config was read from (bookkeeping)
    config_file: Optional[str] = None

    def build_kwargs(self) -> dict:
        """Keyword arguments for :func:`cosmo.models.buildCoarseGrainModel`.

        Passes ``minimize``, ``model`` and ``box_dimension`` (the builder's own
        defaults cover the optional ``frozen_indices`` / ``except_chains`` /
        ``nb_exclusions`` arguments).
        """
        return dict(minimize=self.minimize, model=self.model,
                    box_dimension=self.box_dimension)

    def output_path(self, suffix: str = '') -> str:
        """Path for a generated output file: ``<output_dir>/<outname><suffix>``.

        Examples (defaults): ``output_path('.dcd')`` -> ``traj/traj.dcd``;
        ``output_path('_init.pdb')`` -> ``traj/traj_init.pdb``.

        Built with :class:`pathlib.Path` but returned as ``str`` so it can be
        passed directly to OpenMM/parmed writers.
        """
        return str(Path(self.output_dir) / f'{self.outname}{suffix}')

    def checkpoint_path(self) -> str:
        """Resolved checkpoint path: the explicit ``checkpoint`` option if given,
        otherwise ``<output_dir>/<outname>.chk``.
        """
        return self.checkpoint if self.checkpoint else self.output_path('.chk')

    def prepare_output_dir(self) -> None:
        """Ensure ``output_dir`` exists.

        ``Path.mkdir(parents=True, exist_ok=True)`` creates any missing parents
        and is a no-op if the folder already exists.
        """
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def writes_forcefield(self) -> bool:
        """Whether this model dumps a ``_ff.dat`` force-field file.

        Only the HPS-scale models expose the sigma/epsilon-per-residue format
        that ``dumpForceFieldData`` can write.
        """
        return self.model in ('hps_kr', 'hps_urry', 'hps_ss')

    def make_platform(self):
        """Build the OpenMM ``(platform, properties)`` pair selected by ``device``.

        ``device = 'GPU'`` -> CUDA (mixed precision, device 0); anything else ->
        CPU using ``ppn`` threads.
        """
        if self.device == 'GPU':
            print("Running simulation on GPU (CUDA)")
            return (mm.Platform.getPlatformByName('CUDA'),
                    {'CudaPrecision': 'mixed', 'DeviceIndex': '0'})
        print(f"Running simulation on CPU using {self.ppn} cores")
        return (mm.Platform.getPlatformByName('CPU'), {'Threads': str(self.ppn)})


def read_simulation_config(config_file: str, verbose: bool = True) -> SimulationConfig:
    """Parse a simulation control file into a :class:`SimulationConfig`.

    Parameters
    ----------
    config_file : str
        Path to the control file (e.g. ``md.ini``). Must contain an
        ``[OPTIONS]`` section. Inline comments starting with ``#`` or ``;`` are
        ignored. Underscores in integer counts (e.g. ``10_000``) are allowed.
    verbose : bool, optional (default: True)
        If True, echo each parsed setting to stdout (the original behaviour of
        the inline parser). Set False for a quiet read (e.g. in tests).

    Returns
    -------
    SimulationConfig
        All settings with OpenMM units applied where relevant.

    Notes
    -----
    Behavioural details preserved from the original inline parser
    (:class:`cosmo.dynamics.Dynamics.read_config`):

    - ``ref_t`` / ``tau_t`` get units only when ``tcoupl`` is on.
    - ``box_dimension`` is the raw JSON value (a scalar ``L`` for a cubic box or
      an ``[x, y, z]`` list), and is ``None`` when ``pbc`` is off.
    - ``pcoupl`` requires ``pbc`` (asserted).
    - ``restart`` forces ``minimize`` off.
    """
    def log(msg: str) -> None:
        if verbose:
            print(msg)

    cfg = SimulationConfig(config_file=config_file)

    log('-' * 70)
    log(f"Reading simulation parameters from {config_file} file...")
    config = configparser.ConfigParser(inline_comment_prefixes=("#", ";"))
    config.read(config_file)
    params = config['OPTIONS']

    cfg.md_steps = int(str(params.get('md_steps', cfg.md_steps)).replace('_', ''))
    log(f'  md_steps: {cfg.md_steps}')
    cfg.dt = float(params.get('dt', 0.01)) * unit.picoseconds
    log(f'  dt: {cfg.dt}')
    cfg.nstxout = int(params.get('nstxout', cfg.nstxout))
    log(f'  nstxout: {cfg.nstxout}')
    cfg.nstlog = int(params.get('nstlog', cfg.nstlog))
    log(f'  nstlog: {cfg.nstlog}')
    cfg.nstcomm = int(params.get('nstcomm', cfg.nstcomm))
    log(f'  nstcomm: {cfg.nstcomm}')
    cfg.model = params.get('model', cfg.model)
    log(f'  model: {cfg.model}')

    # Log formatting: precision/width accept an int or the word 'none'/'' to
    # disable (None) -- precision=None keeps full repr precision, width=None
    # disables fixed-width alignment.
    prec_val = params.get('log_precision', None)
    if prec_val is not None:
        cfg.log_precision = None if str(prec_val).strip().lower() in ('none', '') else int(prec_val)
    width_val = params.get('log_width', None)
    if width_val is not None:
        cfg.log_width = None if str(width_val).strip().lower() in ('none', '') else int(width_val)
    log(f'  log columns: precision={cfg.log_precision if cfg.log_precision is not None else "full"}, '
        f'width={cfg.log_width if cfg.log_width is not None else "auto"}')

    cfg.tcoupl = bool(strtobool(str(params.get('tcoupl', cfg.tcoupl))))
    if cfg.tcoupl:
        cfg.ref_t = float(params['ref_t']) * unit.kelvin
        cfg.tau_t = float(params['tau_t']) / unit.picoseconds
        log(f'  tcoupl: on (ref_t={cfg.ref_t}, tau_t={cfg.tau_t})')
    else:
        log('  tcoupl: off')

    cfg.pbc = bool(strtobool(str(params.get('pbc', cfg.pbc))))
    if cfg.pbc:
        cfg.box_dimension = loads(params['box_dimension'])
        log(f'  pbc: on (box_dimension={cfg.box_dimension} nm)')
    else:
        cfg.box_dimension = None
        log('  pbc: off')

    cfg.pcoupl = bool(strtobool(str(params.get('pcoupl', cfg.pcoupl))))
    if cfg.pcoupl:
        assert cfg.pbc, ("Pressure coupling requires box dimensions and periodic "
                         "boundary condition is on")
        cfg.ref_p = float(params['ref_p']) * unit.bar
        cfg.frequency_p = int(params['frequency_p'])
        log(f'  pcoupl: on (ref_p={cfg.ref_p}, frequency={cfg.frequency_p})')
    else:
        log('  pcoupl: off')

    cfg.pdb_file = params['pdb_file']
    log(f'  pdb_file: {cfg.pdb_file}')

    # Output location (topo-style): all files go to <output_dir>/<outname>.*.
    # `protein_code` is the legacy prefix; if given and `outname` is not, it is
    # used as `outname` so existing control files keep their basename.
    cfg.protein_code = params.get('protein_code', cfg.protein_code)
    cfg.output_dir = params.get('output_dir', cfg.output_dir)
    if 'outname' in params:
        cfg.outname = params['outname']
    elif cfg.protein_code:
        cfg.outname = cfg.protein_code
    log(f'  output: {cfg.output_path("")}.* (dir={cfg.output_dir}, name={cfg.outname})')
    cfg.checkpoint = params.get('checkpoint', cfg.checkpoint)

    cfg.device = params.get('device', cfg.device)
    log(f'  device: {cfg.device}')
    if cfg.device == 'CPU':
        cfg.ppn = int(params.get('ppn', cfg.ppn))
        log(f'  threads: {cfg.ppn}')

    cfg.restart = bool(strtobool(str(params.get('restart', cfg.restart))))
    log(f'  restart: {cfg.restart}')
    if cfg.restart:
        # A restart resumes from the checkpoint state, never re-minimizes.
        cfg.minimize = False
    else:
        cfg.minimize = bool(strtobool(str(params.get('minimize', cfg.minimize))))
        log(f'  minimize: {cfg.minimize}')

    log('-' * 70)
    return cfg
