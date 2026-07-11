"""
Run provenance / metadata logging.

Writes a human-readable ``*_runinfo.log`` capturing software versions, hardware,
GPU/CUDA details, and timing for a simulation run. This is a pure side channel:
nothing here changes the simulation, it only records *how* and *where* a run was
produced, which is useful for reproducibility and for debugging performance
differences between machines.

Mirrors the sibling ``topo`` project's ``utils/runinfo.py`` so the two packages
write the same provenance format.

Typical use from a run script::

    import time, cosmo
    start = time.time()
    ...                                   # build model, set up `simulation`
    cosmo.runinfo.write_run_start(
        path, control_file=cfg.config_file, checkpoint_file=checkpoint,
        restart=restart_active, steps_planned=nsteps_remain,
        simulation=simulation, use_gpu=(cfg.device == 'GPU'), ppn=cfg.ppn)
    simulation.step(nsteps_remain)
    cosmo.runinfo.write_run_end(path, simulation=simulation, start_epoch=start)

The log is a simple INI-like text file: a one-line ``# COSMO run info`` banner
followed by aligned ``[run]``, ``[software]``, ``[hardware]``, ``[gpu]`` (GPU runs
only) and ``[result]`` sections. The first four are written before the run so a
crashed run still leaves a "what was attempted" file; ``[result]`` is appended at
the end.
"""
import os
import platform as py_platform
import shutil
import socket
import subprocess
import time
from datetime import datetime

import numpy as np
import openmm as mm
from openmm import unit


def _module_version(module) -> str:
    """Best-effort version string for a Python module."""
    return getattr(module, "__version__", "unknown")


def _safe_platform_property(simulation, key: str) -> str:
    """Best-effort OpenMM platform property lookup."""
    try:
        return simulation.context.getPlatform().getPropertyValue(simulation.context, key)
    except Exception:
        return "not available"


def collect_cuda_metadata(simulation) -> dict:
    """
    Collect CUDA/device metadata when available.

    Combines OpenMM platform properties (device name, driver/runtime version,
    precision, device index) with a best-effort ``nvidia-smi`` query for
    friendlier diagnostics. Every lookup degrades gracefully to a descriptive
    "not available" string rather than raising.
    """
    cuda_info = {}

    # OpenMM platform properties (best effort)
    keys = [
        "CudaDeviceName",
        "CudaDriverVersion",
        "CudaRuntimeVersion",
        "CudaCompiler",
        "CudaPrecision",
        "DeviceIndex",
    ]
    for key in keys:
        value = _safe_platform_property(simulation, key)
        if value != "not available":
            cuda_info[key] = value

    # nvidia-smi query (best effort, friendlier diagnostics)
    nvidia_smi = shutil.which("nvidia-smi")
    if nvidia_smi is None:
        cuda_info["nvidia_smi"] = "not available (nvidia-smi not found on PATH)"
    else:
        try:
            proc = subprocess.run(
                [nvidia_smi, "--query-gpu=name,driver_version", "--format=csv,noheader"],
                capture_output=True,
                text=True,
                timeout=5,
            )
            if proc.returncode == 0:
                smi = proc.stdout.strip()
                cuda_info["nvidia_smi"] = smi if smi else "available but no GPU rows returned"
            elif proc.returncode == 2:
                cuda_info["nvidia_smi"] = (
                    "not available (nvidia-smi returned code 2; "
                    "likely no accessible NVIDIA GPU/driver in this environment)"
                )
            else:
                stderr_msg = (proc.stderr or proc.stdout).strip()
                cuda_info["nvidia_smi"] = (
                    f"not available (nvidia-smi exit code {proc.returncode}: {stderr_msg})"
                )
        except subprocess.TimeoutExpired:
            cuda_info["nvidia_smi"] = "not available (nvidia-smi query timed out)"
        except Exception as e:
            cuda_info["nvidia_smi"] = f"not available (nvidia-smi query failed: {e})"

    return cuda_info


_KEY_W = 16  # key column width for aligned "key : value" lines


def _human_time(seconds: float) -> str:
    """Wall-clock as a single human-friendly figure (s / min / h)."""
    if seconds < 90:
        return f"{seconds:.2f} s"
    if seconds < 5400:
        return f"{seconds / 60:.1f} min"
    return f"{seconds / 3600:.2f} h"


def _context_from_path(path: str) -> str:
    """Short banner context from the output dir, e.g. ``.../L_008/stage_3/x.log`` -> ``L_008/stage_3``."""
    parts = [p for p in os.path.dirname(os.path.abspath(path)).split(os.sep) if p]
    return "/".join(parts[-2:]) if parts else ""


def write_banner(path: str, title: str) -> None:
    """(Over)write the file with the ``# COSMO run info`` banner header."""
    header = "COSMO run info" + (f"  ·  {title}" if title else "")
    with open(path, "w") as f:
        f.write(f"# {header}\n")
        f.write(f"# {'=' * len(header)}\n")


def write_section(path: str, section_name: str, kv_pairs: dict, mode: str = "a") -> None:
    """Append one metadata section (``[name]`` header + aligned ``key : value`` lines)."""
    with open(path, mode) as f:
        f.write(f"\n[{section_name}]\n")
        for k, v in kv_pairs.items():
            f.write(f"{k:<{_KEY_W}}: {v}\n")


def write_run_start(path: str, *, control_file, checkpoint_file, restart,
                    steps_planned, simulation, use_gpu, ppn,
                    coord_source=None, vel_source=None, title=None,
                    append=False, section_label=None) -> None:
    """
    Write the banner + ``[run]``/``[software]``/``[hardware]`` (and ``[gpu]``)
    provenance sections, overwriting any existing file.

    ``append`` / ``section_label`` support a **folded, multi-phase** run-info file
    (one file per residue, one section group per stage -- used by the CSP runner). When
    ``append`` is True the banner is not re-written and the run-invariant
    ``[software]``/``[hardware]``/``[gpu]`` sections are skipped (written by the first
    phase); only the ``[<section_label>]`` run section is appended. ``section_label``
    names it (default ``"run"``), so successive stages read as ``[run: stage 2 ...]``.

    Records timing, input/checkpoint paths, planned step count, the source of the
    initial coordinates/velocities, package versions (Python/NumPy/ParmEd/OpenMM),
    hardware (host, OS, CPU, selected platform), and — for GPU runs — a friendly
    accelerator line plus a detailed ``[gpu]`` section.

    Parameters
    ----------
    path : str
        Output run-info log path (e.g. ``"traj/traj_runinfo.log"``).
    control_file : str
        Path to the simulation control file used for this run.
    checkpoint_file : str
        Resolved checkpoint path for the run.
    restart : bool
        Whether the run actually restarted from a checkpoint.
    steps_planned : int
        Number of integration steps this invocation will run.
    simulation : openmm.app.Simulation
        The constructed simulation (used to read the selected platform / GPU info).
    use_gpu : bool
        Whether the run is using the CUDA platform.
    ppn : int
        Requested CPU threads (informational).
    coord_source : str, optional
        Human-readable description of where the initial coordinates came from
        (checkpoint, init_position, or pdb_file).
    vel_source : str, optional
        Human-readable description of where the initial velocities came from
        (checkpoint or a Boltzmann distribution).
    title : str, optional
        Short banner label (e.g. ``"L=8 · stage 3"``). When omitted, a context is
        derived from the output directory (e.g. ``"L_008/stage_3"``).
    """
    try:
        import parmed as pmd
        parmed_version = _module_version(pmd)
    except Exception:
        parmed_version = "not installed"

    if not append:
        write_banner(path, title if title is not None else _context_from_path(path))

    run = {
        "started": datetime.now().astimezone().isoformat(timespec="seconds").replace("T", " "),
        "restart": "yes" if restart else "no",
        "steps_planned": steps_planned,
        "control_file": control_file if control_file is not None else "(none)",
        "checkpoint_file": checkpoint_file if checkpoint_file is not None else "(none)",
        "init_coords": coord_source if coord_source is not None else "unknown",
        "init_velocities": vel_source if vel_source is not None else "unknown",
    }
    write_section(path, section_label or "run", run)

    # Run-invariant provenance: written once (skipped for appended stages).
    if append:
        return

    software = {
        "python": py_platform.python_version(),
        "numpy": _module_version(np),
        "parmed": parmed_version,
        "openmm": _module_version(mm),
    }
    write_section(path, "software", software)

    hardware = {
        "host": socket.gethostname(),
        "os": py_platform.platform(),
        "platform": simulation.context.getPlatform().getName(),
        "cpu": f"{py_platform.machine()} · {os.cpu_count()} cores · {ppn} threads requested",
    }
    cuda = collect_cuda_metadata(simulation) if use_gpu else None
    if cuda:
        hardware["accelerator"] = (f"{cuda.get('CudaDeviceName', 'GPU')} · "
                                   f"{cuda.get('CudaPrecision', '?')} precision")
    write_section(path, "hardware", hardware)

    if cuda:
        write_section(path, "gpu", cuda)


def write_run_end(path: str, *, simulation, start_epoch: float,
                  final_structure=None, section_label=None) -> None:
    """
    Append the ``[result]`` section (status + wall-clock timing + final state).

    ``section_label`` renames the section for a folded multi-phase file (e.g.
    ``"result: stage 2 translocation"``); defaults to ``"result"``.

    Parameters
    ----------
    path : str
        Output run-info log path (same as passed to :func:`write_run_start`).
    simulation : openmm.app.Simulation
        The simulation, queried for final step count and simulated time.
    start_epoch : float
        ``time.time()`` value captured at the start of the run, used to compute
        elapsed wall-clock time.
    final_structure : str, optional
        Path to the PDB of the last conformation written for this run (usable as
        ``init_position`` for a follow-up run).
    """
    elapsed = time.time() - start_epoch
    state = simulation.context.getState()
    info = {
        "status": "completed",
        "wall_clock": _human_time(elapsed),
        "steps_run": state.getStepCount(),
        "sim_time": f"{state.getTime().value_in_unit(unit.nanosecond):.4g} ns",
        "ended": datetime.now().astimezone().isoformat(timespec="seconds").replace("T", " "),
    }
    if final_structure is not None:
        info["final_structure"] = final_structure
    write_section(path, section_label or "result", info)
