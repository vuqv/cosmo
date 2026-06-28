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

The log is a simple INI-like text file with ``[run_start]``, ``[cuda_metadata]``
(GPU runs only), and ``[run_end]`` sections.
"""
import os
import platform as py_platform
import shutil
import socket
import subprocess
import sys
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


def write_section(path: str, section_name: str, kv_pairs: dict, mode: str = "a") -> None:
    """Write one metadata section (``[name]`` header + ``key: value`` lines)."""
    with open(path, mode) as f:
        f.write(f"\n[{section_name}]\n")
        for k, v in kv_pairs.items():
            f.write(f"{k}: {v}\n")


def write_run_start(path: str, *, control_file, checkpoint_file, restart,
                    steps_planned, simulation, use_gpu, ppn,
                    coord_source=None, vel_source=None) -> None:
    """
    Write the ``[run_start]`` provenance section, overwriting any existing file.

    Records timing, input/checkpoint paths, planned step count, the source of the
    initial coordinates/velocities, package versions (Python/NumPy/ParmEd/OpenMM),
    hardware (host, OS, CPU count, selected platform), and — for GPU runs — a
    ``[cuda_metadata]`` section.

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
    """
    try:
        import parmed as pmd
        parmed_version = _module_version(pmd)
    except Exception:
        parmed_version = "not installed"

    info = {
        "start_time_iso": datetime.now().astimezone().isoformat(),
        "control_file": control_file,
        "checkpoint_file": checkpoint_file,
        "restart_mode": restart,
        "steps_planned": steps_planned,
        "initial_coordinates": coord_source if coord_source is not None else "unknown",
        "initial_velocities": vel_source if vel_source is not None else "unknown",
        "python_version": sys.version.replace("\n", " "),
        "numpy_version": _module_version(np),
        "parmed_version": parmed_version,
        "openmm_version": _module_version(mm),
        "hostname": socket.gethostname(),
        "os": py_platform.platform(),
        "machine": py_platform.machine(),
        "processor": py_platform.processor() or "not reported",
        "cpu_count": os.cpu_count(),
        "requested_threads_ppn": ppn,
        "selected_platform": simulation.context.getPlatform().getName(),
        "use_gpu": use_gpu,
    }
    write_section(path, "run_start", info, mode="w")
    if use_gpu:
        write_section(path, "cuda_metadata", collect_cuda_metadata(simulation))


def write_run_end(path: str, *, simulation, start_epoch: float,
                  final_structure=None) -> None:
    """
    Append the ``[run_end]`` provenance section (wall-clock timing + final state).

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
        "end_time_iso": datetime.now().astimezone().isoformat(),
        "elapsed_seconds": f"{elapsed:.2f}",
        "elapsed_hours": f"{elapsed / 3600:.4f}",
        "final_step": state.getStepCount(),
        "final_time_ns": f"{state.getTime().value_in_unit(unit.nanosecond):.6f}",
    }
    if final_structure is not None:
        info["final_structure"] = final_structure
    write_section(path, "run_end", info)
