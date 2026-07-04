#!/usr/bin/env python
# coding: utf-8
"""
Regression benchmark for COSMO force fields.

For each supported model force field, this script builds a coarse-grained model
from a reference structure (``asyn.pdb`` by default), evaluates the potential
energy of the *initial* (un-minimized) coordinates, and decomposes it into the
individual force groups (Harmonic Bond, Yukawa, PairWise/PairWire, etc.).

The recorded energies act as a known-correct reference. After changing the code,
re-run in check mode to confirm that every model still produces the same initial
energies (within tolerance). A mismatch flags an unintended change in the
physics, not just the API.

Determinism notes
-----------------
* The model is built with ``minimize=False`` so the coordinates come straight
  from the input PDB (no platform-dependent minimization).
* Periodic boundary conditions are disabled so nonbonded forces use a fixed
  non-periodic cutoff -- the energy depends only on the input coordinates.
* Energies are evaluated on OpenMM's ``Reference`` platform (double precision)
  for reproducibility across machines and GPUs.

Usage
-----
Generate / refresh the reference (only do this when the current code is known
to be correct)::

    python benchmark_energies.py --generate

Check the current code against the stored reference (default)::

    python benchmark_energies.py

Options::

    --pdb PATH        input structure (default: ./asyn.pdb, alpha-synuclein)
    --reference PATH  reference JSON (default: ./reference_energies.json)
    --models M [M..]  subset of models to test (default: all supported)
    --rtol FLOAT      relative tolerance (default: 1e-4)
    --atol FLOAT      absolute tolerance in kJ/mol (default: 1e-3)
"""

import argparse
import json
import os
import sys
import warnings
from collections import OrderedDict

import openmm as mm
from openmm import unit

try:
    from parmed.exceptions import OpenMMWarning

    warnings.filterwarnings("ignore", category=OpenMMWarning)
except Exception:  # pragma: no cover - parmed always present alongside cosmo
    pass

# Allow running directly from the benchmarks/ directory without installing.
_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO_ROOT = os.path.dirname(_HERE)
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

import cosmo  # noqa: E402

# Models exercised by the benchmark. All of these support a single-chain
# protein such as alpha-synuclein (asyn.pdb).
DEFAULT_MODELS = ["hps_urry", "hps_kr", "hps_ss", "mpipi"]

DEFAULT_PDB = os.path.join(_HERE, "asyn.pdb")
DEFAULT_REFERENCE = os.path.join(_HERE, "reference_energies.json")


def compute_energies(pdb_file, model_name):
    """Build a model and return its initial per-force-group potential energies.

    Parameters
    ----------
    pdb_file : str
        Path to the input structure.
    model_name : str
        One of the supported COSMO models.

    Returns
    -------
    OrderedDict
        Mapping of ``force-group name -> energy (kJ/mol)`` plus a ``Total``
        entry. Order matches the force-group indices in the system.
    """
    # minimize=False and box_dimension=None keep the calculation deterministic:
    # the energy is a pure function of the input coordinates.
    hps_model = cosmo.models.buildCoarseGrainModel(
        pdb_file, minimize=False, model=model_name, box_dimension=None
    )

    platform = mm.Platform.getPlatformByName("Reference")
    integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
    context = mm.Context(hps_model.system, integrator, platform)
    context.setPositions(hps_model.positions)

    energies = OrderedDict()
    for i, name in enumerate(hps_model.forceGroups):
        e = context.getState(getEnergy=True, groups={i}).getPotentialEnergy()
        energies[name] = e.value_in_unit(unit.kilojoules_per_mole)

    total = context.getState(getEnergy=True).getPotentialEnergy()
    energies["Total"] = total.value_in_unit(unit.kilojoules_per_mole)

    # Release the context/integrator before building the next model.
    del context, integrator
    return energies


def generate(pdb_file, models, reference_path):
    """Compute energies for every model and write them to the reference file."""
    result = OrderedDict()
    result["_meta"] = {
        "pdb_file": os.path.relpath(pdb_file, _HERE),
        "openmm_version": mm.__version__,
        "note": "Initial (un-minimized) potential energy by force group, no PBC, "
        "Reference platform. Regenerate only from known-correct code.",
    }
    for model_name in models:
        print(f"\n=== Generating reference for model: {model_name} ===")
        result[model_name] = compute_energies(pdb_file, model_name)

    with open(reference_path, "w") as fh:
        json.dump(result, fh, indent=2)
    print(f"\nWrote reference energies to {reference_path}")
    _print_table(result, models)


def check(pdb_file, models, reference_path, rtol, atol):
    """Compare current energies against the stored reference.

    Returns
    -------
    bool
        True if all force groups for all models match within tolerance.
    """
    if not os.path.exists(reference_path):
        print(
            f"ERROR: reference file not found: {reference_path}\n"
            f"Run with --generate first (from known-correct code)."
        )
        return False

    with open(reference_path) as fh:
        reference = json.load(fh)

    all_ok = True
    for model_name in models:
        print(f"\n=== Checking model: {model_name} ===")
        if model_name not in reference:
            print(f"  MISSING from reference -- run --generate to add it. FAIL")
            all_ok = False
            continue

        current = compute_energies(pdb_file, model_name)
        ref = reference[model_name]

        # Union of keys so we catch added/removed force groups too.
        keys = list(OrderedDict.fromkeys(list(ref.keys()) + list(current.keys())))
        for key in keys:
            if key not in ref:
                print(f"  [{key}] NEW force group not in reference. FAIL")
                all_ok = False
                continue
            if key not in current:
                print(f"  [{key}] MISSING force group (was in reference). FAIL")
                all_ok = False
                continue

            ref_val = ref[key]
            cur_val = current[key]
            diff = abs(cur_val - ref_val)
            tol = atol + rtol * abs(ref_val)
            status = "ok" if diff <= tol else "FAIL"
            if status == "FAIL":
                all_ok = False
            print(
                f"  [{key:24s}] ref={ref_val:16.6f}  cur={cur_val:16.6f}  "
                f"|diff|={diff:.3e}  tol={tol:.3e}  {status}"
            )

    print("\n" + ("ALL BENCHMARKS PASSED" if all_ok else "BENCHMARK MISMATCH DETECTED"))
    return all_ok


def _print_table(result, models):
    print("\nSummary (kJ/mol):")
    for model_name in models:
        print(f"\n  {model_name}")
        for key, val in result[model_name].items():
            print(f"    {key:24s} {val:16.6f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--generate", action="store_true",
                        help="generate/refresh the reference file (use only on known-correct code)")
    parser.add_argument("--pdb", default=DEFAULT_PDB, help="input structure PDB")
    parser.add_argument("--reference", default=DEFAULT_REFERENCE, help="reference JSON path")
    parser.add_argument("--models", nargs="+", default=DEFAULT_MODELS,
                        help="subset of models to test")
    parser.add_argument("--rtol", type=float, default=1e-4, help="relative tolerance")
    parser.add_argument("--atol", type=float, default=1e-3,
                        help="absolute tolerance in kJ/mol")
    args = parser.parse_args()

    print(f"OpenMM version: {mm.__version__}")
    print(f"Input structure: {args.pdb}")

    if args.generate:
        generate(args.pdb, args.models, args.reference)
        return 0

    ok = check(args.pdb, args.models, args.reference, args.rtol, args.atol)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
