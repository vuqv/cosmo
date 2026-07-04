# COSMO — Project Instructions

COSMO is a coarse-grained (CA/P) simulation engine for intrinsically disordered
proteins and related biomolecules, built on **OpenMM**. These rules are specific
to this repository and combine with the user's global instructions.

## Sibling package: `topo` is the architectural reference (read this first)

COSMO has a **sibling package, `topo`** (`~/work/src/topo`, also reachable as
`../topo` from this repo), written and maintained by the same person. They are two
halves of one design, split only by the physics they model:

- **`topo`** → coarse-grained MD for **globular / folded proteins** (structure-based,
  Gō-like model).
- **`cosmo`** (this repo) → coarse-grained MD for **intrinsically disordered
  proteins** and related biomolecules (RNA/DNA), using hydropathy-scale (HPS) and
  `mpipi` force fields.

`topo` is the **most recently revised** package and the design the maintainer is
happy with. **Treat `topo` as the source of truth for architecture** — file/module
layout, class and function names, signatures and defaults, the config/engine/runner
split, logging format, CLI, I/O layout, packaging, docs scaffolding, and code style.
The explicit goal is to keep the two packages **structurally identical** so they can
be maintained in parallel.

**The rule — when something is (or could be) shared between the two, prioritize
`topo`:**

1. **Look at `topo` first.** Before adding or refactoring anything cross-cutting in
   cosmo (runner, `SimulationConfig`/`read_simulation_config`, `engine`, reporter,
   `mdrun` CLI, output layout, packaging, docs), open the corresponding file in
   `topo` and mirror it: same module/function names, signatures, defaults, control
   flow, and comment style.
2. **Adopt `topo`'s names even when cosmo's legacy name differs** — e.g.
   `buildCoarseGrainModel` (the old `buildHPSModel` alias has been removed), the
   fixed-width reporter (`precision`/`width`, two-space separator),
   `output_dir`/`outname` I/O, the `build_system` / `setup_simulation` /
   `attach_reporters` / `finalize_simulation` engine API.
3. **Diverge only for genuine physics differences** between IDP and folded proteins.
   Cosmo keeps its own force fields and topology handling: HPS/`mpipi` potentials,
   CA **and** P beads (nucleic-acid support), flexible harmonic bonds. Do **not**
   import topo's structure-based physics (rigid `AllBonds` constraints, native Gō
   contacts, STRIDE, `domain.yaml` scaling) into cosmo, and do **not** strip cosmo's
   IDP-specific physics just to look like topo.
4. **When they drift on something that should be shared, cosmo is the side that
   changes.** Surface the gap rather than letting the two diverge.

If you can't tell whether a difference is "shared architecture" (→ follow `topo`) or
"genuine physics" (→ keep cosmo's), ask before proceeding.

## Architecture

The package is organized around three core classes in `cosmo/core/`:

- `geometry` — geometric parameters from input structures.
- `models` — predefined model builders. `models.buildCoarseGrainModel(...)` is the main
  entry point and is exposed as `cosmo.models.buildCoarseGrainModel` (matching the
  sibling `topo` project's name). The legacy `models.buildHPSModel` alias has been
  removed; call `buildCoarseGrainModel`.
- `system` — the workhorse: defines/modifies forces and creates the OpenMM
  `System`. Forces are registered in the `forceGroups` `OrderedDict`
  (name → force); each force gets a unique force-group index in insertion order.

The simulation runner layer (mirrors the sibling `topo` project's design):

- `cosmo/utils/config.py` — `SimulationConfig` dataclass + `read_simulation_config()`.
  The single place `md.ini` is parsed; carries helpers (`build_kwargs`,
  `output_path`, `checkpoint_path`, `make_platform`). Exposed as
  `cosmo.read_simulation_config` / `cosmo.SimulationConfig`.
- `cosmo/engine.py` — reusable `build_system` / `setup_simulation` /
  `attach_reporters` / `finalize_simulation`. `setup_simulation` also handles
  `init_position` (fresh-start coords) and writes the `_runinfo.log` provenance
  header; `finalize_simulation` closes it.
- `cosmo/mdrun/` — the thin CLI orchestration (`cosmo-mdrun` / `python -m cosmo.mdrun`).

Supporting modules: `cosmo/parameters/model_parameters.py` (add new models here),
`cosmo/reporter/` (`cosmoReporter` — the fixed-width `.log` StateDataReporter that
mirrors topo's, with optional per-force-group energy columns),
`cosmo/utils/runinfo.py` (provenance `<outname>_runinfo.log`, mirrors topo's;
exposed as `cosmo.runinfo`), `cosmo/utils/`.

> There is no `Dynamics` class or `cosmo-simulation.py` anymore — that monolithic
> path was replaced by config + engine + mdrun. Don't reintroduce inline `md.ini`
> parsing in scripts; use `cosmo.read_simulation_config`.

## Supported models

`hps_urry` (default), `hps_kr`, `hps_ss` (adds Gaussian angle + torsion bonded
terms), `mpipi` (Wang–Frenkel instead of Ashbaugh–Hatch / LJ). Define new models
in `cosmo/parameters/model_parameters.py`.

## Running a simulation

Simulations are driven by a config file, not hardcoded params. Canonical runner:

```bash
cosmo-mdrun -f md.ini            # console command (after `pip install -e .`)
python -m cosmo.mdrun -f md.ini  # module form, no install needed
```

`md.ini` (`[OPTIONS]` section) controls model, integrator, T/P coupling, PBC,
I/O, and device (CPU/GPU). See `tutorials/01_single_chain_quickstart/md.ini`. The per-tutorial
`run_simulation.py` files are thin wrappers around `cosmo.mdrun.mdrun`.

## Regression benchmark (run before and after non-trivial changes)

`benchmarks/` guards the physics of every force field. It records the initial
(un-minimized) potential energy of `benchmarks/asyn.pdb`, broken
down by force group, for all four models.

```bash
python benchmarks/benchmark_energies.py            # check vs reference (exit 1 on mismatch)
python benchmarks/benchmark_energies.py --generate # refresh reference — ONLY on known-correct code
```

Rules:
- Run the check after any change touching `cosmo/core/system.py`,
  `cosmo/core/models.py`, or `cosmo/parameters/`.
- The check is deterministic by design (`minimize=False`, no PBC, OpenMM
  `Reference` platform) and must reproduce the reference with `|diff| = 0`.
- A mismatch is a regression **unless** you intentionally changed that force
  field; in that case validate the new physics, then `--generate` to update
  `benchmarks/reference_energies.json` and commit it in the same change.

## Data safety

- Never overwrite raw inputs (PDBs, `md.ini`, example data). Write outputs to new
  files/dirs.
- Treat `benchmarks/reference_energies.json` as committed ground truth; only
  regenerate it deliberately from known-correct code.

## Conventions

- Match the surrounding code's style (NumPy-style docstrings, explicit OpenMM
  `unit` usage, `print`-based progress logging in the build pipeline).
- Prefer config-driven behavior over new hardcoded constants.
- Keep force-group names consistent between `system.addSystemForces` and the
  reporter — the benchmark and reporter key off these strings.

### Software design
Keep cosmo's architecture aligned with the `topo` sibling package as closely as
possible. For any shared functionality or behavior, prioritize `topo`'s style — see
**"Sibling package: `topo` is the architectural reference"** at the top of this file
for the full rule.

### Claude behavior
After each task, always review the relevant information files (PLAN, TASK, README, ...) and update them if relevant.