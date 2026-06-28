# COSMO Tutorials

Hands-on, ready-to-run tutorials for **COSMO** — a coarse-grained (one bead per
residue) simulation engine for **intrinsically disordered proteins** and related
biomolecules, built on [OpenMM](https://openmm.org/).

Each subfolder is **one self-contained example**: it ships the input files you
need, and its `README.md` walks you through the run step by step and explains the
concepts involved. Work through them in order.

| # | Tutorial | What you learn |
|---|----------|----------------|
| 1 | [Single-chain quickstart](./01_single_chain_quickstart/) | The minimal workflow: a config file, one PDB, run an MD simulation, read the outputs. |
| 2 | [Models & force fields](./02_models_and_forcefields/) | The four force fields (`hps_urry`, `hps_kr`, `hps_ss`, `mpipi`), when to use each, and seeing the difference in per-force-group energies. |
| 3 | [Periodic box, temperature & pressure](./03_pbc_temperature_pressure/) | Running in a periodic box; the NVT and NPT (barostat) ensembles. |
| 4 | [Restart & outputs](./04_restart_and_outputs/) | Continuing a run from a checkpoint, and a tour of every output file. |
| 5 | [Slab simulation of phase separation (LLPS)](./05_slab_llps/) | The slab method for condensates: NPT-compress, elongate, NVT, and read coexisting densities off ρ(z) to map a phase diagram. |
| 6 | [Protein–RNA complex](./06_protein_rna_complex/) | Multi-component systems: assembling a mixed protein + RNA input and the nucleic-acid-capable models. |

The **ready-to-run files** for each tutorial (PDB, `md.ini`, `run_simulation.py`)
live in the matching folder here under `tutorials/`.

---

## What is COSMO? (the 1-minute version)

COSMO turns a biomolecule structure into a **one-bead-per-residue** coarse-grained
model and simulates it in OpenMM. For a protein each residue becomes a single bead
at the alpha-carbon (Cα); for a nucleic acid, at the phosphate (P). Unlike a
structure-based (Gō-like) model, COSMO is a **transferable hydropathy model**: the
input conformation does *not* define an energy minimum. Each bead carries a size,
a "stickiness", and a charge, and the force field has these terms:

- **Bonds** — harmonic springs along the chain (plus optional **angle/torsion**
  terms in the `hps_ss` model).
- **Electrostatics** — Debye–Hückel (Yukawa) screened Coulomb between charged
  beads.
- **Non-bonded** — the heart of the model: an Ashbaugh–Hatch (LJ-based) or
  Wang–Frenkel (`mpipi`) potential whose attractive depth is set by per-residue
  hydropathy parameters.

Because nothing pins the chain to a folded state, COSMO is ideal for **disordered
proteins**: single-chain dimensions, and — with many chains in a box —
**liquid–liquid phase separation (LLPS)** and biomolecular condensates.

## How you run it

Every simulation is driven by a plain-text config file (conventionally `md.ini`)
and launched with any of:

```bash
cosmo-mdrun -f md.ini            # console command (after `pip install -e .`)
python -m cosmo.mdrun -f md.ini  # module form, no install needed
python run_simulation.py -f md.ini   # same thing, from inside a tutorial folder
```

The runner reads `md.ini`, builds the coarse-grained model from your PDB
(`cosmo.models.buildCoarseGrainModel`), and runs Langevin dynamics. It lives in the
package as `cosmo.mdrun`; each tutorial folder keeps a tiny `run_simulation.py`
that just calls it (`from cosmo.mdrun import mdrun`), so every example stays
self-contained while the runner has a single canonical implementation.

A full reference for `md.ini` options lives in
[Simulation control options](https://vuqv.github.io/cosmo/usage/simulation_control.html).

---

## Prerequisites

1. **Python environment with COSMO + OpenMM installed.** From the repo root the
   `cosmo` package must be importable:
   ```bash
   python -c "import cosmo, openmm; print('OpenMM', openmm.__version__)"
   ```
   See the top-level [README](../README.md) for installation (conda + OpenMM ≥ 7.7,
   then `pip install -e .`).
2. **(Optional) A GPU.** The tutorials default to `device = CPU` so they run
   anywhere; switch to `device = GPU` in `md.ini` if you have a CUDA device. The
   slab tutorial (5) is large and really wants a GPU.

## Conventions used in the tutorials

- Settings are deliberately **small and fast** (`md_steps = 5000`, `device = CPU`)
  so each example finishes in seconds. They are demos, not production runs — for
  real science you would increase `md_steps` to millions and use a GPU.
- Each run writes all its files to one self-contained folder,
  `<output_dir>/<outname>.*` (default `traj/`, e.g. `traj/asyn.dcd`,
  `traj/asyn.log`, `traj/asyn.psf`). These are **not** committed — you generate
  them by running the tutorial. A re-run overwrites them, so change `outname` /
  `output_dir` (or copy the folder aside) to keep a run.
