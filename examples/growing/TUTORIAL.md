# Growing Nascent Chain Simulation — Tutorial

This tutorial explains how to run **iterative nascent chain growth** simulations with COSMO: the ribosome is kept fixed, and the nascent polypeptide is extended one amino acid at a time, with short MD runs after each addition. The workflow is intended for studying translation and co-translational folding.

---

## Table of contents

1. [Overview](#overview)
2. [Prerequisites](#prerequisites)
3. [Input files](#input-files)
4. [Configuration (`md.ini`)](#configuration-mdini)
5. [Running the simulation](#running-the-simulation)
6. [Logic and code flow](#logic-and-code-flow)
7. [Output files](#output-files)
8. [Optional: cropping the ribosome](#optional-cropping-the-ribosome)
9. [References and parameters](#references-and-parameters)

---

## Overview

The script `grow_nascent.py`:

1. **Reads the target sequence** from a nascent-chain PDB file (`nascent.pdb`).
2. **Starts from a ribosome PDB** (`rib.pdb`), which can be a full or cropped structure.
3. **For each amino acid** in order:
   - Appends one new CA bead to the current complex (at a fixed “exit tunnel” position).
   - Builds a coarse-grained HPS model of ribosome + nascent chain.
   - Runs a short MD segment (number of steps derived from codon translation time and `scale_factor`).
   - Writes the final frame and trajectory; the final structure becomes the input for the next iteration.

So the nascent chain is **grown step-by-step** with relaxation after each new residue, while the ribosome is **frozen** and excluded from internal bonding and from nonbonded ribosome–ribosome interactions.

---

## Prerequisites

- **Python** with COSMO and its dependencies (OpenMM, ParmEd, etc.).
- **Input PDBs**: ribosome structure and a nascent-chain PDB that defines the sequence (see below).
- **Config file**: `md.ini` in the same directory (or path given by `-f`).

From the project root you can run the example from `examples/growing/`:

```bash
cd examples/growing
python grow_nascent.py -f md.ini --rib rib.pdb --nascent nascent.pdb
```

---

## Input files

### 1. Ribosome PDB (`rib.pdb`)

- All-atom or coarse-grained structure of the ribosome (or a cropped region).
- Chain IDs must be set; ribosomal chains are automatically detected and excluded from **bonding** and from **nonbonded** ribosome–ribosome interactions.
- All ribosome atoms are **frozen** during the simulation.

To use a smaller region and speed up the run, you can first crop a full ribosome PDB with the COSMO utility (see [Optional: cropping the ribosome](#optional-cropping-the-ribosome)).

### 2. Nascent chain PDB (`nascent.pdb`)

- Defines the **amino acid sequence** to be grown.
- The script only uses **residue names and numbers** from ATOM/HETATM lines (e.g. columns 17–20 and 22–26); coordinates are not used for the sequence.
- Residues should be in the order in which they are added (N- to C-terminus).
- Typically this file is a CA-only or minimal PDB with one line per residue in the desired order.

### 3. Config file (`md.ini`)

- INI-style config read by `grow_nascent.py` (see [Configuration](#configuration-mdini)).

---

## Configuration (`md.ini`)

Relevant options used by the growing simulation:

| Option | Meaning | Example |
|--------|---------|--------|
| `model` | HPS model | `hps_kr`, `hps_urry`, `hps_ss` |
| `rib_structure` | Ribosome PDB path | `rib.pdb` |
| `nascent_structure` | Nascent sequence PDB path | `nascent.pdb` |
| `organism` | Codon translation time table | `ecoli`, `yeast` |
| `nascent_chain_id` | Chain ID for the grown chain | `8` |
| `scale_factor` | Time scaling (real time → simulation time) | `433129300` |
| `translation_type` | Codon choice per residue | `fast` (fastest codon) or `slow` (slowest) |
| `dt` | Integration timestep (ps) | `0.01` |
| `ref_t` | Reference temperature (K) | `310` |
| `tau_t` | Temperature coupling time (1/ps) | `0.01` |
| `pbc` | Periodic boundary conditions | `no` or `yes` |
| `box_dimension` | Box size (nm), used if `pbc=yes` | `30` or `[30, 30, 60]` |
| `device` | Compute device | `CPU` or `GPU` |
| `ppn` | CPU threads (if `device=CPU`) | `4` |
| `nstxout` | Steps between trajectory/checkpoint writes | `100` |
| `nstlog` | Steps between log lines | `100` |

**Important for growth:**

- **`scale_factor`** and **`organism`** + **`translation_type`** determine how many MD steps are run per residue: the script converts each residue’s codon translation time (from the chosen organism/codon set) to simulation time and then to steps using `dt` and `scale_factor`. So the “length” of each iteration is **codon-dependent**, not a fixed step count.
- **`nascent_chain_id`** must match the chain ID used when appending new residues (default `'8'` in the code).

---

## Running the simulation

### Basic command

```bash
python grow_nascent.py -f md.ini --rib rib.pdb --nascent nascent.pdb
```

### Command-line arguments

| Argument | Short | Default | Description |
|----------|--------|---------|-------------|
| `--config` | `-f` | *(required)* | Config file path (e.g. `md.ini`) |
| `--rib` | | `rib.pdb` | Ribosome PDB file |
| `--nascent` | | `nascent.pdb` | Nascent chain PDB (defines sequence) |
| `--start` | | `1` | First residue index to grow (1-based) |
| `--end` | | *(all)* | Last residue index to grow (1-based, inclusive) |
| `--steps` | | `1000` | *(Currently unused)* Steps per iteration are computed from codon time and `scale_factor`. |

**Examples:**

```bash
# Full sequence from md.ini paths
python grow_nascent.py -f md.ini

# Custom PDB paths
python grow_nascent.py -f md.ini --rib my_rib.pdb --nascent my_nascent.pdb

# Grow only residues 1–50
python grow_nascent.py -f md.ini --start 1 --end 50

# Grow only residues 20 to end (chain will have only 20, 21, ...; no 1-19)
python grow_nascent.py -f md.ini --start 20
```

Outputs are written under the `traj/` directory (see [Output files](#output-files)). **Note:** The run always starts from the ribosome PDB. With `--start 20`, the first iteration appends residue 20 to the ribosome (no residues 1-19). True resume from a previous run would require modifying the script to load the last `complex_l{N}_final.pdb` as the initial complex.

---

## Logic and code flow

### High-level pipeline

1. **Parse config** (`md.ini`) and **command-line** arguments.
2. **Parse nascent sequence** from the nascent PDB via `parse_nascent_sequence()` (residue name and number per ATOM/HETATM line, in order).
3. **Generate mRNA (codon) sequence** with `generate_mRNA_sequence()` using `organism` and `translation_type`; compute **steps per residue** from codon translation times and `scale_factor`.
4. **Analyze ribosome**: build a temporary COSMO model of the ribosome only to get:
   - **Ribosome atom indices** (all frozen),
   - **Ribosome chain IDs** (excluded from bonding),
   - **Nonbonded exclusions** for all ribosome–ribosome pairs.
5. **Loop** over residues from `--start` to `--end`:
   - **Append** one CA bead to the current complex with `append_amino_acid_to_complex()` at position (0, 0, 0) Å (exit tunnel).
   - **Run one iteration** with `run_single_iteration()`:
     - Build HPS model with `cosmo.models.buildHPSModel()` (frozen ribosome, `except_chains`, `nb_exclusions`).
     - Add **position restraint** on the **newest** nascent bead at (0.38, 0, 0) nm.
     - Add **planar wall** at *x* = 0.38 nm for all **other** nascent beads (so they stay on the “translated” side).
     - Integrate with Langevin dynamics; write DCD, log, checkpoint, and final PDB.
   - Set **current complex** to the iteration’s final PDB for the next residue.

### Key functions

- **`append_amino_acid_to_complex(input_pdb, output_pdb, res_name, res_num, chain_id='8', position=(0,0,0))`**  
  Appends a single CA atom (and TER/END) to an existing PDB. Used to add the next residue at the tunnel exit.

- **`run_single_iteration(complex_pdb, output_prefix, frozen_indices, except_chains, nb_exclusions, sim_params, md_steps)`**  
  Builds the HPS system for the current complex, adds restraint and wall forces, runs MD for `md_steps` steps, and writes all outputs under `traj_dir` with the given `output_prefix`.

- **`parse_nascent_sequence(nascent_pdb)`** (in `cosmo.utils`)  
  Returns a list of `(residue_name, residue_number)` from the nascent PDB.

- **`generate_mRNA_sequence(aa_three, translation_type, organism)`** (in `cosmo.utils`)  
  Returns mRNA string and total translation time (ms); used to choose codons and to compute steps per residue.

### Physics / restraints

- **Ribosome**: All ribosome atoms are fixed; ribosome–ribosome nonbonded interactions are excluded to avoid artifacts.
- **Newest residue**: Harmonic position restraint at (0.38, 0, 0) nm to mimic the exit tunnel.
- **Older nascent residues**: Planar wall at *x* = 0.38 nm (penalize *x* &lt; 0.38 nm) so the already-translated part stays on one side of the tunnel.

---

## Output files

All per-iteration outputs are written under **`traj/`** (or the directory set in the code).

For iteration *N* (residue *N*), the prefix is `complex_l{N}`. Typical files:

| File | Description |
|------|-------------|
| `traj/complex_l{N}.pdb` | Initial structure for this iteration (ribosome + nascent up to residue *N*). |
| `traj/complex_l{N}_final.pdb` | Final frame after MD; used as input for iteration *N*+1. |
| `traj/complex_l{N}.dcd` | Trajectory for this iteration. |
| `traj/complex_l{N}.chk` | OpenMM checkpoint. |
| `traj/complex_l{N}.log` | Log (step, time, energies, temperature, etc.). |
| `traj/complex_l{N}.psf` | Topology (PSF) for this iteration. |

The **final structure** of the last iteration is the end point of the growing simulation.

---

## Optional: cropping the ribosome

To reduce system size, you can crop a full ribosome PDB to a region around the exit tunnel before running the growing simulation:

```bash
python -m cosmo.utils.crop_ribosome rib_full.pdb rib.pdb
```

Optional arguments: **radius** (Å, default 40) and **x_threshold** (Å, default 60). The filter keeps atoms with *x* ≥ 0; for 0 ≤ *x* &lt; *x_threshold* it keeps only atoms within **radius** of the *x*-axis (cylinder); for *x* ≥ *x_threshold* it keeps all. Then use the cropped `rib.pdb` in `grow_nascent.py` as above.

---

## References and parameters

- **Codon translation times**: The script uses organism-specific median codon translation times (e.g. from `cosmo.utils.ctf_utils`) to assign simulation time per residue. `translation_type` selects the codon (fast or slow) for each amino acid.
- **Scale factor**: Converts real translation time (e.g. in ms) to simulation time; together with `dt` it sets the number of steps per residue.
- **HPS model**: The ribosome–nascent complex is coarse-grained with the chosen HPS model (e.g. `hps_kr`); see COSMO documentation for model options and box/PBC settings.

For more on the COSMO package and HPS models, see the main COSMO documentation and the references cited there.
