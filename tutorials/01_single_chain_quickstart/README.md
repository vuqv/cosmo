# Tutorial 1 — Single-chain quickstart

**Goal:** run your first COSMO simulation end-to-end and understand the inputs and
outputs. We simulate a single intrinsically disordered protein —
α-synuclein (`asyn`, 140 residues) — with the simplest possible configuration.

**Time:** a few seconds on a CPU.

---

## Files in this folder

| File | Role |
|------|------|
| `asyn.pdb` | Input structure (all-atom PDB; COSMO keeps one bead per residue). |
| `md.ini` | Simulation configuration (steps, temperature, I/O, hardware). |
| `run_simulation.py` | Thin runner wrapper (reads `md.ini`, builds the model, runs MD). |

## Background concepts

- **Coarse-graining.** COSMO reads your all-atom PDB but keeps **one bead per
  residue** — the alpha-carbon (Cα) for proteins, the phosphate (P) for nucleic
  acids. A 140-residue protein → 140 particles.
- **Hydropathy-scale model.** The default `hps_urry` is a *hydropathy scale*
  (HPS) model: each residue carries a size (σ) and a "stickiness" (λ) from the
  Urry hydrophobicity scale. Non-bonded interactions use an Ashbaugh–Hatch
  potential plus Debye–Hückel electrostatics between charged residues. Unlike a
  structure-based model, **the input conformation does not define the energy
  minimum** — the chain is free to sample whatever ensemble the force field
  favors, which is exactly what you want for a disordered protein.
- **No periodic box.** A single isolated chain is run with `pbc = no`: there is
  no box and no solvent, so the chain explores its conformational ensemble in
  open space.

> **Want the full theory?** Every energy term (bonds, Ashbaugh–Hatch /
> Wang–Frenkel non-bonded, Debye–Hückel electrostatics, and the optional bonded
> angle/torsion terms), with constants and parameter sources, is documented in
> the [model reference](https://vuqv.github.io/cosmo/). The four force fields are
> compared hands-on in **Tutorial 2**.

## Step-by-step

### 1. Check your environment
From this folder:
```bash
python -c "import cosmo, openmm; print('OK', openmm.__version__)"
```
This must succeed (see the [tutorials overview](../README.md) for setup).

### 2. Look at `md.ini`
Open `md.ini`. The important lines:
```ini
md_steps = 5000      # how long to run (short, for a demo)
model = hps_urry     # the force field (Tutorial 2 compares the four)
ref_t = 300          # temperature in Kelvin
pbc = no             # no periodic box (single chain, no solvent)
pdb_file = asyn.pdb
output_dir = traj    # all outputs go to this folder ...
outname = asyn       # ... named asyn.* (so: traj/asyn.dcd, traj/asyn.log, ...)
device = CPU         # runs anywhere; switch to GPU if you have CUDA
minimize = yes       # relax the input coordinates before stepping
```

### 3. Run it
```bash
python run_simulation.py -f md.ini
```
(equivalently `python -m cosmo.mdrun -f md.ini`, or `cosmo-mdrun -f md.ini` after
`pip install -e .`). COSMO echoes every parsed setting, builds the model (number
of chains, each force term added in order, force-field dump), minimizes, then
steps the dynamics. It ends with `Finished in … seconds`.

### 4. Inspect the outputs
All generated files land in one run folder, `output_dir` (here `traj/`, created
automatically), named `<outname>.*` (here `asyn.*`):

| File | What it is |
|------|------------|
| `traj/asyn.log` | Fixed-width, space-aligned energy/temperature log (one line every `nstlog` steps). |
| `traj/asyn.dcd` | Trajectory (coordinates every `nstxout` steps) — open with VMD/MDAnalysis. |
| `traj/asyn.chk` | Binary checkpoint (positions + velocities) for restarting (Tutorial 4). |
| `traj/asyn.psf` | Topology of the coarse-grained model (load alongside the DCD in analysis tools). |
| `traj/asyn_init.pdb` | The coarse-grained starting structure (one bead per residue). |
| `traj/asyn_final.pdb` | Last conformation; reuse it to seed a follow-up run. |
| `traj/asyn_ff.dat` | Per-residue σ/ε/charge force-field table (HPS models only). |
| `traj/asyn_runinfo.log` | Run provenance: software versions, hardware, GPU, timing. |

Peek at the log:
```bash
head traj/asyn.log
```
The columns are step, time (ps), potential / kinetic / total energy (kJ/mol),
temperature (K), speed, and remaining time. A stable temperature near 300 K and
a non-exploding potential energy mean the run is healthy.

## Try next

- Open the trajectory in VMD: `vmd traj/asyn.psf traj/asyn.dcd`. Watch the
  disordered chain expand and collapse — there is no single folded state.
- Bump `md_steps` to `100000` for a longer ensemble, and compute the radius of
  gyration (R_g) over the trajectory with MDAnalysis.
- Move on to **Tutorial 2** to see how the four force fields differ on this same
  chain.

> Tip: a re-run overwrites the `traj/asyn.*` files. To keep a run, copy the
> folder aside or change `outname`/`output_dir` (e.g. `outname = asyn_T300`, or
> `output_dir = runs/T300`).
