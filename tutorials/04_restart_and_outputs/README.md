# Tutorial 4 — Restart & outputs

**Goal:** continue a simulation from its checkpoint (the normal way to extend a
run or recover from a crash), and get a complete tour of **every file COSMO
writes** and how to read it.

**Time:** a few seconds on a CPU.

---

## Files in this folder

| File | Role |
|------|------|
| `asyn.pdb` | The Tutorial-1 chain. |
| `md.ini` | First leg: runs 2000 steps. |
| `md_restart.ini` | Continuation: `restart = yes`, total target 4000 steps. |
| `run_simulation.py` | Thin runner wrapper. |

## How restart works

COSMO checkpoints positions **and** velocities to `<output_dir>/<outname>.chk`
(here `traj/asyn.chk`) every `nstxout` steps. On a restart the runner:

1. loads the checkpoint and reads the step count already completed
   (`getStepCount()` — this needs OpenMM ≥ 7.7),
2. interprets `md_steps` as the **total target** and runs only the
   *remaining* steps,
3. **appends** new frames to the existing `.dcd` and `.log` (it does not start
   them over).

So `restart = yes` with `md_steps = 4000`, after a first leg of 2000 steps, runs
2000 more. `restart = yes` also forces `minimize = no` — the coordinates come
from the checkpoint, not the PDB.

## Step-by-step

### 1. Run the first leg
```bash
python run_simulation.py -f md.ini
```
This runs 2000 steps and writes `traj/asyn.log`, `traj/asyn.dcd`, `traj/asyn.chk`,
plus the provenance files. Note the final step number at the bottom of
`traj/asyn.log`.

### 2. Continue from the checkpoint
```bash
python run_simulation.py -f md_restart.ini
```
The build log shows it loading `traj/asyn.chk`, reporting 2000 steps already done,
and running the remaining 2000. When it finishes, `traj/asyn.log` and
`traj/asyn.dcd` have **grown** — count the lines / frames to confirm they were
appended, not overwritten:
```bash
wc -l traj/asyn.log     # ~twice the first-leg line count (minus the header)
```

## Every output file, explained

With `output_dir = traj` and `outname = asyn`, a run writes:

| File | When | What it is / how to use it |
|------|------|----------------------------|
| `traj/asyn.log` | during run, every `nstlog` steps | Fixed-width, space-aligned table (`#`-commented header): step, time (ps), potential/kinetic/total energy (kJ/mol), temperature (K), speed (ns/day), remaining time. Parse with `cosmo.reporter.readOpenMMReporterFile`; plot it to check stability. |
| `traj/asyn.dcd` | during run, every `nstxout` steps | Binary trajectory (coordinates only). Load with `traj/asyn.psf` in VMD/MDAnalysis. |
| `traj/asyn.chk` | during run, every `nstxout` steps | Binary checkpoint (positions + velocities). The input to a restart. |
| `traj/asyn.psf` | at build time | CHARMM-style topology of the coarse-grained model (beads, bonds, charges). Pair it with the DCD for analysis. |
| `traj/asyn_init.pdb` | at build time | The coarse-grained **starting** structure (one bead per residue). |
| `traj/asyn_final.pdb` | at finalize | The **last** conformation. Reuse it as the input PDB to seed a fresh run. |
| `traj/asyn_ff.dat` | at build time (HPS models only) | Per-residue σ / ε / charge table actually used by the force field — handy for verifying parameters. |
| `traj/asyn_runinfo.log` | at setup + finalize | Run provenance (INI-style): software versions, hardware, GPU/CUDA, planned/final steps, and wall-clock timing. Pure side channel for reproducibility. |

> **Analysis quickstart (MDAnalysis):**
> ```python
> import MDAnalysis as mda
> u = mda.Universe("traj/asyn.psf", "traj/asyn.dcd")
> print(len(u.atoms), "beads,", len(u.trajectory), "frames")
> ```

## Try next

- Restart a *third* time by bumping `md_steps` again (e.g. 6000) — the trajectory
  keeps growing.
- Seed a brand-new run from a previous run's last frame: set
  `init_position = traj/asyn_final.pdb` (keep `pdb_file` for the topology,
  `restart = no`). Unlike a restart, this starts a fresh step clock with new
  Boltzmann velocities but from those coordinates.
- Move to **Tutorial 5** for a multi-chain slab simulation of phase separation.
