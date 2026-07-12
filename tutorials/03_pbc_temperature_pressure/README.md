# Tutorial 3 — Periodic box, temperature & pressure

**Goal:** move from an isolated chain to a **periodic box**, and learn the two
ensembles COSMO supports out of the box: **NVT** (fixed volume) and **NPT**
(fixed pressure, with a barostat that resizes the box). This is the foundation
for any multi-chain / condensate run (Tutorials 5–6).

**Time:** a few seconds each on a CPU.

---

## Files in this folder

| File | Role |
|------|------|
| `asyn.pdb` | The Tutorial-1 chain again. |
| `md_nvt.ini` | NVT: temperature coupling on, pressure coupling off. |
| `md_npt.ini` | NPT: temperature **and** pressure coupling on. |
| `run_simulation.py` | Thin runner wrapper. |

## Background concepts

- **Periodic boundary conditions (`pbc = yes`).** The chain sits in a box that
  tiles space; a bead leaving one face re-enters the opposite face. Set the size
  with `box_dimension` — a scalar `L` for a cubic box, or `[x, y, z]` for a
  rectangular one (nanometers). PBC is **required** for pressure coupling and for
  any multi-chain density study.
- **NVT (canonical).** `tcoupl = yes`, `pcoupl = no`. The Langevin thermostat
  holds temperature at `ref_t`; the box volume is fixed. This is the default
  ensemble for sampling at a chosen concentration.
- **NPT (isobaric).** `tcoupl = yes`, `pcoupl = yes`. COSMO adds a **Monte Carlo
  barostat** automatically (you do **not** write any extra code) that rescales
  the box to hold the pressure at `ref_p`, attempting a move every
  `frequency_p` steps. Use it to relax a system to a sensible density before a
  production run. **`pcoupl = yes` requires `pbc = yes`** — the runner asserts
  this.

## Step-by-step

### 1. Run NVT
```bash
cosmo-mdrun -f md_nvt.ini                # == python -m cosmo.mdrun -f md_nvt.ini
python run_simulation.py -f md_nvt.ini   # equivalent in-folder shim
```
Outputs go to `traj/asyn_nvt.*`. The box stays 30 nm on a side the whole run.

### 2. Run NPT
```bash
cosmo-mdrun -f md_npt.ini                # == python -m cosmo.mdrun -f md_npt.ini
python run_simulation.py -f md_npt.ini   # equivalent in-folder shim
```
Outputs go to `traj/asyn_npt.*`. In the build log you'll see the barostat being
added (it is **not** present in the NVT run). Over the run the box dimensions
drift as the barostat equilibrates the volume toward `ref_p = 1 bar`.

### 3. Compare
Both runs write the standard `*.log` (see Tutorial 1 for the columns). The NPT
log reflects a system whose volume is being adjusted; the NVT log is at fixed
volume. For a single dilute chain the energetic difference is small — the point
here is the **mechanism**: NPT is how you find a reasonable box/density, after
which you typically switch to NVT for production (exactly the slab recipe in
Tutorial 5).

## Try next

- Shrink the box (`box_dimension = 15`) and watch how PBC + a small box start to
  matter as the chain "sees" its periodic image.
- Use a rectangular box, e.g. `box_dimension = [30, 30, 60]`, the elongated
  geometry used for slab/LLPS simulations (Tutorial 5).
- Move to **Tutorial 4** to restart a run from its checkpoint and tour every
  output file in detail.
