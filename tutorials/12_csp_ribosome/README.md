# Tutorial 12 — Continuous Synthesis on a coarse-grained ribosome (`cosmo-csp`)

**Goal:** synthesize a nascent chain **residue by residue on an explicit
coarse-grained ribosome**, with the elongation cycle split into O'Brien's **three
kinetic sub-stages**. This is the explicit-ribosome half of the Continuous Synthesis
Protocol ported to COSMO's IDP force field (`cosmo.csp`, mirroring `topo.csp`). Its
sibling — the same protocol through an **analytic cylindrical tunnel** — is
[Tutorial 11](../11_csp_cylinder/).

Unlike topo's structure-based `csp` (native contacts, STRIDE, rigid `AllBonds`),
COSMO's port is **sequence-based**: no STRIDE, no `domain.yaml`, no native-contact
map — the nascent chain is a transferable HPS/`mpipi` IDP that grows N→C at the
peptidyl-transferase center (PTC). Everything is **standalone COSMO**: the ribosome
is a truncated CG bead model (`ribosome_trunc.pdb`), no external CHARMM files.

> **This is a proof of concept.** The settings here check that the machinery
> **builds, times, runs and writes outputs** end to end in seconds on CPU. They are
> deliberately tiny (short `L` range, clamped `max_steps_per_stage`) and are **not**
> a physical/scientific validation.

**Prerequisites:** the coarse-grained workflow of
[Tutorial 1](../01_single_chain_quickstart/), the force fields of
[Tutorial 2](../02_models_and_forcefields/), and the explicit-ribosome tutorials
[7–8, 10](../07_translation/) (this is their `cosmo.csp` successor).

Nascent chain: **α-synuclein** (`asyn.pdb`, 140 residues, from Tutorial 1).

## The model

The ribosome enters as a **rigid (mass-0) truncated CG bead model**
(`ribosome_trunc.pdb`, topo P/R/BR rRNA rep, 4576 beads): fixed scenery that
supplies the excluded-volume/electrostatic wall and the **A-/P-site anchors**. Each
elongation cycle is split into O'Brien's three kinetic sub-stages —
**peptidyl-transfer → translocation → tRNA-binding** — with the new residue's
C-terminus restraint switching **A→P** across them, so a residue is delivered to the
A site and ratcheted into the P site. The PTC geometry is always optimized: the new
residue is seeded one peptide bond (0.380 nm) from the C-terminus, EV-clear.

## Files

| File | Role |
|------|------|
| `asyn.pdb` | Full-length nascent protein (only the first `L` residues exist at each length). |
| `ribosome_trunc.pdb` | Truncated CG ribosome (topo P/R/BR rep, 4576 beads) — rigid scenery + A-/P-anchors. |
| `csp.ini` | 3-stage explicit-ribosome config (`model = mpipi`, uniform codon timing, tunnel wall on). |

## Run it

```bash
cd tutorials/12_csp_ribosome

python -m cosmo.csp -f csp.ini      # -> synth_out_csp/L_<L>/stage_<1,2,3>/
#   installed console script:  cosmo-csp -f csp.ini
```

The config defaults to `device = CPU` so it runs anywhere; switch to `device = GPU`
for a real run. Kinetics here use **uniform codon timing** (`codon_times = 0.05`, no
mRNA needed) — the 3-stage split and the seconds→steps mapping still run. To drive it
from a real transcript instead, add `mrna = mrna.txt` (as in
[Tutorial 11](../11_csp_cylinder/)) and drop `codon_times`. `max_steps_per_stage`
clamps each sub-stage's (otherwise ~10⁶-step) dwell down to a traceable test size —
raise it, extend `L_max`, and switch to GPU for production.

### Stitch a movie

```bash
python -m cosmo.csp.movie -o synth_out_csp --ribosome ribosome_trunc.pdb
#   installed console script:  cosmo-csp-movie -o synth_out_csp --ribosome ribosome_trunc.pdb
```

`--ribosome` overlays the static ribosome so you see the chain growing inside it.
Writes a fixed-width VMD-playable movie (`movie.dcd` / `movie.psf` / `movie.tcl`).

## What it produces

- `synth_out_csp/L_001/ … L_008/`, each with `stage_1/2/3/` sub-runs — the
  A→P restraint switch across the three O'Brien sub-stages.
- **Nascent-only trajectories**: each `stage_*/traj.psf`/`.dcd` contains only the
  `L` nascent beads — the 4576-bead ribosome is excluded from the trajectory and
  re-attached as static scenery by the movie stitcher's `--ribosome` overlay.
- `synth_out_csp/dwell_times.dat` — per-residue dwell-time table.

Generated `synth_out_*/` directories are **git-ignored** (bulky
trajectories/checkpoints). Delete them and re-run to reproduce.

## Where this sits in the series

Tutorials 7–10 explored *how to confine the nascent chain*. `cosmo.csp` consolidates
those into one package with codon-resolved O'Brien kinetics and post-elongation
phases. This tutorial is its **explicit coarse-grained ribosome** mode; the
**analytic-tunnel** mode is [Tutorial 11](../11_csp_cylinder/).
