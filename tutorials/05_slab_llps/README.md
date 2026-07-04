# Tutorial 5 — Slab simulation of phase separation (LLPS)

**Goal:** measure whether — and at what temperature — an IDP undergoes
**liquid–liquid phase separation (LLPS)** using the standard **slab method**: pack
many copies of a chain into an elongated box so the system separates into a
**dense** (condensate) and a **dilute** (dispersed) phase coexisting across a flat
interface. This is the workhorse setup of the IDP condensate field.

**Time:** the demo `md_npt.ini` leg is short, but slab systems are large
(thousands of beads, many chains) — **a GPU is strongly recommended** for real
runs. On CPU the demo is slow; treat it as a mechanism walk-through.

---

## Files in this folder

| File | Role |
|------|------|
| `NON.pdb` | A multi-chain starting system (many copies of the IDP in one box). |
| `md_npt.ini` | Stage 1: NPT compression to a sensible density. |
| `run_simulation.py` | Thin runner wrapper. |

## The slab method, in three stages

You do not know the right box size a priori, so you build it in stages
(staged build protocol):

1. **Start big, then compress (NPT).** Put all chains in a generously large cubic
   box and run **NPT** (`pcoupl = yes`) so the Monte Carlo barostat shrinks the
   box to a reasonable, condensed density. This is `md_npt.ini`.
2. **Elongate one axis.** Take the equilibrated *x = y* box edge from stage 1,
   keep `x` and `y`, and **extend `z`** several-fold (e.g. `[L, L, 3L]`). The slab
   of protein now occupies the middle of a long box with empty space on either
   side along `z`.
3. **Hold volume, let it separate (NVT).** Turn the barostat **off**
   (`pcoupl = no`, fixed box) and run a long **NVT** simulation. The chains
   redistribute into a dense slab coexisting with a dilute vapor; the interfaces
   are perpendicular to `z`.

## Step-by-step

### 1. Stage 1 — NPT compression
```bash
python run_simulation.py -f md_npt.ini
```
Watch the box dimensions in the build/run log shrink as the barostat drives the
system toward `ref_p = 1 bar` at `ref_t = 310 K`. Outputs go to `traj/NON.*`.

### 2. Stage 2 — build the elongated box (manual)
Read the final box edge from stage 1 and make a stage-2 config that keeps `x`,`y`
and stretches `z`, e.g.:
```ini
pbc = yes
box_dimension = [L, L, 3L]   ; substitute the equilibrated L from stage 1
pcoupl = no                  ; NVT from here on
restart = yes                ; continue from traj/NON.chk
```
(You start stage 2 from `traj/NON_final.pdb` / `traj/NON.chk`.)

### 3. Stage 3 — production NVT & read off the phase diagram
Run a long NVT simulation, then compute the **density profile along `z`**,
ρ(z). Phase separation shows up as a flat-topped high-density plateau (the
condensate) bracketed by a near-zero baseline (the dilute phase). The two plateau
values are the **coexisting densities** at that temperature.

Repeat the whole procedure across a ladder of `ref_t` values: as you raise the
temperature the dense and dilute densities approach each other and finally merge
at the **critical temperature** — that curve *is* the model's phase diagram.

## Background: reading phase behavior off a slab

- **Why a slab (elongated box)?** The flat interface minimizes curvature/finite-
  size artifacts, so the two coexisting densities are well defined and converge
  quickly compared with a droplet geometry.
- **Concentration & temperature are the knobs.** A single slab at one
  temperature gives one tie-line (two densities). Scanning temperature traces the
  binodal; the system phase-separates below the critical point and stays mixed
  above it.
- **Model choice matters.** The driving force for separation comes straight from
  the hydropathy parameters (Tutorial 2). `mpipi` is tuned for near-quantitative
  phase behavior; the `hps` family is the classic choice. Compare against the
  literature for whichever you pick.

## Try next

- Run stage 1 on a **GPU** (`device = GPU`) for a realistic system size.
- Build the stage-2/stage-3 configs and compute ρ(z) with MDAnalysis (histogram
  bead positions along `z`, average over the trajectory).
- Move to **Tutorial 6** for a multi-component condensate: protein **+ RNA**.
