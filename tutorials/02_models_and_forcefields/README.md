# Tutorial 2 — Models & force fields

**Goal:** understand the four force fields COSMO ships, when to use each, and how
to see — quantitatively — what changes when you swap one for another. We run the
same α-synuclein chain from Tutorial 1 under different models.

**Time:** a few seconds per model on a CPU.

---

## Files in this folder

| File | Role |
|------|------|
| `asyn.pdb` | Same input chain as Tutorial 1. |
| `md.ini` | Configuration with a single line — `model = …` — to change. |
| `run_simulation.py` | Thin runner wrapper. |

## The four models

| Model | Non-bonded form | Hydropathy scale | Adds | Use it for |
|-------|-----------------|------------------|------|------------|
| `hps_urry` | Ashbaugh–Hatch (LJ 12-6) | Urry | — | **Default.** General-purpose IDP / LLPS. |
| `hps_kr` | Ashbaugh–Hatch (LJ 12-6) | Kapcha–Rossy | RNA & phospho-protein parameters | Mixed protein + nucleic-acid systems (see Tutorial 6). |
| `hps_ss` | Ashbaugh–Hatch (LJ 12-6) | Urry | Gaussian **angle** + **torsion** bonded terms | Chains where local backbone stiffness matters. |
| `mpipi` | **Wang–Frenkel** short-range | Mpipi parameter set | — | Near-quantitative phase behavior (Joseph et al. 2021). |

All four share the same backbone: harmonic bonds along the chain and Debye–Hückel
(Yukawa) electrostatics between charged beads. They differ in the **non-bonded
potential** (Ashbaugh–Hatch vs Wang–Frenkel), the **per-residue parameters**
(σ, λ, charge), and — for `hps_ss` — the presence of explicit **angle/torsion**
terms. Define additional models in `cosmo/parameters/model_parameters.py`.

## Step-by-step

### 1. Run the default model
```bash
cosmo-mdrun -f md.ini                # model = hps_urry  (== python -m cosmo.mdrun -f md.ini)
python run_simulation.py -f md.ini   # equivalent in-folder shim
```
You get the usual `traj/asyn.*` outputs (see Tutorial 1).

### 2. Swap the model and re-run
Edit the one line in `md.ini`:
```ini
model = mpipi          # then hps_kr, then hps_ss
```
and (to avoid overwriting) change `outname` to match, e.g.
`outname = asyn_mpipi`. Re-run after each change. Watch the build log: each
model prints a **different set of force terms** as it assembles the system —
`hps_ss` adds angle/torsion groups; `mpipi` reports a Wang–Frenkel force where the
others report Ashbaugh–Hatch.

### 3. See the numbers: per-force-group energies
The most informative way to compare models is the repository's regression
benchmark, which decomposes the **initial potential energy by force group** for
all four models on this exact chain:
```bash
python ../../benchmarks/benchmark_energies.py
```
The per-force-group table makes the differences concrete — e.g. the non-bonded
term changes name and magnitude between `hps_urry` and `mpipi`, and `hps_ss`
gains non-zero angle/torsion rows the others don't have. This is also the script
that guards the physics: any code change that alters these numbers is flagged as
a regression.

## Background: why the choice matters

The model *is* the physics. The hydropathy scale sets how "sticky" each residue
is, which controls chain compaction and the driving force for phase separation;
the non-bonded functional form sets the shape of the attractive well. Two models
can agree on a single chain's size yet disagree sharply on its phase behavior, so
pick the model your scientific question (and the literature you compare against)
calls for — then keep it fixed across a study.

## Try next

- Compute R_g for `asyn` under `hps_urry` vs `mpipi` and compare — the scales are
  calibrated differently.
- Use `hps_ss` and inspect whether the added angle/torsion terms stiffen the
  local backbone in the trajectory.
- Move to **Tutorial 3** to put a chain in a periodic box and run NVT/NPT.
