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

### Co-translational synthesis series (7–12)

Tutorials 7–12 are an **advanced, research-oriented series** on
**co-translational protein synthesis**: instead of folding a full-length chain in
bulk, the nascent chain is grown **N→C, one residue at a time**, at the ribosome's
peptidyl-transferase center (PTC) and extruded down the exit tunnel. Each length is
built, seeded from the previous one, restrained at the new C-terminus, minimized,
run, and saved (`synth_out/L_<L>/`), so the final structure of one length seeds the
next. They share the elongation machinery (`cosmo.translation`,
`cosmo-elongate` / `cosmo-elongate-movie`) and differ only in **how the ribosome /
exit tunnel is modeled** — the open research question of the series.

| # | Tutorial | Ribosome / tunnel model | Outcome |
|---|----------|-------------------------|---------|
| 7 | [Translation (explicit ribosome)](./07_translation/) | v1: tRNA **anchors only**; v2: **rigid explicit ribosome** (topo 3/4-bead P/R/BR rRNA) with a **soft** excluded-volume wall (`ribo_eps`/`ribo_rna_radii` hand-tuned so the chain can thread the tight CG bore). | Chain **collapses** into a globule at the PTC. |
| 8 | [Translation, COSMO RNA + full HPS](./08_translation_cosmo_rna/) | Explicit ribosome with **cosmo 1-bead/nucleotide** rRNA and the **full Ashbaugh–Hatch HPS** potential (a hard, physically-parameterized wall; `model = hps_kr`, no fudge factors). | **Solves it** — the hard, full-size wall gives real radial confinement; the chain stays extended (Rg ≈ 3.6 nm at L=140) and extrudes to ~12 nm in +x. |
| 9 | [Translation through an analytic tunnel](./09_translation_cylinder/) | **No ribosome beads** — an analytic **cylindrical bore through an infinite wall** ("hole in a wall") gives radial confinement; fast, never jams (`cylinder.py`). | Chain stays **extended** and **threads out** the exit; an optional **ejection** phase clears the tunnel into the cytosol. |
| 10 | [Translation kinetics (frozen ribosome)](./10_translation_kinetics/) | **Full frozen explicit CG ribosome** with a PTC at the origin and a forward "piston" wall; port of the legacy `growing_reference/`. Constant steps/residue now; **per-codon kinetics** is a planned follow-up. | Grows + extrudes through the real ribosome; full **ejection** depends on chain length vs. tunnel depth. |
| 11 | [Continuous Synthesis — analytic tunnel (`cosmo-cylinder`)](./11_csp_cylinder/) | **`cosmo.csp`, no ribosome beads** — an analytic cylindrical bore confines the chain radially; one MD segment per residue, per-codon O'Brien kinetics (mrna + *E. coli* table). | Codon-resolved synthesis + `ejection`/`dissociation` free runs; fast, never jams. |
| 12 | [Continuous Synthesis — coarse-grained ribosome (`cosmo-csp`)](./12_csp_ribosome/) | **`cosmo.csp`, explicit rigid CG ribosome** (topo P/R/BR rep, 4576 beads) as scenery + A-/P-anchors; each cycle split into O'Brien's 3 kinetic sub-stages (peptidyl-transfer / translocation / tRNA-binding). | Codon-resolved synthesis on a real CG ribosome; nascent-only trajectories, ribosome overlaid in the movie. |

The arc of the series: tutorial 7's **soft**, hand-tuned excluded-volume wall lets
the nascent IDP collapse at the PTC (Rg plateaus ~2 nm, extrusion stalls at ~7 nm).
Tutorial 8 **fixes this** by building the ribosome wall from the real, consistent HPS
pair potential (hard, full-size beads) — that supplies the radial confinement the
soft wall lacks, so the chain stays extended (Rg ≈ 3.6 nm at L=140) and extrudes to
~12 nm. Tutorials 9 and 10 reach the same goal with **alternative confinement
geometries** — an analytic radial tunnel (9) and a frozen explicit ribosome with a
forward ratchet (10) — and add an optional **post-elongation** phase (`ejection` to
release and diffuse out, or `stallation` to stay threaded). Tutorials 11 and 12 then
**consolidate** the two confinement paths into one package, `cosmo.csp` (mirroring
the sibling `topo.csp`), with codon-resolved O'Brien kinetics — the canonical
synthesis runner going forward, split by geometry: **11** is the analytic tunnel
(`cosmo-cylinder`) and **12** is the explicit coarse-grained ribosome (`cosmo-csp`).

The **ready-to-run files** for each tutorial (PDB, `md.ini` / `elongate.ini`,
`run_simulation.py`) live in the matching folder here under `tutorials/`. The
translation tutorials are launched differently from the `cosmo-mdrun` runs above:
7 and 8 use `cosmo-elongate -f elongate.ini`, 9 uses `python cylinder.py -f
elongate.ini` (analytic tunnel), 10 uses `python grow_nascent.py -f md.ini`
(frozen ribosome), 11 uses `cosmo-cylinder -f cylinder.ini` (analytic tunnel), and
12 uses `cosmo-csp -f csp.ini` (explicit ribosome). See each folder's `README.md`
for the exact command and options.

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
