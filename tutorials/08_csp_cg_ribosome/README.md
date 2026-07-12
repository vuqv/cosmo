# Tutorial 8 — Continuous Synthesis on a coarse-grained *E. coli* ribosome (`cosmo-csp`)

**Goal:** synthesize a nascent chain **residue by residue on an explicit
coarse-grained ribosome**, with the elongation cycle split into **three
kinetic sub-stages**. This is the explicit-ribosome half of the Continuous Synthesis
Protocol ported to COSMO's IDP force field (`cosmo.csp`, mirroring `topo.csp`). Its
sibling — the same protocol through an **analytic cylindrical tunnel** — is
[Tutorial 7](../07_csp_cylinder/).

Unlike topo's structure-based `csp` (native contacts, STRIDE, rigid `AllBonds`),
COSMO's port is **sequence-based**: no STRIDE, no `domain.yaml`, no native-contact
map — the nascent chain is a transferable HPS/`mpipi` IDP that grows N→C at the
peptidyl-transferase center (PTC). Everything is **standalone COSMO**: the ribosome
is a truncated CG bead model, no external CHARMM files.

> **This is the tutorial-scaled version of a production run.** It uses the *same*
> E. coli ribosome, mRNA and codon kinetics as the production configuration
> in [`sandbox/Ecoli/`](../../sandbox/Ecoli/) (`csp_val.ini`), but with a short `L`
> range, clamped `max_steps_per_stage`, and CPU so it **builds, times, runs and
> writes outputs** end to end in seconds. It is a mechanics demo, not a physical
> validation — see the [production values](#going-to-production) below.

**Prerequisites:** the coarse-grained workflow of
[Tutorial 1](../01_single_chain_quickstart/) and the force fields of
[Tutorial 2](../02_models_and_forcefields/). The analytic-tunnel sibling is
[Tutorial 7](../07_csp_cylinder/).

Nascent chain: **α-synuclein** (`asyn.pdb`, 140 residues, from Tutorial 1).

## The ribosome

The ribosome is the **E. coli 50S subunit, PDB [4V9D](https://www.rcsb.org/structure/4V9D)**
(with the [5JTE](https://www.rcsb.org/structure/5JTE) A-site tRNA grafted in),
coarse-grained to the **3/4-bead P/R/BR rRNA representation** (4576 beads) and
truncated to the region around the exit tunnel:
`4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb`. It enters as **rigid (mass-0) scenery**:
fixed beads that supply the ribosome↔nascent **12-10-6 excluded volume + Yukawa
electrostatics** wall and the **A-/P-site anchors**. This is the exact same structure
used by the production `sandbox/Ecoli/` run.

Each elongation cycle is split into three kinetic sub-stages —
**peptidyl-transfer → translocation → tRNA-binding** — with the new residue's
C-terminus restraint switching **A→P** across them, so a residue is delivered to the
A site and ratcheted into the P site. The PTC geometry is always optimized: the new
residue is seeded one peptide bond (0.380 nm) from the C-terminus, EV-clear.

## Files

| File | Role |
|------|------|
| `asyn.pdb` | Full-length nascent protein (only the first `L` residues exist at each length). |
| `4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb` | E. coli 4V9D 50S ribosome, CG P/R/BR rep (4576 beads) — rigid scenery + A-/P-anchors. |
| `mrna.txt` | Demo mRNA for α-synuclein (140 sense codons) — drives the per-codon kinetics. |
| `trans_times.txt` | E. coli 310 K codon → mean translation-time table (seconds), keyed by codon. |
| `csp.ini` | 3-stage explicit-ribosome config (`model = hps_urry`, real E. coli codon kinetics, tunnel wall on). |

## Run it

```bash
cd tutorials/08_csp_cg_ribosome

python -m cosmo.csp -f csp.ini      # -> synth_out_csp/L_<L>/  (per-stage traj_s<1,2,3>.dcd)
#   installed console script:  cosmo-csp -f csp.ini
```

The config defaults to `device = CPU` so it runs anywhere. Kinetics use the **real
E. coli codon table** (`trans_times.txt`) driven by the transcript (`mrna.txt`): each
residue's mean time is looked up by codon, split into the 3 sub-stages, and mapped to
steps via `scale_factor`. `max_steps_per_stage` clamps each sub-stage's (otherwise
~10⁶-step) dwell down to a traceable test size.

### Stitch a movie

```bash
cosmo-csp-movie -o synth_out_csp --ribosome 4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb
#   or: python -m cosmo.csp.movie -o synth_out_csp --ribosome 4v9d_50S_PtR_5jte_AtR_model_cg_trunc.pdb
cd synth_out_csp && vmd -e movie.tcl        # run from synth_out_csp/: movie.tcl loads its files by basename
```

`--ribosome` overlays the static ribosome (copied to `movie_ribosome.pdb`) so you see the
chain growing inside it. Writes a fixed-width VMD-playable movie (`movie.dcd` /
`movie.psf` / `movie.tcl`). For the full movie reference (both synthesis models, all
options), see [Visualizing the synthesis process](../usage/synthesis_visualization.md).

## What it produces

- `synth_out_csp/L_001/ … L_008/`, each with per-stage `traj_s1/2/3.dcd` — the
  A→P restraint switch across the three kinetic sub-stages.
- **Nascent-only trajectories**: the shared `L_<L>/traj.psf` + per-stage `traj_s*.dcd` contain only the
  `L` nascent beads — the 4576-bead ribosome is excluded from the trajectory and
  re-attached as static scenery by the movie stitcher's `--ribosome` overlay.
- `synth_out_csp/dwell_times.dat` — per-residue dwell-time table.

Generated `synth_out_*/` directories are **git-ignored** (bulky
trajectories/checkpoints). Delete them and re-run to reproduce.

## Going to production

The full production run this tutorial is scaled down from lives in
[`sandbox/Ecoli/`](../../sandbox/Ecoli/) (`csp_val.ini`). To reproduce it, start from
this `csp.ini` and:

- **Full length** — remove `L_max` (synthesize all 140 residues).
- **Un-clamp the dwell** — remove `max_steps_per_stage` / `min_steps_per_stage` so each
  stage's step count comes from the real dwell-time calculation.
- **GPU** — `device = GPU`, `ppn = 4`.
- **Post-synthesis** — `ejection_steps = 10_000_000` (let the finished chain diffuse
  out of the exit tunnel).

The same nascent chain + kinetics are run against different organisms' ribosomes in
`sandbox/` (Ecoli / Human / Yeast / Ncrassa) — only the `ribosome` PDB changes — to
compare co-translational behaviour across ribosome structures.

## Where this sits in the series

The two co-translational tutorials split `cosmo.csp` by confinement geometry:
[Tutorial 7](../07_csp_cylinder/) is the **analytic cylindrical tunnel** (no ribosome
beads; fast, never jams), and this one is the **explicit coarse-grained ribosome**
with 3-stage per-codon kinetics.
