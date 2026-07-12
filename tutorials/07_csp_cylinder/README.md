# Tutorial 7 — Continuous Synthesis through an analytic tunnel (`cosmo-cylinder`)

**Goal:** synthesize a nascent chain **residue by residue with codon-resolved
kinetics**, confining it with an **analytic cylindrical tunnel** instead of explicit
ribosome beads. This is the tunnel-geometry half of the Continuous Synthesis
Protocol ported to COSMO's IDP force field (`cosmo.csp`, mirroring `topo.csp`). Its
sibling — the same protocol on an **explicit coarse-grained ribosome** — is
[Tutorial 8](../08_csp_cg_ribosome/).

Unlike topo's structure-based `csp` (native contacts, STRIDE, rigid `AllBonds`),
COSMO's port is **sequence-based**: no STRIDE, no `domain.yaml`, no native-contact
map — the nascent chain is a transferable HPS/`mpipi` IDP that grows N→C at the
peptidyl-transferase center (PTC) and extrudes down the analytic bore.

> **This is a proof of concept.** The settings here check that the machinery
> **builds, times, runs and writes outputs** end to end in seconds on CPU. They are
> deliberately tiny (short `L` range, clamped `max_steps_per_stage`) and are **not**
> a physical/scientific validation.

**Prerequisites:** the coarse-grained workflow of
[Tutorial 1](../01_single_chain_quickstart/) and the force fields of
[Tutorial 2](../02_models_and_forcefields/). This tutorial confines the nascent chain
with an analytic **cylindrical bore** (no explicit ribosome beads) — the
ribosome-free confinement path of `cosmo.csp`.

Nascent chain: **α-synuclein** (`asyn.pdb`, 140 residues, from Tutorial 1).

## The model

**No ribosome beads.** An analytic cylindrical bore of radius `tunnel_radius`
through an otherwise infinite wall (a "hole in a wall") supplies the radial
confinement; the chain extrudes along `+x`. Because the wall is a smooth potential
rather than discrete beads, it is **fast and never jams**. Each residue is grown,
seeded from the previous length, restrained at the new C-terminus, minimized, and
run for **one MD segment** whose length is set by the codon's translation time.

## Files

| File | Role |
|------|------|
| `asyn.pdb` | Full-length nascent protein (only the first `L` residues exist at each length). |
| `mrna.txt` | Demo mRNA (140 cycled sense codons + `UAA` stop) — times each residue. |
| `cylinder.ini` | Analytic-tunnel config (`model = hps_kr`, per-codon mRNA + the *E. coli* 310 K table from `assets/csp/codon_dwell_times/`, then an `ejection` free run). |

## Run it

```bash
cd tutorials/07_csp_cylinder

cosmo-cylinder -f cylinder.ini                  # -> synth_out_cyl/L_<L>/
#   or as a module:  python -m cosmo.csp.cylinder -f cylinder.ini
```

The config defaults to `device = CPU` so it runs anywhere; switch to `device = GPU`
for a real run. Kinetics are **per-codon**: `mrna.txt` gives the codon sequence and
`codon_times` points at the *E. coli* 310 K translation-time table under
`assets/csp/codon_dwell_times/`. `max_steps_per_stage` clamps each residue's dwell
(otherwise ~10⁶ steps) down to a traceable test size — raise it, extend `L_max`, and
switch to GPU for production.

### Stitch a movie

Draw the analytic exit tunnel into the movie by passing the same `cylinder.ini`, then
open it in VMD *from inside* the run folder (`movie.tcl` loads `movie.{psf,dcd}` by
basename, so it only resolves from there):

```bash
cosmo-csp-movie -o synth_out_cyl --tunnel cylinder.ini      # or: python -m cosmo.csp.movie -o synth_out_cyl --tunnel cylinder.ini
cd synth_out_cyl && vmd -e movie.tcl
```

Writes a fixed-width VMD-playable movie (`movie.dcd` / `movie.psf` / `movie.tcl`) that
plays the growing chain threading the (blue) bore, past the red PTC cap and the grey
exit wall. Omit `--tunnel` for a plain movie with no tunnel. For the full movie
reference (both synthesis models, all options), see
[Visualizing the synthesis process](../usage/synthesis_visualization.md).

## What it produces

- `synth_out_cyl/L_005/ … L_010/` — one folder per chain length, each with the MD
  trajectory and the nascent structure at that length.
- `synth_out_cyl/ejection/` — the post-synthesis free run: the C-terminus restraint
  is dropped (`ejection_steps`) so the finished chain diffuses out of the bore. Set to
  `0` to skip.
- `synth_out_cyl/dwell_times.dat` — per-residue codon + dwell-time table; the
  per-codon variation is visible here (e.g. `AUU` ≈ 2 ms vs. `CGU` ≈ 72 ms).

Generated `synth_out_*/` directories are **git-ignored** (bulky
trajectories/checkpoints). Delete them and re-run to reproduce.

## Where this sits in the series

Tutorials 7–10 explored *how to confine the nascent chain*. `cosmo.csp` consolidates
those into one package with codon-resolved kinetics and post-elongation
phases. This tutorial is its **analytic-tunnel** mode; the **explicit
coarse-grained ribosome** mode is [Tutorial 8](../08_csp_cg_ribosome/).
