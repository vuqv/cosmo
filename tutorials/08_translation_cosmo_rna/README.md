# Co-translational synthesis with the full-HPS COSMO RNA ribosome

Same elongation run as [`07_translation`](../07_translation/), with two changes that
go together:

1. The rRNA is coarse-grained with **cosmo's own model — one bead per nucleotide**
   (`rna_model = cosmo`) instead of the topo 3/4-bead P/R/BR representation.
2. The ribosome interacts with the nascent chain through the **full Ashbaugh–Hatch
   HPS potential** (`ribo_full_hps = yes`) — a hard, full-size wall with the *real*
   per-bead σ and hydropathy (ε = 0.8368 kJ/mol) — instead of the soft, separate
   pure-repulsion excluded volume used in tutorial 07.

Because the ribosome now uses the real HPS pair potential, the force field must
parameterise **both protein and RNA**, so this tutorial uses **`model = hps_kr`**
(`hps_urry` has no RNA parameters).

## Why full HPS

The soft excluded-volume wall in tutorial 07 needs two hand-tuned fudge factors
(`ribo_eps`, `ribo_rna_radii`) chosen so the chain can thread the tight CG bore.
Full HPS removes them: the rRNA beads simply use their hps_kr collision radii and
hydropathy, so the ribosome↔nascent interaction is the *same physics* as the
nascent↔nascent interaction — a single, consistent force field for the whole
system. (`ribo_eps` / `ribo_rna_radii` are **ignored** when `ribo_full_hps = yes`.)

## Run

```bash
cd tutorials/08_translation_cosmo_rna
cosmo-elongate -f elongate.ini
```

The tether, tunnel wall, and nascent-only output match tutorial 07; the differences
are the rRNA representation and the full-HPS ribosome. The anchors are read from the
single `P` bead of the P-/A-site tRNA residue 76 (the cosmo rep has no separate
ribose `R` bead).

## Visualize

```bash
cosmo-elongate-movie -o synth_out \
    --ribosome ../../cosmo/translation/structures/4v9d_50S_PtR_5jte_AtR_model_cg_cosmo_trunc.pdb
vmd -e synth_out/movie.tcl
```

## What to expect (and what it tells us)

The full-HPS hard wall **does drive extrusion** — this is the change that fixes the
collapse seen with tutorial 07's soft wall. Growing asyn to full length (L = 140),
the nascent chain stays extended and **threads out along +x** as it grows, instead of
balling up at the PTC (final structures, `synth_out/L_<L>/traj_final.pdb`):

| L | Rg (nm) | end-to-end (nm) | leading-edge x (nm) |
|---|---------|-----------------|---------------------|
| 20  | 1.0 | 2.7  | 3.9  |
| 40  | 1.4 | 4.1  | 5.3  |
| 80  | 2.3 | 8.2  | 9.3  |
| 140 | 3.6 | 11.5 | 11.8 |

By full length Rg ≈ 3.6 nm — essentially the bulk expanded-IDP value for 140-residue
asyn, **not** a collapsed globule (~1–2 nm) — and the leading edge has marched out to
~12 nm. Contrast tutorial 07's **soft** excluded-volume wall: on the same system Rg
plateaus near 2.1 nm and the leading edge stalls at ~7 nm after L ≈ 100, i.e. the
chain collapses and stops extruding.

The lesson: a ribosome wall built from the **real, consistent HPS pair potential**
(hard, full-size beads, same physics as the nascent↔nascent interaction) supplies the
radial confinement a hand-tuned soft wall does not — and that is enough to keep a
nascent IDP extended and extrude it. Tutorials 09 and 10 reach the same goal with
**alternative confinement geometries**:

- [`09_translation_cylinder`](../09_translation_cylinder/) — an analytic tunnel
  (bore + exit wall, no explicit ribosome beads) that confines the chain radially and
  lets the finished protein eject into the cytosol.
- [`10_translation_kinetics`](../10_translation_kinetics/) — growth at the PTC inside
  a frozen explicit ribosome with a tight forward wall (a piston) that ratchets the
  chain out.
