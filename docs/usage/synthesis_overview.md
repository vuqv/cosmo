# Protein synthesis: overview

cosmo can grow a chain **vectorially, N-terminus first**, out of a ribosome while
the nascent chain samples its conformations *as it emerges* — reproducing
co-translational synthesis of an intrinsically disordered protein under
codon-resolved kinetics. This page is the map: the biology being modelled, the two
ribosome models cosmo offers, and which to use. The detailed references are linked at
the end.

## The biology in one minute

A ribosome synthesizes a protein one residue at a time, N→C, threading the growing
("nascent") chain out through the ~80 Å **exit tunnel**, where it begins to sample
conformations before the full sequence exists. How long the ribosome dwells on each
codon sets how much time each segment has to explore before the next residue arrives;
slow (rare) codons act as translational pauses that can change the emerging ensemble.
For a disordered chain the readout is not a native fold but the conformational
statistics — compaction, extension, contacts — of the chain as it extrudes.

Each residue is added through one **elongation cycle** — (1) codon-dependent
aminoacyl-tRNA **decoding** at the A site (usually rate-limiting), (2) **peptidyl
transfer** at the peptidyl-transferase centre (PTC), (3) **translocation** (A→P→E).
cosmo times each residue from its mRNA codon and, in the explicit model, splits the
cycle into these three sub-stages.

## Two ribosome models

Both grow the same sequence-based **IDP** (HPS / `mpipi`) one-bead-per-residue (Cα)
nascent chain — the same force field cosmo uses for equilibrium runs, see
{doc}`../tutorials/02_models` — and both use the same codon kinetics. They differ only
in **how the exit tunnel is represented**:

| | Explicit CG ribosome | Analytic cylinder |
|---|---|---|
| **Runner** | `cosmo-csp` ({doc}`continuous_synthesis`) | `cosmo-cylinder` ({doc}`cylinder_synthesis`) |
| **Tunnel** | a real truncated, coarse-grained large subunit (~4,600 rigid beads; *E. coli* 4V9D) | an analytic cylindrical bore through an infinite wall |
| **Interactions** | ribosome↔nascent **excluded volume + electrostatics**; tRNA tether; A→P translocation | **geometric confinement only** |
| **Needs a structure?** | yes — a truncated CG ribosome PDB ({doc}`ribosome_preparation`) | no |
| **Per residue** | three MD sub-stages (transfer / translocation / wait) | one MD segment |
| **Cost / robustness** | heavier; realistic tunnel wall | fast; never jams |

Both add optional **post-synthesis phases** once the chain reaches full length, run in
order: a `stall` hold at the PTC (restraint still on) then an `ejection` free run
(restraint released) — see {doc}`synthesis_control`.

## Which to use

- **Cylinder** — for fast exploration of how *tunnel geometry + codon kinetics* shape
  the disordered ensemble as it extrudes, or when you have no ribosome structure.
  Simplest starting point.
- **Explicit** — when the **tunnel-wall charge, tunnel shape** (constriction site /
  vestibule), or **translocation-coupled forces** matter to your question.

```{warning}
**The two models are comparable only in the *mean* per-residue dwell time, not in
confinement chemistry.** The cylinder omits the ribosome's electrostatics and surface
excluded volume and has a uniform straight bore. Do not compare conformational
observables (contact patterns, compaction-vs-length, radius of gyration) across the two
models without accounting for those missing terms. See {doc}`cylinder_synthesis`
§"The model: analytic exit tunnel".
```

## Where to go next

* {doc}`ribosome_preparation` — get or build the truncated CG ribosome the explicit
  model needs. cosmo bundles the *E. coli* tutorial ribosome; for other organisms (or
  to build your own) use the sibling `topo` package's scripted preparation pipeline.
* {doc}`codon_dwell_times` — the per-codon dwell-time tables that set the timing
  (*E. coli*, *S. cerevisiae*, *H. sapiens*, *N. crassa* ship with cosmo), and the
  `fastest`/`slowest`/`median` synonymous-codon mRNA.
* {doc}`continuous_synthesis` — the explicit-ribosome runner: the RNC force field,
  the tRNA tether/PTC anchors, the 3-stage mechanics, the codon→MD-step kinetics.
* {doc}`cylinder_synthesis` — the analytic-tunnel runner.
* {doc}`synthesis_control` — the `csp.ini` / `cylinder.ini` parameter reference.
* {doc}`../tutorials/07_csp_cylinder` and
  {doc}`../tutorials/08_csp_cg_ribosome` — hands-on, ready-to-run examples.
