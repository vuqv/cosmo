# `cosmo.translation` — task list

Progress tracker for the co-translational synthesis subsystem. Design and full
rationale: [`cosmo/translation/PLAN.md`](../../cosmo/translation/PLAN.md).
`[x]` + ~~strikethrough~~ = done; `[ ]` = remaining.

## Done — build step v1 (nascent-only elongation) + scaffold

- [x] ~~Scaffold the `cosmo/translation/` package (`__init__`, `__main__`, CLI).~~
- [x] ~~`ElongateConfig` + `read_elongate_config` (INI `[OPTIONS]` parser),
  wrapping the base run options; `cosmo-elongate` console entry point.~~
- [x] ~~`elongate.py` v1 runner: per-length build on the **first `L` residues**
  (sequence-based — no STRIDE / contact map), cold-start layout, A→P new-residue
  placement, C-terminus harmonic restraint to the P-anchor, minimize → run →
  seed next length. Reuses `cosmo.engine` (setup/reporters/finalize).~~
- [x] ~~`read_anchor` — P-/A-anchors (PtR/AtR residue-76 `R` bead) from the
  truncated CG ribosome.~~
- [x] ~~`init_position` as-is coordinate seeding every length (no shift; absolute
  tunnel frame preserved).~~
- [x] ~~`truncate_ribosome.py` — per-residue tunnel truncation (port of topo's;
  supersedes the old per-atom `cosmo/utils/crop_ribosome.py`).~~
- [x] ~~`make_movie.py` + `cosmo-elongate-movie` — stitch the per-length DCDs into
  one fixed-width VMD movie (`movie.psf`/`.dcd`/`.tcl`) that grows the chain N→C.~~
- [x] ~~Example/tutorial: `tutorials/07_translation/` (self-contained `asyn.pdb` +
  `elongate.ini` + README).~~
- [x] ~~Validated v1 on α-synuclein (`L0=5`): chain grows one residue/step,
  C-terminus lands on the P-anchor to ~0.002–0.005 nm, per-length outputs written,
  movie stitches (160 frames).~~

## Done — phase 3 (ribosome coarse-graining)

- [x] ~~`cg_ribosome.py`: all-atom → CG ribosome, **both `rna_model` reps** —
  `topo` (3/4-bead P/R/BR; reproduces 14,662 beads: 3363 CA / 3163 P / 3163 R /
  4973 BR) and `cosmo` (1 bead/nucleotide). Preserves `segID`.~~
- [x] ~~Copied the all-atom source `4v9d_50S_PtR_5jte_AtR_model.pdb` into
  `cosmo/translation/structures/`; generated the truncated CG ribosomes (topo →
  4576 beads, matching O'Brien; cosmo → 1706).~~
- [x] ~~Local P/R/BR bead table (radii 0.71 nm; P −1e) in `ribosome.py`.~~

## Done — phase 4 (build step v2: rigid ribosome scenery)

- [x] ~~`ribosome.py`: append the truncated ribosome as mass-0 beads at `L..N-1`
  (verified frozen: max drift 0.0 nm).~~
- [x] ~~Ribosome↔nascent **excluded volume** (`CustomNonbondedForce`, `ε·(σ/r)¹²`,
  `{nascent}×{ribosome}`) + **Yukawa** extended over ribosome charges; shared
  exclusion list across all `CustomNonbondedForce`s.~~
- [x] ~~`mpipi` dummy in-range `id` + `{nascent}×{nascent}` group on the WF force;
  HPS short-range force restricted to `{nascent}×{nascent}` with dummy beads.~~
- [x] ~~Auto `ptc_offset` (0.4 nm in v2) so the C-terminus clears the P-anchor bead;
  nascent-only output (DCD/PSF/final), full system in the checkpoint.~~
- [x] ~~Wired `rigid_ribosome = yes` through `run_elongation`; validated stable runs
  under **both** `rna_model` (topo & cosmo) on α-synuclein.~~
- [x] ~~`make_movie.py --ribosome` loads the static ribosome as VMD scenery.~~

## Done — phase 5a (v2 extras: tether + wall)

- [x] ~~tRNA tether (`trna_tether`, default on): harmonic bond `CA(L)–tRNA` +
  double-Gaussian orienting angle `CA(L-1)–CA(L)–tRNA` (91.7°/130°, O'Brien basins);
  dedicated angle force so it works for every nascent model. Validated: C-terminus
  held at 0.477 nm from the tRNA bead.~~
- [x] ~~One-sided planar tunnel wall (`tunnel_wall`, default on):
  `U = k·min(x−x0,0)²`; `x0` auto = P-anchor x + tether bond length. Validated:
  nascent chain stays at x ≥ x0 and extrudes +x.~~
- [x] ~~Both wired through `run_elongation` / `elongate.ini` and validated under
  both `rna_model` reps (anchor bead R for topo, P for cosmo).~~
- [x] ~~Equilibrium placement: new residues seeded at the tether/restraint **rest
  position** (not ~1 nm away at the A-site), so the bond starts unstretched — drops
  the per-step temperature spikes (peaks ~780 K → ~480 K).~~
- [x] ~~Diagnosed the "chain penetrates the ribosome" report: the excluded volume is
  on but O'Brien's eps is ~4500× < kT (deliberately soft), because the topo rRNA
  radius (0.71 nm) is large vs cosmo's CG tunnel bore (~0.4–0.6 nm) — a hard wall
  jams. Made `ribo_eps` and `ribo_rna_radii` tunable; `ribo_rna_radii=0.40` +
  `ribo_eps=0.2` confines the chain (overlap 0.22 → 0.04 nm) without jamming. Tutorial
  uses these; defaults stay O'Brien-soft.~~

## Remaining

### Phase 5b — post-elongation + analysis
- [ ] Post-elongation phase: ejection vs. stallation (`make_movie` already appends
  these segments if present).
- [ ] Analysis (Rg / contact fraction vs. length) + a fuller tutorial page.
