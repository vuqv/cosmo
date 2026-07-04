# TASKS — tutorial 10: co-translational synthesis (port of `growing_reference`)

> Tracking checklist for porting `tutorials/growing_reference/` to the current
> cosmo codebase. See `PLAN.md` for the full design. Check items off as completed.

## Key decisions (locked)

- **Schedule for now:** constant **1000 steps per residue** (quick to test).
  Realistic per-codon translation kinetics is **deferred** — added in Phase 5.
- Per-iteration setup: **reuse `cosmo.engine`** (standard log + `_runinfo.log`).
- Output layout: **`L_<L>/traj.*`** so `cosmo-elongate-movie` works directly.
- Scope: **cosmo only** (hps_kr); keep the legacy frozen-full-ribosome + origin-PTC
  method (do not re-architect onto `elongate.py`).
- **Inputs are simulation-ready:** `nascent.pdb` is **CG (CA-only)**, `rib.pdb` is
  **already cropped** → use as-is; no nascent coarse-graining, no `crop_ribosome` step.
- **`scale_factor` omitted** from the test `md.ini` (no real-time→sim-time mapping
  with the constant schedule); it returns in Phase 5 with the codon scheduler.
- **Ribosome–ribosome exclusion via interaction groups** (post-build, as in
  `append_ribosome`) — **not** the reference's O(n²) `nb_exclusions` list.
- **Trajectory writes the nascent chain only** — nascent-only PSF + DCD; the frozen
  ribosome is never written (simulation + per-iteration final PDB still full complex).

---

## Phase 0 — scaffold
- [x] Create `tutorials/10_translation_kinetics/`.
- [x] Copy validated inputs as-is: `nascent.pdb` (CG, CA-only), `rib.pdb` (cropped).
- [x] Draft `md.ini` (model `hps_kr`, ecoli, fast, `n_steps = 1000`, frozen ribosome;
      **no `scale_factor`**).

## Phase 1 — port the core build (fix the 2 broken APIs)
- [x] New module `grow_nascent.py` (consolidate the 3 legacy variants into one).
- [x] `append_amino_acid_to_complex(...)` — keep as-is.
- [x] `buildHPSModel(...)` → `cosmo.models.buildCoarseGrainModel(...)` (keep
      `frozen_indices` / `except_chains`; **drop `nb_exclusions`**).
- [x] `from distutils.util import strtobool` → `from cosmo.utils.config import strtobool`.
- [x] Ribosome index/chain extraction via `cosmo.system().coarseGrainingStructure()`.
- [x] **Interaction groups** post-build: add `{nascent}×{nascent}` and
      `{nascent}×{ribosome}` to `cgModel.ashbaugh_HatchForce` and `cgModel.yukawaForce`
      (omit `{ribosome}×{ribosome}`) — replaces the O(n²) `nb_exclusions`.

## Phase 2 — PTC mechanics (kept from reference)
- [x] Harmonic restraint of the newest residue at `(0.38, 0, 0)` nm (`k = 1e4`).
- [x] One-sided planar wall `x > 0.38 nm` on all earlier nascent beads (`k = 1e4`).

## Phase 3 — modernized run loop
- [x] **Constant 1000 steps/residue** schedule (codon kinetics deferred to Phase 5).
- [x] Per-iteration run via `cosmo.engine` (`setup_simulation` / `attach_reporters`
      / `finalize_simulation`).
- [x] **Nascent-only PSF + DCD** output (`_82` variant pattern / `NascentDCDReporter`);
      frozen ribosome never written to the trajectory.
- [x] Output to `L_<L>/traj.*` layout.
- [x] Small config dataclass + `read_config()` for the translation keys
      (`organism`, `nascent_chain_id`, `rib_structure`, `nascent_structure`, ...).

## Phase 4 — validate
- [x] Smoke run: first ~5 residues, 1000 steps each, GPU.
- [x] Confirm build OK, frozen ribosome stays put, chain grows at the PTC
      (newest ≈ 0.38 nm; earlier beads walled to x > 0.38 nm).
- [x] Confirm `cosmo-elongate-movie -o <outdir>` stitches the outputs.
- [x] Write `README.md` (run command, visualization; mention `crop_ribosome` only
      as the optional upstream prep that produced `rib.pdb` — not a tutorial step).

## Phase 5 — realistic translation kinetics (LATER, not now)
- [ ] Replace the constant schedule with the per-codon scheduler:
      `nsteps = codon_time_ms * 1e9 / (dt_ps * scale_factor)` via
      `ctf_utils.codon_translation_time[organism][codon]['Median_time']`
      (`parse_nascent_sequence` + `generate_mRNA_sequence`).
- [ ] Add `scale_factor` / `translation_type` handling + per-residue codon log.
- [ ] Re-validate; document fast-vs-slow-codon pausing behaviour.

## Phase 4b — ejection (added on request)
- [x] `post_elongation` / `post_elongation_steps` in `GrowConfig` + `md.ini`.
- [x] `restrain_newest=False` path: release the PTC restraint, wall **all** nascent
      beads forward; ribosome stays frozen (only +x escape).
- [x] Post-growth `ejection` / `stallation` phase → `<outdir>/<phase>/` (nascent-only).
- [x] Validate: released chain is mobile + wall holds; movie appends ejection frames.
      Observation: a 40-mer diffuses but does not fully clear the explicit ribosome
      in 300k steps (length/geometry-limited) — documented in README.

## Optional / follow-up
- [ ] (none currently — interaction groups are in the main plan, Phase 1.)
