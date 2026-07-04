# Tutorial 10 — co-translational synthesis by iterative growth

This tutorial grows a nascent chain **one residue at a time** inside a **frozen
explicit ribosome**, the method from the legacy `growing_reference/` ported to the
current cosmo codebase. It differs from tutorials 07–09 in two ways: the ribosome
is the **full frozen CG structure** (not anchors or an analytic tunnel), and the
growth point is a simple **PTC at the origin** rather than a tunnel.

> This version runs a **constant** number of MD steps per residue (quick to test).
> Realistic **per-codon translation kinetics** — variable steps from measured codon
> translation times — is a planned follow-up (see [`PLAN.md`](PLAN.md) Phase 5).

## Method

For each residue `L = 1 … N`:

1. **Append** a new CA bead (chain `8`) at the origin to the current complex
   (ribosome + chain so far) → `L_<L>/complex.pdb`.
2. **Build** the whole complex with `cosmo.models.buildCoarseGrainModel(...,
   frozen_indices=<ribosome>, except_chains=<ribosome chains>)`: the ribosome is
   **frozen** (mass-0) and excluded from bonding; the nascent chain interacts
   normally. Ribosome–ribosome nonbonded is removed with **interaction groups**
   (efficient; not an O(n²) exclusion list).
3. **PTC mechanics:** the newest residue is harmonically restrained at
   `(0.38, 0, 0)` nm; all earlier nascent beads are kept at `x > 0.38 nm` by a
   one-sided planar wall — so the chain extrudes forward (+x).
4. **Run** `md_steps` of Langevin MD (reusing `cosmo.engine`), seeding the next
   residue from the final structure.

The **trajectory is nascent-only** (`L_<L>/traj.psf` + `traj.dcd` contain just the
growing chain); the simulation and `traj_final.pdb` keep the whole complex.

## Inputs (used as-is)

- `nascent.pdb` — the **CG (CA-only)** nascent chain (sequence source, 140 res).
- `rib.pdb` — the **already-cropped** CG ribosome (frozen scenery).
- `md.ini` — control file (`[OPTIONS]`): model `hps_kr`, `md_steps = 1000`,
  `nascent_chain_id = 8`, ecoli/fast metadata (consumed by the Phase 5 codon
  scheduler; just logged for now). `scale_factor` is intentionally absent here.

> `rib.pdb` was cropped upstream from a full ribosome with
> `python -m cosmo.utils.crop_ribosome rib_full.pdb rib.pdb` — not a step of this
> tutorial (the cropped file ships ready to use).

## Run

```bash
cd tutorials/10_translation_kinetics
python grow_nascent.py -f md.ini                 # all 140 residues
python grow_nascent.py -f md.ini --end 20        # first 20 (quicker)
python grow_nascent.py -f md.ini --end 5 --steps 500 -o test_out   # smoke test
```

Per-length outputs go to `synth_out/L_<L>/` (`traj.psf`, `traj.dcd`, `traj.log`,
`traj.chk`, `traj_final.pdb`, `traj_runinfo.log`).

## Post-elongation: ejection

Once the chain reaches full length, an optional **post-elongation** phase continues
from the finished complex (written to `<outdir>/<phase>/`):

- **`ejection`** — releases the PTC restraint so the completed chain is free to
  diffuse. The ribosome stays frozen and the forward wall stays on (`x > 0.38 nm`),
  so the chain can only diffuse **forward (+x), out of the ribosome**. Tests whether
  the nascent chain clears the exit.
- **`stallation`** — keeps the PTC restraint, so the chain stays pinned at the PTC.

```ini
post_elongation       = ejection
post_elongation_steps = 300_000   # use a LONG run; 0 = skip
```

The ejection frames are appended after the growth frames in the movie (below), so
you watch synthesis flow straight into release.

> **Note on what to expect.** Whether the chain fully ejects depends on its length
> vs the ribosome. The explicit frozen ribosome (~11 nm deep here) is a denser,
> longer obstacle than tutorial 09's short analytic tunnel — a short chain stays
> buried and diffuses *inside* without clearing in a few ×10⁵ steps (in a test, a
> 40-mer's COM drifted ~+2 nm and its leading edge just reached the exit, but no
> bead left). To see full ejection, grow to (near) full length so the N-terminus is
> already outside the exit, and/or run more steps.

## Visualize

The nascent-only `L_<L>/traj.*` layout works directly with the shipped movie tool.

**Chain only** — stitch the per-length trajectories into one movie of the chain
growing N→C:

```bash
cosmo-elongate-movie -o synth_out
vmd -e synth_out/movie.tcl
```

**Chain inside the ribosome** — pass `--ribosome` to also load the frozen ribosome
as static scenery, so you see the chain grow and extrude *through* it:

```bash
cosmo-elongate-movie -o synth_out --ribosome rib.pdb
vmd -e synth_out/movie.tcl
```

This copies `rib.pdb` next to the movie (`movie_ribosome.pdb`) and adds it to the
generated `movie.tcl` as a separate, static molecule (drawn as points). The
alignment is automatic: the growth runs use the coordinates **as-is** (no
recentering), so the nascent trajectory and `rib.pdb` share the **same absolute
frame** — the chain appears at the PTC (origin) inside the ribosome and threads out
toward +x. (In VMD you can then restyle the ribosome molecule — e.g. a transparent
surface or VDW — via Graphics → Representations.)
