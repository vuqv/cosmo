# Co-translational synthesis through an analytic exit tunnel

Tutorials [`07`](../07_translation/) and [`08`](../08_translation_cosmo_rna/) build
the ribosome from **explicit beads** and let the nascent chain thread the
coarse-grained exit tunnel. With tutorial 07's **soft** excluded-volume wall the
chain **fails to extrude** — at L ≈ 40 it is a collapsed globule (Rg ≈ 1.2 nm) parked
at the PTC, and by full length Rg plateaus ~2 nm with the leading edge stalled near
7 nm. Tutorial 08's **full-HPS hard wall fixes that** (the chain stays extended and
extrudes to ~12 nm), but it still requires building and equilibrating an explicit CG
ribosome.

This tutorial offers a lighter, purely **analytic** alternative — no ribosome beads
at all: a cylindrical **bore**
of radius `r` along the X-axis drilled through an **infinite wall** at `x_exit` (a
"hole in an infinite wall"). There are **no ribosome beads** — the System is the
nascent chain only, so it is fast and never jams. The bore confines the in-tunnel
segment radially, so the chain stays **extended** and threads out the exit rather
than collapsing. See [`PLAN.md`](PLAN.md) for the full design and the exact force.

## Why a separate script (`cylinder.py`, not `cosmo-elongate`)

The analytic tunnel is a different *physics of confinement* than the explicit-bead
ribosome, so it lives here as a tutorial variant and does **not** modify the
shipped `cosmo.translation` package. It **reuses** that package's tested,
unchanged machinery (per-length build / seed / restrain / output — the "v1 output
path") and adds only the one new force, `add_tunnel_cylinder`, plus a nascent-only
elongation loop.

## The model (forbidden region `S`)

```
   d                              (cytosol: free, any d)
   ^   |##### solid ribosome S #####|
 r |···|············ bore ··········|··············>  allowed past exit
   +---|----------------------------|----------------> x
     x_lo (PTC)                  x_exit
       |##### solid ribosome S #####|
                              ^ infinite exit-face wall (d > r)
```

A bead is penalised by its **penetration depth into `S`** (everything outside the
bore up to the exit face, plus the closed PTC end), escaping via whichever face is
nearer — the bore wall (radial inward push → keeps the chain extended) or the exit
face (`+x` push → a cytosol bead can re-enter the tunnel **only through the bore**,
never off-axis). The 90° mouth corner is rounded by a fillet (radius `rho`) so the
potential is continuous and the MD is stable. The C-terminus is seeded and
position-restrained on the tunnel axis at the PTC `(x_lo, 0, 0)`; new residues are
seeded there.

## Run

```bash
cd tutorials/09_translation_cylinder
python cylinder.py -f elongate.ini
```

All parameters live in [`elongate.ini`](elongate.ini) (`[OPTIONS]` section). The
defaults (PLAN.md §4): bore radius **0.9 nm**, length **10.0 nm** (`x_lo=0`,
`x_exit=10`), axis on X, mouth fillet **0.2 nm**, wall stiffness **8368
kJ/mol/nm²**, `model = hps_urry`, `L0 = 5`.

## Post-elongation: ejection

Once the chain reaches its final length, an optional **post-elongation** phase
continues the same system from the finished structure (written to
`<outdir>/<phase>/`):

- **`ejection`** — releases the C-terminus restraint, so the completed protein is
  free to diffuse. The analytic tunnel stays on (bore + closed PTC end + exit
  wall), so the only way out is +x through the exit. This tests whether the
  nascent protein **diffuses out of the tunnel into the cytosol**.
- **`stallation`** — keeps the restraint, so the chain stays threaded/stalled.

Set in [`elongate.ini`](elongate.ini):

```ini
post_elongation       = ejection
post_elongation_steps = 300_000   # use a LONG run so the protein can clear the tunnel
```

In a test run (elongate to L=40, then 300k ejection steps), the protein **fully
clears the tunnel**: it starts collapsed inside (`x ≈ 0–4 nm`), and once released
the entire chain crosses the exit (`min x > 10 nm` by ~20% into the run), the
center-of-mass moves from ~2 → ~14 nm, and no bead re-enters off-axis (the exit
wall blocks re-entry except through the bore) — the designed behaviour.

## Visualize / validate

Use the tutorial's own movie tool — it stitches the per-length trajectories like
`cosmo-elongate-movie` **and** draws the analytic tunnel (bore tube, closed PTC
cap, and the infinite exit-face wall as an annulus whose hole is the bore),
reading the geometry from the same `elongate.ini` so the drawing matches the run:

```bash
python make_movie_cylinder.py -o synth_out -f elongate.ini
vmd -e synth_out/movie.tcl
```

You then see the chain thread the (blue, transparent) bore, the red PTC end cap it
grows away from, and the grey exit wall it must emerge through. (The plain
`cosmo-elongate-movie -o synth_out` still works — it just omits the tunnel.)

Success = the **N-terminus emerges past `x_exit` (10 nm)** into the cytosol once the
chain exceeds ~26 extended residues — which tutorial 07's soft wall never achieved
(tutorial 08's full-HPS wall does extrude to ~12 nm, but with explicit ribosome beads
rather than this bead-free analytic tunnel). Check
(PLAN.md §6): the leading-edge x climbs toward and past `x_exit`; the in-tunnel
segment stays extended (end-to-end ≫ the 2.3–4.3 nm globule); the run is stable at
~300 K; and there are no off-axis beads at `x < x_exit, d > r` (no re-entry artifact).
