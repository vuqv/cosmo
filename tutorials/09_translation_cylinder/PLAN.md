# Plan — tutorial 09: ribosome as an analytic **tunnel (hole in an infinite wall)**

> Status: design **agreed**, not yet implemented. A new build-step variant ("v3")
> for `cosmo.translation`. Implementation-ready: §3 gives the exact force
> expression, §5 the exact API/wiring. See `cosmo/translation/PLAN.md` for v1/v2.

## 1. Motivation

Tutorials 07 (topo 3/4-bead rRNA) and 08 (cosmo 1-bead rRNA) both fail to extrude
the nascent IDP: at L=40 (contour 14.8 nm) the chain is a **collapsed globule**
(end-to-end 2.3–4.3 nm, Rg ≈ 1–1.3 nm) parked at the PTC; nothing reaches the exit.
The RNA bead model is **not** the cause. The causes are: (a) the only longitudinal
restraint (`tunnel_wall`) is a planar `x ≥ x0` wall with **no radial confinement**,
so the chain balls up in y–z; (b) a hard explicit-bead wall jams (the CG bore is
narrower than the bead contact distance). A real exit tunnel extrudes the chain
because it is a **narrow channel** the chain cannot collapse inside. Tutorial 09
models that directly with an analytic tunnel, replacing the explicit ribosome beads.

## 2. Concept — a hole in an infinite wall

The exit tunnel is a cylindrical **bore** of radius `r` along the X-axis, drilled
through an **infinite wall at x = x_exit** (the ribosome's outer exit face). The
"solid ribosome" (forbidden region) is

```
  S = { x < x_exit  AND  d > r }   ∪   { x < x_lo }
      d = sqrt((y - y0)^2 + (z - z0)^2)      # radial distance from the tunnel axis
```

everything outside the bore up to the exit face, plus the closed PTC end. **Allowed**
= inside the bore (`d < r`, `x_lo ≤ x ≤ x_exit`) ∪ the whole cytosol (`x > x_exit`,
any `d`).

Agreed behaviour:
- In-tunnel beads are confined to the bore → the chain stays **extended** and
  threads out (fixes the 07/08 collapse).
- A cytosol bead (`x > x_exit`) may **re-enter the tunnel only through the bore**
  (`d < r`); an off-axis re-entry (`d > r`) hits the exit-face wall and is pushed
  back out (`+x` restoring force). *(This was the key confirmation.)*
- No explicit ribosome beads (fast; no clashes/jamming). The C-terminus is held on
  the tunnel axis at the PTC; new residues are seeded there.

```
   d                              (cytosol: free, any d)
   ^   |##### solid ribosome S #####|
 r |···|············ bore ··········|··············>  allowed past exit
   +---|----------------------------|----------------> x
     x_lo (PTC)                  x_exit
       |##### solid ribosome S #####|
                              ^ infinite exit-face wall (d > r)
```

## 3. Potential (exact)

One `CustomExternalForce` per nascent bead. Penalise the **penetration depth into
S** — escape via whichever face is nearer, the bore wall (`d − r`) or the exit face
(`x_exit − x`) — plus the closed PTC end. The 90° inner corner at the mouth
`(x_exit, r)` is rounded by a fillet of radius `rho`.

```
  U   = k * max(0, pen)^2  +  k * min(0, x - x_lo)^2
  pen = (rounded) min( x_exit - x , d - r )        # 0 outside S; > 0 inside S
```

Exact OpenMM energy string (globals `k, xlo, xexit, r, y0, z0, rho`; no per-particle
params):

```
"k*max(0, pen)^2 + k*min(0, x - xlo)^2;
 pen = select(corner, pround, psharp);
 corner = step(rho - qx)*step(rho - qd);
 pround = rho - sqrt((rho - qx)^2 + (rho - qd)^2);
 psharp = min(qx, qd);
 qx = xexit - x;
 qd = sqrt((y - y0)^2 + (z - z0)^2) - r"
```

Why this form (checked by hand):
- **Continuous** (U → 0 at both faces; no `step` energy discontinuity) → stable.
- `d` just > r, `x` ≪ x_exit  → `pen = d−r`  → **radial inward** force (confine/thread).
- `x` just < x_exit, `d` ≫ r → `pen = x_exit−x` → **+x** force (exit face blocks re-entry).
- `d < r` (in bore) **or** `x > x_exit` (cytosol) → `qd ≤ 0` or `qx ≤ 0` → `pen ≤ 0`
  → `max(0,pen)=0` → **U = 0**. Hence a cytosol bead enters the tunnel **only through
  the bore**.
- Mouth corner (`qx,qd < rho`): `pround` rounds the 90° edge with fillet radius
  `rho` (the sharp tip is cut to the allowed side: `pround < 0 → max(0,·)=0`).
- `k*min(0, x - x_lo)^2` = the closed PTC end (`x ≥ x_lo`).

The "infinite wall" is a **stiff finite harmonic wall** (large `k`) — MD needs a
smooth/finite potential (agreed).

## 4. Parameters (agreed defaults)

Geometry from the X-aligned truncated ribosome (PTC ≈ x=0; the *E. coli* tunnel is
~80–100 Å long, ~10–20 Å wide).

| key | symbol | default | note |
|---|---|---|---|
| `tunnel_radius` | `r` | **0.9 nm** (Ø 1.8 nm) | bore radius; ~3 CG beads wide — roomy but blocks the globule (Rg ≈ 1–1.3 nm in 07/08). |
| `tunnel_length` | `L_tun` | **10.0 nm** | tunnel length; `x_exit = x_lo + tunnel_length` (= 10.0). Needs ~`L_tun/0.38` ≈ 26 extended residues to span before the N-terminus emerges. |
| `tunnel_x_lo` | `x_lo` | **0.0 nm** | PTC / closed end; new residues seeded on-axis here. |
| `tunnel_x_exit` | `x_exit` | **derived = x_lo + tunnel_length** | exit plane = the infinite wall; not set directly. |
| `tunnel_k` | `k` | **8368 kJ/mol/nm²** | radial + exit-face + PTC-end stiffness (= 20 kcal/mol/Å²); raise for a harder wall. |
| `tunnel_center` | `(y0,z0)` | **(0.0, 0.0) nm** | tunnel axis = X-axis. |
| `tunnel_mouth_round` | `rho` | **0.2 nm** | fillet radius rounding the mouth corner. |

**Geometry summary:** a bore of **radius 0.9 nm**, **length 10.0 nm** along +x from
the PTC (`x_lo=0`) to the exit/infinite-wall (`x_exit=10`), axis on X, mouth rounded
with a 0.2 nm fillet.

## 5. API / wiring (implementation spec)

**`cosmo/translation/ribosome.py`** — add (sibling of `add_tunnel_wall`):
```
def add_tunnel_cylinder(system, nascent_indices, r_nm, x_lo_nm, x_exit_nm,
                        k=TUNNEL_WALL_K, y0_nm=0.0, z0_nm=0.0,
                        mouth_round_nm=0.2) -> mm.Force
```
builds the §3 `CustomExternalForce` over `nascent_indices`, adds it to `system`,
returns it.

**`cosmo/translation/elongate.py`**:
- `ElongationParams`: `ribosome_model: str = "beads"` (`"beads"` = v1/v2 explicit;
  `"cylinder"` = this tunnel), plus `tunnel_radius_nm=0.9`, `tunnel_length_nm=10.0`,
  `tunnel_x_lo_nm=0.0`, `tunnel_center_nm=(0.0,0.0)`, `tunnel_k=TUNNEL_WALL_K`,
  `tunnel_mouth_round_nm=0.2`.
- `read_elongate_config`: parse `ribosome_model`, `tunnel_radius`, `tunnel_length`,
  `tunnel_x_lo`, `tunnel_center` (`"y0,z0"`), `tunnel_k`, `tunnel_mouth_round`. In
  cylinder mode the `ribosome` PDB is **optional** (geometry comes from the params).
- `run_elongation` / `run_length`, cylinder mode:
  - `ribo = None` (no explicit beads / no `append_ribosome`); System is
    nascent-only → reuse the **v1 output path** (seed.pdb via `init_position`,
    `dumpTopology`, `finalize_simulation`; no nascent-only-DCD machinery).
  - `x_exit = x_lo + tunnel_length`; `cterm_seed = (x_lo, 0, 0)` (PTC on the axis).
  - Seed the C-terminus at `cterm_seed` (cold start lays the chain on-axis from
    there; each new residue is seeded at `cterm_seed`, reusing the equilibrium-
    placement logic already added for v2).
  - Hold the C-terminus with a **position restraint to `cterm_seed`** (no tRNA
    tether — no bead in cylinder mode).
  - Add the tunnel once per length:
    `add_tunnel_cylinder(cgModel.system, range(L), r=tunnel_radius, x_lo, x_exit,
    k=tunnel_k, y0,z0=tunnel_center, mouth_round=tunnel_mouth_round)`.
  - Do **not** add the separate planar `tunnel_wall` (the cylinder includes the
    PTC-end term).

**Tutorial files** `tutorials/09_translation_cylinder/`: `elongate.ini`
(`ribosome_model = cylinder`, the §4 params, `model = hps_urry`, `L0 = 5`),
`asyn.pdb` (copy), `README.md`.

## 6. Validation

Re-run the tutorial-08 A/B metric on `asyn` (L0=5 → ~40, n_steps ≥ 2000) vs 07/08:
- **Leading-edge x vs L** should climb toward and past `x_exit` (10 nm) once the
  chain exceeds ~26 residues — 07/08 stalled at ~3–5 nm.
- **In-tunnel segment extended** (end-to-end ≫ the 2.3–4.3 nm globule); residues
  past `x_exit` may fold (folding-outside-the-tunnel signal).
- **Stable, ~300 K** (smooth analytic wall; no bead clashes).
- **No re-entry artifact**: no off-axis beads at `x < x_exit, d > r`.
- Movie (`cosmo-elongate-movie`) to confirm threading + emergence.

Success = the N-terminus emerges past `x_exit` into the cytosol — which neither 07
nor 08 achieved.

## 7. Agreed decisions (resolved)

- **Cylinder only** — no explicit ribosome beads; System is nascent-only + the
  analytic tunnel. The ribosome PDB is optional in this mode.
- **Hole-in-an-infinite-wall** geometry (§2): bore + infinite exit-face wall at
  `x_exit`; cytosol free; **re-entry only through the bore** (`d < r`).
- **On-axis PTC** — cylinder centered on the X-axis; C-terminus seeded and
  position-restrained on-axis at the PTC `(x_lo, 0, 0)`; no tRNA tether.
- **Min-escape potential** (§3): continuous, with the **+x exit-face restoring
  force** that blocks off-axis re-entry (supersedes the earlier hard-step and
  tanh-gate ideas).
- **Stiff finite harmonic wall** (the "infinite" wall is a large-`k` harmonic).
- **Rounded mouth** — 90° corner filleted with radius `rho` (0.2 nm).
- **Purely steric** — geometry only; no ribosome electrostatics.

### Later / validation sweeps (not blockers)
- Sweep `r` ∈ {0.6, 0.9, 1.1} and `tunnel_length` ∈ {5.8 (bead-crop), 10.0}.
- Straight cylinder is the first approximation (the real tunnel is slightly bent).
