# Synthesis through an analytic tunnel (cylinder model)

The **cylinder** runner is the analytic-tunnel variant of protein synthesis: it replaces
the explicit coarse-grained ribosome of the
{doc}`coarse-grained ribosome model <continuous_synthesis>` with a **cylindrical bore
drilled through an infinite wall** (a "hole in an infinite wall"). There are **no ribosome
beads** — the simulated system is the nascent chain only — so it is fast, never jams on
ribosome excluded volume, and the analytic tunnel keeps the in-tunnel segment extended so
the chain threads out the exit.

- **CLI:** `cosmo-cylinder -f cylinder.ini` (or `python -m cosmo.csp.cylinder -f cylinder.ini`)
- **Worked example:** `tutorials/07_csp_cylinder/` (α-synuclein); a larger production
  configuration lives in `sandbox/validate/cylinder.ini`.
- **Module:** `cosmo.csp.cylinder` — a parallel of the explicit-bead `cosmo.csp.protocol`.
  It **reuses** the shared engine `cosmo.csp.core` (length model, seed / restrain / output
  path) and the timing core `cosmo.csp.kinetics`, adding only the one analytic tunnel force
  (`add_tunnel_cylinder`) and a nascent-only synthesis loop.

```{tip}
Same **codon kinetics** as the coarse-grained-ribosome runner (`cosmo-csp`) — each
residue's MD length comes from its codon dwell time — but the cylinder runs **one MD
segment per residue** rather than three sub-stages, because the analytic tunnel has no A→P
translocation to model. For the full kinetics derivation (codon → seconds → integration
steps) see the {doc}`coarse-grained ribosome page <continuous_synthesis>`; this page
describes the tunnel model and the `cylinder.ini` options.
```

```{note}
**Why the cylinder exists.** The explicit-bead exit tunnel in the truncated CG ribosome is
geometrically tight, and a disordered chain tends to jam or ball up at the PTC rather than
thread it. The analytic bore is a clean idealization that reliably extrudes an IDP chain —
purely steric, no electrostatics, no bead clashes.
```

---

## Quick start

```bash
cd tutorials/07_csp_cylinder
cosmo-cylinder -f cylinder.ini         # -> synth_out_cyl/

cosmo-csp-movie -o synth_out_cyl       # stitch per-length trajectories into a movie
vmd -e synth_out_cyl/movie.tcl
```

`cosmo-cylinder` writes, per residue `L`, a standalone trajectory under `<outdir>/L_<L>/`,
optional post-synthesis free runs (`ejection/` then `dissociation/`), and a per-residue
dwell-time log `<outdir>/dwell_times.dat`.

---

## The model: analytic exit tunnel

The tunnel is a cylindrical **bore** of radius `r` along the X-axis, drilled through an
**infinite wall** from the closed PTC end (`x_lo`) to the exit face
(`x_exit = x_lo + tunnel_length`):

```text
   d                              (cytosol: free, any d)
   ^   |##### solid ribosome S #####|
 r |···|············ bore ··········|··············>  allowed past exit
   +---|----------------------------|----------------> x
     x_lo (PTC)                  x_exit
       |##### solid ribosome S #####|
                              ^ infinite exit-face wall (d > r)
```

A single `CustomExternalForce` over every nascent bead penalises its **penetration depth
into the solid region `S`**:

```text
S   = { x < x_exit AND d > r } ∪ { x < x_lo },   d = |(y,z) − (y0,z0)|
U   = k·max(0, pen)²  +  k·min(0, x − x_lo)²
pen = (rounded) min( x_exit − x , d − r )        # 0 outside S; > 0 inside S
```

A bead escapes `S` via whichever face is nearer — the **bore wall** (`d − r`, a radial
inward push that keeps the in-tunnel chain extended) or the **exit face** (`x_exit − x`, a
+x push), so a cytosol bead can only re-enter through the bore, never off-axis. The 90°
inner corner at the mouth is rounded by a fillet of radius `rho = tunnel_mouth_round` so the
potential stays continuous and the MD stays stable.

The C-terminus is **position-restrained on the tunnel axis at the PTC** `(x_lo, y0, z0)`
(stiffness `restraint_k`). Each new residue is *seeded* one equilibrium peptide bond
**deeper** than that rest point — toward the closed PTC end — so the new `L-1↔L` bond
starts at its equilibrium length instead of collapsed onto the previous C-terminus; the
restraint and the closed-end wall then lift the new residue back onto the PTC while the
older chain ratchets forward. (This is the analytic-tunnel analogue of the explicit
protocol's A/P-site offset, without the ribosome, and is what keeps a rigid `AllBonds`
build stable here too.) There is no A/P tRNA tether and no translocation switch — the
chain simply extrudes forward as it grows.

```{note}
**How it differs from the {doc}`coarse-grained ribosome model <continuous_synthesis>`.**
(1) No `ribosome` PDB — the tunnel is analytic, its geometry set by the `tunnel_*` keys.
(2) **One MD segment per residue** (no peptidyl-transfer / translocation / tRNA-binding
sub-stages), so `time_stage_1` / `time_stage_2` are inherited but **unused** — the whole
codon dwell is a single segment. (3) The post-synthesis free runs are `ejection_steps`
then `dissociation_steps` (same keywords as the explicit-ribosome runner, both drop the
C-terminus restraint). (4) The explicit-ribosome
knobs (`trna_tether`, `tunnel_wall`) and the always-on PTC-geometry optimization do not apply.
(5) Purely steric — no ribosome electrostatics.
```

---

## Configuration reference (`cylinder.ini`)

Read by `cosmo.csp.cylinder.read_cylinder_config`. **Units are OpenMM defaults.** The
kinetics and MD keys share the semantics of the {doc}`synthesis_control` page; the tunnel
geometry and post-elongation phase are specific to this runner.

Example `cylinder.ini`:

```ini
[OPTIONS]
; --- inputs (no `ribosome` PDB: the tunnel is analytic) ---
pdb_file = asyn.pdb          ; native PDB (all-atom OR Cα-only CG; only Cα + resnames read)
model    = hps_kr            ; nascent IDP force field

; --- length schedule ---
L0    = 5
L_max = 10                   ; blank -> full residue count

; --- kinetics (same O'Brien codon timing as cosmo-csp) ---
mrna         = mrna.txt      ; one codon per residue (required for per-codon timing);
                             ; or "fastest"/"slowest"/"median" to auto-build a synonymous-codon mRNA
codon_times  = ../../assets/csp/codon_dwell_times/ecoli/ecoli_codon_dwell_times_310K.txt  ; table path (required; or a number of s = uniform)
scale_factor = 4331293
random_seed  = 1
max_steps_per_stage = 60     ; TEST CLAMP (delete for production)
min_steps_per_stage = 30

; --- mechanics / integrator ---
restraint_k = 83680          ; C-terminus -> PTC restraint (kJ/mol/nm^2)
minimize    = yes
dt = 0.01 ; ref_t = 300 ; tau_t = 0.01 ; nstout = 20

; --- analytic exit tunnel ---
tunnel_radius      = 0.9     ; bore radius r (nm); ~3 CG beads wide
tunnel_length      = 10.0    ; bore length (nm); x_exit = x_lo + length
tunnel_x_lo        = 0.0     ; PTC / closed end (nm); C-terminus seeded on-axis here
tunnel_center      = 0.0, 0.0 ; tunnel axis (y0, z0) (nm)
tunnel_k           = 8368    ; wall stiffness (kJ/mol/nm^2 = 20 kcal/mol/A^2)
tunnel_mouth_round = 0.2     ; mouth-corner fillet radius rho (nm)

; --- post-synthesis free runs (both release the C-terminus restraint) ---
ejection_steps     = 40   ; release restraint -> chain diffuses out the exit; 0 -> skip
dissociation_steps = 40   ; second free run -> chain drifts fully clear; 0 -> skip

; --- hardware / output ---
device = CPU ; ppn = 1 ; outdir = synth_out_cyl
```

### Analytic tunnel geometry

| Key | Default | Meaning |
|-----|---------|---------|
| `tunnel_radius` | `0.9` | Bore radius `r` (nm); ~3 CG beads wide. |
| `tunnel_length` | `10.0` | Bore length (nm). `x_exit = tunnel_x_lo + tunnel_length` (derived). |
| `tunnel_x_lo` | `0.0` | PTC / closed end (nm); the C-terminus is seeded on-axis here. |
| `tunnel_center` | `0.0, 0.0` | Tunnel axis `(y0, z0)` (nm). The axis runs along +x. |
| `tunnel_k` | `8368` | Wall stiffness (kJ/mol/nm² = 20 kcal/mol/Å²). |
| `tunnel_mouth_round` | `0.2` | Mouth-corner fillet radius `rho` (nm). |

### Post-elongation

| Key | Default | Meaning |
|-----|---------|---------|
| `ejection_steps` | `0` | Free run with the C-terminus restraint released so the finished chain diffuses out the exit (+x); `0` = skip. Use a **long** run so the chain can clear the tunnel. |
| `dissociation_steps` | `0` | Second free run (restraint still off) so the chain drifts fully clear of the exit; `0` = skip. |

### Shared kinetics & MD keys

These behave exactly as on the {doc}`synthesis_control` page (inherited from the shared
`RunParams`): `model` (default `hps_kr`; any IDP model works), `scale_factor`,
`codon_times`, `random_seed`, `ribosome_traffic` / `initiation_rate` (traffic
correction), `max_steps_per_stage` / `min_steps_per_stage` (testing-only), `constraints`
(`None` flexible default; the equilibrium-bond seed above also keeps `AllBonds` stable),
`restraint_k`, `minimize`, `dt` / `ref_t` / `tau_t` / `nstout`, `device` / `ppn`, and
`resume` (`auto`/`yes`/`no`; CLI `--fresh` forces `no` — see {doc}`synthesis_resume`).

```{warning}
`time_stage_1` / `time_stage_2` are accepted (inherited) but have **no effect** in the
cylinder runner: with a single MD segment per residue, each residue's step count comes
from its **whole** codon dwell, not a three-way split.
```

---

## Outputs

```text
<outdir>/
├── L_<L>/                  # one folder per residue L (single MD segment)
│   ├── traj.dcd            # (nascent-only) trajectory for that length
│   ├── traj_final.pdb      # last conformation (seeds the next residue)
│   └── traj.log, traj.psf, ...
├── ejection/               # free run, restraint off (if ejection_steps > 0)
├── dissociation/           # second free run (if dissociation_steps > 0)
├── dwell_times.dat         # per-residue: codon, sampled dwell (s), ns, integration steps
└── progress.log            # append-only DONE/RUNNING resume status
```

**Movie.** `cosmo-csp-movie -o <outdir>` stitches the per-length trajectories (it
auto-detects the flat `L_<L>/` layout used here vs. the 3-stage layout of `cosmo-csp`).

**Resume.** Like `cosmo-csp`, an interrupted cylinder run continues from the last
completed residue when re-invoked (`resume = auto`, on by default) — the schedule is
re-read from `dwell_times.dat` and the seed reloaded from the last `L_<L>/traj_final.pdb`,
tracked by `progress.log`. See {doc}`synthesis_resume`.

---

## See also

- {doc}`continuous_synthesis` — the explicit coarse-grained ribosome runner (`cosmo-csp`):
  the same codon kinetics with three MD sub-stages per residue and an explicit-bead tunnel.
- {doc}`synthesis_control` — the shared kinetics / MD control options in full.
- {doc}`../cosmo.csp` — the API reference.
