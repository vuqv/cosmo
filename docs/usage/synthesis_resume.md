# Resuming long synthesis runs

A production {doc}`continuous synthesis <continuous_synthesis>` run synthesizes a full
chain one residue at a time, each as three MD sub-stages, at production step counts
(~10⁵–10⁶ steps/stage). A full chain is therefore **hours to days** of wall time — longer
than a typical scheduler slot, and exposed to node failures and pre-emption. **cosmo's CSP
runs are resumable:** re-invoking `cosmo-csp` on an interrupted output directory continues
from the last completed residue instead of restarting from `L0`.

This page covers *operating* a long run: how resume works, the on-disk artifacts, the
config knob, and an HPC requeue pattern. For the physics and the full option reference see
{doc}`continuous_synthesis` and {doc}`synthesis_control`.

```{note}
Resume works for **both** synthesis runners: `cosmo-csp` (explicit ribosome, three
sub-stages per residue) and `cosmo-cylinder` (analytic tunnel, one MD segment per residue).
The `resume` knob, `progress.log`, and the up-front cost report are identical; the only
differences are internal — the cylinder's `dwell_times.dat` carries one step count per
residue (no `#PTC` header, since its tunnel geometry is cheap and deterministic) and its
per-residue finals live in a flat `L_<L>/traj_final.pdb` layout. Everything below applies
to both; examples use `cosmo-csp` unless noted.
```

---

## TL;DR

```ini
# csp.ini  — resume is on by default (nothing to add)
resume = auto        # auto (default) | yes | no
```

```bash
cosmo-csp -f csp.ini              # launch; interrupted at any point...
cosmo-csp -f csp.ini              # ...re-run the SAME command -> continues where it left off
cosmo-csp -f csp.ini --fresh      # ignore any prior run and start over
```

A launched run prints its **exact total planned MD-step count before the first MD step**,
so you can size a scheduler wall-clock request:

```text
[schedule] 140 residues, 23,451,000 planned MD steps (~234.5 ns simulated at dt=0.01 ps; ...)
```

---

## How it works

The insight is that **almost nothing crosses a residue boundary.** Reading the driver
loop `for L in range(L0, L_max+1)`, only two pieces of state carry from one residue to the
next, and both are already on disk:

| State | Recovery on resume |
|---|---|
| the previous residue's final coordinates (the seed for residue `L`) | **reloaded** from `L_<L-1>/traj_final.pdb` |
| the per-residue kinetic schedule `(steps₁, steps₂, steps₃)` | **re-read** from the persisted schedule (no RNG redraw) |

Everything else the loop consumes (ribosome load, PTC anchors, the codon kinetic lists) is
cheap and **recomputed** from the config at startup — there is no STRIDE or native-contact
precompute in cosmo's sequence-based IDP model. The one expensive-but-deterministic
quantity — the optimized PTC restraint geometry — is **persisted at first launch and
re-read on resume** (see below), never re-solved.

So there is **no heavyweight simulation checkpoint**. The resume unit is the **residue**,
and all resume needs to record is *how far the run got*.

### The two on-disk artifacts

At the output root, alongside the per-residue `L_<L>/` directories:

`dwell_times.dat` — the **immutable plan**
: The full per-residue 3-stage schedule, drawn **once** from the seeded generator before
  the main loop, plus a machine-readable `#PTC` header carrying the restraint geometry
  (the A-/P-site target points and the tunnel-wall plane) at full float precision. On
  resume this file is **re-read, never rewritten** — the schedule is *not* redrawn and the
  PTC geometry is *not* re-solved, so both are pinned identical to the residues already on
  disk.

`progress.log` — the **mutable status**
: A human-readable, append-only record of how far the run got:

  ```text
  # cosmo csp progress log -- schema 1
  L_001 RUNNING
  L_001 DONE
  ...
  L_041 DONE
  L_042 RUNNING       <- crash here: L_042 is the partial unit dropped on resume
  ```

  A residue is marked `RUNNING` when it starts and `DONE` only after **all three stages**
  finish. The **`DONE` line is the commit point**: a short append is effectively atomic on
  POSIX, so every crash is recoverable — die before `DONE` and that unit is `RUNNING`
  (dropped and redone); die after `DONE` and resume picks up at the next residue. The
  post-synthesis `ejection` / `dissociation` phases appear as their own units.

On resume, the driver takes the last status per unit, drops the (at most one) `RUNNING`
unit **and its directory**, reloads the last `DONE` residue's final structure as the seed,
and continues from `max(DONE) + 1`. Redone work is bounded to a **single residue** (≤ 3
stages) — negligible against a run measured in hours to days.

```{tip}
Both files are plain text — `cat progress.log` to see exactly how far a run got, or
`tail -f` it to watch a live run advance.
```

---

## The `resume` knob

Set it in the `[OPTIONS]` section of `csp.ini` / `cylinder.ini`, or override on the CLI.

| Value | Behavior |
|---|---|
| `auto` *(default)* | Resume **iff** a `progress.log` is present under `outdir`; otherwise start fresh. |
| `yes` | Require a resumable run — **error** if no `progress.log` is found. |
| `no` | Always start fresh (the pre-resume behavior). |

```bash
cosmo-csp -f csp.ini --no-resume    # or --fresh; forces resume = no for this invocation
```

Because `auto` is the default, **re-running the same command on an existing output
directory resumes** rather than silently overwriting it — the intended behavior for a
requeued job. To deliberately redo a directory from scratch, use `--fresh` (or point at a
new `outdir`).

---

## What resume guarantees (and what it does not)

**Guaranteed — the kinetic schedule is identical to an uninterrupted run.** The schedule
is materialized once and re-read, so `dwell_times.dat` is byte-identical across the
interruption, and the resumed conformation continues from the last completed residue's
final structure (0.001 nm PDB precision, far below the model's thermal noise).

**Not guaranteed — bit-reproducible MD micro-trajectories.** The per-stage dynamics are
stochastic and *not* bit-reproducible across a process boundary: the Langevin thermostat
and velocity initialization draw from OpenMM's own unseeded platform RNG. The intended and
achievable guarantee is the **schedule + conformational continuity**, not a bit-exact
trajectory.

Two consequences worth knowing:

- **No config fingerprint.** Because the schedule is *loaded*, not redrawn, a changed
  `random_seed` or retuned kinetic knob on resume cannot corrupt the tail — the RNG is
  never consulted again. (A genuinely different *force* config on resume is a self-inflicted
  footgun the tool does not police; keep the config stable across a requeue.)
- **Extending a run is a fresh run.** The schedule is fixed at first launch to cover exactly
  `L0..L_max`. Resuming with a **larger `L_max`** is rejected with an actionable error
  (`persisted schedule covers L=1..10 but this run asks for L=1..20 … Extending a run is a
  fresh run`). Decide the final length before launching.

```{warning}
**Presence guard.** `progress.log` records intent; the disk records reality. Before
resuming, CSP verifies that **every** length `L0..last_done` still has its
`traj_final.pdb` on disk. If one is missing (a directory deleted, or lost to a scratch
purge), the resume **aborts** naming the offending length rather than silently continuing
and leaving a permanent hole. Re-run fresh, or restore the missing length.
```

```{note}
**Concurrent invocations.** Two `cosmo-csp` processes on the same `outdir` would race
`progress.log`; don't do that. Resume assumes one writer at a time.
```

---

## Worked example

A fast, self-contained demonstration (small step clamp, CPU) — synthesize `L = 1..8`,
interrupt after residue 4, and resume to completion:

```bash
# 1. Launch. Kill it (Ctrl-C / scheduler pre-emption) after ~residue 4.
cosmo-csp -f csp.ini
# progress.log now ends:  L_004 DONE / L_005 RUNNING   (residue 5 was in flight)

# 2. Re-run the SAME command — resume is automatic.
cosmo-csp -f csp.ini
```

The resume run opens by reporting where it picks up:

```text
[resume] 4 length(s) complete on disk; continuing from L=5.
Done. Synthesized 1 -> 8.
```

`dwell_times.dat` is unchanged (re-read, not rewritten), any partial `L_005/` directory is
dropped and cleanly rebuilt from its persisted schedule row, and the run reaches `L_max`.

---

## HPC requeue pattern

Because a re-run resumes automatically, a self-requeuing SLURM job survives wall-clock
limits with no bookkeeping — the same `cosmo-csp` line runs each slot and advances the run:

```bash
#!/bin/bash
#SBATCH --job-name=csp_synth
#SBATCH --time=48:00:00
#SBATCH --signal=B:USR1@120          # warn the script 120 s before the wall-clock kill
#SBATCH --gres=gpu:1

# Re-queue this same script if we hit the time limit, then continue running.
trap 'scontrol requeue "$SLURM_JOB_ID"' USR1

# resume = auto (default): first slot starts fresh, every requeue continues.
cosmo-csp -f csp.ini &
wait
```

Each slot: the schedule + PTC geometry are read back from `dwell_times.dat`, the last
completed residue's final is reloaded, and synthesis continues. Point `outdir` at
persistent (not node-local) storage so the partial output survives between slots.

---

## See also

- {doc}`continuous_synthesis` — the CSP runner, physics, and full output layout.
- {doc}`synthesis_control` — the canonical per-key `csp.ini` reference (incl. `resume`).
- `cosmo.csp.resume` — the module implementing the schedule/progress persistence and the
  resume actions.
