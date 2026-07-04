# COSMO regression benchmarks

A lightweight regression test that guards the physics of each force field.
It builds a coarse-grained model from a reference structure
(`benchmarks/asyn.pdb`, alpha-synuclein), evaluates the potential
energy of the **initial, un-minimized** coordinates, and decomposes it into the
individual force groups. The recorded values become a known-correct reference;
after changing the code you re-run in *check* mode to confirm the energies are
unchanged.

Models covered: `hps_urry`, `hps_kr`, `hps_ss`, `mpipi`.

## Why this is deterministic

The energy is evaluated as a pure function of the input coordinates:

- `minimize=False` — coordinates come straight from the PDB.
- no periodic boundary conditions — nonbonded forces use a fixed non-periodic
  cutoff, so no box geometry enters.
- OpenMM `Reference` platform (double precision) — reproducible across machines
  and GPUs.

As a result, check mode reproduces the reference to the bit (`|diff| = 0`).

## Usage

Check the current code against the stored reference (the normal case):

```bash
python benchmarks/benchmark_energies.py
```

Exit code is `0` if every force group of every model matches within tolerance,
`1` otherwise — so it can be dropped into CI.

Regenerate the reference — **only when the current code is known to be
correct** (e.g. after an intentional, validated change to the force field):

```bash
python benchmarks/benchmark_energies.py --generate
```

Useful options:

| Option | Default | Meaning |
|--------|---------|---------|
| `--pdb PATH` | `benchmarks/asyn.pdb` | input structure |
| `--reference PATH` | `benchmarks/reference_energies.json` | reference file |
| `--models M [M ...]` | all four | subset to test |
| `--rtol FLOAT` | `1e-4` | relative tolerance |
| `--atol FLOAT` | `1e-3` | absolute tolerance (kJ/mol) |

## Files

- `benchmark_energies.py` — the benchmark (generate / check modes).
- `reference_energies.json` — committed reference values, with a `_meta` block
  recording the structure and OpenMM version used to generate them.

## Interpreting a failure

A mismatch means the initial energy of at least one force group changed. That is
expected and fine *if* you intentionally changed that force field — regenerate
the reference. If you did **not** intend to change the physics, the diff is
pointing at an unintended regression: the per-group breakdown tells you which
force term moved.
