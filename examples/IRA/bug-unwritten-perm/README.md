# Unwritten permutation structures (PR 16)

Two 38-point synthetic structures that make `ira_unify` return `ierr = 0`
with an unwritten `p_perm` entry. Each is a compact cubic core plus
isolated points at different distances (A: 5 points at about 6; B: 1
point at about 52).

## Files

- `ira_bug_structures.xyz` — saved two-frame xyz (structure A, then B)
- `reproducer.py` — builds the measured pair (A seed 1 / 5@6, B seed 2 /
  1@80) through IRA's own `ira_mod` and prints the permutation diagnostics

The script is the generator for the numbers below. It does no caller-side
index arithmetic.

## Run

From the repository root, after a build that produces `libira.so`:

    IRA_LIB=/path/to/libira.so PYTHONPATH=interface \
      python examples/IRA/bug-unwritten-perm/reproducer.py

To write the generated pair:

    python examples/IRA/bug-unwritten-perm/reproducer.py --write-xyz \
      examples/IRA/bug-unwritten-perm/generated.xyz

Unpatched `libira_match` returns success, `hd = 9.500803`, 37 distinct
perm entries of 38, index 37 never assigned. With the permutation check
on this branch the same call reports that the permutation is not a
bijection.

A consumer-side off-by-one (subtracting 1 from a permutation that is
already C-style) looks the same from outside and is not this bug.
