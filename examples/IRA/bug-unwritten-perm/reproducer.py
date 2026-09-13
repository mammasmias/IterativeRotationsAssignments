"""Build the 38-point pair that leaves an unwritten p_perm entry.

Uses IRA's own ira_mod (no caller-side index arithmetic). From the
repository root, after the shared library is built:

    IRA_LIB=/path/to/libira.so PYTHONPATH=interface \\
      python examples/IRA/bug-unwritten-perm/reproducer.py

Optional: write the generated pair (A seed 1 / 5@6, B seed 2 / 1@80):

    python examples/IRA/bug-unwritten-perm/reproducer.py --write-xyz \\
      examples/IRA/bug-unwritten-perm/generated.xyz
"""
import argparse
import os
import sys

import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.abspath(os.path.join(_HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(_REPO, "interface"))
import ira_mod


def core(k, spacing=1.1):
    pts, i = [], 0
    side = int(np.ceil(k ** (1 / 3))) + 1
    for a in range(side):
        for b in range(side):
            for c in range(side):
                if i < k:
                    pts.append([a * spacing, b * spacing, c * spacing])
                    i += 1
    return np.array(pts[:k], dtype=np.float64)


def frag(n, n_loose, spread, seed):
    rng = np.random.default_rng(seed)
    X = [core(n - n_loose)]
    for _ in range(n_loose):
        v = rng.normal(size=3)
        v /= np.linalg.norm(v)
        X.append((v * rng.uniform(0.6 * spread, spread))[None, :])
    Y = np.vstack(X)[:n]
    return np.ascontiguousarray(Y - Y.mean(axis=0))


def write_xyz(path, frames):
    with open(path, "w") as fh:
        for name, coords in frames:
            fh.write(f"{len(coords)}\n{name}\n")
            for row in coords:
                fh.write(
                    "H {: .17e} {: .17e} {: .17e}\n".format(row[0], row[1], row[2])
                )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write-xyz",
        metavar="PATH",
        help="write the generated A/B pair as a two-frame xyz",
    )
    args = parser.parse_args()

    n = 38
    A = frag(n, 5, 6.0, 1)
    B = frag(n, 1, 80.0, 2)
    for name, X in (("A", A), ("B", B)):
        r = np.linalg.norm(X - X.mean(axis=0), axis=1)
        print(f"{name}: max radius {r.max():8.3f}")

    if args.write_xyz:
        write_xyz(args.write_xyz, (("structure A", A), ("structure B", B)))
        print(f"wrote {args.write_xyz}")

    shlib = os.environ.get("IRA_LIB")
    ira = ira_mod.IRA(shlib=shlib) if shlib else ira_mod.IRA()
    typ = np.ones(n, dtype=np.int32)
    rot, trans, perm, hd = ira.match(n, typ, A, n, typ, B, 1.8)
    p = np.asarray(perm).ravel()
    print(f"cerr: none raised (success)   hd {float(hd):.6f}")
    print(f"permutation length {len(p)}  min {p.min()}  max {p.max()}")
    print(f"distinct entries {len(set(p.tolist()))} of {n}")
    missing = sorted(set(range(n)) - set(p.tolist()))
    dupes = sorted({v for v in p.tolist() if list(p).count(v) > 1})
    print(f"missing indices {missing}")
    print(f"duplicated indices {dupes}")
    print("A permutation that is not a bijection was returned as success.")


if __name__ == "__main__":
    main()
