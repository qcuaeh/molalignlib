#!/usr/bin/env python3
"""
Run Molecule.confo_rmsd over every pose pair listed in a reference CSV and
write the results to a CSV file with the same format (pose1,pose2,rmsd with
six decimals) and the same row order as the reference.

The pair list and its order are taken from the reference CSV. Its rows are
not sorted (pairs and molecules are interleaved), so this is the only way
to reproduce the order exactly.

The new values are also compared with the reference and a summary is
printed (rows that differ by more than --tol, worst cases).

Usage
-----
    python test_conformsd.py reference.csv -o results.csv
    python test_conformsd.py reference.csv -d path/to/poses -o out.csv -j 8

Pose identifiers in the CSV are file paths relative to the data directory
(--data-dir, the current directory by default), including the extension (e.g.
mol_0001/pose_3.sdf); the file read is <data_dir>/<pose>, and its format is
taken from that extension.
"""

import argparse
import csv
import math
import sys
import time
from multiprocessing import Pool
from pathlib import Path

from molalignlib import Molecule

# --------------------------------------------------------------------------- #
# Pair list                                                                   #
# --------------------------------------------------------------------------- #

def pairs_from_reference(path):
    """Return [(pose1, pose2, ref_rmsd), ...] in the file order."""
    pairs = []
    with open(path, newline="") as f:
        reader = csv.reader(f)
        header = next(reader)
        if [h.strip() for h in header[:3]] != ["pose1", "pose2", "rmsd"]:
            raise SystemExit("Unexpected header in {}: {}".format(path, header))
        for row in reader:
            if not row:
                continue
            try:
                ref = float(row[2])
            except (IndexError, ValueError):
                ref = math.nan
            pairs.append((row[0], row[1], ref))
    return pairs


# --------------------------------------------------------------------------- #
# RMSD calculation (runs in worker processes)                                 #
# --------------------------------------------------------------------------- #

_CFG = {}
_CACHE = {}


def _init_worker(cfg):
    _CFG.clear()
    _CFG.update(cfg)
    _CACHE.clear()


def _load(pose_id):
    mol = _CACHE.get(pose_id)
    if mol is None:
        # Keep only the poses of the current molecule in the cache.
        prefix = pose_id.split("/")[0] + "/"
        for key in [k for k in _CACHE if not k.startswith(prefix)]:
            del _CACHE[key]
        path = Path(_CFG["data_dir"]) / pose_id
        mol = Molecule.from_file(str(path))
        _CACHE[pose_id] = mol
    return mol


def _compute(task):
    index, pose1, pose2 = task
    t0 = time.perf_counter()
    try:
        mol0 = _load(pose1)
        mol1 = _load(pose2)
        result = mol0.confo_rmsd(mol1, **_CFG["options"])
        rmsd, error = result.rmsd, None
    except Exception as exc:  # report and continue with the next pair
        rmsd, error = math.nan, "{}: {}".format(type(exc).__name__, exc)
    return index, rmsd, error, time.perf_counter() - t0


# --------------------------------------------------------------------------- #
# Main                                                                        #
# --------------------------------------------------------------------------- #

def parse_args():
    ap = argparse.ArgumentParser(
        description="Test confo_rmsd on the pose pairs of a reference CSV.")
    ap.add_argument("reference",
                    help="reference CSV giving the pairs and their order")
    ap.add_argument("-d", "--data-dir", default=".",
                    help="directory the pose paths in the CSV are relative "
                         "to (default: current directory)")
    ap.add_argument("-o", "--output", default="rmsds_confo.csv",
                    help="output CSV (default: rmsds_confo.csv)")
    ap.add_argument("-j", "--jobs", type=int, default=1,
                    help="number of worker processes (default: 1)")
    ap.add_argument("--tol", type=float, default=1e-4,
                    help="tolerance in Å for the comparison with the "
                         "reference (default: 1e-4)")
    ap.add_argument("--limit", type=int, default=None,
                    help="only run the first N pairs (for quick tests)")
    # confo_rmsd options
    ap.add_argument("--heavy-only", action="store_true")
    ap.add_argument("--mass-weight", action="store_true")
    ap.add_argument("--mirror", action="store_true")
    ap.add_argument("--use-bond-type", action="store_true")
    ap.add_argument("--rebond", action="store_true")
    ap.add_argument("--min-freq", type=int, default=100)
    ap.add_argument("--max-trials", type=int, default=10000)
    ap.add_argument("--random", action="store_true",
                    help="seed the random search from the clock")
    return ap.parse_args()


def main():
    args = parse_args()

    pairs = pairs_from_reference(args.reference)
    if args.limit is not None:
        pairs = pairs[:args.limit]
    if not pairs:
        raise SystemExit("No pairs to compute.")

    cfg = {
        "data_dir": args.data_dir,
        "options": dict(
            align=True,
            remap=True,
            heavy_only=args.heavy_only,
            mass_weight=args.mass_weight,
            mirror=args.mirror,
            use_bond_type=args.use_bond_type,
            rebond=args.rebond,
            random=args.random,
            min_freq=args.min_freq,
            max_trials=args.max_trials,
        ),
    }

    n = len(pairs)
    tasks = [(i, p1, p2) for i, (p1, p2, _) in enumerate(pairs)]
    rmsds = [math.nan] * n
    times = [0.0] * n
    errors = []

    print("Computing {} pairs with {} process(es)...".format(n, args.jobs),
          file=sys.stderr)
    t_start = time.perf_counter()

    if args.jobs > 1:
        pool = Pool(args.jobs, initializer=_init_worker, initargs=(cfg,))
        # Tasks are in reference order, so chunks mostly stay within one
        # molecule and the per-process cache is reused.
        results = pool.imap_unordered(_compute, tasks, chunksize=16)
    else:
        pool = None
        _init_worker(cfg)
        results = map(_compute, tasks)

    try:
        for done, (i, rmsd, error, dt) in enumerate(results, start=1):
            rmsds[i] = rmsd
            times[i] = dt
            if error:
                errors.append((i, error))
                print("ERROR {} {}: {}".format(pairs[i][0], pairs[i][1], error),
                      file=sys.stderr)
            if done % 500 == 0 or done == n:
                print("  {}/{} pairs, {:.1f} s".format(
                    done, n, time.perf_counter() - t_start), file=sys.stderr)
    finally:
        if pool is not None:
            pool.close()
            pool.join()

    # Output in the same format and order as the reference.
    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f, lineterminator="\n")
        writer.writerow(["pose1", "pose2", "rmsd"])
        for (p1, p2, _), rmsd in zip(pairs, rmsds):
            writer.writerow([p1, p2, "{:.6f}".format(rmsd)])

    wall = time.perf_counter() - t_start
    print("\nWrote {} rows to {}".format(n, args.output))
    print("Wall time {:.1f} s, CPU time in confo_rmsd {:.1f} s, "
          "slowest pair {:.2f} s".format(wall, sum(times), max(times)))
    print("Errors: {}".format(len(errors)))

    # Comparison with the reference.
    diffs = [(abs(new - ref), new - ref, i)
             for i, ((_, _, ref), new) in enumerate(zip(pairs, rmsds))
             if not (math.isnan(new) or math.isnan(ref))]
    if diffs:
        above = [d for d in diffs if d[0] > args.tol]
        lower = sum(1 for d in above if d[1] < 0)
        print("\nComparison with {} (tol {:g} Å):".format(
            args.reference, args.tol))
        print("  compared pairs: {}".format(len(diffs)))
        print("  equal within tol: {}".format(len(diffs) - len(above)))
        print("  differ: {} ({} lower, {} higher than reference)".format(
            len(above), lower, len(above) - lower))
        print("  max |diff|: {:.6f} Å".format(max(d[0] for d in diffs)))
        for absd, d, i in sorted(above, reverse=True)[:10]:
            p1, p2, ref = pairs[i]
            print("    {} {}  ref {:.6f}  new {:.6f}  diff {:+.6f}".format(
                p1, p2, ref, rmsds[i], d))

    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
