#!/usr/bin/env python3
import statistics
import sys


def load(path):
    with open(path) as handle:
        return [float(x.strip()) for x in handle if x.strip()]


def pct(vals, p):
    if not vals:
        return 0.0
    vals = sorted(vals)
    k = (len(vals) - 1) * p / 100.0
    f = int(k)
    c = min(f + 1, len(vals) - 1)
    if f == c:
        return vals[f]
    return vals[f] + (vals[c] - vals[f]) * (k - f)


def top3_or_zeros(vals):
    if not vals:
        return [0.0, 0.0, 0.0]
    top3 = sorted(vals, reverse=True)[:3]
    while len(top3) < 3:
        top3.append(0.0)
    return top3


def print_compare_summary(compared, matched, failed, skipped, total, abs_tol, rel_tol):
    print(f"compare_total={total}")
    print(f"compare_compared={compared}")
    print(f"compare_matched={matched}")
    print(f"compare_failed={failed}")
    print(f"compare_skipped={skipped}")
    print(f"pf_abs_tol={abs_tol}")
    print(f"pf_rel_tol={rel_tol}")


def print_abs_stats(abs_vals):
    top3 = top3_or_zeros(abs_vals)
    print(f"pf_abs_top3_max={top3[0]},{top3[1]},{top3[2]}")
    if not abs_vals:
        print("pf_abs_median=0")
        print("pf_abs_p90=0")
        print("pf_abs_p95=0")
        return
    print(f"pf_abs_median={statistics.median(abs_vals)}")
    print(f"pf_abs_p90={pct(abs_vals, 90)}")
    print(f"pf_abs_p95={pct(abs_vals, 95)}")


def print_rel_stats(rel_vals):
    top3 = top3_or_zeros(rel_vals)
    print(f"pf_rel_top3_max={top3[0]},{top3[1]},{top3[2]}")
    if not rel_vals:
        print("pf_rel_median=0")
        print("pf_rel_p90=0")
        print("pf_rel_p95=0")
        return
    print(f"pf_rel_median={statistics.median(rel_vals)}")
    print(f"pf_rel_p90={pct(rel_vals, 90)}")
    print(f"pf_rel_p95={pct(rel_vals, 95)}")


def parse_args(argv):
    if len(argv) != 10:
        raise SystemExit(
            "usage: pf_alignment_stats.py <abs_path> <rel_path> <compared> <min_req> <abs_tol> <rel_tol> <matched> <failed> <skipped> <total>"
        )
    abs_path, rel_path, compared, min_req, abs_tol, rel_tol, matched, failed, skipped, total = argv
    return {
        "abs_path": abs_path,
        "rel_path": rel_path,
        "compared": int(compared),
        "min_req": int(min_req),
        "abs_tol": abs_tol,
        "rel_tol": rel_tol,
        "matched": int(matched),
        "failed": int(failed),
        "skipped": int(skipped),
        "total": int(total),
    }


def check_gates(compared, min_req, failed):
    if compared < min_req:
        raise SystemExit(f"error: compared {compared} < required {min_req}")
    if failed != 0:
        raise SystemExit("error: pf compare failed")


def main(argv):
    args = parse_args(argv)
    abs_vals = load(args["abs_path"])
    rel_vals = load(args["rel_path"])
    print_compare_summary(
        args["compared"],
        args["matched"],
        args["failed"],
        args["skipped"],
        args["total"],
        args["abs_tol"],
        args["rel_tol"],
    )
    print_abs_stats(abs_vals)
    print_rel_stats(rel_vals)
    check_gates(args["compared"], args["min_req"], args["failed"])


if __name__ == "__main__":
    main(sys.argv[1:])
