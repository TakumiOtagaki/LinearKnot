#!/usr/bin/env python3
"""Value snapshots for cost-only refactors of the beam kernel.

    regression_snapshot.py snapshot OUT_DIR     # run the CLIs, store outputs
    regression_snapshot.py compare OLD NEW      # compare two snapshots

A snapshot runs linearknot fold (MFE union structure, bpp), linearknot access
(ensemble free energy, accessibility) and linearknot profile (contexts) on a fixed
input set for several beam sizes, and records elapsed time.  `compare`
checks every number with a relative / absolute tolerance and every MFE
structure for equality, and prints the time ratio.  A refactor that only
changes cost must pass `compare` against the snapshot taken before it.
"""

import argparse
import json
import math
import os
import random
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BUILD = os.path.join(ROOT, "build")
MAKE_SCAFFOLD = os.path.join(ROOT, "scripts", "scaffold", "make_scaffold.py")
BEAMS = [10, 50, 1000]
ACCESS = "1,3,6"


def inputs():
    """Fixed random sequences with dense and sparse MEA scaffolds, plus witnesses."""
    records = []
    for n in (60, 120, 200, 300):
        for seed in range(3):
            rng = random.Random(4099 * n + seed)
            records.append(("n%d_s%d" % (n, seed), "".join(rng.choice("ACGU") for _ in range(n))))
    fasta = "".join(">%s\n%s\n" % r for r in records)
    out = []
    for threshold, tag in ((0.5, "dense"), (0.9, "sparse")):
        text = subprocess.run([sys.executable, MAKE_SCAFFOLD, "--threshold", str(threshold),
                               "--min-helix", "2"], input=fasta, capture_output=True,
                              text=True, check=True).stdout.split("\n")
        for a in range(0, len(text) - 2, 3):
            if text[a].startswith(">"):
                out.append((text[a][1:] + "_" + tag, text[a + 1], text[a + 2]))
    out.append(("wmb_witness", "GCGCGCGCGCGCGC", "(.(......)...)"))
    out.append(("vpr_guard", "UAUCUUAGGCAUUGAAAACGUCUG", "(.........)..(.......).."))
    out.append(("kissing", "GGGAAACCCAGGGAAACCCAAUCUUAGGCAUUGAAAACGUCUG",
                "(((...)))...........(.........)..(.......)."))
    return out


def run(cmd):
    started = time.perf_counter()
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=ROOT)
    elapsed = time.perf_counter() - started
    if result.returncode != 0:
        raise RuntimeError(" ".join(cmd[:3]) + ": " + result.stderr.strip())
    return result.stdout, elapsed


def numbers(text):
    out = []
    for token in text.replace(",", " ").replace(";", " ").replace("(", " ").replace(")", " ").split():
        try:
            out.append(float(token))
        except ValueError:
            pass
    return out


def snapshot(out_dir):
    os.makedirs(out_dir, exist_ok=True)
    data = {}
    for name, seq, g in inputs():
        for beam in BEAMS:
            key = "%s_b%d" % (name, beam)
            common = ["--sequence", seq, "--scaffold", g, "-b", str(beam)]
            bpp_path = os.path.join(out_dir, key + ".bpp")
            mfe, t1 = run([os.path.join(BUILD, "linearknot"), "fold"] + common + ["--bpp", bpp_path])
            acc, t2 = run([os.path.join(BUILD, "linearknot"), "access"] + common +
                          ["--access-len", ACCESS, "--probability"])
            capr, t3 = run([os.path.join(BUILD, "linearknot"), "profile"] + common)
            acc = "\n".join(l for l in acc.splitlines() if not l.startswith("# elapsed"))
            data[key] = {
                "mfe_structure": mfe.splitlines()[3].split()[0],
                "mfe_numbers": numbers(mfe.splitlines()[3]),
                "bpp": numbers(open(bpp_path).read().split("\n", 1)[1]),
                "raccess": numbers(acc),
                "capr": numbers(capr),
                "seconds": t1 + t2 + t3,
            }
            print(key, "%.2fs" % (t1 + t2 + t3), flush=True)
    with open(os.path.join(out_dir, "snapshot.json"), "w") as handle:
        json.dump(data, handle)


def close(a, b, rel=1e-9, abs_tol=1e-11):
    return abs(a - b) <= max(abs_tol, rel * max(abs(a), abs(b)))


def compare(old_dir, new_dir):
    old = json.load(open(os.path.join(old_dir, "snapshot.json")))
    new = json.load(open(os.path.join(new_dir, "snapshot.json")))
    failures = 0
    t_old = t_new = 0.0
    for key in sorted(old):
        a, b = old[key], new.get(key)
        if b is None:
            print("MISSING", key)
            failures += 1
            continue
        t_old += a["seconds"]
        t_new += b["seconds"]
        if a["mfe_structure"] != b["mfe_structure"]:
            print("MFE structure differs", key)
            failures += 1
        for field in ("mfe_numbers", "bpp", "raccess", "capr"):
            x, y = a[field], b[field]
            if len(x) != len(y):
                print("%s length differs %s: %d vs %d" % (field, key, len(x), len(y)))
                failures += 1
                continue
            worst = max((abs(p - q) for p, q in zip(x, y)), default=0.0)
            if not all(close(p, q) for p, q in zip(x, y)):
                print("%s differs %s (max abs %.3e)" % (field, key, worst))
                failures += 1
    print("cases=%d failures=%d time old=%.1fs new=%.1fs ratio=%.3f" %
          (len(old), failures, t_old, t_new, t_new / t_old if t_old else math.nan))
    return failures == 0


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="mode", required=True)
    s = sub.add_parser("snapshot")
    s.add_argument("out_dir")
    c = sub.add_parser("compare")
    c.add_argument("old_dir")
    c.add_argument("new_dir")
    args = parser.parse_args()
    if args.mode == "snapshot":
        snapshot(args.out_dir)
    else:
        sys.exit(0 if compare(args.old_dir, args.new_dir) else 1)


if __name__ == "__main__":
    main()
