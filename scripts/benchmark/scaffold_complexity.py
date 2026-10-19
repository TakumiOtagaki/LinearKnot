#!/usr/bin/env python3
"""How the beam kernel's cost depends on the scaffold G.

For random sequences and scaffolds of varied shape (make_scaffold.py with
several thresholds / minimum helix lengths), record features of G and the
beam inside cost measured by scfg2_online_beam_work_profile_tool: elapsed
time and the deterministic work counters (summed probes / visits).

G features
  n, g_pairs, helices (maximal stacked runs), hairpins, branch_points
  (pairs or the exterior with >= 2 pair children), max_depth,
  unpaired_inside (unpaired bases enclosed by a G pair), unpaired_outside,
  loops (G loops holding an unpaired base), crossable_pairs: pairs of
  unpaired bases in different G loops.  Any pair joining two different
  loops crosses some G pair, so this is the VP candidate space.

    scripts/benchmark/scaffold_complexity.py --out results/x.tsv
"""

import argparse
import os
import random
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
TOOL = os.path.join(ROOT, "build", "scfg2_online_beam_work_profile_tool")
MAKE_SCAFFOLD = os.path.join(ROOT, "scripts", "scaffold", "make_scaffold.py")
PARAMS = os.path.join(ROOT, "params", "rna_Turner04.par")


def features(g):
    n = len(g)
    partner = [0] * (n + 2)
    stack = []
    for pos, char in enumerate(g, start=1):
        if char == "(":
            stack.append(pos)
        elif char == ")":
            i = stack.pop()
            partner[i], partner[pos] = pos, i
    pairs = [(i, partner[i]) for i in range(1, n + 1) if partner[i] > i]
    helices = sum(1 for i, j in pairs if partner[i - 1] != j + 1 or i == 1)
    # Loop id of each unpaired base: the innermost enclosing pair (0 = exterior).
    loop_of = [0] * (n + 1)
    children = {0: 0}
    depth, max_depth, enclosing = 0, 0, [0]
    for pos in range(1, n + 1):
        if partner[pos] > pos:
            children[enclosing[-1]] = children.get(enclosing[-1], 0) + 1
            enclosing.append(pos)
            children.setdefault(pos, 0)
            depth += 1
            max_depth = max(max_depth, depth)
        elif partner[pos] and partner[pos] < pos:
            enclosing.pop()
            depth -= 1
        else:
            loop_of[pos] = enclosing[-1]
    hairpins = sum(1 for i, _ in pairs if children.get(i, 0) == 0)
    branch_points = sum(1 for count in children.values() if count >= 2)
    sizes = {}
    unpaired_inside = unpaired_outside = 0
    for pos in range(1, n + 1):
        if partner[pos]:
            continue
        sizes[loop_of[pos]] = sizes.get(loop_of[pos], 0) + 1
        if loop_of[pos]:
            unpaired_inside += 1
        else:
            unpaired_outside += 1
    unpaired = unpaired_inside + unpaired_outside
    crossable = unpaired * (unpaired - 1) // 2 - sum(s * (s - 1) // 2 for s in sizes.values())
    return {
        "n": n, "g_pairs": len(pairs), "helices": helices, "hairpins": hairpins,
        "branch_points": branch_points, "max_depth": max_depth,
        "unpaired_inside": unpaired_inside, "unpaired_outside": unpaired_outside,
        "loops": len(sizes), "crossable_pairs": crossable,
    }


def work(seq, g, beam):
    result = subprocess.run([TOOL, seq, g, str(beam), "inside", PARAMS, "1"],
                            capture_output=True, text=True, cwd=ROOT, check=True)
    lines = result.stdout.splitlines()
    header, values = lines[0].split("\t"), lines[1].split("\t")
    row = dict(zip(header, values))
    probes = sum(int(float(v)) for k, v in row.items()
                 if k.endswith("_probes") or k.endswith("_visits"))
    family_visits = 0
    for line in lines[3:]:
        fields = line.split("\t")
        if len(fields) == 3 and fields[1].isdigit():
            family_visits += int(fields[1])
    return float(row["elapsed_mean_ms"]) / 1000.0, probes, family_visits


def job(task):
    name, seq, g, beam, threshold, min_helix = task
    seconds, probes, visits = work(seq, g, beam)
    row = {"name": name, "threshold": threshold, "min_helix": min_helix, "beam": beam,
           "seconds": round(seconds, 4), "work_probes": probes, "item_visits": visits}
    row.update(features(g))
    return row


def scaffolds(records, threshold, min_helix):
    fasta = "".join(">%s\n%s\n" % (name, seq) for name, seq in records)
    out = subprocess.run([sys.executable, MAKE_SCAFFOLD, "--threshold", str(threshold),
                          "--min-helix", str(min_helix)],
                         input=fasta, capture_output=True, text=True, check=True).stdout
    lines = out.split("\n")
    return {lines[a][1:]: (lines[a + 1], lines[a + 2])
            for a in range(0, len(lines) - 2, 3) if lines[a].startswith(">")}


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", required=True)
    parser.add_argument("--lengths", default="200,400,800")
    parser.add_argument("--count", type=int, default=15)
    parser.add_argument("--thresholds", default="0.3,0.6,0.9")
    parser.add_argument("--min-helix", default="1,3")
    parser.add_argument("--beam", type=int, default=100)
    parser.add_argument("--jobs", type=int, default=8)
    args = parser.parse_args()
    tasks = []
    for n in [int(v) for v in args.lengths.split(",")]:
        records = []
        for s in range(args.count):
            rng = random.Random(7919 * n + s)
            records.append(("n%d_s%d" % (n, s), "".join(rng.choice("ACGU") for _ in range(n))))
        for threshold in [float(v) for v in args.thresholds.split(",")]:
            for min_helix in [int(v) for v in args.min_helix.split(",")]:
                scaf = scaffolds(records, threshold, min_helix)
                for name, _ in records:
                    seq, g = scaf[name]
                    tasks.append((name, seq, g, args.beam, threshold, min_helix))
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    columns = None
    with open(args.out, "w") as out, ProcessPoolExecutor(args.jobs) as pool:
        for row in pool.map(job, tasks):
            if columns is None:
                columns = list(row)
                out.write("\t".join(columns) + "\n")
            out.write("\t".join(str(row[c]) for c in columns) + "\n")
            out.flush()


if __name__ == "__main__":
    main()
