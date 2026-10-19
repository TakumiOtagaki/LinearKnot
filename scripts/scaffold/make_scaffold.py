#!/usr/bin/env python3
"""Build a pseudoknot-free scaffold G for LinearKnot from LinearPartition.

G = the LinearPartition MEA structure, keeping only pairs whose base-pair
probability is at least --threshold, then dropping helices shorter than
--min-helix (stacked runs of kept pairs).  The result is pseudoknot-free and
canonical by construction, so it is a valid scaffold.

Input: FASTA (or one sequence per line).  Output: records of ">name",
sequence, G, the input format of LinearKnot / linearknot profile /
linearknot access / linearknot sample.

    scripts/scaffold/build_linearpartition.sh
    scripts/scaffold/make_scaffold.py --threshold 0.5 seqs.fa > scaffolds.txt
"""

import argparse
import os
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_LP = os.path.join(ROOT, "build", "third_party", "LinearPartition", "linearpartition")


def read_sequences(handle):
    name, chunks, count = None, [], 0
    for raw in handle:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            if chunks:
                yield name, "".join(chunks)
            name, chunks = line[1:].strip(), []
            continue
        if name is None:
            count += 1
            yield "seq%d" % count, line
            continue
        chunks.append(line)
    if chunks:
        yield name, "".join(chunks)


def pairs_of(structure):
    stack, pairs = [], {}
    for pos, char in enumerate(structure, start=1):
        if char == "(":
            stack.append(pos)
        elif char == ")":
            i = stack.pop()
            pairs[i] = pos
    return pairs


def run_linearpartition(binary, sequence, beam, gamma):
    with tempfile.TemporaryDirectory() as tmp:
        bpp_path = os.path.join(tmp, "bpp.txt")
        cmd = [binary, "-V", "-M", "--gamma", str(gamma), "-b", str(beam), "-r", bpp_path]
        result = subprocess.run(cmd, input=sequence + "\n", capture_output=True,
                                text=True, check=True)
        lines = [line.strip() for line in result.stdout.splitlines() if line.strip()]
        mea = next(line for line in reversed(lines)
                   if len(line) == len(sequence) and set(line) <= set("().")
                   and line != sequence)
        bpp = {}
        with open(bpp_path) as handle:
            for line in handle:
                fields = line.split()
                if len(fields) == 3:
                    bpp[(int(fields[0]), int(fields[1]))] = float(fields[2])
    return mea, bpp


def scaffold(sequence, mea, bpp, threshold, min_helix):
    kept = {i: j for i, j in pairs_of(mea).items() if bpp.get((i, j), 0.0) >= threshold}
    if min_helix > 1:
        # Helix = maximal run of stacked kept pairs (i,j), (i+1,j-1), ...
        filtered = {}
        for i, j in kept.items():
            if kept.get(i - 1) == j + 1:
                continue  # not the outermost pair of its helix
            run = [(i, j)]
            while kept.get(run[-1][0] + 1) == run[-1][1] - 1:
                run.append((run[-1][0] + 1, run[-1][1] - 1))
            if len(run) >= min_helix:
                filtered.update(run)
        kept = filtered
    out = ["."] * len(sequence)
    for i, j in kept.items():
        out[i - 1], out[j - 1] = "(", ")"
    return "".join(out)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", nargs="?", default="-", help="FASTA file (default: stdin)")
    parser.add_argument("--threshold", type=float, default=0.5,
                        help="minimum base-pair probability of a kept MEA pair (default 0.5)")
    parser.add_argument("--min-helix", type=int, default=1,
                        help="drop helices with fewer stacked pairs (default 1: keep all)")
    parser.add_argument("--gamma", type=float, default=3.0, help="MEA gamma (default 3.0)")
    parser.add_argument("--beam", type=int, default=100, help="LinearPartition beam (default 100)")
    parser.add_argument("--linearpartition", default=DEFAULT_LP,
                        help="path to the linearpartition wrapper")
    args = parser.parse_args()
    if not os.path.exists(args.linearpartition):
        sys.exit("linearpartition not found; run scripts/scaffold/build_linearpartition.sh")
    handle = sys.stdin if args.input == "-" else open(args.input)
    for name, raw in read_sequences(handle):
        sequence = raw.upper().replace("T", "U")
        mea, bpp = run_linearpartition(args.linearpartition, sequence, args.beam, args.gamma)
        g = scaffold(sequence, mea, bpp, args.threshold, args.min_helix)
        print(">%s" % name)
        print(sequence)
        print(g)


if __name__ == "__main__":
    main()
