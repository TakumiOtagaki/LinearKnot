#!/usr/bin/env python3
"""Benchmarks for linearknot access / linearknot profile.

Modes (each writes one TSV to --out):
  scaling   time and peak memory of linearknot access against length and beam
  accuracy  beam vs exact: max / mean absolute error of contexts and
            accessibility
  pkfree    empty scaffold vs LinearRaccess exact (-beam=0, turner2004)

Random sequences are uniform over ACGU; scaffolds come from
scripts/scaffold/make_scaffold.py (LinearPartition MEA pairs with probability
>= --threshold, helices of at least --min-helix pairs).
"""

import argparse
import os
import random
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
BUILD = os.path.join(ROOT, "build")
MAKE_SCAFFOLD = os.path.join(ROOT, "scripts", "scaffold", "make_scaffold.py")
CONTEXTS = ["Bulge", "Exterior", "Hairpin", "Internal", "Multiloop", "PseudoLoop", "Stem"]


def random_sequences(lengths, seeds):
    for n in lengths:
        for seed in seeds:
            rng = random.Random(1000003 * n + seed)
            yield "n%d_s%d" % (n, seed), "".join(rng.choice("ACGU") for _ in range(n))


def scaffolds(records, threshold, min_helix):
    fasta = "".join(">%s\n%s\n" % (name, seq) for name, seq in records)
    result = subprocess.run(
        [sys.executable, MAKE_SCAFFOLD, "--threshold", str(threshold),
         "--min-helix", str(min_helix)],
        input=fasta, capture_output=True, text=True, check=True)
    lines = result.stdout.split("\n")
    out = {}
    for at in range(0, len(lines) - 2, 3):
        if lines[at].startswith(">"):
            out[lines[at][1:]] = (lines[at + 1], lines[at + 2])
    return out


def timed(cmd):
    """Run cmd under /usr/bin/time; return (stdout, seconds, peak_kb)."""
    with tempfile.NamedTemporaryFile("r", suffix=".time") as stamp:
        result = subprocess.run(["/usr/bin/time", "-f", "%e %M", "-o", stamp.name] + cmd,
                                capture_output=True, text=True, cwd=ROOT)
        if result.returncode != 0:
            raise RuntimeError("%s failed: %s" % (cmd[0], result.stderr.strip()))
        seconds, kb = stamp.read().split()[-2:]
    return result.stdout, float(seconds), int(kb)


def parse_access(text):
    out = {}
    for line in text.splitlines():
        if not line[:1].isdigit():
            continue
        pos, _, rest = line.partition("\t")
        for item in rest.split(";"):
            if "," in item:
                w, v = item.split(",")
                out[(int(pos), int(w))] = float(v)
    return out


def parse_capr(text):
    out = {}
    for line in text.splitlines():
        fields = line.split()
        if fields and fields[0] in CONTEXTS:
            out[fields[0]] = [float(v) for v in fields[1:]]
    return out


def raccess(seq, g, lengths, beam=None, extra=()):
    cmd = [os.path.join(BUILD, "linearknot"), "access", "--sequence", seq, "--scaffold", g,
           "--access-len", ",".join(map(str, lengths)), "--probability"]
    cmd += ["--exact"] if beam is None else ["-b", str(beam)]
    return timed(cmd + list(extra))


def capr(seq, g, beam=None):
    cmd = [os.path.join(BUILD, "linearknot"), "profile", "--sequence", seq, "--scaffold", g]
    cmd += ["--exact"] if beam is None else ["-b", str(beam)]
    return timed(cmd)


def run_scaling(args, out):
    records = list(random_sequences(args.lengths, args.seeds))
    scaf = scaffolds(records, args.threshold, args.min_helix)
    out.write("name\tn\tg_pairs\tbeam\tseconds\tpeak_kb\n")
    for name, _ in records:
        seq, g = scaf[name]
        for beam in args.beams:
            try:
                _, seconds, kb = raccess(seq, g, args.access_len, beam)
                out.write("%s\t%d\t%d\t%d\t%.3f\t%d\n" % (name, len(seq), g.count("("), beam, seconds, kb))
            except RuntimeError as error:
                # e.g. "beam retained no structure": record and continue.
                out.write("%s\t%d\t%d\t%d\tNA\tNA\n" % (name, len(seq), g.count("("), beam))
                print(error, file=sys.stderr)
            out.flush()


def run_accuracy(args, out):
    records = list(random_sequences(args.lengths, args.seeds))
    scaf = scaffolds(records, args.threshold, args.min_helix)
    out.write("name\tn\tbeam\tquantity\tmax_abs_err\tmean_abs_err\n")
    for name, _ in records:
        seq, g = scaf[name]
        exact_access = parse_access(raccess(seq, g, args.access_len)[0])
        exact_capr = parse_capr(capr(seq, g)[0])
        for beam in args.beams:
            got = parse_access(raccess(seq, g, args.access_len, beam)[0])
            for w in args.access_len:
                errs = [abs(got[k] - exact_access[k]) for k in exact_access if k[1] == w]
                out.write("%s\t%d\t%d\taccess_w%d\t%.3e\t%.3e\n" % (
                    name, len(seq), beam, w, max(errs), sum(errs) / len(errs)))
            got_capr = parse_capr(capr(seq, g, beam)[0])
            for c in CONTEXTS:
                errs = [abs(a - b) for a, b in zip(got_capr[c], exact_capr[c])]
                out.write("%s\t%d\t%d\tcontext_%s\t%.3e\t%.3e\n" % (
                    name, len(seq), beam, c, max(errs), sum(errs) / len(errs)))
            out.flush()


def run_pkfree(args, out):
    if not args.linracc or not os.path.exists(args.linracc):
        sys.exit("--linracc PATH to the LinearRaccess LinRacc binary is required")
    out.write("name\tn\tw\tmax_abs_diff\tmean_abs_diff\n")
    for name, seq in random_sequences(args.lengths, args.seeds):
        g = "." * len(seq)
        ours = parse_access(raccess(seq, g, args.access_len)[0])
        with tempfile.TemporaryDirectory() as tmp:
            fa = os.path.join(tmp, "in.fa")
            res = os.path.join(tmp, "out.txt")
            with open(fa, "w") as handle:
                handle.write(">%s\n%s\n" % (name, seq))
            subprocess.run([args.linracc, "-seqfile=" + fa, "-outfile=" + res,
                            "-access_len=" + ",".join(map(str, args.access_len)),
                            "-beam=0", "-probabilities", "-energy=turner2004"],
                           capture_output=True, check=True)
            theirs = parse_access(open(res).read())
        for w in args.access_len:
            keys = [k for k in ours if k[1] == w and k in theirs]
            diffs = [abs(ours[k] - theirs[k]) for k in keys]
            out.write("%s\t%d\t%d\t%.3e\t%.3e\n" % (name, len(seq), w, max(diffs),
                                                    sum(diffs) / len(diffs)))
        out.flush()


def ints(text):
    return [int(v) for v in re.split(r"[ ,]+", text.strip()) if v]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("mode", choices=["scaling", "accuracy", "pkfree"])
    parser.add_argument("--out", required=True)
    parser.add_argument("--lengths", type=ints, default=[100, 200, 400])
    parser.add_argument("--seeds", type=ints, default=[1, 2])
    parser.add_argument("--beams", type=ints, default=[20, 50, 100])
    parser.add_argument("--access-len", type=ints, default=[1, 4, 8])
    parser.add_argument("--threshold", type=float, default=0.5)
    parser.add_argument("--min-helix", type=int, default=2)
    parser.add_argument("--linracc", default="")
    args = parser.parse_args()
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as out:
        {"scaling": run_scaling, "accuracy": run_accuracy, "pkfree": run_pkfree}[args.mode](args, out)


if __name__ == "__main__":
    main()
