#!/usr/bin/env python3
import math
import statistics
import sys


def to_float(value):
    try:
        return float(value)
    except Exception:
        return None


def cmd_abs_diff(a, b):
    fa = to_float(a)
    fb = to_float(b)
    if fa is None or fb is None:
        raise SystemExit(1)
    print(abs(fa - fb))


def cmd_scale(a, b):
    fa = to_float(a)
    fb = to_float(b)
    if fa is None or fb is None:
        raise SystemExit(1)
    print(max(abs(fa), abs(fb)))


def cmd_rel_diff(abs_diff, scale):
    fa = to_float(abs_diff)
    fs = to_float(scale)
    if fa is None or fs is None:
        raise SystemExit(1)
    if fs == 0.0:
        print(0.0)
        return
    print(fa / fs)


def load_values(path):
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


def cmd_stats(path):
    vals = load_values(path)
    if not vals:
        print("0 0 0")
        return
    print(f"{statistics.median(vals)} {pct(vals, 90)} {pct(vals, 95)}")


def cmd_top3(path):
    vals = load_values(path)
    top3 = sorted(vals, reverse=True)[:3]
    while len(top3) < 3:
        top3.append(0.0)
    print(f"{top3[0]},{top3[1]},{top3[2]}")


def cmd_finite_rate(finite_count, compared):
    fc = int(finite_count)
    comp = int(compared)
    if comp == 0:
        print("0")
        return
    print((fc * 100.0) / comp)


def cmd_check_finite(a, b):
    fa = to_float(a)
    fb = to_float(b)
    if fa is None or fb is None:
        raise SystemExit(1)
    if math.isfinite(fa) and math.isfinite(fb):
        return
    raise SystemExit(1)


def cmd_check_nonzero(value):
    v = to_float(value)
    if v is None:
        raise SystemExit(1)
    if abs(v) > 1e-12:
        return
    raise SystemExit(1)


def cmd_check_within_tol(abs_diff, scale, abs_tol, rel_tol):
    d = to_float(abs_diff)
    s = to_float(scale)
    at = to_float(abs_tol)
    rt = to_float(rel_tol)
    if None in (d, s, at, rt):
        raise SystemExit(1)
    ok = (d <= at) or (s > 0.0 and d <= rt * s)
    raise SystemExit(0 if ok else 1)


def cmd_check_rate_ge(rate_str, limit):
    rate = to_float(str(rate_str).rstrip("%"))
    lim = to_float(limit)
    if rate is None or lim is None:
        raise SystemExit(1)
    raise SystemExit(0 if rate >= lim else 1)


def cmd_check_le(value, limit):
    val = to_float(value)
    lim = to_float(limit)
    if val is None or lim is None:
        raise SystemExit(1)
    raise SystemExit(0 if val <= lim else 1)


def cmd_check_thresholds(*values):
    if len(values) != 12:
        raise SystemExit(1)
    keys = ["abs_median", "abs_p90", "abs_p95", "rel_median", "rel_p90", "rel_p95"]
    actual_vals = values[:6]
    limit_vals = values[6:]
    actual = dict(zip(keys, actual_vals))
    limits = dict(zip(keys, limit_vals))
    violations = []
    for key in keys:
        limit_raw = limits.get(key, "")
        if limit_raw in ("", None):
            continue
        limit = to_float(limit_raw)
        if limit is None or not math.isfinite(limit):
            continue
        val = to_float(actual.get(key, "nan"))
        if val is None or not math.isfinite(val):
            continue
        if val > limit:
            violations.append(f"{key}={val}>{limit}")
    if violations:
        print(";".join(violations))
        raise SystemExit(1)


def command_table():
    return {
        "abs_diff": (cmd_abs_diff, 2),
        "scale": (cmd_scale, 2),
        "rel_diff": (cmd_rel_diff, 2),
        "stats": (cmd_stats, 1),
        "top3": (cmd_top3, 1),
        "finite_rate": (cmd_finite_rate, 2),
        "check_finite": (cmd_check_finite, 2),
        "check_nonzero": (cmd_check_nonzero, 1),
        "check_within_tol": (cmd_check_within_tol, 4),
        "check_rate_ge": (cmd_check_rate_ge, 2),
        "check_le": (cmd_check_le, 2),
        "check_thresholds": (cmd_check_thresholds, 12),
    }


def main(argv):
    if len(argv) < 1:
        raise SystemExit("usage: mfe_alignment_math.py <command> [args...]")
    cmd = argv[0]
    args = argv[1:]
    table = command_table()
    entry = table.get(cmd)
    if not entry or len(args) != entry[1]:
        raise SystemExit(
            "usage: mfe_alignment_math.py <abs_diff|scale|rel_diff|stats|top3|finite_rate|check_finite|check_nonzero|check_within_tol|check_rate_ge|check_le|check_thresholds> ..."
        )
    entry[0](*args)


if __name__ == "__main__":
    main(sys.argv[1:])
