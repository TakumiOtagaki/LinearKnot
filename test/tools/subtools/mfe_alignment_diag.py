#!/usr/bin/env python3
import collections
import csv
import sys


def normalize_reason_prefix(reason):
    if reason in ("", "none", None):
        return "none"
    for marker in ("_i_", "_state_", "_rule_"):
        idx = reason.find(marker)
        if idx > 0:
            return reason[:idx]
    return reason


def read_rows(path):
    with open(path, encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def row_value(row, key, legacy_key=None):
    value = row.get(key)
    if value not in (None, ""):
        return value
    if legacy_key is not None:
        return row.get(legacy_key)
    return value


def print_mismatch_summary(rows):
    mismatches = [r for r in rows if r.get("mismatch") == "1"]
    print(f"diag_total={len(rows)}")
    print(f"diag_mismatched={len(mismatches)}")
    by_family = collections.Counter(r.get("topology_family") or "none" for r in mismatches)
    for fam, cnt in sorted(by_family.items()):
        print(f"diag_mismatch_family_{fam}={cnt}")
    by_rule = collections.Counter(
        (row_value(r, "decision_path_rule", "dp_rule") or "none") for r in mismatches
    )
    for rule, cnt in by_rule.most_common(5):
        print(f"diag_mismatch_rule_{rule}={cnt}")


def print_real_reason_summaries(rows):
    by_debug_reason = collections.Counter(
        (row_value(r, "debug_score_reason", "real_score_reason") or "none") for r in rows
    )
    for reason, cnt in by_debug_reason.most_common(5):
        print(f"diag_debug_reason_{reason}={cnt}")
    by_density2_dp_reason = collections.Counter(
        (
            row_value(
                r,
                "debug_score_comparison_path_reason",
                "debug_score_density2_dpdecisions_reason",
            )
            or row_value(
                r,
                "real_score_comparison_path_reason",
                "real_score_density2_dpdecisions_reason",
            )
            or "none"
        )
        for r in rows
    )
    for reason, cnt in by_density2_dp_reason.most_common(5):
        print(f"diag_debug_density2_dp_reason_{reason}={cnt}")
    by_density2_dp_reason_prefix = collections.Counter(
        normalize_reason_prefix(
            row_value(
                r,
                "debug_score_comparison_path_reason",
                "debug_score_density2_dpdecisions_reason",
            )
            or row_value(
                r,
                "real_score_comparison_path_reason",
                "real_score_density2_dpdecisions_reason",
            )
            or "none"
        )
        for r in rows
    )
    for reason, cnt in by_density2_dp_reason_prefix.most_common(5):
        print(f"diag_debug_density2_dp_reason_prefix_{reason}={cnt}")
    lane_a_pkfree = sum(
        1
        for r in rows
        if (r.get("topology_family") == "pk_free")
        and (
            row_value(r, "debug_score_reason", "real_score_reason") or ""
        ).startswith("real_no_positive_weight_state_V_rule_V_INTERNAL")
    )
    print(f"diag_lane_a_pkfree_no_positive_v_internal={lane_a_pkfree}")
    by_debug_mode = collections.Counter(
        (row_value(r, "debug_score_mode", "real_score_mode") or "none") for r in rows
    )
    for mode, cnt in by_debug_mode.most_common(5):
        print(f"diag_debug_mode_{mode}={cnt}")


def print_decision_path_score_summary(rows):
    decision_path_scores = []
    for r in rows:
        s = (row_value(r, "decision_path_score", "dp_score") or "").strip()
        if not s:
            continue
        try:
            v = float(s)
        except ValueError:
            continue
        if v == v:
            decision_path_scores.append(v)
    if not decision_path_scores:
        return
    decision_path_scores_sorted = sorted(decision_path_scores)
    m = len(decision_path_scores_sorted) // 2
    if len(decision_path_scores_sorted) % 2 == 1:
        decision_path_med = decision_path_scores_sorted[m]
    else:
        decision_path_med = 0.5 * (
            decision_path_scores_sorted[m - 1] + decision_path_scores_sorted[m]
        )
    near_zero = sum(1 for v in decision_path_scores if abs(v) <= 1e-9)
    print(f"diag_decision_path_score_total={len(decision_path_scores)}")
    print(f"diag_decision_path_score_median={decision_path_med}")
    print(f"diag_decision_path_score_near_zero={near_zero}")


def print_shadow_summary(rows):
    shadow_rows = [
        r
        for r in rows
        if (row_value(r, "debug_score_energy", "real_score_energy") or "") not in ("", "nan", "NaN")
    ]
    if not shadow_rows:
        return
    api_abs = []
    real_abs = []
    real_better = 0
    by_family = collections.defaultdict(lambda: {"total": 0, "real_better": 0})
    for r in shadow_rows:
        try:
            legacy = float(r.get("legacy_mfe_energy") or "nan")
            api = float(row_value(r, "union_energy", "current_union_energy", "scfg_rescoring_energy") or "nan")
            debug = float(row_value(r, "debug_score_energy", "real_score_energy") or "nan")
        except ValueError:
            continue
        if not (legacy == legacy and api == api and debug == debug):
            continue
        api_d = abs(api - legacy)
        real_d = abs(debug - legacy)
        api_abs.append(api_d)
        real_abs.append(real_d)
        fam = r.get("topology_family") or "none"
        by_family[fam]["total"] += 1
        if real_d < api_d:
            real_better += 1
            by_family[fam]["real_better"] += 1
    if not (api_abs and real_abs):
        return
    api_abs.sort()
    real_abs.sort()
    mid = len(api_abs) // 2
    if len(api_abs) % 2 == 1:
        api_med = api_abs[mid]
        real_med = real_abs[mid]
    else:
        api_med = 0.5 * (api_abs[mid - 1] + api_abs[mid])
        real_med = 0.5 * (real_abs[mid - 1] + real_abs[mid])
    print(f"diag_debug_shadow_total={len(api_abs)}")
    print(f"diag_debug_shadow_api_abs_median={api_med}")
    print(f"diag_debug_shadow_debug_abs_median={real_med}")
    print(f"diag_debug_shadow_debug_better_than_api={real_better}")
    for fam in sorted(by_family.keys()):
        total = by_family[fam]["total"]
        better = by_family[fam]["real_better"]
        print(f"diag_debug_shadow_family_{fam}_total={total}")
        print(f"diag_debug_shadow_family_{fam}_debug_better={better}")


def main(argv):
    if len(argv) != 1:
        raise SystemExit("usage: mfe_alignment_diag.py <diag_csv>")
    rows = read_rows(argv[0])
    print_mismatch_summary(rows)
    print_real_reason_summaries(rows)
    print_decision_path_score_summary(rows)
    print_shadow_summary(rows)


if __name__ == "__main__":
    main(sys.argv[1:])
