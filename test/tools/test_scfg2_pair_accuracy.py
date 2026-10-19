#!/usr/bin/env python3
"""Unit tests for the SCFG2 pair-set evaluator."""

from __future__ import annotations

import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent))

from scfg2_pair_accuracy import (  # noqa: E402
    evaluate_pair_sets,
    parse_dot_bracket,
)


def main() -> int:
    assert parse_dot_bracket("(....)") == frozenset({(1, 6)})
    assert parse_dot_bracket("[....]") == frozenset({(1, 6)})
    assert parse_dot_bracket("{....}") == frozenset({(1, 6)})
    assert parse_dot_bracket("<....>") == frozenset({(1, 6)})
    assert parse_dot_bracket("x_X...") == frozenset()

    reference = parse_dot_bracket("(.[..).]")
    assert reference == frozenset({(1, 6), (3, 8)})
    perfect = evaluate_pair_sets(reference, reference)
    assert perfect.all_pairs.f1 == 1.0
    assert perfect.crossing.f1 == 1.0
    assert perfect.noncrossing.f1 == 1.0
    assert perfect.exact

    one_pair = evaluate_pair_sets({(1, 6)}, reference)
    assert one_pair.all_pairs.true_positive == 1
    assert one_pair.all_pairs.false_positive == 0
    assert one_pair.all_pairs.false_negative == 1
    assert one_pair.all_pairs.ppv == 1.0
    assert one_pair.all_pairs.sensitivity == 0.5
    assert one_pair.all_pairs.f1 == 2.0 / 3.0
    # The coordinate is a reference crossing pair, but it is not crossing in
    # the incomplete prediction because its partner stem was missed.
    assert one_pair.crossing.true_positive == 0
    assert one_pair.reference_crossing_coordinate_recovered == 1

    empty = evaluate_pair_sets(set(), set())
    assert empty.all_pairs.ppv == 1.0
    assert empty.all_pairs.sensitivity == 1.0
    assert empty.all_pairs.f1 == 1.0
    assert empty.exact

    no_prediction = evaluate_pair_sets(set(), {(1, 6)})
    assert no_prediction.all_pairs.ppv == 0.0
    assert no_prediction.all_pairs.sensitivity == 0.0
    assert no_prediction.all_pairs.f1 == 0.0

    false_prediction = evaluate_pair_sets({(1, 6)}, set())
    assert false_prediction.all_pairs.ppv == 0.0
    assert false_prediction.all_pairs.sensitivity == 1.0
    assert false_prediction.all_pairs.f1 == 0.0

    for invalid in ("(", ")", "A"):
        try:
            parse_dot_bracket(invalid)
        except ValueError:
            pass
        else:
            raise AssertionError(f"accepted invalid structure: {invalid}")

    print("scfg2_pair_accuracy=ok")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
