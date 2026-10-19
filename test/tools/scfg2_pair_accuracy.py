#!/usr/bin/env python3
"""Pair-set accuracy metrics for pseudoknotted dot-bracket structures."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Iterable, Sequence


Pair = tuple[int, int]
OPEN_TO_CLOSE = {"(": ")", "[": "]", "{": "}", "<": ">"}
CLOSE_TO_OPEN = {close: open_ for open_, close in OPEN_TO_CLOSE.items()}
UNPAIRED = {".", "_", "x", "X"}


def parse_dot_bracket(structure: str) -> frozenset[Pair]:
    stacks: dict[str, list[int]] = {symbol: [] for symbol in OPEN_TO_CLOSE}
    pairs: set[Pair] = set()
    partners: set[int] = set()
    for position, symbol in enumerate(structure, start=1):
        if symbol in OPEN_TO_CLOSE:
            stacks[symbol].append(position)
            continue
        if symbol in CLOSE_TO_OPEN:
            open_symbol = CLOSE_TO_OPEN[symbol]
            if not stacks[open_symbol]:
                raise ValueError(
                    f"unbalanced close {symbol!r} at position {position}")
            left = stacks[open_symbol].pop()
            if left in partners or position in partners:
                raise ValueError("one endpoint occurs in multiple pairs")
            partners.update((left, position))
            pairs.add((left, position))
            continue
        if symbol not in UNPAIRED:
            raise ValueError(
                f"unsupported dot-bracket symbol {symbol!r}"
                f" at position {position}")
    unfinished = [symbol for symbol, stack in stacks.items() if stack]
    if unfinished:
        raise ValueError(
            f"unbalanced open symbols: {','.join(sorted(unfinished))}")
    return frozenset(pairs)


def crossing_pairs(pairs: Iterable[Pair]) -> frozenset[Pair]:
    ordered = tuple(sorted(pairs))
    crossing: set[Pair] = set()
    for index, (i, j) in enumerate(ordered):
        for k, l in ordered[index + 1:]:
            if (i < k < j < l) or (k < i < l < j):
                crossing.update(((i, j), (k, l)))
    return frozenset(crossing)


@dataclass(frozen=True)
class PairScore:
    true_positive: int
    false_positive: int
    false_negative: int

    @property
    def predicted(self) -> int:
        return self.true_positive + self.false_positive

    @property
    def reference(self) -> int:
        return self.true_positive + self.false_negative

    @property
    def ppv(self) -> float:
        if self.predicted == 0:
            return 1.0 if self.reference == 0 else 0.0
        return self.true_positive / self.predicted

    @property
    def sensitivity(self) -> float:
        if self.reference == 0:
            return 1.0
        return self.true_positive / self.reference

    @property
    def f1(self) -> float:
        denominator = (
            2 * self.true_positive + self.false_positive + self.false_negative)
        return 1.0 if denominator == 0 else 2 * self.true_positive / denominator


def score_sets(predicted: Iterable[Pair],
               reference: Iterable[Pair]) -> PairScore:
    predicted_set = frozenset(predicted)
    reference_set = frozenset(reference)
    return PairScore(
        true_positive=len(predicted_set & reference_set),
        false_positive=len(predicted_set - reference_set),
        false_negative=len(reference_set - predicted_set),
    )


@dataclass(frozen=True)
class StructureAccuracy:
    all_pairs: PairScore
    crossing: PairScore
    noncrossing: PairScore
    reference_crossing_coordinate_recovered: int
    reference_crossing_count: int
    exact: bool


def evaluate_pair_sets(predicted: Iterable[Pair],
                       reference: Iterable[Pair]) -> StructureAccuracy:
    predicted_set = frozenset(predicted)
    reference_set = frozenset(reference)
    predicted_crossing = crossing_pairs(predicted_set)
    reference_crossing = crossing_pairs(reference_set)
    return StructureAccuracy(
        all_pairs=score_sets(predicted_set, reference_set),
        crossing=score_sets(predicted_crossing, reference_crossing),
        noncrossing=score_sets(
            predicted_set - predicted_crossing,
            reference_set - reference_crossing,
        ),
        reference_crossing_coordinate_recovered=len(
            predicted_set & reference_crossing),
        reference_crossing_count=len(reference_crossing),
        exact=predicted_set == reference_set,
    )


def score_fields(prefix: str, score: PairScore) -> dict[str, object]:
    return {
        f"{prefix}_tp": score.true_positive,
        f"{prefix}_fp": score.false_positive,
        f"{prefix}_fn": score.false_negative,
        f"{prefix}_ppv": score.ppv,
        f"{prefix}_sensitivity": score.sensitivity,
        f"{prefix}_f1": score.f1,
    }


def accuracy_fields(accuracy: StructureAccuracy) -> dict[str, object]:
    return {
        **score_fields("all", accuracy.all_pairs),
        **score_fields("crossing", accuracy.crossing),
        **score_fields("noncrossing", accuracy.noncrossing),
        "reference_crossing_coordinate_recovered":
            accuracy.reference_crossing_coordinate_recovered,
        "reference_crossing_count": accuracy.reference_crossing_count,
        "exact": int(accuracy.exact),
    }


def main(argv: Sequence[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction", required=True)
    parser.add_argument("--reference", required=True)
    args = parser.parse_args(argv)
    if len(args.prediction) != len(args.reference):
        parser.error("prediction and reference lengths differ")
    accuracy = evaluate_pair_sets(
        parse_dot_bracket(args.prediction),
        parse_dot_bracket(args.reference),
    )
    fields = accuracy_fields(accuracy)
    print("\t".join(fields))
    print("\t".join(str(value) for value in fields.values()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(__import__("sys").argv[1:]))
