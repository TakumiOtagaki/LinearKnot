#!/usr/bin/env python3
"""Generate reproducible random sequence/structure pairs for density-2 tests.

Structure grammar:
  S -> S '.' | S '(' S ')' | epsilon
"""

from __future__ import annotations

import argparse
import random
import sys
from typing import List, Optional, Sequence, Tuple


ALPHABET: Tuple[str, ...] = ("A", "U", "G", "C")
ALLOWED_PAIRS: Tuple[Tuple[str, str], ...] = (
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),
)


class GenerationError(RuntimeError):
    pass


def generate_structure(rng: random.Random, length: int, turn: int) -> str:
    if length < 0:
        raise ValueError("length must be >= 0")
    if length == 0:
        return ""

    def build(n: int) -> str:
        if n == 0:
            return ""

        options: List[Tuple[str, Optional[int]]] = []
        # S -> S '.'  (left-recursive form reinterpreted with bounded size)
        if n >= 1:
            options.append(("dot", None))
        # S -> S '(' S ')' : choose inner span size k, remaining prefix uses n-k-2.
        # Enforce TURN: inner length must be >= turn to satisfy j - i - 1 >= turn.
        if n >= 2:
            min_inner = max(turn, 0)
            for inner in range(min_inner, n - 1):
                options.append(("pair", inner))

        if not options:
            return ""

        kind, inner = options[rng.randrange(len(options))]
        if kind == "dot":
            return build(n - 1) + "."

        assert inner is not None
        left_len = n - inner - 2
        if left_len < 0:
            raise GenerationError("invalid grammar split")
        return build(left_len) + "(" + build(inner) + ")"

    return build(length)


def assign_sequence(rng: random.Random, structure: str) -> str:
    seq = ["?"] * len(structure)
    stack: List[int] = []

    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise GenerationError("unbalanced structure while assigning sequence")
            j = stack.pop()
            b1, b2 = ALLOWED_PAIRS[rng.randrange(len(ALLOWED_PAIRS))]
            seq[j] = b1
            seq[i] = b2
        elif ch == ".":
            seq[i] = ALPHABET[rng.randrange(len(ALPHABET))]
        else:
            raise GenerationError(f"invalid structure character: {ch!r}")

    if stack:
        raise GenerationError("unbalanced structure after assignment")

    if any(ch == "?" for ch in seq):
        raise GenerationError("internal assignment error")

    return "".join(seq)


def validate_pairing(seq: str, structure: str) -> None:
    if len(seq) != len(structure):
        raise GenerationError("sequence and structure length mismatch")

    stack: List[int] = []
    allowed = set(ALLOWED_PAIRS)

    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            if not stack:
                raise GenerationError("unbalanced structure")
            j = stack.pop()
            pair = (seq[j], seq[i])
            if pair not in allowed:
                raise GenerationError(f"disallowed pair at ({j},{i}): {pair}")
        elif ch != ".":
            raise GenerationError(f"invalid structure character: {ch!r}")

    if stack:
        raise GenerationError("unbalanced structure")

    for ch in seq:
        if ch not in ALPHABET:
            raise GenerationError(f"invalid sequence character: {ch!r}")


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate random (seq, structure) pairs")
    parser.add_argument("--seed", type=int, default=1, help="PRNG seed")
    parser.add_argument("--count", type=int, default=10, help="number of cases")
    parser.add_argument("--length", type=int, default=24, help="sequence/structure length")
    parser.add_argument("--turn", type=int, default=3, help="minimum unpaired bases between paired positions")
    parser.add_argument(
        "--seq-only",
        action="store_true",
        help="emit random sequences only (structure field left empty)",
    )
    parser.add_argument(
        "--prefix",
        default="rand",
        help="case id prefix (case ids are <prefix>_0001, ...)",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)

    if args.count <= 0:
        print("error: --count must be > 0", file=sys.stderr)
        return 2
    if args.length < 0:
        print("error: --length must be >= 0", file=sys.stderr)
        return 2

    rng = random.Random(args.seed)

    print("case_id\tseq\tG")
    for i in range(1, args.count + 1):
        if args.seq_only:
            seq = "".join(rng.choice(ALPHABET) for _ in range(args.length))
            print(f"{args.prefix}_{i:04d}\t{seq}\t")
            continue
        structure = generate_structure(rng, args.length, args.turn)
        seq = assign_sequence(rng, structure)
        validate_pairing(seq, structure)
        print(f"{args.prefix}_{i:04d}\t{seq}\t{structure}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
