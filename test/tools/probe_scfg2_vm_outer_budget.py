#!/usr/bin/env python3
"""Adversarial index-budget experiment, not an RNA semantic/runtime oracle.

Each right column has <= b WM records, but a fixed left endpoint has O(n)
positive records. Check a monotone prefix top-k index against a full reference
and count index work without relying on sparsity of positive WM values.
The exhaustive sorted reference is validation-only and excluded from work.
"""

import argparse
import csv
import hashlib
import json
import shutil
from dataclasses import dataclass
from pathlib import Path


@dataclass
class Work:
    build_comparisons: int = 0
    build_moves: int = 0
    rollback_moves: int = 0
    cursor_steps: int = 0
    membership_marks: int = 0
    probes: int = 0

    def total(self):
        return sum(vars(self).values())


class PrefixTop:
    """Immutable entries in increasing end order; query upper bounds decrease.

    Building the top-k array costs O(k) per input. One constant-size undo
    record per input suffices to roll back a prefix in O(k) per removed input.
    All queries together cost O(C*k + M*k), space O(C+k).
    """

    def __init__(self, entries, width, work):
        self.entries = entries  # (end, score); no sorting inside the index
        self.width = width
        self.work = work
        self.top = []  # entry indices, descending score / ascending end
        self.undo = []
        self.cursor = len(entries)
        self.last_upper = float("inf")
        if any(a[0] >= b[0] for a, b in zip(entries, entries[1:])):
            raise ValueError("entries must arrive in increasing endpoint order")
        for idx in range(len(entries)):
            pos = 0
            key = self.key(idx)
            while pos < len(self.top):
                work.build_comparisons += 1
                if key < self.key(self.top[pos]):
                    break
                pos += 1
            if pos >= width:
                self.undo.append((-1, None))
                continue
            work.build_moves += len(self.top) - pos + 1
            self.top.insert(pos, idx)
            evicted = self.top.pop() if len(self.top) > width else None
            self.undo.append((pos, evicted))

    def key(self, idx):
        end, score = self.entries[idx]
        return -score, end

    def query(self, upper):
        if upper > self.last_upper:
            raise ValueError("prefix upper bounds must not increase")
        self.last_upper = upper
        while self.cursor and self.entries[self.cursor - 1][0] > upper:
            self.cursor -= 1
            self.work.cursor_steps += 1
            pos, evicted = self.undo[self.cursor]
            if pos < 0:
                continue
            if self.top[pos] != self.cursor:
                raise AssertionError("rollback does not undo the last insertion")
            self.work.rollback_moves += len(self.top) - pos
            self.top.pop(pos)
            if evicted is not None:
                self.top.append(evicted)
                self.work.rollback_moves += 1
        return self.top


def run(n, b):
    work = Work()
    # Populate by increasing right column. No per-left O(log n) sorting.
    left = [[] for _ in range(n + 2)]
    max_right_wm = 0
    for end in range(2, n):
        starts = range(2, min(b + 2, end + 1))
        max_right_wm = max(max_right_wm, len(starts))
        for start in starts:
            # Deliberately dense and nonmonotone weights, with ties.
            left[start].append((end, (31 * start + 17 * end) % 101))
    # Two right families (WMv/WMp), each with <= b records per endpoint.
    right = [list(range(2, min(b + 2, end - 3))) for end in range(n + 1)]
    max_right = max(map(len, right))
    width = b + max_right
    indices = [PrefixTop(xs, width, work) for xs in left]
    # Coordinate stamps give deterministic O(1) membership without clearing
    # an n-entry array for each parent or assuming hash-table behavior.
    right_seen = [0] * (n + 2)
    entries = sum(map(len, left))
    queries = checks = old_enumerations = max_left_degree = 0
    new_right_edges = preserved_right_edges = left_edges = 0
    # Shared family/span traversal makes j decrease for each fixed i/start.
    for span in range(n, 6, -1):
        for i in range(1, min(b + 1, n - span + 2)):
            j = i + span - 1
            prefix = indices[i + 1]
            candidates = prefix.query(j - 5)  # k <= j-turn-1, end=k-1
            queries += 1
            max_left_degree = max(max_left_degree, prefix.cursor)
            old_enumerations += 2 * prefix.cursor
            for k in right[j - 1]:
                right_seen[k] = queries
                work.membership_marks += 1
            emitted = []
            for family in range(2):
                novel = []
                for k in right[j - 1]:
                    work.probes += 1
                    if not i + 2 <= k <= j - 4:
                        continue
                    left_edges += 1  # known right supplies a left-hole edge
                    # Preserve every positive-positive right-target edge.
                    if i + 1 <= k - 1 <= n - 1:
                        preserved_right_edges += 1
                for idx in candidates:
                    work.probes += 1
                    end, _ = prefix.entries[idx]
                    k = end + 1
                    if right_seen[k] == queries:
                        continue
                    novel.append(idx)
                    if len(novel) == b:
                        break
                new_right_edges += len(novel)
                emitted.append(novel)
            # Independent full-history reference: this is NOT algorithm work.
            expected = sorted(
                (idx for idx in range(prefix.cursor)
                 if right_seen[prefix.entries[idx][0] + 1] != queries),
                key=prefix.key)[:b]
            if emitted != [expected, expected]:
                raise AssertionError("bounded novel edges differ from prefix top-b")
            if prefix.top != sorted(range(prefix.cursor), key=prefix.key)[:width]:
                raise AssertionError("prefix ranking differs from exhaustive reference")
            checks += 1
    if max_right_wm > b or max_right > b:
        raise AssertionError("fixture violates a right-column beam")
    if entries > n * b or queries > n * b:
        raise AssertionError("fixture violates the chart/parent budget")
    # Structural bounds apply for arbitrary scores/density, not their average.
    if work.build_comparisons > entries * width:
        raise AssertionError("build comparison budget exceeded")
    if work.build_moves > entries * (width + 1):
        raise AssertionError("build move budget exceeded")
    if work.rollback_moves > entries * (width + 1) or work.cursor_steps > entries:
        raise AssertionError("rollback budget exceeded")
    if work.probes > 2 * queries * (max_right + width):
        raise AssertionError("split-probe budget exceeded")
    if work.membership_marks > queries * max_right:
        raise AssertionError("membership budget exceeded")
    bound = entries * (3 * width + 3) + queries * (3 * max_right + 2 * width)
    if work.total() > bound:
        raise AssertionError("total structural budget exceeded")
    return dict(n=n, b=b, wm_entries=entries, vm_parents=queries,
                max_right_degree=max_right_wm, max_left_degree=max_left_degree,
                uncapped_left_enumerations=old_enumerations,
                capped_probes=work.probes, index_work=work.total(),
                structural_work_bound=bound, index_work_over_nb2=work.total() / (n * b * b),
                preserved_right_edges=preserved_right_edges, left_hole_edges=left_edges,
                capped_novel_right_edges=new_right_edges, reference_checks=checks,
                **vars(work))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sizes", default="64,128,256,512")
    parser.add_argument("--widths", default="4,20")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sizes = list(map(int, args.sizes.split(',')))
    widths = list(map(int, args.widths.split(',')))
    if min(sizes) < 8 or min(widths) < 1:
        parser.error("sizes >= 8 and widths >= 1 required")
    results = []
    for b in widths:
        for n in sizes:
            row = run(n, b)
            results.append(row)
            print(f"n={n} b={b}: uncapped={row['uncapped_left_enumerations']} "
                  f"capped={row['capped_probes']} index_work={row['index_work']} "
                  f"work/(nb²)={row['index_work_over_nb2']:.3f}", flush=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "counts.tsv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(results)
    (args.output_dir / "scope.json").write_text(json.dumps({
        "scope": "adversarial endpoint-index experiment; no RNA semantic or whole-outer claim",
        "sizes": sizes, "widths": widths,
        "validation": "all query prefixes compared with a full sorted reference",
        "work": "index build, rollback, cursor movement, coordinate marks and bounded split probes; validation excluded",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }, indent=2) + '\n')
    shutil.copyfile(__file__, args.output_dir / Path(__file__).name)


if __name__ == "__main__":
    main()
