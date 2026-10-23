"""Independent test reference for the PseudoLoop definition introduced in 0.2.0.

Uses only the completed union pair map: no rule IDs, scaffold provenance,
energies or probabilities. Run this file to check the hand-specified fixtures.
See documents/pseudoloop_context_change.md for scope and the loop convention.
"""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import re


def pairs_from_db(db):
    stacks = {x: [] for x in '([{<'}
    close = dict(zip(')]}>', '([{<'))
    pairs = set()
    for i, x in enumerate(db, 1):
        if x in stacks:
            stacks[x].append(i)
        elif x in close:
            if not stacks[close[x]]:
                raise ValueError('unbalanced structure')
            pairs.add((stacks[close[x]].pop(), i))
        elif x != '.':
            raise ValueError(x)
    if any(stacks.values()):
        raise ValueError('unbalanced structure')
    return pairs


def classify(n, pairs):
    pairs = tuple(sorted(tuple(pair) for pair in pairs))
    mate = [0] * (n + 1)
    pair_at = {}
    for i, j in pairs:
        if not (1 <= i < j <= n) or mate[i] or mate[j]:
            raise ValueError('invalid pair map or repeated endpoint')
        mate[i], mate[j] = j, i
        pair_at[i] = pair_at[j] = (i, j)
    def crosses(a, b):
        i, j = a
        k, l = b
        return i < k < j < l or k < i < l < j
    crossing = {p: frozenset(q for q in pairs if crosses(p, q)) for p in pairs}

    @lru_cache(None)
    def closed_right(lo):
        hi = mate[lo]
        if hi <= lo:
            return None
        pos = lo
        while pos <= hi:
            if mate[pos]:
                if mate[pos] < lo:
                    return None
                hi = max(hi, mate[pos])
            pos += 1
        return hi

    @lru_cache(None)
    def closed_left(hi):
        lo = mate[hi]
        if not 0 < lo < hi:
            return None
        pos = hi
        while pos >= lo:
            if mate[pos]:
                if mate[pos] > hi:
                    return None
                lo = min(lo, mate[pos])
            pos -= 1
        return lo

    def branches(lo, hi):
        out = []
        pos = lo
        while pos <= hi:
            if not mate[pos]:
                pos += 1
                continue
            end = closed_right(pos)
            if end is None or end > hi:
                raise ValueError(('nonclosed child in ordinary-loop side', lo, hi, pos))
            out.append((pos, end))
            pos = end + 1
        return out

    labels = []
    explanations = []
    for x in range(1, n + 1):
        if mate[x]:
            labels.append('S')
            explanations.append(None)
            continue
        left = x - 1
        while left:
            if not mate[left]:
                left -= 1
                continue
            start = closed_left(left)
            if start is None:
                break
            left = start - 1
        right = x + 1
        while right <= n:
            if not mate[right]:
                right += 1
                continue
            end = closed_right(right)
            if end is None:
                break
            right = end + 1
        explanation = {'left_stop': left, 'right_stop': right}
        if left == 0 and right == n + 1:
            context = 'E'
        elif left == 0 or right == n + 1:
            raise ValueError(('one-sided exposed boundary', x, left, right))
        else:
            a, b = pair_at[left], pair_at[right]
            explanation.update(left_pair=a, right_pair=b)
            if a == b:
                child = branches(a[0] + 1, a[1] - 1)
                explanation['children'] = child
                if not child:
                    context = 'H'
                elif len(child) >= 2 or crossing[pair_at[child[0][0]]]:
                    # A child whose exterior exposes crossing bands supplies
                    # multiple stem attachments. An enclosing noncrossing pair
                    # shields any deeper pseudoknot and supplies only one.
                    context = 'M'
                else:
                    k, l = child[0]
                    context = 'I' if k > a[0] + 1 and l < a[1] - 1 else 'B'
            else:
                outer, inner = sorted((a, b))
                same_band = (outer[0] < inner[0] < inner[1] < outer[1]
                             and crossing[outer] and crossing[outer] == crossing[inner])
                if not same_band:
                    context = 'P'
                else:
                    lchild = branches(outer[0] + 1, inner[0] - 1)
                    rchild = branches(inner[1] + 1, outer[1] - 1)
                    explanation.update(left_children=lchild, right_children=rchild)
                    if lchild or rchild:
                        context = 'M'
                    else:
                        context = 'I' if inner[0] > outer[0] + 1 and inner[1] < outer[1] - 1 else 'B'
        labels.append(context)
        explanations.append(explanation)
    return ''.join(labels), explanations


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fixtures', type=Path,
                        default=Path(__file__).resolve().parents[1] /
                        'fixtures/pseudoloop_context_change.json')
    args = parser.parse_args()
    cases = json.loads(args.fixtures.read_text())['cases']
    windows = 0
    for case in cases:
        labels, _ = classify(case['length'], case['pairs'])
        if labels != case['expected']:
            raise AssertionError((case['name'], labels, case['expected']))
        for run in filter(None, re.split('S+', labels)):
            if len(set(run)) != 1:
                raise AssertionError(('mixed fully unpaired run', case['name'], run))
            windows += len(run) * (len(run) + 1) // 2
    print(json.dumps({'geometry_fixtures': len(cases),
                      'fully_unpaired_windows': windows,
                      'status': 'ok; geometry reference'}))


if __name__ == '__main__':
    main()
