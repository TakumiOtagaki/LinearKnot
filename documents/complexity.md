# Computational complexity

Bounds refer to logical live storage, not allocator RSS. Let n be sequence
length, b the uniform finite beam width, r the number of requested window
length entries, and K the number of emitted structures. Repeated window
lengths in the C++ API count as separate output entries.

## Finite-beam production calls

The following bounds assume 1 <= b <= n, fixed internal-loop, hairpin and
branching-ancestor limits, standard scaffold forcing, word-RAM arithmetic,
and expected constant-time hash lookup. They include scaffold preparation,
forcing, pilot and outer preparation, main DP, and the listed consumer.

| Operation | Time | Peak logical space |
| --- | --- | --- |
| Numeric inside and regenerated outside/BPP | O(nb^2 log(b+2)) | O(nb) |
| Structural loop profile | O(nb^2 log(b+2)) | O(nb) |
| Accessibility for r requested lengths | O(nb^2 log(b+2) + rnb^2 + nr) | O(nb + nr) |
| MFE, MEA and gamma-centroid, including structure output | O(nb^2 log(b+2)) | O(nb) |

For fixed r, accessibility has the same asymptotic bound as inside.
For b >= 2, log(b+2) can be written as Theta(log b). Fixed beam width and
fixed local limits make these bounds linear in n. They are not claims for
exact/KeepAll calls, arbitrary forced sets, unbounded geometry limits, or
every member of the toolset.

## Storage and operation budgets

For scaffold G, let m be its pair count, Delta its maximum number of direct
pair children, and E_d the number of endpoint-to-branching-ancestor entries
kept under ancestor limit d. Then m and Delta are O(n), and E_d <= nd.

| Component | Budget |
| --- | --- |
| Scaffold geometry and RMQ preprocessing | O(n) time and space |
| Capped endpoint/ancestor metadata | O(n+E_d) space |
| Per-child summaries and bounded registries | O(n+mb) space |
| Public numeric states and retained private carriers | O(nb) space |
| Current-column join scratch | O(Delta b + b^2) space |
| Candidate generation and bounded selection | O(nb^2 log(b+2)) total time |

Incoming prefix candidates are bounded by b within each activation group;
nested frames each retain at most b candidates. Deduplication is charged to
retained WI states. Public selection retains a bounded prefix instead of
sorting all historical candidates. Parent groups use bounded sorts and
ordered merges. Consumed histories and column scratch are released.
Thus O(nb + E_d + Delta b + b^2) = O(nb) under the stated assumptions.

The numeric chart's optional lookup index adds O(n+L) space for L live
primary keys: O(n) row headers and geometrically sized open-addressed rows.
Rows shrink after erasure and release their slots when empty, so old
candidate peaks do not accumulate into a dense n-by-n allocation. The
canonical sparse map remains the source of presence and iteration order.
Queued V sources store only their rule and split coordinates. Bounded
completion maxima use an indexed heap with O(log(b+2)) updates; querying
immutable sorted alpha rows costs O(log(b+2)). These changes preserve the
budgets above. An optional sorted diagnostic view has a separate
O(b log(b+2)) observation cost per changed pool.

For standard forcing width beta = max(1000,b) and pilot width p = ceil(b^q),
q in {1/2,2/3}, the sequential phase budgets are

```text
T = O(n beta^2 log(beta+2) + n p^2 log(p+2) + n b^2 log(b+2))
S = O(n beta + beta^2 + np + nb + b^2).
```

The forcing floor is constant, beta = O(b), and p <= b. The guide itself
uses O(np^2 log(p+2)) time and O(np) saved scores. Regenerated outside
reconstructs products from retained states. Pair export uses coordinate
counting sorts. Profiles use difference arrays, chain prefix log products
and monotone hinted insertion, avoiding a scan proportional to window
length or an additional log n factor on those production paths.

MFE stores selected choices. Beam posterior decoding replays the
alpha-ranked execution with a maximum additive reward and stores selected
decompositions. Silent unary steps are aliased; losing alternatives are
released. Selected-DAG output and pair assembly stay within the same
budgets. Keeping a complete alternative forest is a different contract.

## Exact numeric outside and profiles

Exact BPP and profiles store numeric states and regenerate deductions.
Their primary numeric states occupy fourteen dense triangular planes;
presence flags distinguish an absent item from an explicitly stored zero.
Four-coordinate BE states retain sparse keys. Initializing these planes
costs O(n^2) time and space. Numeric marginal collection traverses them
directly without materializing a second hash map. Finite-beam charts remain
sparse, so they do not acquire this quadratic allocation.

With fixed local-loop model size, peak logical space is O(n^2 + nr).
This bound does not include the separate traceback storage used by exact
MFE, posterior decoding, sampling or subopt.

Let S be exact schedule size, D the number of generated deductions, and
T_exact the complete exact forward cost, including geometry and generation.
S = O(n^2), including nested scaffold BE pairs. With expected constant-time
sparse BE lookup, the bounds are

```text
BPP:      O(T_exact + D + S + n) time
Profiles: O(T_exact + (r+1)D + rS + nr) time.
```

The finite-beam theorem does not establish an exact whole-call cubic-time
bound. Fixed local-loop size does not by itself bound every geometry scan.
Skipping evaluation after a zero child and omitting VM without an admissible
V consumer reduce work without changing the nonzero derivations. They do
not establish a smaller worst-case bound on deduction generation.

## Sampling and suboptimal enumeration

Sampling retains a derivation representation for repeated traceback:
a factorized forest in finite-beam mode and an arithmetic LogSumExp tape
in exact mode. Subopt uses a factorized forest. Their storage is not
covered by the numeric O(nb) theorem. Account for representation construction,
sampling choices, traceback work and output separately. Finite-beam sampling
scans term alternatives at each visited state; it does not provide a
constant-time choice guarantee independent of that state's alternatives.

For factorized forests, let R = S+T+F count stored states, terms and
child-factor occurrences. For an exact sampling tape, count its stored
arithmetic nodes instead.

Ordinary CLI and `_each` APIs deliver one completed structure at a time,
requiring O(n) completed-output scratch. Collecting APIs additionally retain
O(nK) output. Sampling `--unique` stores O(nU) structure/count entries for
U distinct sampled structures.

For lazy subopt, let H count inserted heap tuples, J their rank-vector work,
and Lout the state/factor occurrences traversed for emitted roots. A
conservative expected-hash time bound is
O(R + J + H log(H+2) + Lout + nK). Space includes the forest, stored rank
vectors, candidate/cache/seen entries and up to O(nK) ambiguity keys.
Streaming removes the completed-output vector; it does not remove those
internal caches or the K dependence. Output alone costs Omega(nK) time.

## Scope and measurement

Scaffold proposal and evaluation of multiple scaffolds are separate costs.
Caller-retained outputs and optional diagnostics also need their own budgets.
Measured wall time or RSS is evidence about particular runs, not a proof of
an asymptotic bound. Hash collision worst cases are outside the expected
dictionary model used above.

See [validation](validation.md), the [C++ API](../src/scfg2/api/README.md),
and the [shared executor](../src/scfg2/README.md) for implementation entry points.
