# Public C++ APIs

See [computational complexity](../../../documents/complexity.md) for time,
space, geometry, forcing and output assumptions. Under fixed local caps,
standard forcing, 1 <= b <= n and expected constant-time dictionaries,
finite-beam numeric inside/outside, loop profiles, MFE and MEA/centroid use
O(nb^2 log(b+2)) time and O(nb) peak logical space. Accessibility for r
requested window entries adds O(rnb^2+nr) time and O(nr) space.
Exact calls, sampling/subopt traceback storage, caller-retained output and variable
output enumeration have separate contracts.

## Structural profiles and accessibility

Include `scfg2/api/profiles.hh` and link `LinearKnot::scfg2_beam` for both exact
and beam profile calls:

```cpp
using namespace cparty::scfg2;
const auto profile = api::beam_structural_profiles(sequence, scaffold_g,
                                                  {1, 3, 5}, {}, true);
const auto exact = api::exact_structural_profiles(sequence, scaffold_g,
    {1, 3, 5}, WmbClosurePolicy::RequireWeakClosure, true);
if (profile.ok) {
  // Window index 1 requests length 3; starts are one-based.
  const double total = profile.accessibility[1][2];
  const auto &by_loop = profile.accessibility_by_loop[1][2];
  // by_loop: Hairpin, Bulge, Internal, Multi, Exterior, PseudoLoop.
  // Sum over these six entries equals total.
}
```

The final `by_loop` argument defaults to false; the additional result vector
is then empty. When enabled, `accessibility_by_loop[window][start][context]`
is the joint probability that the **entire window** is unpaired in that
context. It is neither a conditional type distribution nor a product of
per-base marginals. Context order is `accessibility_loop_context_names`;
for length one, entry `c` equals `contexts[start][c+1]`. Entries beyond the
last valid start (and oversized windows) are zero. Requested lengths retain
their input order, including duplicates. The opt-in six channels use the
same difference-array and chain pass as total accessibility, adding O(nr)
space without changing the existing profile time bound for r lengths.

`linearknot access --by-loop` exposes the same values in an extra column,
ordered `E,H,B,I,M,P` to follow LinearRaccess, with PseudoLoop added as `P`.
The existing total column follows `--probability`; the extra column always
contains joint probabilities. `--byloop` and `-byloop` are aliases.

## Linear/beam API

From an external CMake project, link the beam target:

```cmake
add_subdirectory(/path/to/LinearKnot LinearKnot-build)
target_link_libraries(my_folder PRIVATE LinearKnot::scfg2_beam)
```

Include `scfg2/api/beam.hh`. The simplest LinearKnot Viterbi prediction is:

```cpp
#include "scfg2/api/beam.hh"

using cparty::scfg2::api::BeamSession;

BeamSession session(sequence, scaffold_g, 100);
const auto best = session.viterbi();

if (!best.ok) {
  // Invalid input or internal/numerical failure: inspect best.error.
} else if (!best.reachable) {
  // This finite beam disconnected the root; increase the beam size.
} else {
  // best.generated_structure is G'.
  // best.union_structure is R = G union G'.
  // best.energy is available even when exp(best.log_weight) underflows.
}
```

`BeamSession::viterbi()` runs the shared online deduction kernel in LogMax,
using the configured public-state ranking. Inside uses LogSumExp. These main
semirings can retain different beams, even with the same completion guide.

Beam MFE stores compact selected rule/split descriptors for traceback, with
full records for exceptional deductions. It still retains the LogMax chart
and selected choices, so its constant memory cost can exceed inside alone.
See [the MFE storage contract](../../../documents/complexity.md).

The configurable form is:

```cpp
cparty::scfg2::api::BeamConfig config;
config.beam_size = 100;
config.max_hairpin_unpaired = 30;
config.max_case3_ancestor_depth = 30;
config.ranking = cparty::scfg2::api::BeamRanking::Outer; // default
config.outer_pilot_exponent = cparty::scfg2::api::OuterPilotExponent::TwoThirds;

BeamSession session(sequence, scaffold_g, config);
```

### Outer and baseline ranking

The finite-beam default is **Outer**. A LogSumExp pilot of width
`p=ceil(b^(2/3))` builds bounded max completion-context scores `h_p(s)`.
Public main candidates use `log alpha(s) + log h_p(s)` when a guide value is
available, and keep their previous priority when it is absent. The candidate
need not have appeared in the pilot. `OuterPilotExponent::Half` selects
`p=ceil(sqrt(b))`; the additional guide frontier is fixed at `4p`.
Private partial-product selectors and WMB's raw internal view keep their
existing selection. The guide affects retention, not thermodynamic weights.

Set `config.ranking = BeamRanking::Baseline` to reproduce the previous
ranking: the W-prefix priority for V/WMB's root view, raw priorities for
other public families. This is independent of WMB closure/model policy.
KeepAll (`beam_size == size_t::max()`) skips guide construction because no
public ranking decision can change retention.

The seven public `linearknot` subcommands (`fold`, `access`, `profile`,
`sample`, `subopt`, `mea` and `centroid`) share:

```sh
build/linearknot fold -i scaffolds.txt --beam 100                 # outer, q=2/3
build/linearknot fold -i scaffolds.txt --ranking baseline         # prior ranking
build/linearknot fold -i scaffolds.txt --outer-pilot-exponent 0.5  # smaller pilot
build/linearknot fold -i scaffolds.txt --exact                    # no pruning
```

`--ranking outer|baseline` defaults to outer. The pilot option accepts only
`0.5` and `2/3` and requires outer ranking. Explicit ranking/pilot options
with `--exact` are errors; an ordinary `--exact` call remains valid. Raccess
metadata and MEA/Centroid/Subopt `--stats` report the selected ranking.

A serial BeamSession caches immutable numeric guide scores and scaffold-only
forced keys, so inside, outside, Viterbi and target-probability calls reuse
one guide. Configure the thermodynamic model before constructing the session
and keep it fixed for its lifetime. Outside reconstructs WMB consumer views
and BE prefix FIFO ordering with the same guide. Sampling and subopt
build their factorized forest from the guided LogSumExp execution.
Posterior decoding replays the same alpha-ranked execution with selected traces.

The guide adds O(np) saved scores under the bounded guide assumptions.
Its conditional construction bound is O(np² log(p+2)), contained in O(nb²)
for either fixed exponent; this is not a whole-pipeline proof including
scheduling, forced exceptions and outputs. Benchmark the complete pilot,
outer and main computation. Reported elapsed times include any preparation
performed by that call; CLI wall time includes preparation even when Viterbi
constructs the guide first. Subsequent session calls reuse it. `BeamWorkStats`
counts the main DP, excluding pilot/guide work. Standalone sampling/profile/decoding
calls include their own preparation in total elapsed time.
Outer is a heuristic; retained mass and speed need evaluation together,
and neither improvement for every input nor equal-time superiority is promised.
See [ranking model and scope](../../../documents/model.md#beam-and-exact-computation).

The two geometry limits are independent of beam width. The defaults define
the finite candidate policy used by the benchmarks.
`max_case3_ancestor_depth` counts only
branching scaffold ancestors (pairs or the exterior with at least two pair
children), because only those can be the common ancestor of a CASE3
crossing pair's endpoints. Stacked and interior-loop pairs are not counted,
so the limit is exact whenever no unpaired position has more than that many
branching ancestors. Every result's `work.case3_truncated_endpoints` counts
the positions that exceeded it; zero means the limit removed no CASE3
candidate from that run. Report it with benchmark results. Set either field to
`std::numeric_limits<std::size_t>::max()` only for an explicit unlimited
geometry experiment; finite-width resource bounds assume bounded geometry
limits as specified in the complexity contract.

The same session exposes log-domain inside:

```cpp
const auto z = session.inside();
```

`z.log_partition_weight` is canonical. `z.partition_weight` is only a
convenience projection and may underflow to zero or overflow to infinity.
`z.ensemble_free_energy` is computed directly from the log weight.
`z.elapsed_seconds` and `z.work` report timing and the same production-path
work counters exposed by Viterbi.

The retained-DAG outside and base-pair marginals are:

```cpp
const auto outside = session.outside();

for (const auto& pair : outside.pair_probabilities) {
  // pair coordinates are 1-based.
  // pair.round_probability: V/Round ownership.
  // pair.square_probability: VP-family/Square ownership.
  // pair.union_probability: probability that the union contains this pair.
}
// outside.unpaired_probabilities is also 1-based; index zero is unused.
```

Fixed scaffold pairs have union and Round probability one. Generated Round and
Square probabilities come from the analytic beta recursion on immutable
retained numeric inside states. Outside regenerates the admitted deductions
using the shared grammar/weights; it does not independently prune a second
beam. Consequently
`outside.log_partition_weight` and the same-session inside log partition refer
to the same retained derivations.

Beam outside saves public alphas, one rule mask per surviving public item,
compact geometric/rank indices and the genuine VPR/BTAIL private cells. It
saves no branch forest or arithmetic tape. Public pruning commits before
downstream consumption, so a removed public item has no outside contribution.
For finite public/private widths bounded by b, the reverse representation uses
O(nb + b²) space (O(nb) for b <= n); forced-state budgets and requested results
are additional. Outer's O(np) numeric guide is separate. Whole-process memory
also includes inside scheduling and scaffold-forcing preparation; see the
complexity contract for the O(nb) whole numeric-call bound under standard
forcing and b <= n. Arbitrary forced sets are outside that claim.
`reverse_records`, `reverse_terms`, `reverse_factors` and
`reverse_storage_bytes` describe the backend. The retained backend counts
public/private cells in `reverse_records`; terms, factors and tape nodes are
zero. Its hash-allocation byte estimate is unavailable, reported as zero;
use process RSS to measure memory. `elapsed_seconds` includes forward, beta,
marginal collection and any uncached guide preparation. Beam CapR/Raccess also
regenerate outside and stream emissions into output-sized range differences.
Sampling/Subopt keep their forest; beam MEA/centroid use numeric regenerated
posteriors and selected max traces.
`ExactSession::outside()` and exact CapR/Raccess also use numeric alpha/beta
and regenerate the original exact deductions. They retain the exact span
schedule and boundary child aliases, with no beam, pilot or arithmetic tape.
Exact profiles stream whole-run emissions into difference arrays and keep
only state chain steps for joint-unpaired windows. Their logical peak space
is O(n²+nr), with r requested window-length entries, including duplicates
(r=0 for BPP/CapR), under the
fixed local-loop model. Exact runtime still includes its exhaustive candidate
generation; the finite-beam time theorem does not cover it. The exact CLI
uses inside alone for ensemble free energy when BPP is not requested.
Tape BPP/profile references remain in tests; exact MEA/centroid, sampling and
subopt keep their separate graph/forest contracts. See the
[exact migration and proof](../../../documents/complexity.md#exact-numeric-outside-and-profiles).

For an inside width-saturation ladder:

```sh
build/scfg2_beam_inside_saturation_tool \
  SEQUENCE SCAFFOLD_OR_empty 20,50,100,200,400 \
  params/rna_Turner04.par noexact
```

Use `exact` instead of `noexact` only where the exact calculation is practical.
Its canonical log partition need not have a representable raw projection.
The tool reports consecutive
`delta_log_z`, free-energy change, monotonicity, timing, and work counters.

For a base-pair/unpaired probability saturation ladder:

```sh
build/scfg2_beam_outside_saturation_tool \
  SEQUENCE SCAFFOLD_OR_empty 50,100,200,400 \
  params/rna_Turner04.par
```

The tool runs each width independently and compares every pair probability
and every unpaired probability both with the preceding width and with the
largest requested width. Pair means use all `n(n-1)/2` coordinates, treating
absent sparse entries as zero. It also reports elapsed time and analytic
reverse-record count; the largest width is a convergence reference, not an
exact oracle.

Finite target probability is fail-closed:

```cpp
const auto p = session.probability(target_gprime);

if (!p.ok) {
  // Invalid input or numerical failure.
} else if (!p.denominator_reachable) {
  // The inside-ranked finite beam disconnected the root.
} else if (!p.target_included) {
  // The target was pruned from that denominator DAG. No quotient is returned.
} else {
  // p.probability <= 1 on one shared inside-ranked retained DAG.
}
```

This probability run is intentionally different from independently dividing
the LogMax Viterbi result by a separately pruned inside result. Independent
semirings can retain different DAGs, for which the quotient need not be a
probability. `target_free_energy` is the free energy of all retained
derivations producing the target pair map; it is not an unambiguity claim.

The current predictor consumes `(sequence, scaffold G)`. It predicts the
additional pair map `G'` and the union `R=G union G'`. Supplying an all-dot
scaffold is accepted, but this API does not claim to solve the separate
sequence-only scaffold-proposal problem. A proposal stage such as the optional
HotKnots sidecar may be composed before `BeamSession`.

Public configuration and result fields use standard-library types, shared
public pair types and policy enums. The session's private cache has an opaque
implementation type; callers do not consume internal replay/frontier objects.

## Streaming structures

`api/sampling.hh` exports `exact_sample_each` / `beam_sample_each`;
`api/subopt.hh` exports `exact_subopt_each` / `beam_subopt_each`.
They deliver one structure synchronously to a callback and return metadata
without a completed-output vector:

```cpp
#include "scfg2/api/sampling.hh"
auto result = cparty::scfg2::api::beam_sample_each(sequence, scaffold_g, K, seed,
    [&](auto &&sample) {
      output << sample.union_structure << '\n';
      return true;  // false makes this the last delivered structure
    }, beam_config);
if (!result.ok) throw std::runtime_error(result.error);
```

For subopt, use `beam_subopt_each(sequence, scaffold_g, consumer, config, beam)`
(or `exact_subopt_each(sequence, scaffold_g, consumer, config, policy)`).
The argument can be moved into caller-owned storage, or is released after
the callback returns. Callbacks run in draw/rank order and may stop after
any structure; no following draw/rank is requested. The runtime's serial
model-configuration restriction still applies. A pull `next()` iterator is
not exposed by this synchronous interface.

`SampleStreamResult::samples_emitted` / `SuboptStreamResult::structures_emitted`
count callback invocations, including a callback that returns false or throws.
A successful consumer stop sets `stopped_by_consumer=true` (sampling) or
`stop_reason=Consumer` (subopt). Standard and nonstandard callback exceptions
return `ok=false` with an error. Already delivered results cannot be revoked;
inspect the final summary. Elapsed time includes callback work.

The existing collecting APIs use the same generator and retain their error
rollback behavior. Their O(nK) completed-result storage is opt-in by choosing
those APIs. Sampling's streaming API retains only one completed structure
in addition to its existing forest/traceback storage. Subopt also removes
its completed-output vector, but retains rank/heap caches and exact physical
pair-map keys for ambiguity detection; those keys alone may use O(nK).
Streaming subopt is consequently not a K-independent-memory algorithm.

`linearknot sample` and `linearknot subopt` use streaming delivery for ordinary
output. Sampling `--unique` keeps U distinct strings/counts, using O(nU)
aggregation storage and sorting frequencies after sampling. A late error may
leave a previously written prefix; stderr and a nonzero exit mark failure.
See [storage budgets and validation](../../../documents/complexity.md#sampling-and-suboptimal-enumeration).

## Ranked suboptimal structures

`api/subopt.hh` exports `exact_subopt(sequence, G, SuboptConfig, policy)` and
`beam_subopt(sequence, G, SuboptConfig, BeamConfig)`. Both are linked through
`LinearKnot::scfg2_beam` and use the same lazy k-best consumer and local weights.

```cpp
#include "scfg2/api/subopt.hh"
cparty::scfg2::api::SuboptConfig config;
config.top_k = 20;
config.delta_energy = 3.0;  // optional kcal/mol, relative to this forest's MFE
auto result = cparty::scfg2::api::exact_subopt(sequence, scaffold_g, config);
```

K must be positive. A supplied delta must be finite and nonnegative. Results
contain up to K `SuboptStructure` entries in nondecreasing energy order, with
canonical log weight, energy, ordered generated pairs and generated/union
strings. Equal-energy ordering is deterministic within a forest; exact and
beam do not promise the same tie order. No extra heuristic pruning is used
by the enumerator. Fewer than K results are valid when support or band runs
out. `stop_reason` is `TopK`, `DeltaEnergy` or `Exhausted`; `TopK` means the
requested prefix was produced, without an extra exhaustion lookahead.

Exact ranks the full shared grammar at fixed G. Beam ranks the **inside-ranked
retained forest**, as sampling/MEA do; its best energy can differ from the
separately LogMax-pruned `BeamSession::viterbi()` result. The optional band is
inclusive relative to this same forest's optimum, with `1e-9` kcal/mol
rounding tolerance. It cannot recover a derivation removed during forward
pruning or by the beam hairpin/CASE3 caps.

Unambiguity is a grammar requirement. Public wrapper aliases may emit the
same pair within one derivation; a linear partner map coalesces those
observations. Two complete derivations yielding the same physical generated
pair map are an error: a collecting call clears its structures and returns
`ok=false`; a streaming call reports the already-delivered prefix count.
It does not silently deduplicate and continue. Explicit `LegacyCompatible`
closure can therefore fail on the known duplicate witness.

```bash
build/linearknot subopt -i scaffolds.txt --top-k 20 --stats
build/linearknot subopt -i scaffolds.txt --exact -k 20 -e 3
```

Output is name, sequence, scaffold, then `union_structure (energy)` lines.
`--stats` reports ranking scope, MFE, stop reason and forest/candidate/cache
counts. This consumer stores the branch forest and lazy candidate lists;
its memory contract differs from regenerated outside. DP construction,
enumeration, derivation traversal and the `O(nK)` output are separate costs.
See [the complexity and witness note](../../../documents/complexity.md#sampling-and-suboptimal-enumeration).

## MEA decoding

Link `LinearKnot::scfg2_beam` and include `scfg2/api/mea.hh`:

```cpp
auto prediction = cparty::scfg2::api::beam_mea(sequence, scaffold_g, 1.0, config);
// Or exact_mea(sequence, scaffold_g, gamma) for an uncapped reference.
```

`expected_accuracy` maximizes `2*gamma*sum_pairs p(i,j) + sum_unpaired q(i)`
for the completed union `R = G union G'`. Fixed scaffold pairs contribute
their constant `2*gamma` reward and stay in the output. Gamma defaults to one
and must be finite and nonnegative. Results expose `generated_pairs`,
`scaffold_structure`, `generated_structure` and `union_structure`.

Beam posterior and output support come from the alpha-ranked inside execution.
Numeric regenerated outside obtains the marginals without saving an edge
forest. A second execution of the same shared grammar/executor carries alpha
and an additive posterior maximum together. All pruning and representative
choices still use alpha, so the admitted support is reproduced; the posterior
score never selects the beam. Only the best decomposition is retained, with
first forward alternative winning exact ties. The baseline is
`sum_i q(i) + 2*gamma*|G|`; each generated pair adds
`2*gamma*p(i,j)-q(i)-q(j)`. Negative connector rewards remain available.
Thermodynamic weights determine the ensemble/support; they are not added
to the MEA objective. Exact decoding retains the frozen-forest reference
consumer. The serial-call restriction of the runtime adapter still applies.

`ok=true, reachable=false` means the retained root was lost. Other errors,
including invalid gamma or posterior bounds, have `ok=false`.
For beam calls, `decoding_states` is the peak number of live selected-DAG
nodes and `decoding_terms/factors` count the chosen root DAG visited at output.
These counts depend on the objective; they no longer describe every forest
alternative. Exact calls retain the former whole-forest counts.
`decoding_seconds` includes reward assembly, alpha/max replay and output;
`elapsed_seconds` also includes numeric inside/outside and ranking preparation.
Under fixed local caps, standard forcing, `1 <= b <= n` and expected constant
time dictionary operations, beam MFE/MEA/centroid take `O(nb² log(b+2))` time
and `O(nb)` peak logical space. Exact/unlimited-cap calls are excluded.
See [the consumer storage contract](../../../documents/complexity.md).

The CLI reads the same sequence/scaffold records as the other consumers:

```sh
build/linearknot mea -i scaffolds.txt -b 100 --gamma 1
build/linearknot mea --sequence AAAGGGAAAGGGCCCAAACCCAAA \
  --scaffold '...(((......))).........' --gamma 1000
```

Output includes the MEA union and its `EA` score. `--exact` selects the exact
reference, `-o FILE` writes records, and `--stats` prints timing/work to stderr.
The scaffold's LinearPartition MEA preprocessing is separate from this PK
decoder. See `documents/model.md` for
the legacy source audit, exhaustive oracle and length-scaling evidence.

## Hamada gamma-centroid decoding

Include `scfg2/api/centroid.hh` and link `LinearKnot::scfg2_beam`:

```cpp
auto prediction = cparty::scfg2::api::beam_centroid(sequence, scaffold_g, 1.0, config);
// Or exact_centroid(sequence, scaffold_g, gamma).
```

This maximizes the pair-level gain `E[gamma*TP + TN]`, using generated-pair
reward `(gamma+1)*p(i,j)-1`. It shares the MEA forest, band geometry,
pair-owner mapping and max-plus consumer. Gamma affects decoding only.
Fixed G is preserved and contributes gamma per pair. There is no independent
probability-threshold projection; admissibility comes from the shared grammar.

`CentroidResult` and `MeaResult` alias `PosteriorDecodingResult`.
`expected_accuracy` is the full expected gain of the called objective;
`expected_true_positive_pairs` and `expected_true_negative_pairs` are exposed.
TN ranges over every unordered position pair, including impossible pairs.
Its constant is assembled from sparse posterior mass in O(pair-count+n),
without a dense pair scan. The finite/exact, gamma, serial-call, error and
work/timing contracts of MEA apply.

```sh
build/linearknot centroid -i scaffolds.txt -b 100 --gamma 1 --stats
```

Output includes `EG`, `ETP` and `ETN`, with `EG=gamma*ETP+ETN`.
`linearknot mea` retains the legacy base-level objective and its `EA` output.
See `documents/model.md` for the
gain definition, different-output witness and exhaustive/length evidence.

## Hotspot scaffold proposal API

The optional sequence-only scaffold proposer is a separate lightweight
library:

```cmake
target_link_libraries(my_folder PRIVATE LinearKnot::scfg2_hotspot_proposal)
```

```cpp
#include "scfg2/proposal/hotspot_wavefront.hh"

cparty::scfg2::proposal::WavefrontHotspotConfig config;
config.anchor_cap = 8;
config.beam_width = 16;
config.passes = 8;
config.max_results = 24;
config.diversity_mode =
    cparty::scfg2::proposal::HotspotDiversityMode::
        DeterministicPassDiversification;
config.diversity_priority_jitter_amplitude = 0.25;

const auto proposed =
    cparty::scfg2::proposal::propose_hotspots_wavefront(sequence, config);
```

Every successful result contains bounded, deduplicated, pseudoknot-free
candidate scaffolds plus an empty-scaffold baseline. `heuristic_score` ranks
proposal paths only; it is not a CParty energy. To predict a completed
structure, run the same `BeamSession` configuration for every proposed `G`
and select the reachable result with minimum returned CParty energy.
Equal completed maps `R=G union G'` must retain their distinct split
provenance because the current legacy-aligned energy is not split-invariant.

The diversity mode is a proposal-search option, not part of scoring.
`DeterministicPassDiversification` varies bounded anchor representatives and
intermediate beam-retention priority across passes while remaining
reproducible. It does not alter exported `heuristic_score` or final CParty
energy. Set `diversity_mode = HotspotDiversityMode::None` for a single
unperturbed search policy; additional identical passes are then redundant.
Increasing `max_results` returns more completed survivors but cannot recover
an alignment discarded at an earlier beam boundary.

The high-level repeated-Viterbi API performs that policy directly:

```cmake
target_link_libraries(my_folder PRIVATE LinearKnot::scfg2_hotspot_viterbi)
```

```cpp
#include "scfg2/api/hotspot_viterbi.hh"

const auto prediction =
    cparty::scfg2::api::predict_hotspot_viterbi(
        sequence, proposal_config, beam_config);

if (!prediction.ok) {
  // Proposal failure or an incomplete Viterbi comparison.
} else if (!prediction.reachable) {
  // Every proposed scaffold lost the finite-beam root.
} else {
  const auto& winner =
      prediction.structures[prediction.selected_structure];
  const auto& selected_split =
      prediction.records[winner.selected_record];
}
```

Every proposal is run once, serially, with the exact same `BeamConfig`.
`records` retains successful, unreachable, and provenance-bearing candidate
runs. Reachable records are grouped by the completed base-pair map, ignoring
whether its rendered pairs use round or square brackets. Each group exposes
all source record indices, minimum/maximum energy, energy spread, and the
selected split. A failed Viterbi call makes the whole comparison fail closed;
an ordinary unreachable finite-beam root remains a valid record.

`select_minimum_energy_records(...)` exposes the same deterministic grouping
and selection policy for consumers with a different scaffold proposer.

`proposal_elapsed_seconds`, `viterbi_elapsed_seconds`, and
`total_elapsed_seconds` report stage timing. Each
`ScaffoldViterbiRecord::elapsed_seconds` reports its individual Viterbi time.
Each `BeamViterbiResult::work` also reports aggregate item/frontier counts,
pk-free join and extension probes, `V` parent probes, public `VP/WMBP/CASE3`
candidate counts, and per-nonterminal counts from the exact LogMax prediction
run. These are observations of the production path rather than a separate
profiling semiring.
For command-line measurement, build and run:

```sh
build/scfg2_hotspot_viterbi_tool \
  SEQUENCE params/rna_Turner04.par \
  BEAM_WIDTH ANCHOR_CAP HOTSPOT_WIDTH PASSES MAX_RESULTS \
  deterministic
```

`MAX_RESULTS` includes the mandatory empty-scaffold baseline. Thus `1` means
baseline only, while `2` is the smallest hotspot-versus-baseline comparison.

Cluster consumers that need candidate-level process parallelism use the
observation-only tools:

```text
build/scfg2_hotspot_proposal_tool
build/scfg2_fixed_scaffold_viterbi_tool
test/tools/run_scfg2_parallel_hotspot_viterbi.py
```

The proposer is invoked once. Each fixed scaffold is then evaluated in an
isolated worker process, after which the coordinator applies the same
pair-map grouping and deterministic minimum-energy selection as
`select_minimum_energy_records()`. This is process parallelism rather than
thread parallelism because thermodynamic parameter loading under `BeamSession`
uses Vienna's process-global state. It does not introduce a second DP kernel.

## Completed-WMB model policy

All exact and beam calls, including CapR/Raccess, sampling and posterior
decoding, default to `WmbClosurePolicy::RequireWeakClosure`. Completed
`WMB(i,j)` requires a weakly closed interval in fixed scaffold G. Partial
`WMBP/WMBW` carriers retain their existing domains. Local weights are unchanged.

To reproduce the previous ensemble, select the policy explicitly:

```cpp
using cparty::scfg2::api::WmbClosurePolicy;
cparty::scfg2::api::ExactSession historical(
    sequence, scaffold_g, WmbClosurePolicy::LegacyCompatible);
config.wmb_closure_policy = WmbClosurePolicy::LegacyCompatible;
cparty::scfg2::api::BeamSession historical_beam(sequence, scaffold_g, config);
```

Exact free functions accept a final optional policy argument, after gamma,
seed, target or accessibility lengths where applicable. Beam functions use
`BeamConfig::wmb_closure_policy`. Results carry `wmb_closure_policy`;
`wmb_closure_policy_name(...)` returns a reproducible name. Session inside
caches cannot mix policies. Existing CLIs use the corrected default.

This correction changes the conditional ensemble, not just its evaluator.
Unguarded legacy parity is scoped to `LegacyCompatible`; independent corrected
oracle evidence and remaining Psm/split-invariance limitations are recorded in
`documents/model.md`.
No general unambiguity or full SI weight-equivalence theorem is claimed.

As of 2026-10-11 the official oracle is upstream `800e95a` plus local
0001/0002/0006/0007, without completed-WMB closure. It compares exact inside
log Z and PF free energy in compatibility mode; BPP and legacy-MFE equality
are not part of that contract. Pin, build identity and reproduction are in
`test/oracle/legacy_cparty/README.md`.

## Exact API

From an external CMake project, add this repository and link its public target:

```cmake
add_subdirectory(/path/to/LinearKnot LinearKnot-build)
target_link_libraries(my_inverse_folder PRIVATE LinearKnot::scfg2_exact)
```

The target exports its C++17 language requirement and public include path.

Include `scfg2/api/exact.hh` and use
`cparty::scfg2::api::ExactSession` for external exact-folding consumers.

```cpp
#include "scfg2/api/exact.hh"

cparty::scfg2::api::ExactSession session(sequence, scaffold_g);

auto z = session.inside();
auto best = session.viterbi();
auto energy = session.evaluate_energy(target_gprime);
auto p = session.probability(target_gprime);
```

`ExactSession::inside()` is cached.  Repeated energy/probability evaluations
for the same `(S,G)` reuse `Z`.

Exact and beam production use unscaled log local weights from the same PF
parameter source. Exact `InsideResult::log_partition_weight`,
`ViterbiResult::log_weight`, `EnergyEvaluationResult::log_weight` and
`ProbabilityResult::{log_partition_weight,log_target_weight}` are
canonical. Raw `partition_weight`, `weight` and `target_weight` fields are
convenience projections that may be zero or infinity even for a successful,
reachable call. Use `ok` and `reachable` for status. Free energies use the
finite log value directly; exact `log_probability` is
`log_target_weight-log_partition_weight`. The known best-derivation versus
same-output ambiguity contract is unchanged. Production does not build
length-dependent PF scale tables or apply a PF scale correction.

Pair-list overloads use 1-based positions:

```cpp
using cparty::scfg2::api::BasePair;
auto p = session.probability(std::vector<BasePair>{{5, 20}, {6, 19}});
```

Each input pair accepts either endpoint order and is normalized to `i < j`.
Self-pairs, coordinates outside the sequence, and repeated endpoints are
rejected. Output pair lists always use `i < j`.

`scaffold_g` and a dot-bracket target are strings of the same length as the
sequence. The scaffold must be noncrossing and encoded with balanced `()`;
`[]`, `{}`, `<>` and unknown symbols are rejected even when those brackets
would describe a noncrossing pair map. Use `.` or `_` for free bases and
`x` or `X` for forced unpaired bases. The target may use `()`, `[]`, `{}`, or
`<>`; all its pairs are treated as the generated `G'` map relative to the
fixed scaffold supplied to the session. A full union target containing
scaffold pairs is also accepted; fixed scaffold pairs are removed by the
evaluator boundary.

`ViterbiResult::generated_pairs` and `generated_structure` contain every
selected pair not already fixed by scaffold `G`. This ownership classification
is independent of the grammar lane: `V` and `BE` are Round grammar pairs and
`VP` is Square, but a newly selected Round `V` pair still belongs to the public
generated `G'` map. In `generated_structure`, generated Round pairs retain
`()`, while Square `VP` pairs retain `[]`. A Square pair is part of the
density-2 crossing extension relative to fixed scaffold `G`. If `G` is empty,
the public/open `VP` family is inactive and the result is necessarily
pseudoknot-free.
`scaffold_structure` is the original session scaffold and is entirely Round.
`union_structure` is the completed structure `R = G union G'`, formed by
overlaying those two endpoint-disjoint maps.  Thus the three public fields are:

- `scaffold_structure`: fixed `G`;
- `generated_structure`: generated difference `G'`;
- `union_structure`: completed `R = G union G'`.

Check both status fields:

```cpp
if (!p.ok) {
  // Invalid input or numerical failure: inspect p.error.
} else if (!p.reachable) {
  // Fail-closed constrained DP returned zero for this target.
} else {
  // p.probability, p.log_probability, p.energy are available.
}
```

The public header contains only standard-library types; the internal replay
and weight snapshots stay behind the linked library.

Thermodynamic parameter loading uses Vienna's process-global state.
Configure `CPARTY_FIXED_ENERGY_PARAM_FILE` before constructing an exact
or beam session (for example `params/rna_Turner04.par`) and do not call
sessions concurrently. Keep this configuration unchanged throughout a
session's lifetime, including calls that use cached inside results. Each
weight owner loads an explicitly selected file before snapshotting its PF
tables and long-loop parameters; an empty, unreadable or unloadable file
selection returns an error. Files must have the v2.0 header and at least one
recognized parameter section. Valid partial files overlay the current global
model; header-only and unknown-section-only files are rejected before any
model mutation. Section-body parsing remains Vienna's existing implementation.
Without an override, the current Vienna global
parameters are used (Turner2004 in a fresh process). CLIs use these built-in
defaults without requiring a relative `params/` directory.

The weight boundary normalizes sequence letters to uppercase RNA (`T` to
`U`) for both encoded bases and special-loop lookup. Session input getters
retain the supplied spelling. Thread-safe ownership is a later boundary
cleanup; the APIs do not hide the current limitation.

For an installed library, use `cmake --install build --prefix PREFIX`, then
`find_package(LinearKnot CONFIG REQUIRED)` with `PREFIX` in `CMAKE_PREFIX_PATH`.
The installed `LinearKnot::scfg2_exact` and `LinearKnot::scfg2_beam` targets propagate
the bundled static `LinearKnot::RNA` dependency and include path.

See `src/scfg2/README.md` and
`documents/scfg2/energy_evaluator_contract_20260716.md` for semantic guarantees
and non-guarantees.
