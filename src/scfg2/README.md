# Shared folding kernel

`src/scfg2/` contains LinearKnot's shared grammar, deduction generator,
local weights and DP executor. Exact and finite-beam computation use this
surface; swapping LogSumExp and LogMax selects inside or MFE evaluation.
Outside, sampling, ranked structures and posterior decoding are consumers.

For CLI usage, start with the [main README](../../README.md). For library
integration, see the [public C++ API](api/README.md).

## Layers

| Directory | Responsibility |
| --- | --- |
| `grammar/` | Nonterminals, rules, splits and child schemas |
| `deductions/` | Shared geometric admission and deduction generation |
| `weights/` | Local log weights and thermodynamic parameter boundary |
| `engine/`, `replay/`, `semirings/`, `storage/` | Execution, values, choices and stores |
| `beam/` | Retained states, bounded joins, ranking and regenerated consumers |
| `decoding/`, `sampling/` | Posterior rewards, sampling and lazy k-best |
| `proposal/` | Optional bounded hotspot scaffold proposals |
| `api/`, `cli/` | Library and command-line boundaries |

## Execution contract

The shared executor owns schedule traversal, deduction generation,
admissibility, local-weight lookup, child lookup, semiring aggregation and
chart writes. A valid deduction consumes children earlier in the schedule.
Violations report rule/parent/child/split diagnostics, including in Release
builds. There is no separate pseudoknot-free grammar or scoring kernel.

Inside uses LogSumExp and MFE uses LogMax over model log weights. Exact
computation is unpruned. Finite-beam computation retains bounded public and
private candidates. The default outer ranking uses a smaller inside pilot
to estimate completion context; it changes retention, not local weights.

Production BPP and structural profiles use numeric inside with regenerated
outside contributions. MFE retains selected choices. Beam MEA/centroid
replay alpha-ranked execution while retaining selected posterior maxima.
Sampling and subopt retain alternative derivation representations and have
separate output-dependent resource contracts.

Numeric exact inside, MFE, BPP and profiles explicitly select
`ExactExecutionOptions::numeric()`. Primary spans use dense triangles with
separate presence flags; BE retains its full sparse key. A zero child
(`-infinity`, never a finite cutoff) skips local-weight evaluation. Live
products retain local/child arithmetic order and candidate order. VM(i,j)
is omitted only when the shared V admission predicate is false: its sole
canonical consumer is V_VM at the same endpoints. These changes preserve
all other state values and root results; unused VM entries become zero.
WM, VM, WI, WIP and VPR split selections use the actual completed child
values; WMBP selects nonzero VP starts before its geometry scan. These are
ordered subsequences of the shared generator, with no finite-mass cutoff.
For the standard numeric context, dense charts also index nonzero primary
starts and ends with bit rows. Split selection enumerates these rows in
ascending coordinate order, then applies the shared expansion and complete
child checks. Sparse charts, custom contexts and unsupported shapes retain
the full split scan. The optional index falls back to that scan if its
allocation fails; it does not change chart values or presence.
Numeric runs reuse one candidate buffer while retaining complete generation
before evaluation and the original order of rule blocks and alternatives.

Internal exact runners default to the reference policy, preserving full
VM values and zero-contribution observations. Numeric execution explicitly
allows those observations to be omitted. The additional pre-generation
selection option requires a pure admission predicate: calls on proved-zero
branches may be omitted as well. Arithmetic-tape/factorized
consumers, raw carrier observations and finite-beam retention keep their
existing execution policy. Regenerated outside never applies the forward
zero-product shortcut: a zero child can still receive a nonzero adjoint.

Numeric finite-beam inside and stored MFE select a primary lookup index.
The sparse chart remains authoritative for values, presence and iteration;
the index mirrors checkpointed doubles in per-family, per-end rows. BE and
noncanonical keys use the original full-key lookup. Rows shrink after
pruning, and an index allocation failure falls back to the sparse chart.
Scaffold forcing, the pilot and the completed numeric guide inherit the
option; other semirings and derivation representations retain their
original storage policy.

During outer preparation, immutable canonical alpha rows can use a presence
bitmap and prefix popcounts to locate the original sorted record directly.
An indexed row spans no more bitmap words than retained records. Sparse or
noncanonical rows, BE and optional allocation failures retain binary search.
This changes only point lookup; values, candidate order and pruning are unchanged.

Private agendas confined to one end column store canonical parents in
start-indexed slots. Their independent parent lists and per-parent source
vectors retain traversal and deduction order. Other complete keys use a
sparse fallback. Column boundaries release only occupied values and source
buffers; metadata grows geometrically and contains no default graph values.
If metadata allocation fails, that column uses the full-key map. Aggregate
pointers remain stable within the column and expire at its next reset.
Source vectors request an initial capacity of eight entries; failure of this
optional reservation falls back to ordinary vector growth.

The V agenda queues three-byte rule/wing records in first-arrival order and
expands them through the shared child schema before evaluation. Its duplicate
bitmap covers only wing pairs whose total length is at most MAXLOOP.
Standard numeric runs without memory observations can index pending parents
directly in MAXLOOP+2 rolling end columns. Each occupied cell releases its entire source
buffer on consumption or column expiry. Column stamps prevent aliasing;
unsupported events and optional allocation failure use the original map.
Memory observations and generic contexts retain the map representation.
Completion ranking reads the immutable sorted alpha rows and uses an indexed heap for
bounded maxima, reusing a deduction buffer between fully generated parents.
Its canonical map mutation order, score/key cutoff, equal
score predecessor rule and mandatory-state ownership remain unchanged.
Proposal positions use a reusable bit set for ascending deduplication.
Only touched words are sorted and reset; raw proposal budget checks precede
deduplication, and unsupported positions or allocation failure use the
original sort/unique path.
These storage changes do not alter beam widths, geometry limits or retained
states. The diagnostic sorted rank view is generated only when observed.
Private aggregate reads dispatch by rule and use one lookup for admission
and value retrieval. VP's first three rule blocks share one immutable
geometry observation while preserving their original deduction order.

The model, completed-WMB policy and decoding objectives are documented in
[the model guide](../../documents/model.md). Resource assumptions and bounds
are in [complexity](../../documents/complexity.md).

## C++ integration

For an installed package:

```cmake
find_package(LinearKnot CONFIG REQUIRED)
target_link_libraries(my_folder PRIVATE LinearKnot::scfg2_beam)
```

Or add a source checkout:

```cmake
add_subdirectory(/path/to/LinearKnot LinearKnot-build)
target_link_libraries(my_folder PRIVATE LinearKnot::scfg2_beam)
```

```cpp
#include "scfg2/api/beam.hh"

cparty::scfg2::api::BeamSession session(sequence, scaffold_g, 100);
const auto ensemble = session.inside();
const auto best = session.viterbi();
```

The internal `cparty` namespace is retained. Session results expose status,
log weights and completed/generated pair maps. Configure the thermodynamic
model before creating a session and keep it fixed throughout its lifetime;
the runtime adapter supports serial calls. See the
[API contracts](api/README.md) for lifetimes and errors.

## Validation

From the checkout root:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLINEARKNOT_ENABLE_LTO=OFF
test/tools/run_required_scfg2_tests.sh
```

The maintained suite builds with -j8. It includes exact-inside and local
weight witnesses, cross-mode consistency, outside/posterior checks,
sampling/decoding, installed consumers, CLI and resource checks.
See the [validation guide](../../documents/validation.md) and
[test tools](../../test/tools/README.md).

The external semantic reference is pinned, patched CParty inside.
Its [configuration and patches](../../test/oracle/legacy_cparty/README.md)
are separate from the production runtime. Oracle comparison uses explicit
LegacyCompatible closure policy; the production default requires completed-WMB
weak closure. Legacy Viterbi, sampling and BPP are not inside semantic oracles.
