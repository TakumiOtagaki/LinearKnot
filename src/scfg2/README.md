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
