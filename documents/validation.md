# Build and validation

## Native build

Use a C compiler and a C++17 compiler with CMake >=3.15. From the checkout:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLINEARKNOT_ENABLE_LTO=OFF
cmake --build build -j8 --target scfg2_exact_inside_oracle_test scfg2_rules_core_exact_weight_test
build/scfg2_exact_inside_oracle_test
build/scfg2_rules_core_exact_weight_test
test/tools/run_required_scfg2_tests.sh
```

The maintained runner builds with -j8 and checks grammar, weights, frozen
exact-inside witnesses, retained/regenerated outside, profile normalization,
sampling distributions, ranked decoding, exact/beam consistency, CLI behavior,
installed C++ consumers and resource witnesses. Optional benchmark drivers are
not a substitute for this acceptance suite.

## External semantic oracle

The primary oracle is upstream CParty inside at a pinned commit plus the
ordered local corrections in [the oracle configuration](../test/oracle/legacy_cparty/config.json).
Use the [oracle build instructions](../test/oracle/legacy_cparty/README.md)
for renewed random comparison. Building that optional oracle requires network
access; the maintained frozen-witness suite does not require downloading it.

Compare the oracle with explicit `WmbClosurePolicy::LegacyCompatible`.
The production default requires completed-WMB weak closure. Legacy MFE,
sampling or BPP are not interchangeable semantic references for inside.

## Focused regression witnesses

- [Private column storage](../test/scfg2_column_store_test.cc) compares
  indexed source buffers with an independent full-key map, including ordered
  deduction fields, empty-but-present moved buffers, noncanonical keys,
  pointer lifetimes, copy/move, graph ownership, throwing value constructors
  and allocation failures in both metadata and optional source reservations.
- [Beam numeric lookup](../test/scfg2_beam_numeric_lookup_test.cc) compares
  indexed and reference execution under both ranking and closure policies,
  including full retained chart iteration, value bits, choices, ties and
  work counts. Storage tests cover [row growth and deletion](../test/scfg2_primary_lookup_index_test.cc)
  and [chart fallback and ownership](../test/scfg2_numeric_primary_lookup_test.cc),
  including allocation failure after a canonical map update.
- [Completion lookup](../test/scfg2_bounded_outer_alpha_lookup_test.cc) checks
  full-key and absent-key values against the original alpha map, as well as
  bitmap ranks across word boundaries, nonfinite and signed-zero values,
  sparse/mixed-row fallback, large coordinates and index allocation failures.
  It also compares position deduplication against sort/unique across word boundaries,
  out-of-range sentinels and optional allocation failures.
  [Bounded maxima](../test/scfg2_outer_max_pool_test.cc) compare the heap
  against a frozen set implementation after each offer, including map
  iteration/buckets, links, rank ties, eviction and copy/move lifetimes.
  Trusted point mirrors additionally check fixed-family/end and fixed-family/span
  domains, full-key fallback and disabling metadata before ownership moves.
  [Numeric beam frontiers](../test/scfg2_right_boundary_item_frontier_test.cc)
  compare start-indexed storage against the default map, including retained
  callback order, pruned sets, chart iteration/bits and unsupported-key fallback.
  Memory-observed runs keep their original containers; the beam lookup witness
  also compares the unobserved numeric frontier against that path.
  [Scaffold forcing](../test/scfg2_beam_scaffold_forcing_test.cc) compares the
  forced traceback set and iteration order against the full-choice observation
  path; the private standard-context path omits only unreachable choice records,
  preserving all numeric states and finite-beam competition.
  [V agenda](../test/scfg2_v_retained_agenda_test.cc) checks complete expanded
  deductions and arrival order against shared generation, with nonempty
  scaffolds, excluded positions, duplicate sources and column release.
- [Numeric exact execution](../test/scfg2_exact_numeric_execution_test.cc)
  compares reference execution with dense storage, zero-child cutoff and
  unused-VM omission, split selection and candidate-buffer reuse, separately
  and together, under both closure policies.
  It checks state bits, Viterbi choices/ties, outside adjoints, chart presence
  and checkpoint boundaries. A zero-child witness retains its nonzero
  outside derivative; reference observations and custom products retain
  their original contract. Ordered deduction fingerprints are frozen from
  the preceding implementation, including inactive child fields and WMBP
  rule-block order. Selected split streams also match filtering the actual
  expanded children, including finite log masses of `-1e300`.
- [Exact live support](../test/scfg2_exact_live_support_test.cc) compares
  dense support rows with full scans across word boundaries, all primary
  families, nonfinite values, erasure, copy/move and optional allocation
  failures. Selected splits match sparse reference scans under both closure
  policies and numeric semirings, including custom-context fallbacks.
  [Schedule checks](../test/scfg2_schedule_invariant_test.cc) compare lazy
  lexicographic dependency checks with the original complete rank tuples;
  later invalid children still produce diagnostics.
- [Shared BE/WMBP availability](../test/scfg2_wmbp_be_availability_test.cc)
  checks the no-overlap boundary against patched-legacy inside/marker values
  and exact/saturated-beam outside.
- [Log-weight tests](../test/scfg2_log_weight_test.cc) cover multiloop closing
  neighbors, PF/log factor agreement and fixed-structure consumers.
- [Runtime snapshot tests](../test/scfg2_runtime_snapshot_test.cc) cover the
  parameter owner and its independent test-only observations.
- [Constrained output grouping](../test/scfg2_constrained_small_n_grouping_test.cc)
  compares small-input pair-map groups and poisons diagnostic carriers to
  check that they do not supply hidden grammar mass.
- [VP split-range fixture](../test/fixtures/scfg2_vp_split_range_20261009/README.md)
  records the single-maximal-scaffold-pair case where the full split sum is
  required. Current energy goldens are in the
  [inside-oracle test](../test/scfg2_exact_inside_oracle_test.cc).

These are finite witnesses for the documented contracts, not general proofs
of unambiguity, constraint completeness or full published-model equivalence.

## Package and installed-command checks

The Python distribution installs a native CLI. Source/wheel checks enforce
payload scope, command names and bundled notices. The package workflows build
sdists and Linux/macOS wheels; the Conda workflow builds and tests the native
package. A configured build workflow is distinct from a completed platform
validation or a published package.

Use [Python installation](../python/README.md) and
[Conda instructions](../conda_recipe/README.md) for local packaging.
`python/smoke_test.py` exercises all seven installed subcommands outside the
source checkout. Publication remains separate from building artifacts.

## Reproducible examples

The [example gallery](../examples/README.md) records the input, settings,
commands and binary hashes. Its plotter checks context normalization,
agreement with pair probabilities, joint accessibility, opening-energy
conversion and suboptimal energy ordering.
