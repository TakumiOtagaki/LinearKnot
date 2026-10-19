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
