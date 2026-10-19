# Test and benchmark tools

Configure the checkout, then run the maintained acceptance suite:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DLINEARKNOT_ENABLE_LTO=OFF
test/tools/run_required_scfg2_tests.sh
```

The runner builds with -j8 and checks grammar, execution, weights,
frozen exact-inside witnesses, retained/regenerated outside, profiles,
sampling, subopt, posterior decoding, installed C++ APIs and CLI behavior.
It also runs the dependency gate and resource witnesses.

The mandatory exact-inside and local-weight checks can be run separately:

```bash
cmake --build build -j8 --target scfg2_exact_inside_oracle_test scfg2_rules_core_exact_weight_test
build/scfg2_exact_inside_oracle_test
build/scfg2_rules_core_exact_weight_test
```

The optional [external CParty oracle](../oracle/legacy_cparty/README.md)
is pinned and patched. Its disposable checkout is generated under `build/`;
it is not production code or a required download for the frozen suite.

Other scripts provide targeted random probes and benchmark coordination.
Run each script with `--help` for its inputs. Dataset runners validate
input, executable, parameter and configuration provenance before reusing
results. External predictors and datasets are user-supplied.

For resource guarantees and model assumptions, see
[complexity](../../documents/complexity.md) and
[validation](../../documents/validation.md).
