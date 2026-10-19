# WMBP leftward-BE availability witnesses

These 100/130/160-nt fixtures preserve the shared-BE overlap geometry from
the 1600-nt outside-bound failures. A generated GC stem inside a fixed G
subtree makes duplicate consumption numerically visible. Targets are
V(27,69), V(27,99), and V(27,129), respectively.

Before repair, SCFG2 exact counted the target about twice. Patched legacy
inside counts it at most once: its i-descending schedule reads BE with an
opening left of the WMBP parent before that BE is computed. The shared
production geometry now makes that restriction explicit. The maintained
`scfg2_wmbp_be_availability_test` checks legacy inside/marker goldens, maximum
target use one, exact and saturated beam outside, and shared Viterbi parity.
`exact-legacy-order` remains as the earlier diagnostic filter; it is redundant
with the production gate after repair.

```sh
cmake --build build -j8 --target scfg2_outside_pair_audit_tool
build/scfg2_outside_pair_audit_tool test/fixtures/scfg2_wmbp_leftward_be_20261008/stem20.txt exact 27 69
build/scfg2_outside_pair_audit_tool test/fixtures/scfg2_wmbp_leftward_be_20261008/stem20.txt exact-legacy-order 27 69
python3 test/tools/legacy_pair_occurrence_probe.py test/fixtures/scfg2_wmbp_leftward_be_20261008/stem20.txt 27 69 --trace-wmbp 17 98 81
```

The external probe leaves the primary legacy oracle intact, uses all six
oracle patches, disables sampling, and fixes MFE seed 0 / pf_scale=1.
It marks target V visits with x and evaluates Z(x)/Z(1) at 0,0.5,1,2,4.
Its trace observes the zero BE read and the later positive BE write.

Full diagnosis and measurements:
[pair occurrence audit](../../../documents/validation.md).
