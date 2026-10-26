# Scaffold helper defaults matching the paper (0.2.1)

Status: implemented in 0.2.1 (2026-10-26). The comparison below describes
the previous defaults and the paper rule now used by default.

## Problem

Before 0.2.1, `scripts/scaffold/make_scaffold.py` used a different rule from
the paper. Its defaults could produce different scaffolds and therefore
different downstream results. The helper now uses the paper's rule.

| Step | Helper default before 0.2.1 | Paper's experiments |
|---|---|---|
| Starting structure | LinearPartition MEA structure, γ = 3 | none: all pairs from the BPP list |
| Probability threshold | 0.5 | **0.7** |
| Minimum helix length | 1 stacked pair (keep all) | **2 stacked pairs** |
| LinearPartition call | (helper's own) | `linearpartition -V -b 200 --dangles 2 -c 0.7 -r BPP_FILE` |

The paper (Methods, Experimental settings) says: "we construct experimental
scaffolds from LinearPartition pairing probabilities at beam size 200. We
select all pairs with probability at least 0.7, then retain only helices of
at least two stacked pairs."

## The paper's rule (reference implementation)

The rule is `threshold_scaffold(sequence, bpp_path, threshold=0.7,
min_helix=2)` in the experiments repository
(`LinearKnotExperiments/scripts/prepare_biological_inputs.py`, about line 77).
It is called with LinearPartition commit `b450fb3e`, settings
`{"beam": 200, "dangles": 2, "pair_probability_threshold": 0.7, "min_helix": 2}`.

1. Read the LinearPartition BPP file (`i j p`, 1-based, `i<j`).
2. Keep every pair with `p >= threshold`. Each kept pair must be canonical
   (AU, UA, CG, GC, GU, UG). Pairs with p > 0.5 never share a base and never
   cross, so the kept set is a valid secondary structure. The reference
   implementation still checks both and raises an error otherwise.
3. Split the kept pairs into maximal helices: runs (i, j), (i+1, j−1), … of
   stacked pairs. Keep a helix only if it has at least `min_helix` pairs.
4. Write the kept pairs as a dot-bracket string with round brackets.

The `-c 0.7` option only limits which pairs LinearPartition writes. It does
not change the probabilities.

## Change

The helper defaults to threshold 0.7, minimum helix 2, direct BPP selection,
LinearPartition-V beam 200 and dangles 2. The previous behavior is available
as `--from-mea --gamma 3 --threshold 0.5 --min-helix 1 --beam 100`.
Settings and LinearPartition provenance are written in a comment and an
optional JSON sidecar. The README and script docstring describe the rule.

## Test

Check the helper against the experiments repository's saved scaffolds.
Use a few inputs from the bpRNA and RNAcentral cohorts: same sequence and
LinearPartition version should give a byte-identical scaffold. Include:
- one input whose kept pairs contain an isolated pair, which must be dropped;
- one input with no pair above 0.7, which must give an all-dot scaffold.

## Version

This changes the source-checkout helper, not the inference model.
VERSION and release notes are updated to 0.2.1. The paper's existing
v0.2.0 computations used scaffolds made by the experiments repository,
so they do not require recomputation for this helper change.

## Implemented interface

The default is BPP threshold 0.7, minimum helix 2, beam 200, dangles 2.
`--from-mea` selects the previous rule, including its default beam 100;
all settings can be overridden explicitly. `--gamma` requires `--from-mea`.
Incompatible selected pairs are rejected rather than resolved by silently
choosing one. Malformed/non-finite BPP records and selected noncanonical
pairs fail; the probability reader uses the paper's 1e-5 float tolerance.

A `# linearknot_scaffold_metadata=` JSON comment records the actual mode,
threshold, minimum helix, beam, dangles, gamma, LinearPartition commit and
wrapper/native hashes. LinearKnot's record reader ignores the comment.
The commit is null for wrappers whose revision cannot be established;
standalone files do not inherit the enclosing LinearKnot checkout's commit.
`--metadata FILE` writes the same settings with input/scaffold hashes per
record. `--no-metadata-header` restores the plain three-line record output.

The helper remains a source-checkout script. The LinearPartition wrapper
requires `python-gflags` in its Python environment. No experimental energy
model or accuracy claim is introduced by this change.


## Validation evidence

The controlled-wrapper test covers direct BPP selection, the exact 0.7
boundary, stacked-helix filtering, empty output, conflicting pairs,
noncanonical pairs, batch input, metadata and the legacy MEA route.

With the paper's pinned LinearPartition commit
`b450fb3e63189073b68d385589035f992080aa3a`, fresh CLI runs match all six
selected saved scaffolds: four bpRNA inputs (52, 55, 75 and 132 nt) and two
RNAcentral inputs (4,793 and 5,482 nt). The scaffold data records are byte
identical after omitting the new metadata comment. The 75-nt input loses an
isolated threshold pair; the 55-nt input has no pair at probability >=0.7.
The previous helper at v0.2.0 and the new `--from-mea` mode also give
byte-identical records on two real inputs with the comment disabled.

The saved-cohort comparison is optional and needs the paper's data checkout:

```bash
python3 test/tools/compare_make_scaffold_saved.py \
  --linearpartition /path/to/LinearPartition/linearpartition \
  --output build/scaffold-saved-comparison
```

Use `--python /path/to/venv/bin/python` if the wrapper's `python-gflags`
dependency is installed in a dedicated environment.
