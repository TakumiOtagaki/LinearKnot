# Experimental CONTRAfold / DP09 backend

This optional backend uses fixed CONTRAfold-v2 log scores for ordinary RNA
loops and the retained pseudoknot parameter set named DP09 in CParty
(`rna_DirksPierce09.par`). It reuses LinearKnot's exact deduction generator, schedule and regenerated
outside pass. A scaffold is a fixed, noncrossing set of base pairs that all
conditional structures must contain. Base-pair probability (BPP) means the
marginal probability of a pair in this conditional structure distribution.

This is a research interface, separate from the stable `linearknot` command
and Python API. Enabling the build option does not change their backend.
The release contains original coefficients, not fitted research coefficients,
RNet2 outputs or neural-network weights. It provides **exact conditional**
BPP and parameter derivatives; it does not claim linear running time and
does not include scaffold sampling or averaging over multiple scaffolds.

## Build and check

From the repository root, with CMake, a C/C++17 compiler and Python 3:

```sh
cmake -S . -B build/cf -DLINEARKNOT_BUILD_CONTRAFOLD=ON -DLINEARKNOT_ENABLE_LTO=OFF
cmake --build build/cf -j4 --target linearknot_contrafold_probe linearknot_contrafold_score
ctest --test-dir build/cf -R '^contrafold_exact_backend$' --output-on-failure
```

The option defaults to `OFF`. Both executables are under
`build/cf/experimental/contrafold/`. The test needs only the Python standard
library. It enumerates all allowed noncrossing structures of five short
sequences and checks partition functions and BPP against an independent
structure-tree scorer, then checks conditional scaffold pairs, crossing
pairs, incident probability sums, twelve-coordinate BPP derivatives and
structure-likelihood derivatives by finite differences. These checks are
implementation checks, not a prediction-accuracy benchmark.

## Line-oriented protocol

The probe accepts one request per line and returns one JSON object per line.
Positions are one-based. Sequences use uppercase RNA bases; scaffold lengths
must match. Ordinary generated pairs are AU, UA, CG, GC, GU or UG and have
at least three intervening bases. A minimal invocation is:

```sh
build/cf/experimental/contrafold/linearknot_contrafold_probe <<'REQUESTS'
M metadata
R example GGGGGAACCAACCC ((.....)).....
E example 0 0 0 0 0 0 0 0 0 0 0 0
L example (([[[..))..]]]
D example
REQUESTS
```

`R id sequence scaffold` registers a case; `D id` releases it. `M id` returns
backend and unit metadata. `E id d0 ... d11` evaluates BPP, log partition
function `logZ`, and expected parameter features; `B` evaluates BPP without
the feature means. A parameter feature is `-RT * d(log local weight)/d(displacement)`;
`expected_feature` reports its conditional expectation. Additive energy
coordinates have count-like features, whereas multiplier coordinates carry
energy-equivalent features; these are not all occurrence probabilities.
`RT` is the gas constant times absolute temperature in kcal/mol.
`logZ` excludes scaffold-constant scores; fixed scaffold
pairs still have BPP one. Different scaffolds' `logZ` values must not be used
as their mixing probabilities without accounting for the model definition.

The twelve displacements, in order, are:

| Index | Name | Displacement unit |
| --- | --- | --- |
| 0 | exterior | kcal/mol |
| 1 | in_multiloop | kcal/mol |
| 2 | in_pseudoloop | kcal/mol |
| 3 | band | kcal/mol |
| 4 | unpaired | kcal/mol |
| 5 | nested_region | kcal/mol |
| 6 | multiloop_branch | ignored; ordinary CONTRAfold branch is fixed |
| 7 | band_multiloop_closing | kcal/mol |
| 8 | band_multiloop_branch | kcal/mol |
| 9 | band_multiloop_unpaired | kcal/mol |
| 10 | stack_multiplier | dimensionless additive displacement |
| 11 | internal_multiplier | dimensionless additive displacement |

A multiloop is a loop with multiple helical branches. Ordinary multiloop
branches, including pseudoknot-as-branch boundaries, use the fixed
CONTRAfold `multi_paired = -0.9253883752` log score. Consequently coordinate
6 has zero derivative. The DP09 stack and internal-loop multipliers default
to **0.89 and 0.74**; these are multipliers, not exponents. They can be moved
through coordinates 10 and 11. Scalars in the core parameter object are
stored as `double` to preserve fractional displacements, with unchanged
integer-valued defaults. The adapter currently requires 37 degrees Celsius.
CONTRAfold values are dimensionless log potentials; the hybrid parameters
should not be interpreted as independently measured physical free energies.

After `E` or `B`, `V id version count (i j adjoint)...` returns the gradient
of `sum(adjoint_ij * BPP_ij)` with respect to the twelve displacements. This
is a vector–Jacobian product, a way to differentiate a scalar loss through
the pair probabilities. Use the `version` from that case's latest evaluation.
Fixed scaffold pair adjoints contribute zero gradient. `L id structure`,
after `E`, returns the full-structure log likelihood and its gradient when
the structure belongs to the grammar and contains the scaffold; unreachable
structures have `reachable: false` and `log_probability: null`. Its gradient
field is not meaningful in that unreachable case. Dot-bracket notation uses
matching bracket types to represent crossing pairs. `T id index step`
evaluates an independently displaced snapshot's `logZ` for diagnostics.
Errors are reported as JSON objects containing `error`; clients must check
for them rather than treating them as valid numeric results.

The companion score executable reads `sequence noncrossing_structure`
records and emits one ordinary CONTRAfold log score per line. It is an
oracle for valid canonical noncrossing structures, not a general input
validator or a pseudoknot scorer.

## Origin and boundaries

`provenance.json` identifies copied research adapters and their source
commits. `third_party/linearsampling/{utility.h,feature_weight.h}` are
unmodified files from LinearSampling commit
`1492a93f3616bf8ae8ab9177f9bb85a53da4b155`. Their full upstream license is
retained in that directory; the repository's GPL text does not relicense
these files. See the root `THIRD_PARTY_NOTICES.md`.

The parameter-feature view was adapted from the earlier LinearKnot
parameter-learning experiments. Pair terms are added exactly once at the
generated-pair emission boundary. Fixed scaffold pair scores cancel within
a conditional distribution and are omitted from numerator and denominator
together. This distribution is conditional on the supplied scaffold, not
a complete unrestricted pseudoknot ensemble. The default weak-closure and
legacy span-support restrictions are retained. No other research model,
dataset, decoder, scaffold-retention policy or learned correction is part
of this backend snapshot.
