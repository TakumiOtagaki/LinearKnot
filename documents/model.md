# Folding model

## Scaffold-conditioned folding

LinearKnot takes a sequence and a fixed noncrossing scaffold G. It searches
the density-2 hierarchical model for additional pairs G', returning the
completed structure R = G union G'. Scaffold selection is separate from a
single fixed-scaffold call. The CLI takes the supplied scaffold; the optional
C++ hotspot API proposes a bounded set of scaffolds for MFE prediction.

Scaffold pairs remain fixed. An all-dot scaffold is valid and produces
pseudoknot-free folding: crossing extensions require fixed scaffold pairs.
Scores and probabilities are conditional on the chosen scaffold. Changing
the scaffold is not assumed to preserve the energy assigned to a completed
pair map.

## Shared grammar and weights

Exact and finite-beam consumers share the deduction generator and local
weights. Inside uses LogSumExp; minimum-energy prediction uses LogMax.
Posterior calculations, sampling and decoding consume this model.
The beam restricts retained support; it does not define a second energy model.

Local thermodynamic weights use the bundled ViennaRNA PF parameter source
and the CParty pseudoknot parameterization. Production operates directly on
log weights. Free energies are computed from the log values without a
length-dependent numerical scaling factor. Raw weight projections may
underflow or overflow; public C++ callers should use the log fields and
check `ok` and `reachable`.

Completed WMB structures require weak closure relative to G by default.
Open partial band-chain states remain available. The C++
`WmbClosurePolicy::LegacyCompatible` policy selects the historical completed-WMB
domain for controlled comparisons with the external oracle. It is separate
from beam ranking and is not an alternative thermodynamic parameter set.

Long scaffold BE loops remain admissible and use explicit logarithmic
initiation extrapolation beyond the tabulated loop sizes. Local search
windows for newly generated internal loops remain part of the model.

See the [weight-layer contracts](../src/scfg2/weights/README.md) and
[third-party notices](../THIRD_PARTY_NOTICES.md) for parameter provenance.

## Derivation correspondence

The adopted density-2 model uses the completed-region weak-closure condition
from CParty's mathematical decomposition. Under the same boundary and energy
conventions, the shared grammar regroups that decomposition's child products
in a fixed order. Intermediate states retain the original split coordinates
and scaffold anchors; contracting those states recovers the original rule.
Conversely, the fixed order gives one expansion of each original derivation.
Each local energy term occurs once, preserving total energy and Boltzmann
weight before pruning.

The correspondence covers these intermediate combinations:

| Combination | Preserved information and completion condition |
| --- | --- |
| VP three-gap product | Combine the last two WI gaps, then the first; retain the original split tuple and all three factors. |
| BE wings | Combine the inner BE and right WIP, then the left WIP and local factor; retain both scaffold-pair anchors. |
| WMBP band join and WMB completion | Retain the shared BE anchors and parent-boundary condition; require weak closure when the WMB chain completes. Partial WMBP/WMBW chains may remain open. |
| Empty WI gap | Restore the omitted gap as a unit factor, without an extra split choice or unpaired assignment. Empty WIP remains invalid. |

Together with CParty's decomposition theorem, this is the model-level
correspondence used in the paper. It excludes `LegacyCompatible`, which omits
completed-WMB weak closure and can admit duplicate derivations. Beam pruning
restricts the retained support. The [implemented boundary contracts](../src/scfg2/deductions/README.md)
and [finite validation witnesses](validation.md) are separate evidence about
the program; enumeration alone is not a proof for arbitrary inputs.

## Beam and exact computation

The default beam width is 100 with outer ranking. An inside pilot of width
ceil(b^(2/3)) supplies completion-context scores for ranking; exponent 0.5
is also supported. The guide affects selection, not local weights.
`--ranking baseline` selects the baseline priorities.

Finite-beam calls also use local limits for new hairpins and CASE3 branching
scaffold ancestors, both defaulting to 30. A wider beam leaves those limits
in place. `--exact` removes beam pruning and these beam-only limits for the
supplied scaffold. Exact computation still uses the fixed-scaffold model.

See [API ranking options](../src/scfg2/api/README.md#outer-and-baseline-ranking)
and [resource bounds](complexity.md) for the configuration contracts.

## Probabilities and structural contexts

Pair probabilities refer to completed union structures, including fixed
scaffold pairs with probability one. For each base, the unpaired probability
is q(i) = 1 - sum_j p(i,j).

Loop profiles partition each base into Stem, Hairpin, Bulge, Internal,
Multiloop, Exterior and PseudoLoop. Since 0.2.0, PseudoLoop labels unpaired
gaps between crossing bands (direct WI emissions). Multiloop also includes
the wings of multiloops closed by pairs of one band, independently of the
energy parameters used there. Labels apply to directly emitted bases;
a nested child keeps its own context. See the
[geometric definition](pseudoloop_context_change.md). Window accessibility is
the joint probability that the entire requested window is unpaired. The optional
loop channels partition that window probability by context. Opening energy
is -RT log P; a zero probability has infinite opening energy.

## Structure selection

MFE maximizes the model's thermodynamic log weight. Sampling draws according
to Boltzmann probability. Subopt enumerates in nondecreasing thermodynamic
energy within the chosen exact or retained support.

MEA maximizes

```text
2 gamma sum_(i,j in predicted pairs) p(i,j) + sum_(i predicted unpaired) q(i).
```

Gamma-centroid maximizes expected pair-level gain gamma TP + TN. Its
variable pair reward is `(gamma+1) p(i,j) - 1`. Both decoders preserve G
and optimize over grammar-admitted structures. Their posterior objectives
are distinct and are not thermodynamic energies.

## Fixed-target energy and probability

`ExactSession::evaluate_energy()` constrains the shared exact grammar to a
requested generated pair map and selects its maximum-weight derivation.
The boundary removes any fixed scaffold pairs from the target and checks
that traceback returns exactly the requested generated pairs. A mismatch
is returned as unreachable. This uses the same local weights as inside.

`ExactSession::probability()` divides that selected derivation's weight by
the inside partition function. If multiple derivations emit the same pair
map, this numerator is their maximum weight, rather than their summed mass.
The API therefore computes a best-derivation energy. Under the adopted
model's one-to-one correspondence, a reachable target has a single derivation
and the numerator is its structure weight. This interpretation does not
extend to `LegacyCompatible`: the constrained output-grouping witness includes
a compatibility-mode target with two derivations. The API does not sum those
derivations.

Successful reachable evaluation is sound for the requested pair map. The
local constraint predicates and final check have regression coverage, but
no general target-completeness proof: an unreachable result is not a
certificate that no grammar derivation exists. The model-level correspondence
does not assert scaffold-split invariance or cover different boundary and
energy conventions. See the [constrained output-grouping witness](../test/scfg2_constrained_small_n_grouping_test.cc)
and [validation guide](validation.md).

For API result fields, error contracts, and streaming callbacks, see the
[C++ API](../src/scfg2/api/README.md).
