# engine

Own the generic DP executor:

- chart storage
- exact fill order / agenda policy
- state updates via a supplied semiring
- optional backpointer recording
- deduction-tree traceback from recorded choices
- outside propagation over regenerated admitted deductions (`outside.hh`)

This layer should not own rule-specific semantics.

`CompactChoiceMap` stores selected rule/split descriptors and decodes them
through the common deduction expansion and runtime child mapping. A full
deduction is retained only when that descriptor cannot reproduce every field,
including the inactive child slots used by the historical tie key. Reads
return owning values without a decoded-record cache. Public beam MFE uses
this store; the original `ChoiceMap` remains available for detailed traces
and reference comparisons. Both stores feed the same forward executor and
the same chosen-pair output consumer.

Regenerated outside reads immutable inside values in reverse topological order.
The caller supplies the same branch admission, local weights and boundary child
mapping as inside. Each child occurrence receives the product of its siblings,
without division; zero values and repeated child keys are supported. This helper
does not recover old online pre-prune snapshots from a final pruned chart, and
is used by the production beam's retained-state program in
`beam/retained_outside.hh`.

Complete-output probabilities use `GeneratedTargetIndex::matches_local_emissions`:
generated pairs must belong to the requested output, and requested endpoints
cannot be emitted unpaired. This reuses the shared pair and unpaired emission
observers. Together with the complete derivation's per-base partition and an
initial scaffold-compatibility check, it enforces the full union output.
`exactly_partitions` retains its narrower interval-ownership contract for
observations; it is not a general output-conditioning predicate. In particular,
a VP interval can span an ordinary scaffold interior whose generated pairs are
consumed by another branch. The 96-nt regression in the composed-beam test
checks that target mass is finite and that denominator/maximum weights remain
unchanged in both unlimited and finite beam runs.
