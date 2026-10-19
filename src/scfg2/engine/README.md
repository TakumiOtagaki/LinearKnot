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
