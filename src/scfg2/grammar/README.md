# grammar

Own the shared grammar vocabulary:

- nonterminals
- rule ids
- rule specs
- split schema
- child schema
- grammar-level predicates

This layer should not depend on scheduler policy, beam policy, or CLI/runtime
adapters.

`surface_kind.hh` owns `WmbpStepPrefixKind`, used to classify step prefixes in
exact replay observation. It is grammar/storage vocabulary and does not create
a second recurrence or scoring semantics. The former VP/WMB carrier-surface
enums have no remaining consumers and were retired.
