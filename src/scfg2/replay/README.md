# replay

Own derived consumers of the shared kernel:

- selected-decision replay
- score replay
- structure reconstruction helpers
- exact and beam execution adapters using the shared grammar and local weights
- sampling and ranked-derivation consumers

`w_final_exact_inside.hh` provides the snapshot-backed raw/log weight context
and the exact schedule adapter. The context reads public `RuntimeSnapshot`
accessors directly and queries its own immutable scaffold geometry. The
geometry owns the unpaired-run counts; the context keeps no duplicate array
or long-lived reference to those counts.

Execution adapters:

- build or receive a shared `scfg2` schedule
- run the shared semiring/local-weight kernel
- expose its chart, decisions or derivation graph to the requested consumer

Owner-local diagnostic projections live only in `test/` observation helpers;
production snapshots own none. The raw observation context and the separate
`storage::exact_store` remain for independent witnesses. Replay and exported-
structure consumers derive their behavior from the shared grammar; they must
not introduce another folding semantics or restore the retired `src/scfg/`
runtime.
