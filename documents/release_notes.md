# Release notes


## 0.2.1

The source-checkout scaffold helper now uses the paper's construction rule:
all LinearPartition-V BPP pairs with probability at least 0.7, keeping helices
of at least two stacked pairs, with beam 200 and dangles 2. It no longer
starts from an MEA structure by default. `--from-mea` retains the previous
mode and its defaults (gamma 3, threshold 0.5, minimum helix 1, beam 100).

Output comments record the settings and LinearPartition provenance;
`--metadata FILE` also writes a JSON sidecar. LinearKnot accepts the comments.
The inference model and supplied-scaffold CLI behavior are unchanged. The
paper's v0.2.0 computations already used this construction rule externally.

## 0.2.0

PseudoLoop now labels unpaired gaps between crossing bands. Unpaired bases
in band multiloop wings (WIP, VPL, VPR and BE wing emissions) move to
Multiloop. This definition applies to exact and beam profiles and to
`access --by-loop`. A nested child retains its own loop context.

The seven context names and their order are preserved. Energy parameters,
partition functions, pair probabilities, total window accessibility and
structure selection are unchanged. Floating-point regrouping can affect
the final digits of context totals. Context results from 0.1.1 must be
recomputed before comparing them with 0.2.0.

Incremental CMake builds now detect changes to VERSION and refresh the CLI
version and package metadata. Example runners reject an executable whose
version differs from the checkout.

See the [definition and validation](pseudoloop_context_change.md).
