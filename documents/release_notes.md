# Release notes

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
