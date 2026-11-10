# Release notes

## 0.2.3

An optional CONTRAfold-v2 / DP09 research backend now provides exact
scaffold-conditioned base-pair probabilities (BPP), parameter-feature expectations,
BPP vector–Jacobian products and whole-structure log likelihoods. Enable it
with `-DLINEARKNOT_BUILD_CONTRAFOLD=ON`; see the
[build and protocol guide](../experimental/contrafold/README.md).

This interface uses original coefficients and a supplied noncrossing scaffold.
It does not include learned coefficients, neural-network weights, scaffold
sampling or averaging over scaffolds. The stable `linearknot` CLI and Python
API keep their existing backend; the optional probe is a separate executable.
DP09 scalar storage now accepts fractional values; its defaults are unchanged.
The ordinary CONTRAfold multiloop branch is fixed, including pseudoknot branch
boundaries. The default band multipliers remain 0.89 and 0.74.

The snapshot includes independent provenance review and implementation checks
against short complete noncrossing enumerations, finite-difference gradients
and the research adapter. Those checks do not measure prediction accuracy or
establish linear running time. The original adapter snapshot is also preserved
at tag `contrafold-v1`; version and citation metadata are updated for this
release. This backend is currently available from a repository checkout and is
excluded from the explicitly packaged Python source distribution.

## 0.2.2

Whole-structure probability evaluation now handles generated base pairs inside
a scaffold component crossed by a VP state. The former interval-partition
check could reject a valid output and report zero probability because it
assigned an independently generated interior pair to the crossing VP branch.
The corrected evaluator checks local paired and unpaired emissions against
the requested complete output, with an initial scaffold-compatibility check.

This fixes LinearKnot's target-probability consumer; energy parameters,
partition functions, base-pair probabilities and MFE selection are unchanged.
Regression tests cover the local ownership boundary and a 96-nt complete
output in both unlimited and finite beam runs.


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
