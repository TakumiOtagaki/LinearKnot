# PseudoLoop: gaps between crossing bands (0.2.0)

Status: implemented in 0.2.0 (2026-10-23). The specification and preflight
evidence below describe the change from v0.1.1 (`4f97fd3`), with documentation
baseline `a8ace24`. The implementation uses one fixed-rule context map in
`loop_context.hh` for both direct emissions and beam aggregates.
Weights, grammar support and the seven output labels are preserved.

## Problem

The loop profile assigns each unpaired base one of six contexts. The v0.1.1
rule was "PseudoLoop = unpaired bases charged by a pseudoknot penalty"
(`src/scfg2/profiles/loop_context.hh`, comment above the PseudoLoop cases).
That covers two different kinds of position:

1. **Gaps between crossing bands** (WI regions, penalty PUP per base).
   These loops exist only because two bands cross. Example: the three loops
   of an H-type pseudoknot.
2. **Wings of a multiloop closed by two pairs of one band** (WIP, VPL, VPR
   and the unpaired BE wing, penalty `c_p` per base). Geometrically this is
   an ordinary multiloop: an outer and an inner pair of the same band close
   it, and branches sit in its wings. Pairs crossing the enclosing band do
   not end in the wings. A closed child branch can contain its own pseudoknot.

Type 2 is labelled PseudoLoop only because the energy model uses the
pseudoknot multiloop parameters (`ap`, `bp`, `c_p`) there. Readers expect
"PseudoLoop" to mean type 1, so the old label was misleading.

```
 (   ..  (...)  ..  (   ...[[[...   )   ....  )   ...]]]
 i   a   branch  a  k                l   a     j
```

`(i,j)` and `(k,l)` are two pairs of one scaffold band; `[[[ ]]]` crosses
both. The bases marked `a` lie in a three-way multiloop closed by `(i,j)`,
`(k,l)` and the branch. They were labelled PseudoLoop in v0.1.1. The bases in
`...[[[` and `...]]]` that lie between the two bands are type 1.

## Definition since 0.2.0

- **PseudoLoop**: unpaired bases in a gap between crossing bands, that is,
  bases directly emitted by WI rules. This does not include every unpaired
  base inside the interval of a WI state.
- **Multiloop**: unpaired bases of ordinary multiloops, as now, and also
  unpaired bases in the wings of a multiloop closed by two pairs of one band
  (scaffold band via BE, or band of crossing pairs via VP), even though their
  energy uses `c_p`.

Hairpin, Bulge and Internal are unchanged; they already keep their labels
inside a pseudoknot. The number of contexts stays seven, and the output
labels and column order stay the same.

## Geometry reference independent of the grammar

Classification takes the complete union pair map R, not rule IDs, energy
parameters or the scaffold/generated origin of a pair. It labels a completed
structure; posterior probabilities still depend on the chosen scaffold and
its weights. Equality of per-structure labels for the same R does not imply
scaffold-independent energies or probabilities.

An interval is **closed in R** if each pair has either both endpoints inside
or both outside. This is a condition on the complete structure, distinct
from the completed-WMB admission policy's weak closure relative to G.
The smallest closed block beginning at an opening endpoint is obtained by
extending its right boundary until all partners of included bases are inside;
it fails if any partner lies before its left boundary. The leftward operation
is symmetric. A block skipped during a scan must lie wholly on that side of
the query base; an enclosing block containing the query must never be skipped.

For each unpaired base x:

1. Scan left and right, ignoring unpaired bases and skipping the smallest
   closed child blocks. Stop at the first paired endpoint that cannot be
   skipped in that direction.
2. If both scans reach sequence ends, assign Exterior. A one-sided end with
   an exposed endpoint on the other side is an error in this classifier,
   rather than an arbitrary Exterior/PseudoLoop fallback.
3. If both stops belong to the same pair, classify its loop after identifying
   its closed child blocks. With no child, it is Hairpin. With one ordinary
   stem child, it is Bulge or Internal according to whether one or both
   sides have unpaired bases. Multiple stem attachments give Multi.
4. A closed child with crossing bands exposed at its exterior supplies
   multiple stem attachments, even if it occupies one closed interval.
   Therefore a loop around that child is Multi. A noncrossing enclosing
   pair shields a deeper pseudoknot and supplies one attachment instead.
   Counting every closed block as one stem would incorrectly relabel existing
   Multi positions as Bulge/Internal; those changes are outside this proposal.
5. Otherwise, let P and Q be the pairs at the two stops. They are in the same
   band here when one strictly encloses the other and their sets of crossing
   partners are equal and **nonempty**. Equal empty sets do not define a band.
   Examine the two wings between these pairs: a closed child branch in either
   wing gives Multi; without a branch, one or two nonempty unpaired wings
   give Bulge or Internal. Wing contents must be dots or complete closed
   children; an escaping child is a diagnostic, not something to skip.
6. Other exposed-boundary configurations are PseudoLoop.

The [test reference](../test/tools/pseudoloop_context_reference.py) implements
this definition without importing the deduction generator or weight layer.
It is not a folding/admissibility oracle, and arbitrary pair maps outside
the adopted model are not thereby certified as valid model structures.

## Local assignment and nested structures

WIP is directly consumed by BE/VP wings, the VPR right wing, and WIP's own
continuations. VPL and VPR are consumed by the corresponding VP wing rules.
However, there are recursive paths such as
`WI -> WMB -> BE -> WIP`: a completed child pseudoknot can occur inside a WI
gap. Its own loop labels are determined locally, not inherited from WI.
Likewise, a V or WMB branch inside WIP retains its own Hairpin, Bulge,
Internal, Multi or PseudoLoop labels. Only the unpaired runs emitted by the
listed rules change. Whole state intervals must not be recolored.

## Rule mapping

| Rule | Assigned positions | Current | Proposed |
| --- | --- | --- | --- |
| `WI_BASE_SINGLE`, `WI_EXTEND_UNPAIRED` | `j` | PseudoLoop | PseudoLoop |
| `WIP_EXTEND_UNPAIRED` | `j` | PseudoLoop | **Multi** |
| `WIP_BASEPAIR_V`, `WIP_BASEPAIR_WMB` | `[i, k-1]` | PseudoLoop | **Multi** |
| `VPL_SPLIT_VP` | `[i, k-1]` | PseudoLoop | **Multi** |
| `VPR_SPLIT_VP_BASEPAIR` | `[k+1, j]` | PseudoLoop | **Multi** |
| `BE_WIP_BASEPAIR` | `[l+1, j-1]` | PseudoLoop | **Multi** |
| `BE_BASEPAIR_WIP` | `[i+1, k-1]` | PseudoLoop | **Multi** |

The six changing table rows contain **seven rule IDs**, because
`WIP_BASEPAIR_V` and `WIP_BASEPAIR_WMB` share a row.

## Code locations

All three places must agree; the exact and finite-beam profiles are compared
in validation.

1. `src/scfg2/profiles/loop_context.hh`: `fixed_rule_context()` is the single
   source for fixed labels. `unpaired_emission()` retains each rule's run
   boundaries and chain behavior while obtaining its context from that map.
2. `src/scfg2/profiles/beam_profiles.hh`: `aggregated_run_context()` checks
   that the rule is a supported aggregate, then uses the shared mapping.
   Both analytic and taped VPR right-run chains use the same map. Leftward
   start-prefix chains retain their original geometry and ownership.
3. Comments and descriptions: `src/scfg2/api/profiles.hh` (comment at the
   context names), `src/scfg2/cli/profile.cc` header comment,
   `documents/model.md` (Probabilities and structural contexts),
   `README.md` (profile section), `examples/README.md`.

Use one shared mapping for rules with a fixed unpaired context, including
the beam shortcuts. The internal-loop choice between Bulge and Internal
still depends on geometry. The two VPR-chain consumers should obtain their
context from the same fixed-rule mapping. An unknown aggregate rule must be
rejected, not silently classified as either PseudoLoop or Multi.

## Window accessibility and the non-mixing condition

Define the by-loop value as the joint probability that every base of a
window is unpaired and belongs to the specified loop context. Summing the
six values must recover total joint unpaired accessibility.

The geometry reference has a useful structural property. All bases in a
maximal consecutive run of unpaired bases have the same exposed boundaries:
moving within the run only passes dots, so the scans skip the same paired
blocks and reach the same stops. The child branches and band relation are
also unchanged. Thus the entire run has one label, and every fully unpaired
window, being a subset of such a run, has one label. Changing between a WI
gap and a wing context requires passing a paired endpoint.

This establishes non-mixing for the specified geometry classifier. Applying
it to production additionally requires that the rule-emitted labels agree
with that classifier. A fully unpaired window containing different proposed
labels is a **blocking counterexample**: do not silently drop its mass,
normalize the six values, or weaken the total-accessibility definition.
Conversely, changing all labels to Multi because one base is in a wing is
not a valid repair. Existing run/chain ownership must continue to recover
the joint event, not products of per-base marginals.

## What must not change

- Partition function, pair probabilities, unpaired probabilities, MFE,
  samples and decoders. Only the context labels move.
- At every base, Multi(new) + PseudoLoop(new) = Multi(old) + PseudoLoop(old),
  and the other five context probabilities are mathematically unchanged.
  Multi can only gain mass and PseudoLoop can only lose mass. The seven still
  sum to one.
- For by-loop accessibility, the six components still sum to total window
  accessibility. Once non-mixing is established, the Multi + PseudoLoop
  invariant and gain/loss inequalities also hold for each window; the other
  four unpaired loop channels are unchanged.

These are mathematical invariants, not promises of bit-identical probability
arrays. In particular, `ProfileAccumulator::finish()` reconstructs Stem as
one minus the sum of the six unpaired channels. Redistributing their mass
can change the last floating-point bits. Use the existing absolute 1e-9
tolerance on the small probability witnesses, record observed errors, and
investigate larger discrepancies rather than silently increasing tolerance.
Keep discrete structures, sample counts for a fixed seed, and scalar DP
results as separate regression checks with identical settings/build flags.

## Tests and validation

- Freeze v0.1.1 outputs before changing production, with explicit closure,
  ranking, beam and geometry settings and binary/parameter hashes.
- Check the geometry reference against the
  [14 explicit fixtures](../test/fixtures/pseudoloop_context_change.json):
  ordinary loops, an H-type pseudoknot, exposed/shielded nested pseudoknots,
  BE/VP wings on both sides, a WMB branch in a wing, and a wing inside a WI gap.
- Compare proposed rule-emitted labels against the independent pair-map
  classifier on every completed structure/context group in enumeration.
  The existing profile oracle uses `unpaired_emission()` to assign labels,
  so passing that oracle alone does not validate this geometric definition.
  Preserve each discrepancy as a counterexample; do not tune the reference
  to agree without a structural explanation and a regression fixture.
- Exact versus KeepAll beam profile agreement, with beam-only geometry limits
  disabled, including the aggregated-run and VPR-chain paths. A large beam
  with the CLI's default caps is not necessarily the exact candidate set.
- For finite beams, compare analytic retained-outside against the taped
  reference on the same retained calculation. Exercise both rankings and
  both closure policies. Independent finite beams need not equal exact.
- Check the old/new per-base and per-window invariants above on the same
  fixtures, including mirrored chains, cuts, empty runs and boundary windows.
- `test/tools/test_linearknot_access.py`: labels are unchanged; update any
  expected values.
- Regenerate `examples/terc_profile` and `examples/hairpin_exterior` outputs
  and figures.

### Concrete BE and VP witness

One 32-nt target has pairs
`(1,27), (3,7), (10,23), (13,32)`. Use
`G={(1,27),(3,7),(10,23)}` for BE wings or
`G={(3,7),(13,32)}` for VP wings; both targets were found by exact enumeration.
The same complete structure has the same proposed labels under either
scaffold: positions 2, 8–9 and 24–26 are Multi, 4–6 are Hairpin, and
11–12, 14–22 and 28–31 are PseudoLoop. Mirror the fixture to cover the other
wing. The fixture file supplies canonical sequences and complete labels.

The nested-WI witness is 45 nt: a complete 32-nt wing example sits inside
one gap of the 14-nt H-type example. Position 7 is Multi in the child wing,
while position 38 remains PseudoLoop in the surrounding WI gap. This checks
the difference between direct emission and recursive containment.

### Preflight evidence, 2026-10-23

Run the geometry-only fixtures with:

```bash
python3 -B test/tools/pseudoloop_context_reference.py
```

An observation-only extension of the existing exact profile enumerator at
`a8ace24` recorded the current labels, the proposed seven-rule relabeling,
and complete pair maps. It did not change production headers, weights or
the DP. For each closure policy, 27 cases on 26 distinct sequence/scaffold
inputs produced 461 structure/context groups, covering all seven changed
rule IDs. All 10,565 base-instance labels agreed with the geometry reference;
all 19,987 fully unpaired window instances had one proposed label. The six
new BE/VP/nested target structures were explicitly found in those groups.
Both policies had the same label-group counts; this does not imply equal
partition functions or unambiguity in compatibility mode.

The 14 standalone geometry fixtures also passed, including 545 fully
unpaired windows. These are finite preflight checks, not an arbitrary-input
proof that every admitted production derivation matches the reference.
The maintained production regression is now:

```bash
python3 -B test/tools/test_pseudoloop_geometry.py build
```

It compares actual emitted labels under both closure policies with the
independent reference, requires all seven changed rules to occur in completed
derivations, and requires all six concrete target structures to be reached.
The beam profile regression covers both closure policies and rankings at
widths 2, 5 and 20, plus KeepAll without geometry caps, comparing the analytic
and taped reverses and comparing KeepAll with exact.

### Implementation comparison

A frozen pre-change executable and numeric observer were compared with the
new implementation on all six concrete inputs. The numeric comparison uses
108 configurations: both closure policies, exact, and both beam rankings at
widths 2, 5 and 20 plus unlimited KeepAll. Log partition weights and total
window accessibility were bit-identical. Multi + PseudoLoop mass was
preserved at every base and requested window; Multi only gained mass and
PseudoLoop only lost it (absolute tolerance 1e-9). The other five per-base
contexts differed by at most 2.22e-16 from floating-point regrouping.

For the six inputs, exact and beam-5 CLI runs under both rankings gave
162 unchanged comparisons covering MFE, fixed-seed samples and unique
counts, subopt, MEA, centroid, BPP, total accessibility and opening energy.
Output comparisons omit elapsed-time metadata. These finite checks
complement the independent geometry test and the maintained acceptance suite.

## Expected effect on reported results

- The archived paper exact-sample file, SHA-256
  `41ecd697f944ceae1aca785e7c6da927d92ffdc98ebf74dc5e3c0a325987c2b7`,
  contains 1,000 draws. Reclassification of window 123–126 found 871 fully
  unpaired draws: 848 PseudoLoop, 18 Hairpin, 3 Bulge and 2 Multi. The earlier
  count of 850 PseudoLoop was not reproduced from this file. The two Multi
  cases have ordinary enclosing pairs with multiple child stems; they are
  not WI-gap instances. The sampled joint PseudoLoop contribution is
  848/1,000, not the conditional fraction 848/871.
- Recomputing the saved Terc beam-200 example with identical settings gives
  PseudoLoop probability 0.8397007149209372 at base 123, compared with
  0.8397009254807158 before the change. Multi increases from
  0.001058880685586579 to 0.0010590912453651063. The README's approximate
  84% statement still holds; the raw profile, manifest and figures have
  been regenerated. These are beam marginals, distinct from the archived
  exact-sample frequencies above.
- Loop-context comparisons in the paper (Figure 1C context track, context
  TV in Figure 2B and Supplementary tables, by-loop validation table,
  Terc profiles) must be recomputed. The intended change preserves DP
  support, weights and asymptotic resource bounds; identical wall times or
  measured memory peaks are not guaranteed by relabeling alone.

## Version

This changes the meaning of an output column. VERSION is 0.2.0; see
[release notes](release_notes.md). The number, names and order of columns
are unchanged. No new energy model or calculation mode is introduced.
