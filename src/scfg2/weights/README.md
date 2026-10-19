# weights

Own rule-local transition weights for the shared kernel:

- family-scoped context interfaces
- rule-local weight helpers keyed by shared `RuleId`
- thin adapters to current runtime objects only at the boundary

This layer should not own deduction generation, fill order, or broad runtime
state.

## Source and model provenance

`log_loops.hh` is a log-space translation of the bundled ViennaRNA PF
hairpin/internal/stem helpers; its location under SCFG2 does not make that
implementation independent of ViennaRNA. The shared raw/log parameter-table
builder remains in `src/ViennaRNA/params/params.c`. The runtime snapshot and
pair-table wrapper call that backend through the documented local-weight
boundary; they do not supply a second thermodynamic model.

`pseudoknot_parameters.hh` retains the CParty/DP09 model values and band
multipliers. The six top-level `params/*.par` files are inherited data, with
their pinned source and hashes recorded in `THIRD_PARTY_NOTICES.md`. Source
origin, model data and permission to relicense are separate questions. This
cleanup keeps the existing license grants and scientific attribution.

The VP/WMB condition table lives in `deductions/README.md`. Its BE-availability
restriction belongs to deduction support, while the unchanged `PB` factors
belong to this layer; organizing that distinction adds no energy correction.

## Ordinary multiloop closing boundary (2026-10-17)

The closing stem reverses pair orientation to `(j,i)` and uses the two
**inner** neighbors `s[j-1]`, `s[i+1]` in the existing dangles=1/2 branch.
They exist even when `i==1` or `j==n`; exterior endpoint tests do not apply.
Branch/exterior stems still use outer neighbors `s[i-1]`, `s[j+1]` with
their original endpoint guards. Dangles=0 still supplies neither neighbor.
The shared raw/log contexts expose the same corrected factor to exact,
beam, outer and reverse consumers. Do not add a separate scoring kernel.

The old guard was inherited from CParty PF `exp_Mbloop`; oracle patch `0009`
corrects that helper too. `LegacyCompatible` selects the historical WMB
closure domain, not the retired multiloop boundary error. PF smoothing and
LogMax optimization of the PF-derived score remain unchanged; this is not
a migration to Vienna's integer MFE energy. See
`documents/eval_structure/scfg2_multiloop_closing_boundary_20261017.md`.

## Runtime snapshot boundary (2026-10-16)

`RuntimeSnapshot` replaces the inherited `W_final_pf` owner. Its constructor
takes sequence, scaffold and `RuntimeOptions`; it owns log parameters only,
without materialized DP tables. It owns normalized/encoded sequence, Vienna
parameter tables, long-loop extrapolation and pseudoknot model snapshots.
`pseudoknot_parameters.hh` stores the unchanged model data and its serial
configuration boundary. Log PK factors are snapshot-owned;
creating another owner cannot overwrite an earlier owner's factors.
The 0.89/0.74 band multipliers are also copied with that model.

Owner-local diagnostic projections are test support only (2026-10-19).
`RuntimeSnapshot` owns no diagnostic object, prefix, offsets or span vectors;
`RuntimeOptions` selects dangles only. The snapshot/storage and poisoned-carrier
witnesses explicitly own `testing::ObservedRuntimeSnapshot` and its unchanged
29-table observations under `test/`. They are not installed or visible through
the production context. Keep the separate `storage::exact_store`: independent
raw-carrier witnesses still read it, and raw traceback wrappers still reset it.
Public exact/beam charts, tapes, forests and the shared grammar/executor/context
are unchanged. No old folding or sampling implementation is carried by the
snapshot. See `documents/eval_structure/scfg2_diagnostic_storage_cleanup_20261019.md`.

The replay-side weight context reads the snapshot's public scaffold/weight
accessors directly; the former `WFinalAccess` forwarding header is retired.
Unpaired-run counts stay in its immutable scaffold geometry, with local const
views instead of a second copied array. Tests and tools use the same direct
snapshot constructor, preserving their raw/log contexts and semirings. The
retired `ReferenceRuntimeSnapshot` adapter and its ignored constructor flags
are removed; only the two explicit observation witnesses own diagnostic tables.
The independent internal-loop witness projects log factors locally. Retired
`part_func_*`, `sparse_tree` and penalty-global sources must not be restored. See
`documents/eval_structure/scfg2_runtime_boundary_replacement_20261016.md`.

Numerical-scaling removal (2026-10-16): the snapshot is log-only.
`ParameterDomain`, raw PF/PK copies, MFE seeds, scale/length arrays, environment
scaling controls and free-energy scale corrections are removed. Bounded raw
observation views exponentiate the same log factors on demand; they do not own
a second parameter model. Independent raw Vienna comparisons remain in tests.
Physical length penalties and the 0.89/0.74 BE coefficients are unchanged.
The historical inactive `VM_SCALE2` rule ID remains for rule-inventory stability;
its runtime weight is the multiplicative identity, with no rescaling.

Retired WMBP right-open observation bundles and VP/WMB profile accessors no
longer belong to the local-weight interface. They had no current callers and
kept diagnostic-store references in production virtual tables. Rule-local
weights and necessary raw witness projections are unchanged.

A follow-up removes the uncalled `vp_case3_surface_total`,
`vp_case12_surface_total` and `wmb_surface_total` hooks as well. Local-weight
contexts do not read materialized diagnostic totals. Witnesses that require
raw carrier projections continue to use the independent `storage::exact_store`;
owner-local materialized observations are confined to test support.

The active grammar has 15 nonterminals after removal of root-dead `VP_DIRECT`
and `VP_CLOSED`. `WMBP_DIRECT_VP` directly names its canonical `VP` child; no
child identity projection is needed. Local weights and physical pair ownership
are unchanged. Historical 17-state notes describe the earlier baseline.

## Runtime input and parameter boundary (2026-10-11)

`runtime_parameters.hh` honors `CPARTY_FIXED_ENERGY_PARAM_FILE` before the
owner snapshots PF tables or long-loop extrapolation values. Explicit empty,
unreadable or unloadable selections fail. The bundled loader now also rejects
a missing v2.0 header or absence of any recognized parameter section, before
mutating tables or the loaded model name. Valid partial files remain overlays
on the current globals; a complete file is needed to select a model independent
of prior loads. Identical values are valid and are not rejected by a table-hash
comparison. This is envelope validation, not a replacement for Vienna's
section-body parser or a general malformed-body validator.
No pathname cache is independent of Vienna's mutable globals.
Without an override the current global model
is used; a fresh process has the bundled Turner2004 defaults. Configuration
must stay fixed for the lifetime of public sessions, which may cache inside.
Calls remain serial.

`sequence.hh` provides uppercase RNA spelling (`T` to `U`) for the owner's
encoded bases and sequence-dependent special-loop factors. Normalization
does not change grammar admission, local-factor ownership or the caller's
stored input. The owner releases both encoded sequences and its log PF
parameter snapshot through RAII, including constructor-unwinding paths.

## Completed-WMB domain policy (2026-10-10)

`WmbClosurePolicy::RequireWeakClosure` is the default for every consumer.
The shared deduction generator enforces completed-WMB weak closure; the weight
layer does not zero factors to emulate admission. `LegacyCompatible` restores
the previous WMB domain explicitly. Open WMBP/WMBW prefixes and all rule-local
factors, including Psm ownership, are unchanged. Do not impose closure on a
partial carrier or compensate for removed paths with profile emissions.

The official legacy-inside oracle now pins upstream `800e95a` plus local
0001/0002/0006/0007, without closure. Compare its log Z / PF free energy with
explicit `LegacyCompatible`; this oracle choice does not change the default
domain or any local-weight ownership. See `test/oracle/legacy_cparty/README.md`.

## Production value domain (2026-10-10)

All public exact/beam consumers use `WFinalLogWeightContext` and tagged
`LogWeight` values. `AllLocalWeightModel::local_weight_for` is the kernel
entry; the tag prevents an accidental second logarithm. Impossible factors
are `-infinity`, neutral factors are zero, and products add their log factors.
`weight_for()` raw projections and `WFinalExactInsideContext`
remain bounded comparison/observation interfaces only.

The bundled Vienna PF parameter builder has one shared raw/log table builder.
Its log entry evaluates the exponent directly, preserving temperature,
`pf_smooth` truncation, smoothed dangles/mismatches and special-loop tables.
`log_loops.hh` follows PF helper branches; integer MFE helpers are not the
weight source. Vienna's forbidden short-hairpin `INF` sentinel remains impossible, whereas a finite
energy whose exponential underflows still has a finite log weight.

The runtime owner snapshots these log tables and PK parameters. Its
`exp*` fields contain log weights and must never be passed to raw `exp_E_*`
helpers. Access local factors through the context. There are no `scale`,
`expMLbase`, `expcp_pen` or `expPUP_pen` length arrays. Length factors are length times
the corresponding log constant; numerical PF scaling is absent. Model PK
penalties and the 0.89/0.74 BE coefficients remain part of the score.

Forward, retained outside, profiles and forest consumers share this boundary.
See [implementation and evidence](../../../documents/model.md#shared-grammar-and-weights).

Thermodynamic k-best traceback observes an explicit `local_log_weight` on
factorized terms. It is accumulated from these same injected factors during
forward multiplication; it is not obtained by subtracting child inside
values from a large term mass. Child/wrapper edges carry zero local weight.
This observation adds no scoring model and changes no ownership or admission
rule. Inside/outside continue using the existing forward mass fields.

Long internal-loop boundary (2026-10-10):

- `internal_loop.hh` preserves the raw reference weight evaluation. Sizes
  <= `MAXLOOP` delegate unchanged to the existing Vienna helper. Larger
  loops use explicit logarithmic initiation from size 30, the appropriate
  generic/1xn mismatch factors or bulge terminal-AU factors, and capped
  NINIO asymmetry without indexing a size-limited array.
- BE follows fixed scaffold pairs and remains admissible for sizes >30.
  V/VP search windows and all deduction guards remain unchanged.
- The runtime boundary snapshots the unrounded temperature-scaled size-30
  initiation energies, NINIO slope and cap with each owner's PF tables;
  later global parameter loads do not change that owner's long-loop weights.
- `pf_smooth=0` truncates the whole initiation energy and the NINIO slope
  as Vienna's expanded PF parameter tables do. `pf_smooth=1` preserves the
  continuous energies. BE's existing pseudoknot multiplier stays
  outside this raw loop helper. Production `log_loops.hh` expresses the
  same loop/PK ownership in log space; forward and regenerated outside share it.
- This intentional correction replaces inherited undefined behavior.
  Patched legacy oracle `0007` provides an independent implementation;
  the author reported it upstream on 2026-10-10 (issue number not recorded).

Regenerated outside reads the same local weight as forward. Existing runtime
child identity resolution is shared in `engine/runtime_child_item.hh`: the
`WMBP_DIRECT_VP` grammar carrier `VP_DIRECT` reads the public `VP` alpha.
Forward, retention checks, and reverse propagation use that same mapping;
this is an observation of the existing child boundary, not a new weight or
permission for additional same-span/nested mass.

MEA and Hamada gamma-centroid use the same thermodynamic local weights to
build their posterior/support. Their additive posterior rewards live in
`scfg2/decoding/pair_posterior.hh`, rather than a second thermal weight model.
Beam decoding couples unchanged alpha with a selected posterior max trace
through the existing deduction-value transform; exact keeps its frozen forest.
Every pruning/representative decision uses alpha, including private BE/WMBP
and prefix paths. Neither objective may alter inside ranking, support or
cached forward values. Band/lane carriers affect geometry;
both objectives score physical pair marginals without separate band gains.

Exact WMBP step-prefix observation uses `WmbpStepPrefixKind` from
`grammar/surface_kind.hh`. It classifies storage geometry without adding a
scoring semantics. The unused VP/WMB carrier enums are retired; family weight
interfaces must not include the retired `scfg/rules_api.hh`.

Initial scaffolded families:

- `w`
- `v`
- `wi`
- `vm`
- `wmv_wmp`
- `wm`
- `wip`
- `vpl`
- `vpr`
- `vp`
- `wmbw`
- `wmbp`
- `wmb`
- `be`

Status as of 2026-07-01:

- `v`, `vm`, `wm`, `wmv_wmp`
  - local-weight wiring is live in the current pk-free exact surface
- `w`, `wi`, `wip`, `vpl`, `vpr`, `vp`, `wmbw`, `wmbp`, `wmb`, `be`
  - family-local context interfaces and transition helpers now mirror the
    legacy restricted recurrence factors
  - replay-side `W_final_pf` exact root execution now runs through
    `AllLocalWeightModel<WFinalExactInsideContext>`
  - the weight layer is no longer just a `0.0` placeholder or test-only
    scaffold

Operational rule:

- witness discovery is no longer the main blocker for this layer
- implement transition/local weights monotonically from legacy restricted
  recurrence into `scfg2/weights/`
- only go back to witness mining when a concrete family needs a new targeted
  acceptance/debug probe

Geometry/guard rule:

- `square` is only a crossing-geometry carrier relative to the round scaffold
  `G`; it is not a separate scoring semantics
- keep one shared weight/kernel surface and add crossing-aware helpers only for
  families that need them
- when current `scfg2` deduction guards still look "too legacy", prefer that
  temporary duplication over prematurely factoring away parent-side
  admissibility checks before inside-alignment is secure
- shared `WMBP_SPLIT_BE_WMBP_VP` / `WMBP_SPLIT_BE_WMBW_VP` require the BE
  outer opening to be at or right of the parent WMBP opening. Patched legacy's
  i-descending fill reads a leftward BE before computation, giving zero;
  the shared deduction geometry makes that availability boundary explicit
  under other schedules. A cached prefix must come from an admitted source.
  This is a deduction/ownership gate, not an extra weight factor: retain the
  existing local `PB^2` multiplier and the explicit BE child.

Support/local-factor boundary:

- anchor-rich `WMBP` / `WMB` rules may need BE-derived geometry at the weight
  boundary even when the rule-local multiplier stays small
- `WMB_SPLIT_BE_WMBP_WI` currently keeps its local factor at `PB`; anchored BE
  bulk belongs to the support/child side, not to an extra local multiplier
- `WMBP_SPLIT_BE_*_VP` currently uses a provisional replay-local compatibility
  factor inherited from legacy/current; keep that choice isolated so later
  exact-slice audits can narrow it without touching deduction/transition shape
- `WMBP_SPLIT_BE_WI_VP` now carries the immediate child `VP` span into the
  weight layer, but this is geometry only for now
- do not treat that span as permission to read nested `VP_WI_CASE3` ownership
  through a generic `VP` public state; the current conservative proxy remains
  the anchored BE interval until an explicit right-open VP surface exists
- replay/current exact observation helpers may bundle, for
  `WMBP_SPLIT_BE_WI_VP`, the immediate child `VP` surface profile together
  with that anchored-BE proxy
- replay/current exact observation helpers may also expose the child span's
  split-keyed `FutureKeptLeft` exact-carrier items from `exact_store`
- those helpers may further distinguish the subset whose split is itself
  admissible as direct `VP_WI_CASE3` on that observed span
- bundled observation may include aggregate totals for:
  - all same-span `FutureKeptLeft` exact carriers
  - the direct-`VP_WI_CASE3`-admissible subset
- this bundle is for audit/debug and contract design only; it still does not
  authorize the rule to consume same-span or nested `VP_CASE3_LEFT` mass
- replay/current exact contexts may also expose a canonical direct
  `VP_WI_CASE3` right-open candidate profile on a child `VP` span
- that profile may carry:
  - the direct case3 split itself
  - direct left / middle / right child factors on that split
  - right-tail / middle-right public and exact surface totals
  - the matching split-keyed `FutureKeptLeft` exact carrier total
  - a reweighted right-open export candidate derived from the carrier's open
    product and the direct child right-block scale
- replay-side `WMBP_SPLIT_BE_WI_VP` local weight should stay at `PB^2`
- any anchored-BE proxy scale or direct child right-block scale belongs to the
  support/audit side, not to the rule-local multiplier itself, because the
  shared scfg2 deduction already carries anchored-BE support on an explicit
  `BE` child
- the candidate's open product may now be reconstructed directly from the
  canonical split's child factors; `FutureKeptLeft` remains as an audit bridge
- replay/current exact contexts may also expose a minimal
  `WMBP_SPLIT_BE_WI_VP` right-open support contract carrying only:
  - anchored-BE proxy fallback
  - direct case3 availability
  - child right-block scale
  - right-open export candidate total
- that support contract now belongs to the shared `WMBP` weight boundary
  itself, not just to replay-local observation code
- replay-side observation may still bundle richer child audit data, but keep
  that under an audit/observation sub-surface rather than mixing it back into
  the support contract
- when replay bundles that richer audit data, prefer a standalone audit helper
  plus a bundled projection, so support and audit stay independently testable
- apply the same rule to immediate child observation when a family-specific
  child-facing replay helper exists
- keep `FutureKeptLeft` and other child-factor detail outside that support
  contract; those stay in observation/audit surfaces only
- this still does not authorize `WMBP_SPLIT_BE_WI_VP` to consume same-span or
  nested `VP_CASE3_LEFT` mass as support
- replay/current exact contexts may expose `VP_CASE3_*` and `VP_C12_*`
  observation surfaces at the `VP` boundary
- `VP_STACK` and `VP_INTERNAL_LOOP` may share one unary-wrapper
  observation/storage helper once their immediate child span is known
- that shared helper is only a compact replay/weight boundary contract; it
  does not merge the rule identities or erase their different child-geometry
  conditions
- replay/current exact contexts may also expose `WMB_BTAIL / WMB_BSTEP /
  WMB_NONEMPTY / WMB_EMPTY_LOCAL` observation surfaces at the `WMB` boundary
- those surfaces are for measurement and future child-aware contracts; exposing
  them does not by itself authorize new ownership transfer into `WMBP/WMB`
- current exact output grouping must not treat these storage-only surfaces as
  hidden grammar children: their totals are observer-derived provenance and
  are not consumed by the active deduction generator or transition weights
- if a future transition begins consuming a carrier total, that is a semantic
  ownership change and requires output-group expansion plus a targeted
  inside/Viterbi/constrained parity witness
- `scfg2_constrained_small_n_grouping_test` enforces the current boundary by
  poisoning carrier surfaces with a huge value and requiring output groups,
  summed/best weights, and derivation counts to remain exactly invariant

## VM base-prefix factorization

- `VM_SPLIT_WMp_BASE` owns a parent-dependent multiloop closing factor and an
  `expMLbase[left_unpaired]` gap factor
- the gap factor is multiplicative, so the private beam helper may expose
  `vm_base_unpaired_step_weight() == expMLbase[1]`
- the helper recurrence aggregates only the child/gap tail; it must apply
  `split_wm_wmp_weight(i,j)` exactly once after that aggregation
- this is an algebraic factorization of the existing rule family, not a new
  energy term or a transfer of ownership from `VM` to `WMp`
- prefix extension admissibility must use the same adjacent
  `can_pair_left_span` contract as the original split generator; do not infer
  it from `is_unpaired` alone
- finite VM prefix ranking is performed on the child/gap tail before the
  parent closing factor; `vm_base_gap_weight(length)` exposes the exact gap
  table for sparse candidate evaluation
- the selected parent receives `split_wm_wmp_weight(i,j)` exactly once after
  that private-tail beam decision

## WM start-prefix factorization

- `WM_START_V` and `WM_START_WMB` share the same `expMLbase[k-i]` left-gap
  factor and differ only in the branch seed weight
- their private prefix may therefore aggregate the two seeds and extend with
  `extend_unpaired_weight()`, while retaining the original best `(rule,k)` for
  MaxProduct traceback
- this factor replaces the two start families exactly; recursive
  `WM_SPLIT_*` and `WM_EXTEND_UNPAIRED` remain separate deductions

## WIP start-prefix factorization

- `WIP_BASE_V/WMB` and `WIP_BASEPAIR_V/WMB` share the same branch seed
  factors, with the padded variants additionally owning
  `expcp_pen[k-i]`
- the private prefix may extend a seed by
  `wip_extend_unpaired_weight() == expcp_pen[1]`, but only under the original
  adjacent `can_pair_left_span` admissibility contract
- the aggregate retains an original best `(rule,k)` representative for
  MaxProduct traceback; it does not create a new grammar rule or energy term
- `WIP_SPLIT_V/WMB` and `WIP_EXTEND_UNPAIRED` remain ordinary independent
  deductions outside this factor

## VPL start-prefix factorization

- `VPL_SPLIT_VP` owns `expcp_pen[k-i]`; the private prefix exposes the same
  factor as repeated `vpl_split_vp_weight(1)` steps
- the original unsigned `border_b/Bp` gate is evaluated once per retained VP
  seed only because its result is constant on that seed's complete admissible
  left-padding interval
- one aggregate and original best `k` representative replace the deduction
  fan-out; no rule or energy ownership is added

## VPR base-prefix factorization

- `VPR_SPLIT_VP_BASEPAIR` owns `expcp_pen[j-k]`, the right unpaired run
  `[k+1,j]` (SI eq. vi; legacy's `expcp_pen[k-i]` is a bug reported as
  TheCOBRALab/CParty#5 and corrected on 2026-10-07)
- across one admissible column the private aggregate is multiplied by
  `expcp_pen[1]` before the new `VP(i,j-1) * expcp_pen[1]` seed is added
- the private state adds each VP seed with its complete original factor and
  carries the aggregate only across the original adjacent
  `can_pair_right_span` predicate
- `max(border_B,bp)` remains a public-parent admissibility gate and does not
  become a new local multiplier
- the aggregate retains the original best `k` for MaxProduct traceback;
  `VPR_SPLIT_VP_WIP` remains a separate ordinary binary deduction

## WMB private BTAIL factorization

- the complete private tail is
  `sum_k WMBP(i,k) * WI(k+1,t-1)` grouped by the inner BE close `t`
- public WMBW is not interchangeable with this factor because WMBW owns an
  additional support restriction
- `wmb_split_be_wmbp_wi_weight`, currently PB, stays outside BTAIL and is
  applied exactly once when the retained BE frame attaches
- BTAIL introduces no new energy term; it only changes the multiplication and
  summation order of the existing three-child WMB rule

## Regenerated outside weight boundary

`engine/outside.hh` reads the same supplied local weight as inside, exactly
once per deduction. It multiplies that factor by parent beta and sibling inside
values; it does not introduce an outside-specific energy/Boltzmann factor.
BE's six-rule expansion and geometric admission stay on the shared deduction
surface. Retention and boundary child aliases belong to the caller, not weights.
Production beam BPP and CapR/Raccess now use the retained numeric outside
driver. Ordinary start aggregates expand original deductions over retained
seeds, with the same exact split guards and local factors. Genuine VPR carry
and BTAIL cells reverse their shared product recurrence. No outside-specific
weight is introduced. Sampling and MEA/centroid still use the factorized forest.

Public outer ranking uses a separate numeric completion-context guide.
That score selects public states and the matching WMB/BE consumer views;
it never enters a local thermodynamic weight, posterior beta multiplier or
main semiring graph/tape. Pilot values are constants for the ranking decision,
not differentiated factors in the retained main ensemble. Both outer and
baseline ranking consume the same deduction and local-weight surfaces.

VP recursive cases 6–9 use the complete legacy border ranges. A single
maximal scaffold pair does not contract that sum to its endpoint. This is a
deduction/admission correction shared by exact and beam; the VP split local
factor and energy ownership are unchanged.
