# Local weights

This layer owns family-scoped transition contexts and local factors keyed by
shared `RuleId`. It does not own deduction generation, fill order or chart
storage. Exact, beam, reverse and replay consumers use the same weight source.

## Source and model provenance

`log_loops.hh` is a log-space translation of the bundled ViennaRNA PF
hairpin/internal/stem helpers; its location under SCFG2 does not make that
implementation independent of ViennaRNA. The shared raw/log parameter-table
builder remains in [params.c](../../ViennaRNA/params/params.c). The runtime
snapshot calls that backend through the local-weight boundary; it does not
supply a second thermodynamic model.

`pseudoknot_parameters.hh` retains the CParty/DP09 model values and band
multipliers. The six top-level `params/*.par` files are inherited data, with
their pinned source and hashes recorded in the
[third-party notices](../../../THIRD_PARTY_NOTICES.md). Source origin, model
data and permission to relicense are separate questions. The existing license
grants and scientific attribution remain in place.

## Runtime snapshot and parameter ownership

[RuntimeSnapshot](runtime_snapshot.hh) takes sequence, scaffold and
`RuntimeOptions`. It owns normalized/encoded sequence, log Vienna parameters,
long-loop extrapolation values and pseudoknot parameters, including the
0.89/0.74 band multipliers. Constructing another snapshot does not overwrite
an existing owner's local factors. Its immutable scaffold geometry owns the
unpaired-run counts; the replay context reads them without a duplicate array.

The production snapshot owns no materialized DP or diagnostic tables.
Independent observations belong to test support, including
`testing::ObservedRuntimeSnapshot` and `storage::exact_store`. Their totals
must not become hidden grammar children or local-weight inputs. The
[runtime snapshot witness](../../../test/scfg2_runtime_snapshot_test.cc) and
[carrier-poisoning witness](../../../test/scfg2_constrained_small_n_grouping_test.cc)
check those boundaries.

`runtime_parameters.hh` honors `CPARTY_FIXED_ENERGY_PARAM_FILE` before
snapshotting. Empty, unreadable or unloadable selections fail. The loader
requires a v2.0 header and at least one recognized parameter section before
mutating the model. Valid partial files overlay the current globals;
a complete file is needed for a model independent of previous loads.
Section-body parsing remains Vienna's existing implementation. Without an
override, a fresh process starts with the bundled Turner2004 defaults.
Parameter loading is process-global: keep configuration fixed for a public
session's lifetime and use serial calls.

`sequence.hh` normalizes bases to uppercase RNA (`T` to `U`) for encoded
bases and special-loop lookup. Public session input getters retain the
supplied spelling. Encoded sequences and parameter snapshots use RAII,
including constructor-unwinding paths.

## Production value domain

Public exact/beam consumers use `WFinalLogWeightContext` and tagged
`LogWeight` values. `AllLocalWeightModel::local_weight_for` is the kernel
entry; the tag prevents an accidental second logarithm. Impossible factors
are `-infinity`, neutral factors are zero, and products add log factors.
`weight_for()` raw projections and `WFinalExactInsideContext` are bounded
comparison/observation interfaces.

Numeric exact forward execution may skip this pure local-weight lookup when
a canonical child is exactly semiring zero. It preserves the local factor
and multiplication order for every live product; finite log weights are
never thresholded. Reference execution still evaluates zero-child candidates
for observation. Numeric split selection may omit such candidates before
construction, including their pure admission-predicate calls. It uses the
shared child identities and preserves the order of surviving alternatives;
candidate-buffer reuse does not interleave generation and evaluation.
Reverse propagation evaluates the shared local weight even
when the inside product is zero, since a child adjoint can remain nonzero.

The bundled Vienna PF builder has one shared raw/log table builder.
Its log entry evaluates exponents directly, preserving temperature,
`pf_smooth`, dangle/mismatch smoothing and special-loop tables.
`log_loops.hh` follows PF helper branches; integer MFE helpers are not the
weight source. An `INF` short-hairpin sentinel is impossible, whereas a
finite energy whose exponential underflows still has a finite log weight.

The snapshot's `exp*` fields contain log weights and must not be passed to
raw `exp_E_*` helpers. Length factors use length times the corresponding
log constant. There are no numerical PF scale tables or free-energy scale
corrections. Physical length penalties and BE coefficients remain in the
score. The inactive `VM_SCALE2` inventory entry has identity weight and
performs no rescaling.

Thermodynamic k-best terms observe `local_log_weight` accumulated from the
same injected factors during forward multiplication. It is not recovered by
subtracting child inside values from a large term mass. Wrapper/child edges
carry zero local weight. See the [shared model](../../../documents/model.md#shared-grammar-and-weights).

## Loop-weight boundaries

### Ordinary multiloop closure

The closing stem reverses pair orientation to `(j,i)` and uses the two
**inner** neighbors `s[j-1]`, `s[i+1]` for dangles=1/2. They exist even when
`i==1` or `j==n`; exterior endpoint guards do not apply. Branch/exterior
stems use outer neighbors `s[i-1]`, `s[j+1]` with their endpoint guards.
Dangles=0 supplies neither neighbor. Raw/log contexts expose the same
factor to exact, beam, outer and reverse consumers.

Oracle patch `0009` corrects CParty PF `exp_Mbloop` at the same boundary.
`LegacyCompatible` selects the historical completed-WMB domain, not this
retired endpoint error. The [log-weight tests](../../../test/scfg2_log_weight_test.cc)
cover closing-neighbor and fixed-structure cases.

### Long scaffold internal loops

- `internal_loop.hh` preserves the bounded raw reference evaluation. Sizes
  <= `MAXLOOP` use the Vienna helper. Larger loops use explicit logarithmic
  initiation from size 30, generic/1xn mismatches or bulge terminal-AU
  factors, and capped NINIO asymmetry without indexing a size-limited array.
- BE follows fixed scaffold pairs and remains admissible for sizes >30.
  V/VP search windows remain part of deduction support.
- Each snapshot owns the unrounded temperature-scaled size-30 initiation
  energies, NINIO slope and cap. Later global parameter loads do not change
  that owner's long-loop factors.
- `pf_smooth=0` truncates the whole initiation energy and NINIO slope as
  Vienna's expanded PF tables do; `pf_smooth=1` preserves continuous energies.
  BE's pseudoknot multiplier stays outside the loop helper.

Production `log_loops.hh` expresses the same ownership in log space.
Oracle patch `0007` supplies an independent long-loop implementation;
see the [oracle inventory](../../../test/oracle/legacy_cparty/README.md#patch-inventory)
and [internal-loop witness](../../../test/scfg2_internal_loop_weight_test.cc).

## Deduction support and local factors

The [VP/WMB geometry contract](../deductions/README.md#vpwmb-geometry-contract)
belongs to deduction support. It must not be emulated by a second weight
model. `square` carries crossing geometry relative to scaffold G and does
not introduce a separate scoring semantics.

`WmbClosurePolicy::RequireWeakClosure` is the default. Completed WMB admission
requires weak closure; partial WMBP/WMBW prefixes keep their existing domains.
`LegacyCompatible` restores the historical completed-WMB domain explicitly.
Rule-local factors, including Psm ownership, are unchanged. The external
oracle's pin and ordered patches have one maintained source in
[config.json](../../../test/oracle/legacy_cparty/config.json); compare its
inside result with `LegacyCompatible`.

For BE/VP joins, preserve these ownership boundaries:

- `WMBP_SPLIT_BE_WMBP_VP` and `WMBP_SPLIT_BE_WMBW_VP` require
  `BE.outer_open >= WMBP.i`. This explicit geometry gate preserves patched
  legacy fill availability. Their local multiplier remains `PB^2`.
- `WMBP_SPLIT_BE_WI_VP` also adds only `PB^2`; anchored scaffold support is
  carried by the explicit BE child. Its anchor checks do not authorize
  reading same-span/nested VP carrier totals as extra mass.
- `WMBP_DIRECT_VP` names its canonical `VP` child directly and adds `PB`.
- `WMB_SPLIT_BE_WMBP_WI` adds `PB`; the BE and private-tail factors remain
  on the child/support side.

Storage provenance such as `WmbpStepPrefixKind` in
[the surface vocabulary](../grammar/surface_kind.hh) classifies geometry.
It does not add hidden children or a second scoring rule. Changing a
transition to consume diagnostic totals would change semantics and requires
a targeted inside/Viterbi/constrained-output witness.

MEA and Hamada gamma-centroid derive marginals from the same local weights.
Their additive rewards live in [pair_posterior.hh](../decoding/pair_posterior.hh),
outside the thermal weight model. Beam decoding replays unchanged inside alpha with selected
posterior maxima; all pruning and representative choices still use alpha.
Exact decoding uses its frozen forest. Both objectives score physical pair
marginals without a separate band gain. See the
[decoding objectives](../../../documents/model.md#structure-selection).

## Prefix factorizations

The formulas below use ordinary multiplicative weight notation. Production
injects their log factors through the semiring boundary; indexed length
expressions do not imply that the snapshot stores raw length arrays.

### VM base-prefix factorization

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

### WM start-prefix factorization

- `WM_START_V` and `WM_START_WMB` share the same `expMLbase[k-i]` left-gap
  factor and differ only in the branch seed weight
- their private prefix may therefore aggregate the two seeds and extend with
  `extend_unpaired_weight()`, while retaining the original best `(rule,k)` for
  MaxProduct traceback
- this factor replaces the two start families exactly; recursive
  `WM_SPLIT_*` and `WM_EXTEND_UNPAIRED` remain separate deductions

### WIP start-prefix factorization

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

### VPL start-prefix factorization

- `VPL_SPLIT_VP` owns `expcp_pen[k-i]`; the private prefix exposes the same
  factor as repeated `vpl_split_vp_weight(1)` steps
- the original unsigned `border_b/Bp` gate is evaluated once per retained VP
  seed only because its result is constant on that seed's complete admissible
  left-padding interval
- one aggregate and original best `k` representative replace the deduction
  fan-out; no rule or energy ownership is added

### VPR base-prefix factorization

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

### WMB private BTAIL factorization

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
weight is introduced. Beam MEA/centroid replay alpha-ranked execution with
selected posterior maxima; exact decoding keeps a frozen forest. Sampling and
subopt have separate derivation representations. See the
[consumer storage contracts](../../../documents/complexity.md).

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
