# semirings

Own the algebra used by the shared DP engine:

- exact/inside semiring
- viterbi semiring
- optional debug/counting semirings later

Semirings should not enumerate deductions themselves.

## Value-domain contract

The executor and factored beam helpers must not assume that a raw local
Boltzmann weight already has the semiring's representation. Built-in semirings
provide:

- `from_weight(double)`: convert a positive raw local weight into the value
  domain (bounded compatibility/reference use);
- `from_log_weight(double)` and `from_weight(weights::LogWeight)`: inject a
  tagged log local factor directly; production uses this boundary;
- `is_zero(value)`: recognize the additive identity without assuming it is
  numeric `0`;
- `priority(value)`: expose the higher-is-better scalar used only at a beam
  ranking boundary.

SumProduct and MaxProduct preserve raw values and exponentiate log inputs.
LogMax and
LogSumExp are production algebras: both use log weights, `-infinity` as zero,
addition for `times`, and log-domain values as priorities. LogMax uses `max`
for `plus`; LogSumExp uses a stable `max + log1p(exp(...))` implementation.
Both run through the same deduction generator, online executor, and factored
caches as their product-domain counterparts.

`traceback_priority(value)` provides a representation-independent comparison
scale for compressed Viterbi witnesses. It is `log(value)` for product-domain
semirings and identity for LogMax. When mathematically equal terms round
differently across representations, compressed caches compare this priority
with a small tolerance and resolve ties by a stable deduction/split key.

Finite-width beam boundaries themselves retain their ordinary exact priority
ordering. Product and log executions may consequently materialize different
unreachable side charts when boundary candidates are numerically tied. The
cross-domain contract is equality of the root score and the full
root-reachable selected path, not identity of every discarded internal item.

When adding another value representation or a new finite-path owner, require
explicit `Semiring::zero()`, route every local factor through the typed
`from_weight` boundary (or explicit `from_log_weight`),
and use `is_zero` at liveness boundaries.

## Factorized analytic outside

`FactorizedLogSumExpSemiring` evaluates exactly the same log-domain forward
operations and rankings. `checkpoint_value` seals item/private aggregate
boundaries into a `FactorizedInsideGraph`. It is a no-op for numeric and tape
values, so all consumers share the existing deduction/weight executor.

A state stores its actual forward alpha and terms `(c_actual, child IDs)`.
Products flatten their unobserved local factors, and a running sum appends
terms without a LogAdd chain. Reused private aggregates, cache factors and
observed contributions are sealed before their consumers; this preserves
shared identity and separate WM heap/forced evaluation instances. Ranking
products do not create operation nodes; disconnected evaluated aggregate
states can remain with zero posterior. State IDs follow children, allowing
one generic topological analytic beta sweep:

```
beta(root) = 0
beta(child) = logadd(beta(child), beta(parent) + c_actual - alpha(child))
```

The factorized consumer uses the equivalent normalized form for numerical stability:
`term_post = post(parent) * exp(c_actual-alpha(parent))`, accumulated into
each child occurrence. Here `post` is the expected number of state visits
in a derivation. It is a binary presence probability only when that state
cannot occur twice in one derivation. The long-input shared-BE/WMBP overlap
audit led to an explicit shared-geometry guard restoring legacy's BE
availability boundary; see
`documents/eval_structure/scfg2_outside_pair_occurrence_audit_20261008.md`.
The same terms support stochastic traceback and posterior decoding.
Their finite uniform-width forest target is O(n b²); admitted
candidate/aggregate terms are stored instead of re-ranking a second beam.
Public beam outside and CapR/Raccess use retained numeric regeneration
instead, with O(nb + b²) reverse storage; see `api/README.md`. Zeroed chart entries are erased for this
representation; captured old state IDs remain valid. Nonzero VP/WMBP chart
entries and the existing indexes retain their original forward semantics.

The exact path and taped beam/profile helpers remain validation oracles.
The arithmetic and full-kernel oracle tests cover frozen derivatives, shared
aliases, pre-commit captures, BPP, contexts and window accessibility. KeepAll
is an uncapped comparison mode; the O(n b²) target applies to finite uniform
widths and bounded configured local/ancestor limits. Arbitrary externally
forced sets add their own candidate cost.

## Inside alpha with a selected posterior maximum (2026-10-14)

`inside_pair_max.hh` is the beam MEA/centroid replay value. Its alpha follows
LogSumExp exactly; score follows additive max with first-alternative ties.
All `priority`, `traceback_priority`, `better` and ordering decisions use
alpha to reproduce the original public/private support. Do not use the
generic thermal Viterbi choice observer to recover its posterior choices;
its own immutable selected DAG carries those choices.

Pair rewards enter through the existing deduction-value transform from
`decoding/pair_posterior.hh`; local thermodynamic weights always add zero to
the posterior objective. Constants/empty traces need no node, silent unary
checkpoints alias their child, and only the selected product survives a sum.
No permanent allocation pool retains losing alternatives. Node/factor/product
live and peak counters are separate from cumulative arithmetic work.
Last-reference destruction uses an intrusive iterative queue; root output
also visits a DAG rather than expanding a shared occurrence tree.
