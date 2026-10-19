# Finite-beam execution

This layer consumes the shared exact grammar, deductions, weights and
semirings. It manages retained states, factorized joins, ranking and source
regeneration. It does not define an independent energy model.

Public calls use outer ranking by default: a LogSumExp pilot of width
ceil(b^(2/3)) builds bounded completion-context scores. Exponent 0.5 selects
a smaller pilot. Baseline ranking is also available. Missing guide values
fall back to existing priorities; candidate admission does not require a
pilot occurrence.

Production guide preparation saves numeric scores. Reference/observation
helpers may retain more provenance. Main LogSumExp and LogMax calls can
retain different supports. MEA/centroid selection follows the alpha-ranked
inside support, not a beam selected by posterior reward.

Completed states are pruned before downstream consumption. Private prefix
selectors bound partial-product candidates within activation groups.
Source agendas transfer candidate ownership to their consuming phase;
expired column histories and consumed event vectors are released.

Numeric outside reconstructs contributions from retained states using the
same product-rule derivatives and admission/selection contracts. MFE keeps
selected decompositions. Sampling and subopt use a factorized forest of
alternatives, so their storage differs from numeric outside.

`BeamConfig::memory_observer` can observe forcing, pilot, outer and main
without changing selection. Live family counts and cumulative work counters
have different meanings; optional diagnostic storage is accounted separately.

See [ranking configuration](../api/README.md#outer-and-baseline-ranking),
[model scope](../../../documents/model.md), and
[complexity](../../../documents/complexity.md).
