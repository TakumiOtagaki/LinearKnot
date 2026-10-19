#ifndef CPARTY_SCFG2_ENGINE_VITERBI_RESULT_HH
#define CPARTY_SCFG2_ENGINE_VITERBI_RESULT_HH

#include "scfg2/engine/chart.hh"
#include "scfg2/engine/choice_map.hh"
#include "scfg2/engine/traceback.hh"
#include "scfg2/engine/traceback_pairs.hh"
#include "scfg2/engine/traceback_render.hh"

#include <utility>
#include <vector>

namespace cparty::scfg2::engine {

template <class Value>
struct ViterbiRunResult {
  Chart<Value> chart;
  ChoiceMap<Value> choices;
};

template <class Value>
struct ViterbiTracebackResult {
  ViterbiRunResult<Value> run;
  deductions::ItemKey root;
  DeductionTraceNode<Value> trace;
  std::vector<TracebackPair> pairs;
  std::vector<TracebackPair> generated_pairs;
  TracebackLaneStructures lane_structures;
};

template <class Value, class PairabilityContext>
ViterbiTracebackResult<Value> build_viterbi_traceback_result(
    const deductions::ItemKey &root,
    ViterbiRunResult<Value> run,
    const PairabilityContext &ctx) {
  auto trace = build_traceback_tree(root, run.chart, run.choices);
  auto pairs = collect_traceback_pairs(trace, ctx);
  auto generated_pairs = collect_generated_traceback_pairs(trace, ctx);
  auto lanes = render_traceback_pairs(root.j, pairs);
  return {std::move(run),
          root,
          std::move(trace),
          std::move(pairs),
          std::move(generated_pairs),
          std::move(lanes)};
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_VITERBI_RESULT_HH
