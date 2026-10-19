#ifndef CPARTY_SCFG2_SAMPLING_TAPE_SAMPLER_HH
#define CPARTY_SCFG2_SAMPLING_TAPE_SAMPLER_HH

#include "scfg2/engine/traceback_pairs.hh"
#include "scfg2/grammar/nonterminal.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <algorithm>
#include <cmath>
#include <random>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::sampling {

// Stochastic traceback on a logsumexp tape.
//
// The taped inside pass (exact or beam) is the whole derivation forest as an
// arithmetic DAG: each LogAdd node sums alternatives, each Product node
// multiplies the parts of one alternative.  Walking it from the root, picking
// one side of every LogAdd in proportion to its value and both sides of every
// Product, draws a derivation with probability weight / Z on exactly the DAG
// that produced Z, including every beam aggregate and cache.  A generated
// pair (i,j) is in the sample iff the walk visits the V or VP item node of
// (i,j), the same ownership the outside pair marginals use.
class TapeStructureSampler {
 public:
  template <class Chart, class PairabilityContext>
  TapeStructureSampler(const semirings::LogSumExpTape &tape,
                       std::size_t root,
                       const Chart &chart,
                       const PairabilityContext &ctx)
      : tape_(tape), root_(root) {
    for (const auto &[item, value] : chart.values()) {
      if (value.node == semirings::invalid_tape_node) continue;
      const auto nt = item.nonterminal;
      const bool round = nt == grammar::NonTerminal::V;
      const bool square = nt == grammar::NonTerminal::VP;
      if (!round && !square) continue;
      if (ctx.pair_partner(item.i) == item.j) continue;  // pair of G
      pair_nodes_[value.node] = engine::TracebackPair{
          item.i, item.j,
          round ? engine::TracebackPairLane::Round
                : engine::TracebackPairLane::Square};
    }
  }

  // Generated pairs of one sampled structure, sorted by (i, j).
  template <class Rng>
  std::vector<engine::TracebackPair> sample(Rng &rng) const {
    const auto &nodes = tape_.nodes();
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    std::vector<engine::TracebackPair> pairs;
    std::vector<std::size_t> stack{root_};
    while (!stack.empty()) {
      const std::size_t index = stack.back();
      stack.pop_back();
      const auto found = pair_nodes_.find(index);
      if (found != pair_nodes_.end()) {
        pairs.push_back(found->second);
      }
      const auto &node = nodes[index];
      switch (node.operation) {
        case semirings::LogSumExpTapeNode::Operation::Constant:
          break;
        case semirings::LogSumExpTapeNode::Operation::Product:
          stack.push_back(node.lhs);
          stack.push_back(node.rhs);
          break;
        case semirings::LogSumExpTapeNode::Operation::LogAdd: {
          const double lhs = nodes[node.lhs].log_value;
          const double p_lhs =
              std::isfinite(lhs) ? std::exp(lhs - node.log_value) : 0.0;
          stack.push_back(uniform(rng) < p_lhs ? node.lhs : node.rhs);
          break;
        }
      }
    }
    engine::stable_unique_traceback_pairs(&pairs);
    std::sort(pairs.begin(), pairs.end(), [](const auto &a, const auto &b) {
      return a.i != b.i ? a.i < b.i : a.j < b.j;
    });
    return pairs;
  }

 private:
  const semirings::LogSumExpTape &tape_;
  std::size_t root_;
  std::unordered_map<std::size_t, engine::TracebackPair> pair_nodes_;
};

}  // namespace cparty::scfg2::sampling

#endif  // CPARTY_SCFG2_SAMPLING_TAPE_SAMPLER_HH
