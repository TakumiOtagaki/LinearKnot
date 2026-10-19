#ifndef CPARTY_SCFG2_ENGINE_FACTORIZED_PAIR_OWNERS_HH
#define CPARTY_SCFG2_ENGINE_FACTORIZED_PAIR_OWNERS_HH

#include "scfg2/engine/traceback_pairs.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"

#include <unordered_map>

namespace cparty::scfg2::engine {

// Observation only: public V/VP items emit generated pairs; private prefix
// nodes are silent, and fixed scaffold pairs are overlaid by the consumer.
class FactorizedPairOwners {
 public:
  template <class Chart, class Context>
  FactorizedPairOwners(const semirings::FactorizedInsideGraph &graph,
                       const Chart &chart, const Context &ctx) {
    for (const auto &[item, value] : chart.values()) {
      if (value.node >= graph.states().size()) continue;
      const bool round = item.nonterminal == grammar::NonTerminal::V;
      const bool square = item.nonterminal == grammar::NonTerminal::VP;
      if ((!round && !square) || ctx.pair_partner(item.i) == item.j) continue;
      pairs_[value.node] = TracebackPair{item.i, item.j,
          round ? TracebackPairLane::Round : TracebackPairLane::Square};
    }
  }

  const TracebackPair *find(std::uint32_t node) const {
    const auto it = pairs_.find(node);
    return it == pairs_.end() ? nullptr : &it->second;
  }

 private:
  std::unordered_map<std::uint32_t, TracebackPair> pairs_;
};

}  // namespace cparty::scfg2::engine
#endif
