#ifndef CPARTY_SCFG2_SAMPLING_FACTORIZED_SAMPLER_HH
#define CPARTY_SCFG2_SAMPLING_FACTORIZED_SAMPLER_HH

#include "scfg2/engine/factorized_pair_owners.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"

#include <random>

namespace cparty::scfg2::sampling {

// Choose a term with exp(c(t)-alpha(P)), then visit its factors. The same
// item/private-state forest serves the beta recursion and stochastic replay.
class FactorizedStructureSampler {
 public:
  template <class Chart, class Context>
  FactorizedStructureSampler(const semirings::FactorizedInsideGraph &graph,
                              std::uint32_t root, const Chart &chart,
                              const Context &ctx)
      : graph_(graph), root_(root), pairs_(graph, chart, ctx) {}

  template <class Rng>
  std::vector<engine::TracebackPair> sample(Rng &rng) const {
    std::vector<engine::TracebackPair> pairs;
    std::vector<std::uint32_t> pending{root_};
    std::uniform_real_distribution<double> uniform(0.0, 1.0);
    while (!pending.empty()) {
      const auto node = pending.back();
      pending.pop_back();
      const auto *owner = pairs_.find(node);
      if (owner) pairs.push_back(*owner);
      const auto &state = graph_.states()[node];
      auto chosen = state.first_term;
      if (state.term_count > 1) {
        const double draw = uniform(rng);
        double cumulative = 0.0;
        for (std::size_t t = state.first_term;
             t < std::size_t(state.first_term) + state.term_count; ++t) {
          chosen = static_cast<std::uint32_t>(t);
          cumulative += std::exp(graph_.terms()[t].log_value - state.log_value);
          if (draw < cumulative) break;
        }
      }
      const auto &term = graph_.terms()[chosen];
      for (std::size_t f = term.first_factor;
           f < std::size_t(term.first_factor) + term.factor_count; ++f)
        pending.push_back(graph_.factors()[f]);
    }
    engine::stable_unique_traceback_pairs(&pairs);
    std::sort(pairs.begin(), pairs.end(), [](const auto &a, const auto &b) {
      return a.i != b.i ? a.i < b.i : a.j < b.j;
    });
    return pairs;
  }

 private:
  const semirings::FactorizedInsideGraph &graph_;
  std::uint32_t root_;
  engine::FactorizedPairOwners pairs_;
};

}  // namespace cparty::scfg2::sampling

#endif
