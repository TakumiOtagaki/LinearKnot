#ifndef CPARTY_SCFG2_DECODING_INSIDE_PAIR_MAX_HH
#define CPARTY_SCFG2_DECODING_INSIDE_PAIR_MAX_HH

#include "scfg2/decoding/pair_posterior.hh"
#include "scfg2/semirings/inside_pair_max.hh"

#include <unordered_set>

namespace cparty::scfg2::decoding {

// Read a chosen DAG once; count occurrences with saturation instead of
// expanding shared subtrees. Repeated generated features remain errors.
template <class Context>
DecodedPairPosterior decode_inside_pair_max(
    const semirings::InsidePairMaxValue &root,
    const engine::OutsideMarginalSummary &marginal, cand_pos_t n, const Context &ctx,
    double gamma, PairPosteriorObjective objective) {
  using Node = semirings::PairMaxNode;
  const PairPosteriorScoring scoring(marginal, n, gamma, objective);
  const auto sealed = root.checkpointed();
  std::unordered_set<const Node *> visited;
  std::vector<const Node *> order;
  std::vector<std::pair<const Node *, bool>> stack;
  if (sealed.node) stack.emplace_back(sealed.node.get(), false);
  std::size_t factors = 0;
  while (!stack.empty()) {
    const auto [node, ready] = stack.back();
    stack.pop_back();
    if (ready) {
      order.push_back(node);
      continue;
    }
    if (!visited.insert(node).second) continue;
    stack.emplace_back(node, true);
    factors += node->factors.size();
    for (auto it = node->factors.rbegin(); it != node->factors.rend(); ++it)
      stack.emplace_back(it->get(), false);
  }
  std::unordered_map<const Node *, unsigned char> count;
  if (sealed.node) count[sealed.node.get()] = 1;
  std::vector<PosteriorPairFeature> selected;
  for (auto it = order.rbegin(); it != order.rend(); ++it) {
    const auto *node = *it;
    const auto occurrences = count[node];
    if (!occurrences) continue;
    if (node->feature) {
      if (occurrences > 1)
        throw std::runtime_error("posterior decoding derivation repeats a generated-pair state");
      selected.push_back(*scoring.feature(*node->feature, ctx));
    }
    for (const auto &child : node->factors) {
      auto &c = count[child.get()];
      c = static_cast<unsigned char>(std::min<unsigned>(2, c + occurrences));
    }
  }
  return finish_pair_posterior(selected, marginal, n, ctx, gamma, objective,
                               root.score, order.size(), factors);
}

} // namespace cparty::scfg2::decoding
#endif
