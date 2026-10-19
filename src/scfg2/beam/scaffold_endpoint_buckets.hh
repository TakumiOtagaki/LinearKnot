#ifndef CPARTY_SCFG2_BEAM_SCAFFOLD_ENDPOINT_BUCKETS_HH
#define CPARTY_SCFG2_BEAM_SCAFFOLD_ENDPOINT_BUCKETS_HH

#include "scfg2/types.hh"

#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

// A pair child has exactly one immediate scaffold parent. The child opening
// is therefore a global bucket ID; an ancestor-by-child matrix is unnecessary.
class ScaffoldEndpointBuckets {
 public:
  template <class Context>
  static ScaffoldEndpointBuckets build(cand_pos_t n, const Context &ctx) {
    ScaffoldEndpointBuckets out;
    out.parents_.assign(n + 1, -1);
    out.children_.resize(n + 1);
    out.endpoints_.resize(n + 1);
    for (cand_pos_t child = 1; child <= n; ++child) {
      if (ctx.pair_partner(child) <= child) continue;
      const auto parent = ctx.parent_index(child);
      if (parent < 0) continue;
      out.parents_[child] = parent;
      out.children_[parent].push_back(child);
    }
    return out;
  }
  void append(cand_pos_t child, cand_pos_t endpoint) {
    endpoints_[child].push_back(endpoint);
  }
  const std::vector<cand_pos_t> &bucket(cand_pos_t ancestor, cand_pos_t child) const {
    if (child < 0 || static_cast<std::size_t>(child) >= parents_.size() ||
        parents_[child] != ancestor) return empty();
    return endpoints_[child];
  }
  const std::vector<cand_pos_t> &children(cand_pos_t ancestor) const {
    if (ancestor < 0 || static_cast<std::size_t>(ancestor) >= children_.size()) return empty();
    return children_[ancestor];
  }
 private:
  static const std::vector<cand_pos_t> &empty() {
    static const std::vector<cand_pos_t> out;
    return out;
  }
  std::vector<cand_pos_t> parents_;
  std::vector<std::vector<cand_pos_t>> children_, endpoints_;
};
} // namespace cparty::scfg2::beam
#endif
