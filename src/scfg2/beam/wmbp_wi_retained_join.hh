#ifndef CPARTY_SCFG2_BEAM_WMBP_WI_RETAINED_JOIN_HH
#define CPARTY_SCFG2_BEAM_WMBP_WI_RETAINED_JOIN_HH

#include "scfg2/beam/wmbp_wi_prefix_store.hh"
#include "scfg2/deductions/generate_exact_basic.hh"

#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct WmbpWiRetainedJoinStats {
  std::size_t prefix_vp_probes = 0;
  std::size_t direct_probes = 0;
  std::size_t emitted_sources = 0;
};

template <class Semiring, class Context>
class WmbpWiRetainedJoin {
 public:
  void build_column(cand_pos_t j,
                    const std::vector<cand_pos_t> &vp_starts,
                    const WmbpWiPrefixStore<Semiring> &prefixes,
                    Context &ctx) {
    column_ = j;
    stats_ = {};
    parents_.clear();
    deductions_by_parent_.clear();
    for (const cand_pos_t l : vp_starts) {
      const deductions::ItemKey direct_parent{
          grammar::NonTerminal::WMBP, l, j, -1, -1};
      ++stats_.direct_probes;
      if (deductions::exact_wmbp_item_active(direct_parent, ctx)) {
        auto direct = deductions::generate_exact_basic_deductions(
            direct_parent, grammar::RuleId::WMBP_DIRECT_VP, ctx);
        for (auto &source : direct) append(std::move(source));
      }
      for (const auto &prefix : prefixes.retained(l)) {
        ++stats_.prefix_vp_probes;
        const deductions::ItemKey parent{
            grammar::NonTerminal::WMBP, prefix.be.i, j, -1, -1};
        if (!deductions::exact_wmbp_item_active(parent, ctx)) continue;
        const std::vector<cand_pos_t> one_start{l};
        auto family = deductions::
            generate_exact_wmbp_be_vp_family_deductions_for_vp_starts(
                parent, ctx, one_start, false, false, true);
        for (auto &source : family.wi) {
          const bool same_be = deductions::ItemKeyEq{}(
              source.children[0].item, prefix.be);
          const bool same_wi = prefix.has_wi
              ? source.child_count == 3 && deductions::ItemKeyEq{}(
                    source.children[1].item, prefix.wi)
              : source.child_count == 2;
          if (same_be && same_wi) append(std::move(source));
        }
      }
    }
    for (auto &[parent, sources] : deductions_by_parent_) {
      (void)parent;
      std::sort(sources.begin(), sources.end(), [](const auto &lhs,
                                                   const auto &rhs) {
        if (lhs.rule != rhs.rule) {
          return lhs.rule == grammar::RuleId::WMBP_DIRECT_VP;
        }
        return lhs.split.p < rhs.split.p;
      });
    }
    std::sort(parents_.begin(), parents_.end(), [](const auto &lhs,
                                                   const auto &rhs) {
      return lhs.i > rhs.i;
    });
  }

  const std::vector<deductions::ItemKey> &parents() const { return parents_; }
  const std::vector<deductions::Deduction> &deductions_for(
      const deductions::ItemKey &parent) const {
    const auto found = deductions_by_parent_.find(parent);
    if (found == deductions_by_parent_.end()) {
      static const std::vector<deductions::Deduction> empty;
      return empty;
    }
    return found->second;
  }
  cand_pos_t column() const { return column_; }
  const WmbpWiRetainedJoinStats &stats() const { return stats_; }

 private:
  void append(deductions::Deduction source) {
    const auto parent = source.parent;
    auto [found, inserted] = deductions_by_parent_.try_emplace(parent);
    if (inserted) parents_.push_back(parent);
    found->second.push_back(std::move(source));
    ++stats_.emitted_sources;
  }

  cand_pos_t column_ = -1;
  WmbpWiRetainedJoinStats stats_{};
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      deductions_by_parent_;
  std::vector<deductions::ItemKey> parents_;
};

}  // namespace cparty::scfg2::beam

#endif
