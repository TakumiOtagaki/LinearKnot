#ifndef CPARTY_SCFG2_BEAM_WMBP_SHARED_RETAINED_JOIN_HH
#define CPARTY_SCFG2_BEAM_WMBP_SHARED_RETAINED_JOIN_HH

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/beam/column_store.hh"

#include <algorithm>
#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

struct WmbpSharedRetainedJoinStats {
  std::size_t wmbp_left_probes = 0;
  std::size_t wmbw_left_probes = 0;
  std::size_t emitted_sources = 0;
};

// Reverse-discover the two shared-BE WMBP modes.  A retained right VP(l,j)
// and retained left state (i,l-1) uniquely determine the WMBP parent (i,j)
// and the split coordinate l.  The shared exact generator remains responsible
// for all scaffold-border and BE-anchor admissibility.
template <class Context>
class WmbpSharedRetainedJoin {
 public:
  template <class StartsFn>
  void build_column(cand_pos_t j,
                    const std::vector<cand_pos_t> &vp_starts,
                    StartsFn &&retained_starts,
                    Context &ctx) {
    column_ = j;
    stats_ = {};
    parents_.clear();
    deductions_by_parent_.reset_column(j);
    for (const cand_pos_t l : vp_starts) {
      if (l <= 1 || l >= j) continue;
      append_mode(l,
                  grammar::NonTerminal::WMBP,
                  grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP,
                  retained_starts(grammar::NonTerminal::WMBP, l - 1),
                  &stats_.wmbp_left_probes,
                  ctx);
      append_mode(l,
                  grammar::NonTerminal::WMBW,
                  grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP,
                  retained_starts(grammar::NonTerminal::WMBW, l - 1),
                  &stats_.wmbw_left_probes,
                  ctx);
    }
    // Sorting each parent's buffer is independent of the parent traversal
    // order. parents_ lists each created entry once, before its final sort.
    for (const auto &parent : parents_) {
      auto &sources = deductions_by_parent_.find(parent)->second;
      std::sort(sources.begin(), sources.end(), [](const auto &lhs,
                                                   const auto &rhs) {
        if (lhs.rule != rhs.rule) {
          return lhs.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP;
        }
        return lhs.split.k < rhs.split.k;
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
  const WmbpSharedRetainedJoinStats &stats() const { return stats_; }

 private:
  template <class Starts>
  void append_mode(cand_pos_t l,
                   grammar::NonTerminal left_family,
                   grammar::RuleId rule,
                   const Starts &starts,
                   std::size_t *probe_count,
                   Context &ctx) {
    const std::vector<cand_pos_t> one_vp_start{l};
    for (const cand_pos_t i : starts) {
      ++*probe_count;
      const deductions::ItemKey parent{
          grammar::NonTerminal::WMBP, i, column_, -1, -1};
      if (!deductions::exact_wmbp_item_active(parent, ctx)) continue;
      auto family = deductions::
          generate_exact_wmbp_be_vp_family_deductions_for_vp_starts(
              parent,
              ctx,
              one_vp_start,
              left_family == grammar::NonTerminal::WMBP,
              left_family == grammar::NonTerminal::WMBW,
              false);
      auto &sources = left_family == grammar::NonTerminal::WMBP
          ? family.wmbp
          : family.wmbw;
      for (auto &source : sources) {
        if (source.rule != rule || source.child_count != 3 ||
            source.children[1].item.nonterminal != left_family ||
            source.children[1].item.i != i ||
            source.children[1].item.j != l - 1) {
          continue;
        }
        auto [found, inserted] = deductions_by_parent_.try_emplace(parent);
        if (inserted) parents_.push_back(parent);
        found->second.push_back(std::move(source));
        ++stats_.emitted_sources;
      }
    }
  }

  cand_pos_t column_ = -1;
  WmbpSharedRetainedJoinStats stats_{};
  ColumnStore<std::vector<deductions::Deduction>, grammar::NonTerminal::WMBP, 8>
      deductions_by_parent_;
  std::vector<deductions::ItemKey> parents_;
};

}  // namespace cparty::scfg2::beam

#endif
