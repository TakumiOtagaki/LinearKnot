#ifndef CPARTY_SCFG2_BEAM_WMBW_RETAINED_JOIN_AGENDA_HH
#define CPARTY_SCFG2_BEAM_WMBW_RETAINED_JOIN_AGENDA_HH

#include "scfg2/deductions/generate_exact_basic.hh"

#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct WmbwRetainedJoinAgendaStats {
  std::size_t binary_join_probes = 0;
};

class WmbwRetainedJoinAgenda {
 public:
  template <class Context, class WiStartsFn, class WmbpStartsFn>
  void build_column(cand_pos_t j,
                    const Context &ctx,
                    WiStartsFn &&wi_starts,
                    WmbpStartsFn &&wmbp_starts) {
    column_ = j;
    deductions_by_parent_.clear();
    parents_.clear();
    for (const cand_pos_t wi_i : wi_starts(j)) {
      const cand_pos_t wmbp_j = wi_i - 1;
      if (wmbp_j < 1) continue;
      for (const cand_pos_t parent_i : wmbp_starts(wmbp_j)) {
        ++stats_.binary_join_probes;
        const deductions::ItemKey parent{
            grammar::NonTerminal::WMBW, parent_i, j, -1, -1};
        std::vector<deductions::Deduction> generated;
        deductions::append_exact_wmbw_split_for_wmbp_end(
            generated, parent, ctx, wmbp_j);
        if (generated.empty()) continue;
        auto [found, inserted] = deductions_by_parent_.try_emplace(parent);
        if (inserted) parents_.push_back(parent);
        found->second.insert(found->second.end(),
                             generated.begin(), generated.end());
      }
    }
    std::sort(parents_.begin(), parents_.end(), [](const auto &lhs,
                                                   const auto &rhs) {
      return lhs.i > rhs.i;
    });
  }

  std::vector<cand_pos_t> parent_starts() const {
    std::vector<cand_pos_t> starts;
    starts.reserve(parents_.size());
    for (const auto &parent : parents_) starts.push_back(parent.i);
    return starts;
  }

  const std::vector<deductions::Deduction> &deductions_for(
      const deductions::ItemKey &parent) const {
    const auto found = deductions_by_parent_.find(parent);
    if (found == deductions_by_parent_.end()) {
      static const std::vector<deductions::Deduction> empty;
      return empty;
    }
    return found->second;
  }

  // Consume this completed parent's source buffer once. Aggregate values and
  // parent lists stay available until the next column; source views expire.
  std::vector<deductions::Deduction> take_for(
      const deductions::ItemKey &parent) {
    const auto found = deductions_by_parent_.find(parent);
    if (found == deductions_by_parent_.end()) {
      return {};
    }
    return std::move(found->second);
  }

  cand_pos_t column() const { return column_; }
  const WmbwRetainedJoinAgendaStats &stats() const { return stats_; }

 private:
  cand_pos_t column_ = -1;
  WmbwRetainedJoinAgendaStats stats_{};
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      deductions_by_parent_;
  std::vector<deductions::ItemKey> parents_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMBW_RETAINED_JOIN_AGENDA_HH
