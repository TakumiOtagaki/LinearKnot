#ifndef CPARTY_SCFG2_BEAM_BE_CORE_LEFT_PREFIX_AGENDA_HH
#define CPARTY_SCFG2_BEAM_BE_CORE_LEFT_PREFIX_AGENDA_HH

#include "scfg2/beam/be_core_left_rule_index.hh"

#include <algorithm>
#include <cstddef>
#include <tuple>
#include <utility>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct BeCoreLeftPrefixAgendaStats {
  std::size_t max_left_prefixes_per_core = 0;
  std::size_t retained_core_events = 0;
  std::size_t internal_parent_probes = 0;
  std::size_t left_parent_probes = 0;
  std::size_t emitted_sources = 0;
  std::size_t rejected_candidates = 0;
  std::size_t duplicate_sources = 0;
  std::size_t max_parent_group = 0;
  std::size_t unordered_parent_columns = 0;
};

// Only the current column's join results are materialized. Future columns
// retain child identities and prefix events through BeCoreLeftRuleIndex.
class BeCoreLeftPrefixAgenda {
 public:
  template <class Context>
  static BeCoreLeftPrefixAgenda build(cand_pos_t n, const Context &ctx) {
    BeCoreLeftPrefixAgenda out;
    out.index_ = BeCoreLeftRuleIndex::build(n, ctx);
    return out;
  }

  template <class Context>
  void prepare_column(cand_pos_t j, Context &ctx) {
    current_.clear();
    parents_.clear();
    std::vector<std::pair<cand_pos_t, std::vector<deductions::ItemKey>>> groups;
    bool ordered = true;
    last_prepared_ = j;
    const auto probes_before =
        stats_.internal_parent_probes + stats_.left_parent_probes;
    const auto emitted_before = stats_.emitted_sources;
    index_.for_each_column_source(j, ctx, [&](const auto &source) {
      auto [found, inserted] = current_.try_emplace(source.parent);
      found->second.push_back(source);
      if (inserted) {
        const auto child = source.split.k;
        if (groups.empty() || groups.back().first != child) {
          if (!groups.empty() && groups.back().first > child) ordered = false;
          groups.push_back({child, {}});
        }
        groups.back().second.push_back(source.parent);
      }
      ++stats_.emitted_sources;
    }, [&](grammar::RuleId rule) {
      if (rule == grammar::RuleId::BE_INTERNAL_LOOP)
        ++stats_.internal_parent_probes;
      else
        ++stats_.left_parent_probes;
    });
    stats_.rejected_candidates +=
        stats_.internal_parent_probes + stats_.left_parent_probes -
        probes_before - (stats_.emitted_sources - emitted_before);
    const auto before = [](const auto &a, const auto &b) {
      return std::tie(a.i, a.ip, a.jp) > std::tie(b.i, b.ip, b.jp);
    };
    // Readiness gives increasing direct scaffold children k. Their spans
    // are disjoint, and each group contains at most the retained b cores at
    // k. Sort each group, then concatenate in reverse child order.
    for (auto group = groups.rbegin(); group != groups.rend(); ++group) {
      auto &items = group->second;
      stats_.max_parent_group = std::max(stats_.max_parent_group, items.size());
      std::sort(items.begin(), items.end(), before);
      parents_.insert(parents_.end(), items.begin(), items.end());
    }
    if (!ordered) {
      // Preserve the standalone helper for nonmonotone event submission.
      ++stats_.unordered_parent_columns;
      std::sort(parents_.begin(), parents_.end(), before);
    }
  }

  template <class Context>
  void observe_retained_wip_column(const std::vector<cand_pos_t> &starts,
                                   cand_pos_t right, const Context &ctx) {
    index_.observe_retained_wip_column(starts, right, ctx);
    stats_.max_left_prefixes_per_core = index_.max_left_prefixes_per_core();
  }

  template <class Context>
  void observe_committed_be_column(
      const std::vector<deductions::ItemKey> &retained_be, const Context &) {
    index_.observe_committed_be_column(retained_be);
    stats_.retained_core_events = index_.core_count();
  }

  const std::vector<deductions::Deduction> &deductions_for(
      const deductions::ItemKey &parent) const {
    const auto found = current_.find(parent);
    if (found == current_.end()) {
      static const std::vector<deductions::Deduction> empty;
      return empty;
    }
    return found->second;
  }

  std::vector<deductions::ItemKey> parents() const {
    return parents_;
  }

  const BeCoreLeftPrefixAgendaStats &stats() const { return stats_; }
  bool pending_empty() const { return !index_.has_events_after(last_prepared_); }

 private:
  BeCoreLeftRuleIndex index_;
  cand_pos_t last_prepared_ = -1;
  std::unordered_map<deductions::ItemKey, std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash, deductions::ItemKeyEq> current_;
  std::vector<deductions::ItemKey> parents_;
  BeCoreLeftPrefixAgendaStats stats_{};
};

}  // namespace cparty::scfg2::beam

#endif
