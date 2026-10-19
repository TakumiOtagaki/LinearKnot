#ifndef CPARTY_SCFG2_BEAM_VP_LOOP_RETAINED_EVENTS_HH
#define CPARTY_SCFG2_BEAM_VP_LOOP_RETAINED_EVENTS_HH

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/beam/memory_observation.hh"

#include <algorithm>
#include <cstddef>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct VpLoopRetainedEventStats {
  std::size_t parent_probes = 0;
  std::size_t stack_sources = 0;
  std::size_t internal_sources = 0;
};

template <class Context>
class VpLoopRetainedEvents {
 public:
  explicit VpLoopRetainedEvents(
      grammar::NonTerminal carrier = grammar::NonTerminal::VP)
      : carrier_(carrier) {}

  void observe_committed_column(cand_pos_t child_j,
                                const std::vector<cand_pos_t> &child_starts,
                                cand_pos_t n,
                                Context &ctx) {
    for (const cand_pos_t child_i : child_starts) {
      const deductions::ItemKey child{
          carrier_, child_i, child_j, -1, -1};
      // INTERNAL requires left+right < MAXLOOP. The adjacent STACK parent
      // is also in this window; keep the original increasing traversal.
      const cand_pos_t min_parent_i =
          std::max<cand_pos_t>(1, child_i - MAXLOOP);
      for (cand_pos_t parent_i = min_parent_i;
           parent_i < child_i;
           ++parent_i) {
        const cand_pos_t left = child_i - parent_i - 1;
        const cand_pos_t max_parent_j =
            std::min<cand_pos_t>(n, child_j + MAXLOOP - left);
        for (cand_pos_t parent_j = child_j + 1;
             parent_j <= max_parent_j;
             ++parent_j) {
          ++stats_.parent_probes;
          const deductions::ItemKey parent{
              carrier_, parent_i, parent_j, -1, -1};
          const bool active = deductions::exact_vp_public_item_active(parent, ctx);
          if (!active) continue;
          // A STACK's sole child is (parent.i+1,parent.j-1). Other
          // coordinates cannot match this delivered child.
          if (parent_i == child_i - 1 && parent_j == child_j + 1)
            append_matching(parent, grammar::RuleId::VP_STACK, child, ctx);
          append_matching(parent,
                          grammar::RuleId::VP_INTERNAL_LOOP,
                          child,
                          ctx);
        }
      }
    }
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

  std::vector<cand_pos_t> parent_starts(cand_pos_t j) const {
    const auto found = parent_starts_by_column_.find(j);
    return found == parent_starts_by_column_.end()
        ? std::vector<cand_pos_t>{}
        : found->second;
  }

  const VpLoopRetainedEventStats &stats() const { return stats_; }

  // Call after this carrier's entire parent phase. No future column reads
  // these sources: the retained chart carries their completed numeric value.
  // Future STACK/INTERNAL targets remain queued (at most MAXLOOP columns).
  void release_column(cand_pos_t j) {
    const auto found = parent_starts_by_column_.find(j);
    if (found == parent_starts_by_column_.end()) return;
    for (const cand_pos_t i : found->second)
      deductions_by_parent_.erase({carrier_, i, j, -1, -1});
    parent_starts_by_column_.erase(found);
  }

  std::size_t pending_column_count() const { return parent_starts_by_column_.size(); }
  std::size_t pending_parent_count() const { return deductions_by_parent_.size(); }

  void append_memory(std::vector<ContainerMemory> &out,
                     const char *sources_name, const char *starts_name) const {
    out.push_back(hash_vector_memory(sources_name, deductions_by_parent_));
    out.push_back(hash_vector_memory(starts_name, parent_starts_by_column_));
  }

 private:
  void append_matching(const deductions::ItemKey &parent,
                       grammar::RuleId rule,
                       const deductions::ItemKey &child,
                       Context &ctx) {
    std::vector<deductions::Deduction> generated;
    if (rule == grammar::RuleId::VP_INTERNAL_LOOP) {
      // One (parent, child) probe: test the split directly instead of
      // regenerating every internal loop of the parent (O(L^2) per probe).
      if (deductions::exact_vp_internal_loop_split_allows(
              parent, child.i, child.j, ctx)) {
        grammar::RuleSplit split{};
        split.k = child.i;
        split.l = child.j;
        auto source = deductions::expand_basic_rule(parent, rule, split);
        if (deductions::is_valid(source)) generated.push_back(std::move(source));
      }
    } else {
      generated = deductions::generate_exact_basic_deductions(parent, rule, ctx);
    }
    for (auto source : generated) {
      if (source.child_count != 1 ||
          !deductions::ItemKeyEq{}(source.children[0].item, child)) {
        continue;
      }
      auto [found, inserted] = deductions_by_parent_.try_emplace(parent);
      if (inserted) {
        auto &starts = parent_starts_by_column_[parent.j];
        starts.push_back(parent.i);
      }
      const bool duplicate = std::any_of(
          found->second.begin(), found->second.end(), [&](const auto &old) {
            return old.rule == source.rule && old.split.k == source.split.k &&
                   old.split.l == source.split.l;
          });
      if (duplicate) continue;
      found->second.push_back(std::move(source));
      if (rule == grammar::RuleId::VP_STACK) {
        ++stats_.stack_sources;
      } else {
        ++stats_.internal_sources;
      }
    }
  }

  VpLoopRetainedEventStats stats_{};
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      deductions_by_parent_;
  std::unordered_map<cand_pos_t, std::vector<cand_pos_t>>
      parent_starts_by_column_;
  grammar::NonTerminal carrier_ = grammar::NonTerminal::VP;
};

}  // namespace cparty::scfg2::beam

#endif
