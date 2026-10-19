#ifndef CPARTY_SCFG2_BEAM_VP_RECURSIVE_SPLIT_RETAINED_EVENTS_HH
#define CPARTY_SCFG2_BEAM_VP_RECURSIVE_SPLIT_RETAINED_EVENTS_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <algorithm>
#include <array>
#include <cstddef>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

inline constexpr std::size_t vp_recursive_split_rule_index(
    grammar::RuleId rule) {
  switch (rule) {
    case grammar::RuleId::VP_WIP_VP_LEFT: return 0;
    case grammar::RuleId::VP_VP_WIP_RIGHT: return 1;
    case grammar::RuleId::VP_WIP_VPR: return 2;
    case grammar::RuleId::VP_VPL_WIP: return 3;
    default: return 4;
  }
}

struct VpRecursiveSplitRetainedStats {
  std::array<std::size_t, 4> join_probes{};
  std::array<std::size_t, 4> emitted_sources{};
  std::size_t rejected_candidates = 0;
  std::size_t duplicate_sources = 0;
};

// Reverse events for the four binary public-VP split rules.  Left-oriented
// rules join WIP(a,k-1) with VP/VPR(k,j-1); right-oriented rules join
// VP/VPL(a,k) with WIP(k+1,j-1).  With width-b retained child columns each
// rule performs at most b^2 probes for a target column.
class VpRecursiveSplitRetainedEvents {
 public:
  explicit VpRecursiveSplitRetainedEvents(
      grammar::NonTerminal carrier = grammar::NonTerminal::VP)
      : carrier_(carrier) {}

  template <class Context,
            class VpStartsAt,
            class WipStartsAt,
            class VplStartsAt,
            class VprStartsAt>
  void prepare_column(cand_pos_t j,
                      VpStartsAt &&vp_starts_at,
                      WipStartsAt &&wip_starts_at,
                      VplStartsAt &&vpl_starts_at,
                      VprStartsAt &&vpr_starts_at,
                      Context &ctx) {
    sources_.clear();
    std::set<std::tuple<grammar::RuleId, cand_pos_t, cand_pos_t>> seen;
    if (j <= 2) return;

    join_left(grammar::RuleId::VP_WIP_VP_LEFT,
              j,
              vp_starts_at(j - 1),
              wip_starts_at,
              ctx,
              seen);
    join_left(grammar::RuleId::VP_WIP_VPR,
              j,
              vpr_starts_at(j - 1),
              wip_starts_at,
              ctx,
              seen);
    join_right(grammar::RuleId::VP_VP_WIP_RIGHT,
               j,
               wip_starts_at(j - 1),
               vp_starts_at,
               ctx,
               seen);
    join_right(grammar::RuleId::VP_VPL_WIP,
               j,
               wip_starts_at(j - 1),
               vpl_starts_at,
               ctx,
               seen);
  }

  std::vector<deductions::Deduction> take_for(
      const deductions::ItemKey &parent) {
    const auto found = sources_.find(parent);
    if (found == sources_.end()) return {};
    auto out = std::move(found->second);
    sources_.erase(found);
    return out;
  }

  std::vector<cand_pos_t> parent_starts() const {
    std::vector<cand_pos_t> out;
    out.reserve(sources_.size());
    for (const auto &[parent, sources] : sources_) {
      if (!sources.empty()) out.push_back(parent.i);
    }
    return out;
  }

  const VpRecursiveSplitRetainedStats &stats() const { return stats_; }

 private:
  template <class Context, class WipStartsAt>
  void join_left(grammar::RuleId rule,
                 cand_pos_t j,
                 const std::vector<cand_pos_t> &right_starts,
                 WipStartsAt &&wip_starts_at,
                 Context &ctx,
                 std::set<std::tuple<grammar::RuleId,
                                     cand_pos_t,
                                     cand_pos_t>> &seen) {
    for (const cand_pos_t k : right_starts) {
      if (k <= 1) continue;
      for (const cand_pos_t a : wip_starts_at(k - 1)) {
        probe(rule, a - 1, j, k, ctx, seen);
      }
    }
  }

  template <class Context, class LeftStartsAt>
  void join_right(grammar::RuleId rule,
                  cand_pos_t j,
                  const std::vector<cand_pos_t> &wip_starts,
                  LeftStartsAt &&left_starts_at,
                  Context &ctx,
                  std::set<std::tuple<grammar::RuleId,
                                      cand_pos_t,
                                      cand_pos_t>> &seen) {
    for (const cand_pos_t wip_start : wip_starts) {
      const cand_pos_t k = wip_start - 1;
      if (k < 1) continue;
      for (const cand_pos_t a : left_starts_at(k)) {
        probe(rule, a - 1, j, k, ctx, seen);
      }
    }
  }

  template <class Context>
  void probe(grammar::RuleId rule,
             cand_pos_t i,
             cand_pos_t j,
             cand_pos_t k,
             Context &ctx,
             std::set<std::tuple<grammar::RuleId,
                                 cand_pos_t,
                                 cand_pos_t>> &seen) {
    const std::size_t idx = vp_recursive_split_rule_index(rule);
    ++stats_.join_probes[idx];
    if (i < 1 || i >= j) return;
    const deductions::ItemKey parent{
        carrier_, i, j, -1, -1};
    // The historical public-VP event inventory is deliberately broader than
    // the item-level parent gate; the online executor applies that gate when
    // it visits the parent.
    // Test the shared full-range predicate in O(1), without scanning the
    // parent's interval or restricting k to an internal scaffold pair.
    const bool admitted = deductions::exact_vp_recursive_split_allows(parent, rule, k, ctx);
    bool matched = false;
    if (admitted) {
      grammar::RuleSplit split{};
      split.k = k;
      auto source = deductions::expand_basic_rule(parent, rule, split);
      engine::AllLocalWeightModel<Context> weight(&ctx);
      if (deductions::is_valid(source)) {
        engine::require_legacy_span_schedule_support(source);
        if (!weights::is_zero(weight.local_weight_for(source))) {
          matched = true;
          if (!seen.emplace(rule, i, k).second) {
            ++stats_.duplicate_sources;
          } else {
            sources_[parent].push_back(std::move(source));
            ++stats_.emitted_sources[idx];
          }
        }
      }
    }
    if (!matched) ++stats_.rejected_candidates;
  }

  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      sources_;
  VpRecursiveSplitRetainedStats stats_{};
  grammar::NonTerminal carrier_ = grammar::NonTerminal::VP;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_RECURSIVE_SPLIT_RETAINED_EVENTS_HH
