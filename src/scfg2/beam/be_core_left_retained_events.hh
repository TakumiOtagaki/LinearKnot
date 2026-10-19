#ifndef CPARTY_SCFG2_BEAM_BE_CORE_LEFT_RETAINED_EVENTS_HH
#define CPARTY_SCFG2_BEAM_BE_CORE_LEFT_RETAINED_EVENTS_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct BeCoreLeftRetainedStats {
  std::size_t child_events = 0;
  std::size_t ancestor_probes = 0;
  std::size_t max_ancestors_per_child = 0;
  std::array<std::size_t, 2> emitted_sources{};
  std::size_t rejected_candidates = 0;
  std::size_t duplicate_sources = 0;
};

// Exact observation surface for the two recursive BE modes whose last
// semantic child is the anchored BE core. The right-WIP modes are owned by
// the separate BE_RTAIL event surface.
class BeCoreLeftRetainedEvents {
 public:
  template <class Context>
  static BeCoreLeftRetainedEvents build(cand_pos_t n, const Context &ctx) {
    BeCoreLeftRetainedEvents out;
    std::vector<cand_pos_t> open;
    for (cand_pos_t pos = 1; pos <= n; ++pos) {
      const cand_pos_t partner = ctx.pair_partner(pos);
      if (partner > pos) {
        out.ancestors_by_open_[pos] = open;
        open.push_back(pos);
      } else if (partner >= 1 && partner < pos && !open.empty() &&
                 open.back() == partner) {
        open.pop_back();
      }
    }
    return out;
  }

  void prepare_column(cand_pos_t j) {
    current_.clear();
    const auto found = pending_.find(j);
    if (found == pending_.end()) return;
    current_ = std::move(found->second);
    pending_.erase(found);
  }

  template <class Context>
  void observe_committed_column(
      const std::vector<deductions::ItemKey> &retained_be,
      Context &ctx) {
    for (const auto &child : retained_be) {
      if (child.nonterminal != grammar::NonTerminal::BE) continue;
      ++stats_.child_events;
      const auto found = ancestors_by_open_.find(child.i);
      if (found == ancestors_by_open_.end()) continue;
      stats_.max_ancestors_per_child = std::max(
          stats_.max_ancestors_per_child, found->second.size());
      for (const cand_pos_t i : found->second) {
        ++stats_.ancestor_probes;
        const cand_pos_t j = ctx.pair_partner(i);
        if (!(j > child.j)) continue;
        const deductions::ItemKey parent{
            grammar::NonTerminal::BE, i, j, child.ip, child.jp};
        append_matching(parent,
                        grammar::RuleId::BE_INTERNAL_LOOP,
                        child,
                        ctx);
        append_matching(parent,
                        grammar::RuleId::BE_WIP_BASEPAIR,
                        child,
                        ctx);
      }
    }
  }

  const std::vector<deductions::Deduction> &deductions_for(
      const deductions::ItemKey &parent) const {
    const auto found = current_.find(parent);
    if (found == current_.end()) return empty();
    return found->second;
  }

  std::vector<deductions::ItemKey> parents() const {
    std::vector<deductions::ItemKey> out;
    out.reserve(current_.size());
    for (const auto &[parent, sources] : current_) {
      if (!sources.empty()) out.push_back(parent);
    }
    std::sort(out.begin(), out.end(), [](const auto &lhs, const auto &rhs) {
      return std::tie(lhs.i, lhs.ip, lhs.jp) >
          std::tie(rhs.i, rhs.ip, rhs.jp);
    });
    return out;
  }

  const BeCoreLeftRetainedStats &stats() const { return stats_; }
  bool pending_empty() const { return pending_.empty(); }

 private:
  using Sources = std::unordered_map<deductions::ItemKey,
                                     std::vector<deductions::Deduction>,
                                     deductions::ItemKeyHash,
                                     deductions::ItemKeyEq>;

  template <class Context>
  void append_matching(const deductions::ItemKey &parent,
                       grammar::RuleId rule,
                       const deductions::ItemKey &child,
                       Context &ctx) {
    engine::AllLocalWeightModel<Context> weight(&ctx);
    bool matched = false;
    for (auto source : deductions::generate_exact_basic_deductions(
             parent, rule, ctx)) {
      engine::require_legacy_span_schedule_support(source);
      if (weights::is_zero(weight.local_weight_for(source)) || source.child_count == 0 ||
          !deductions::ItemKeyEq{}(
              source.children[rule == grammar::RuleId::BE_WIP_BASEPAIR ? 1
                                                                       : 0]
                  .item,
              child)) {
        continue;
      }
      matched = true;
      auto &sources = pending_[parent.j][parent];
      const bool duplicate = std::any_of(
          sources.begin(), sources.end(), [&](const auto &old) {
            return old.rule == source.rule && old.split.k == source.split.k &&
                old.split.l == source.split.l;
          });
      if (duplicate) {
        ++stats_.duplicate_sources;
        continue;
      }
      sources.push_back(std::move(source));
      ++stats_.emitted_sources[
          rule == grammar::RuleId::BE_INTERNAL_LOOP ? 0 : 1];
    }
    if (!matched) ++stats_.rejected_candidates;
  }

  static const std::vector<deductions::Deduction> &empty() {
    static const std::vector<deductions::Deduction> value;
    return value;
  }

  std::unordered_map<cand_pos_t, std::vector<cand_pos_t>>
      ancestors_by_open_;
  std::map<cand_pos_t, Sources> pending_;
  Sources current_;
  BeCoreLeftRetainedStats stats_{};
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_CORE_LEFT_RETAINED_EVENTS_HH
