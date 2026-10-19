#ifndef CPARTY_SCFG2_BEAM_BE_BASE_STACK_RETAINED_EVENTS_HH
#define CPARTY_SCFG2_BEAM_BE_BASE_STACK_RETAINED_EVENTS_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <optional>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct BeBaseStackRetainedStats {
  std::size_t base_events = 0;
  std::size_t stack_sources = 0;
  std::size_t rejected_candidates = 0;
  std::size_t duplicate_sources = 0;
};

class BeBaseStackRetainedEvents {
 public:
  template <class Context>
  void prepare_column(cand_pos_t j, Context &ctx) {
    current_.clear();
    const auto pending = pending_.find(j);
    if (pending != pending_.end()) {
      current_ = std::move(pending->second);
      pending_.erase(pending);
    }
    const cand_pos_t i = ctx.pair_partner(j);
    if (!(i >= 1 && i < j)) return;
    const deductions::ItemKey parent{
        grammar::NonTerminal::BE, i, j, i, j};
    append_matching(parent,
                    grammar::RuleId::BE_BASE_SAMEPAIR,
                    std::nullopt,
                    ctx,
                    current_);
  }

  template <class Context>
  void observe_committed_column(
      const std::vector<deductions::ItemKey> &retained_be,
      cand_pos_t n,
      Context &ctx) {
    for (const auto &child : retained_be) {
      if (child.nonterminal != grammar::NonTerminal::BE) continue;
      const cand_pos_t i = child.i - 1;
      const cand_pos_t j = child.j + 1;
      if (i < 1 || j > n) continue;
      auto &target = pending_[j];
      append_matching(
          deductions::ItemKey{
              grammar::NonTerminal::BE, i, j, child.ip, child.jp},
          grammar::RuleId::BE_STACK,
          child,
          ctx,
          target);
      if (child.ip == child.i && child.jp == child.j) {
        append_matching(
            deductions::ItemKey{grammar::NonTerminal::BE, i, j, i, j},
            grammar::RuleId::BE_STACK,
            child,
            ctx,
            target);
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

  const BeBaseStackRetainedStats &stats() const { return stats_; }

 private:
  using Sources = std::unordered_map<deductions::ItemKey,
                                     std::vector<deductions::Deduction>,
                                     deductions::ItemKeyHash,
                                     deductions::ItemKeyEq>;

  template <class Context>
  void append_matching(const deductions::ItemKey &parent,
                       grammar::RuleId rule,
                       std::optional<deductions::ItemKey> child,
                       Context &ctx,
                       Sources &target) {
    engine::AllLocalWeightModel<Context> weight(&ctx);
    bool matched = false;
    for (auto source : deductions::generate_exact_basic_deductions(
             parent, rule, ctx)) {
      engine::require_legacy_span_schedule_support(source);
      if (weights::is_zero(weight.local_weight_for(source))) {
        continue;
      }
      if (child.has_value() &&
          (source.child_count != 1 ||
           !deductions::ItemKeyEq{}(source.children[0].item, *child))) {
        continue;
      }
      matched = true;
      auto &sources = target[parent];
      const bool duplicate = std::any_of(
          sources.begin(), sources.end(), [&](const auto &old) {
            return old.rule == source.rule;
          });
      if (duplicate) {
        ++stats_.duplicate_sources;
        continue;
      }
      sources.push_back(std::move(source));
      if (rule == grammar::RuleId::BE_BASE_SAMEPAIR) {
        ++stats_.base_events;
      } else {
        ++stats_.stack_sources;
      }
    }
    if (!matched) ++stats_.rejected_candidates;
  }

  static const std::vector<deductions::Deduction> &empty() {
    static const std::vector<deductions::Deduction> value;
    return value;
  }

  std::map<cand_pos_t, Sources> pending_;
  Sources current_;
  BeBaseStackRetainedStats stats_{};
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_BASE_STACK_RETAINED_EVENTS_HH
