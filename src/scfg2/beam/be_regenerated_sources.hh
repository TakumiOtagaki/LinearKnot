#ifndef CPARTY_SCFG2_BEAM_BE_REGENERATED_SOURCES_HH
#define CPARTY_SCFG2_BEAM_BE_REGENERATED_SOURCES_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/be_core_left_rule_index.hh"
#include "scfg2/beam/be_left_right_frame_events.hh"
#include "scfg2/beam/be_right_shell_source_agenda.hh"
#include "scfg2/beam/be_right_tail_retained_core_index.hh"

namespace cparty::scfg2::beam {

// The final retained chart is sufficient for all six BE rules. Build O(nb)
// child/event indices once, then visit columns in either order. Right-tail
// IDs contain no numerical products or saved deductions. Original deductions are reified
// only for the current column and use the shared guards/local weights.
template <class Semiring>
class BeRegeneratedSources {
 public:
  template <class Chart, class Context>
  static BeRegeneratedSources from_retained_chart(
      cand_pos_t n, const Chart &chart, const Context &ctx) {
    BeRegeneratedSources out;
    out.left_ = BeCoreLeftRuleIndex::template from_retained_chart<Semiring>(n, chart, ctx);
    auto right = BasicBeRightTailRetainedCoreIndex<typename Semiring::value_type>::build(n, ctx);
    // Register every retained core before completing WIP dependencies. The
    // original forward order changes no alpha or registration admissibility.
    for (const auto &[item, value] : chart.values())
      if (item.nonterminal == grammar::NonTerminal::BE && !Semiring::is_zero(value))
        right.observe_retained_be_column({item}, ctx);
    out.right_.resize(n + 1);
    right.for_each_pending_child_tail([&](const auto &child, const auto &key) {
      if (!Semiring::is_zero(chart.get(child))) out.right_[key.j].push_back(key);
    });
    return out;
  }

  template <class Chart, class Context>
  void prepare_column(cand_pos_t j, const Chart &chart, Context &ctx,
                       std::size_t tail_width, std::size_t shell_width) {
    current_.clear();
    left_.for_each_column_source(j, ctx, [&](const auto &source) {
      if (!Semiring::is_zero(chart.get(source.parent)))
        current_[source.parent].push_back(source);
    });
    std::vector<BasicScoredBeRightTail<typename Semiring::value_type>> tails;
    for (const auto &key : right_[j]) {
      const deductions::ItemKey core{grammar::NonTerminal::BE,
          key.core.k, key.core.l, key.core.ip, key.core.jp};
      const deductions::ItemKey child{grammar::NonTerminal::WIP,
          key.core.l + 1, j - 1, -1, -1};
      const auto factor = Semiring::times(chart.get(core), chart.get(child));
      if (!Semiring::is_zero(factor)) tails.push_back({key, factor});
    }
    TopWidthBeRightTailSelection tail_selection(tail_width);
    tail_selection.template commit_column<Semiring>(j, tails);
    const auto frames = score_be_right_shell_frames_from_tail_events<Semiring>(
        tails, chart, ctx, [&](const auto &key) { return tail_selection.allows(key); });
    TopWidthBeRightShellFrameSelection shell_selection(shell_width);
    shell_selection.template commit_column<Semiring>(j, frames);
    BeRightShellSourceAgenda shell;
    shell.prepare_column(frames, shell_selection, ctx);
    for (const auto &parent : shell.parents()) {
      if (Semiring::is_zero(chart.get(parent))) continue;
      const auto &sources = shell.deductions_for(parent);
      current_[parent].insert(current_[parent].end(), sources.begin(), sources.end());
    }
  }

  template <class Chart, class Context>
  std::vector<deductions::Deduction> deductions_for(
      const deductions::ItemKey &parent, const Chart &chart, Context &ctx) const {
    if (Semiring::is_zero(chart.get(parent))) return {};
    auto out = deductions::generate_exact_basic_deductions(
        parent, grammar::RuleId::BE_BASE_SAMEPAIR, ctx);
    const auto stack = deductions::generate_exact_basic_deductions(
        parent, grammar::RuleId::BE_STACK, ctx);
    out.insert(out.end(), stack.begin(), stack.end());
    const auto found = current_.find(parent);
    if (found != current_.end()) out.insert(out.end(), found->second.begin(), found->second.end());
    engine::AllLocalWeightModel<Context> weight(&ctx);
    out.erase(std::remove_if(out.begin(), out.end(), [&](const auto &source) {
      engine::require_legacy_span_schedule_support(source);
      if (weights::is_zero(weight.local_weight_for(source)))
        return true;
      for (std::size_t c = 0; c < source.child_count; ++c)
        if (Semiring::is_zero(chart.get(source.children[c].item))) return true;
      return false;
    }), out.end());
    return out;
  }

 private:
  BeCoreLeftRuleIndex left_;
  std::vector<std::vector<BeRightTailKey>> right_;
  std::unordered_map<deductions::ItemKey, std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash, deductions::ItemKeyEq> current_;
};

}  // namespace cparty::scfg2::beam
#endif
