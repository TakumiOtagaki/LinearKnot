#ifndef CPARTY_SCFG2_BEAM_BE_LEFT_RIGHT_FRAME_EVENTS_HH
#define CPARTY_SCFG2_BEAM_BE_LEFT_RIGHT_FRAME_EVENTS_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/be_left_right_frame_selection.hh"
#include "scfg2/beam/be_right_tail_selection.hh"
#include "scfg2/deductions/expand_basic.hh"

#include <vector>

namespace cparty::scfg2::beam {

template <class Semiring, class Chart, class Context, class TailAllowsFn>
std::vector<BasicScoredBeRightShellFrame<typename Semiring::value_type>>
score_be_right_shell_frames_from_tail_events(
    const std::vector<BasicScoredBeRightTail<
        typename Semiring::value_type>> &tails,
    const Chart &chart,
    Context &ctx,
    TailAllowsFn &&tail_allows) {
  engine::AllLocalWeightModel<Context> weight(&ctx);
  using Value = typename Semiring::value_type;
  using Scored = BasicScoredBeRightShellFrame<Value>;
  std::vector<Scored> out;
  for (const auto &tail_candidate : tails) {
    const auto &tail = tail_candidate.key;
    if (!tail_allows(tail)) continue;
    const cand_pos_t i = ctx.pair_partner(tail.j);
    if (!(i >= 1 && i < tail.core.k)) continue;
    const deductions::ItemKey parent{
        grammar::NonTerminal::BE,
        i,
        tail.j,
        tail.core.ip,
        tail.core.jp};
    grammar::RuleSplit split{};
    split.k = tail.core.k;
    split.l = tail.core.l;
    const auto append_mode = [&](grammar::RuleId rule,
                                 BeRightShellFrameMode mode,
                                 Value factor) {
      const auto deduction = deductions::expand_basic_rule(parent, rule, split);
      const auto local = weight.local_weight_for(deduction);
      if (weights::is_zero(local)) return;
      out.push_back(Scored{
          BeRightShellFrameKey{tail, i, mode},
          Semiring::times(Semiring::from_weight(local), factor)});
    };
    if (ctx.is_empty_region(i, tail.core.k)) {
      append_mode(grammar::RuleId::BE_BASEPAIR_WIP,
                  BeRightShellFrameMode::RightOnly,
                  tail_candidate.factor);
    }
    if (ctx.weakly_closed(i + 1, tail.core.k - 1)) {
      const deductions::ItemKey left_wip{grammar::NonTerminal::WIP,
                                         i + 1,
                                         tail.core.k - 1,
                                         -1,
                                         -1};
      append_mode(grammar::RuleId::BE_WIP_WIP,
                  BeRightShellFrameMode::LeftRight,
                  Semiring::times(
                      tail_candidate.factor, chart.get(left_wip)));
    }
  }
  return out;
}

inline bool same_scored_be_right_shell_frames(
    std::vector<ScoredBeRightShellFrame> lhs,
    std::vector<ScoredBeRightShellFrame> rhs) {
  std::sort(lhs.begin(), lhs.end(), be_right_shell_frame_rank_before);
  std::sort(rhs.begin(), rhs.end(), be_right_shell_frame_rank_before);
  if (lhs.size() != rhs.size()) return false;
  for (std::size_t idx = 0; idx < lhs.size(); ++idx) {
    if (!(lhs[idx].key == rhs[idx].key)) return false;
    const double scale = std::max(
        1e-14,
        std::max(std::fabs(lhs[idx].complete_contribution),
                 std::fabs(rhs[idx].complete_contribution)) *
            1e-12);
    if (std::fabs(lhs[idx].complete_contribution -
                  rhs[idx].complete_contribution) > scale) {
      return false;
    }
  }
  return true;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_LEFT_RIGHT_FRAME_EVENTS_HH
