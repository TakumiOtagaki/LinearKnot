#ifndef CPARTY_SCFG2_BEAM_PREPARE_RANKING_HH
#define CPARTY_SCFG2_BEAM_PREPARE_RANKING_HH

#include "scfg2/beam/pilot_outer.hh"
#include "scfg2/beam/public_item_ranking.hh"

namespace cparty::scfg2::beam {

struct PreparedRanking {
  replay::WFinalVpWmbpOnlineWidths widths;
  PublicItemRanking ranking;
  std::size_t pilot_width = 0;
  std::size_t guide_items = 0;
  double pilot_seconds = 0.0;
  double outer_seconds = 0.0;
};

// The completed guide owns only numeric scores and immutable forced keys.
// It can be reused by different semirings/consumers of one serial session.
// Neither pilot alpha nor h is differentiated by a main graph/tape scope.
template <class Context>
PreparedRanking prepare_ranking(
    Context &ctx, cand_pos_t n, replay::WFinalVpWmbpOnlineWidths widths,
    std::size_t main_width, RankingMode mode,
    PilotExponent exponent = PilotExponent::TwoThirds) {
  using namespace outer_detail;
  PreparedRanking out;
  out.widths = std::move(widths);
  if (!main_width) throw std::invalid_argument("beam size must be positive");
  if (mode != RankingMode::Baseline && mode != RankingMode::Outer)
    throw std::invalid_argument("invalid beam ranking mode");
  (void)pilot_exponent_name(exponent); // validate even when no guide is needed
  if (mode == RankingMode::Baseline) return out;
  constexpr auto keep_all = std::numeric_limits<std::size_t>::max();
  if (main_width == keep_all) {
    for (auto width : {out.widths.vp, out.widths.wmbp, out.widths.be,
                       out.widths.vm, out.widths.v, out.widths.wmv,
                       out.widths.wmb, out.widths.wmp, out.widths.wm,
                       out.widths.wi, out.widths.wip, out.widths.wmbw,
                       out.widths.vpl, out.widths.vpr, out.widths.w})
      if (width != keep_all)
        throw std::invalid_argument("outer needs a finite main width when public states are pruned");
    return out; // KeepAll has no ranking-dependent public selection.
  }
  out.pilot_width = pilot_width_for(main_width, exponent);
  if (out.pilot_width > std::numeric_limits<std::size_t>::max() / 4)
    throw std::invalid_argument("outer frontier width overflows");
  if (out.widths.force_scaffold_derivation && !out.widths.forced_items) {
    out.widths.forced_items = replay::scaffold_only_derivation_items(ctx, n, out.widths);
  }
  auto pilot_widths = replay::WFinalVpWmbpOnlineWidths::uniform(out.pilot_width);
  pilot_widths.max_hairpin_unpaired = out.widths.max_hairpin_unpaired;
  pilot_widths.max_case3_ancestor_depth = out.widths.max_case3_ancestor_depth;
  pilot_widths.forced_items = out.widths.forced_items;
  pilot_widths.force_scaffold_derivation = out.widths.force_scaffold_derivation;
  pilot_widths.memory_observer = out.widths.memory_observer;
  pilot_widths.memory_column_interval = out.widths.memory_column_interval;
  pilot_widths.memory_phase = MemoryPhase::Pilot;
  pilot_widths.numeric_lookup_index = out.widths.numeric_lookup_index;
  const auto pilot_start = Clock::now();
  auto pilot = replay::run_w_final_vp_wmbp_online<Sum>(ctx, n, pilot_widths);
  out.pilot_seconds = seconds(pilot_start);
  const auto outer_start = Clock::now();
  if (out.widths.memory_observer)
    out.widths.memory_observer(chart_memory(MemoryPhase::Outer, "begin",
        out.pilot_width, n, pilot.chart));
  // Finalized scores suffice for ranking. Validate predecessor links as their
  // targets finalize instead of retaining a second global link dictionary.
  auto guide = bounded_outer<Context, true, false>(pilot.chart, Item{NT::W, 1, n, -1, -1}, ctx,
                             out.pilot_width, 4 * out.pilot_width, false, false,
                             out.widths.numeric_lookup_index);
  if (out.widths.memory_observer)
    out.widths.memory_observer(chart_memory(MemoryPhase::Outer, "guide_end",
        out.pilot_width, n, guide.outer.values, true));
  out.outer_seconds = seconds(outer_start);
  out.guide_items = guide.outer.values.values().size();
  // No predecessor links or pilot Values survive into the main computation.
  auto scores = std::make_shared<const Chart>(std::move(guide.outer.values));
  out.ranking = [scores = std::move(scores)](
      const Item &item, double alpha, double existing) {
    return outer_priority(scores->get(item), alpha, existing);
  };
  return out;
}

} // namespace cparty::scfg2::beam
#endif
