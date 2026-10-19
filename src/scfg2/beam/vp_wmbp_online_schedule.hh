#ifndef CPARTY_SCFG2_BEAM_VP_WMBP_ONLINE_SCHEDULE_HH
#define CPARTY_SCFG2_BEAM_VP_WMBP_ONLINE_SCHEDULE_HH

#include "scfg2/beam/wmbp_be_vp_prefix_events.hh"

#include <array>
#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

struct VpWmbpOnlineScheduleStats {
  std::size_t vp_parent_visits = 0;
  std::size_t be_parent_visits = 0;
  std::size_t wmbp_parent_visits = 0;
  std::size_t reverse_parent_vp_probes = 0;
  std::array<std::size_t,
             static_cast<std::size_t>(grammar::NonTerminal::BE) + 1>
      family_parent_visits{};
};

// Whole-column phases after the VP, BE, and WMBP barriers.  This ordering is
// deduction-derived: in particular WI must precede WMBW and WIP must precede
// VPR.  Keeping the family loop outside the start loop is the exact KeepAll
// shape needed before a finite right-boundary frontier can replace the full
// start range.
inline constexpr std::array<grammar::NonTerminal, 12>
    vp_wmbp_online_remaining_phases{
        grammar::NonTerminal::VM,
        grammar::NonTerminal::V,
        grammar::NonTerminal::WMv,
        grammar::NonTerminal::VPL,
        grammar::NonTerminal::WMB,
        grammar::NonTerminal::WMp,
        grammar::NonTerminal::WM,
        grammar::NonTerminal::WIP,
        grammar::NonTerminal::WI,
        grammar::NonTerminal::WMBW,
        grammar::NonTerminal::VPR,
        grammar::NonTerminal::W,
    };

// Observable full phase order for dependency tests. VP, BE, and WMBP are the
// three explicit barriers above vp_wmbp_online_remaining_phases.
inline constexpr std::array<grammar::NonTerminal, 15>
    vp_wmbp_online_column_phases{
        grammar::NonTerminal::VP,
        grammar::NonTerminal::BE,
        grammar::NonTerminal::WMBP,
        grammar::NonTerminal::VM,
        grammar::NonTerminal::V,
        grammar::NonTerminal::WMv,
        grammar::NonTerminal::VPL,
        grammar::NonTerminal::WMB,
        grammar::NonTerminal::WMp,
        grammar::NonTerminal::WM,
        grammar::NonTerminal::WIP,
        grammar::NonTerminal::WI,
        grammar::NonTerminal::WMBW,
        grammar::NonTerminal::VPR,
        grammar::NonTerminal::W,
    };

inline constexpr std::size_t vp_wmbp_online_column_phase_rank(
    grammar::NonTerminal nonterminal) {
  for (std::size_t rank = 0; rank < vp_wmbp_online_column_phases.size();
       ++rank) {
    if (vp_wmbp_online_column_phases[rank] == nonterminal) return rank;
  }
  return vp_wmbp_online_column_phases.size();
}

template <class BeginColumnFn,
          class VpParentStartsFn,
          class BeParentsFn,
          class UpdateItemFn,
          class CommitVpFn,
          class CommitBeFn,
          class DiscoverWmbpParentsFn,
          class CommitWmbpFn,
          class RemainingPhaseStartsFn,
          class CommitRemainingPhaseFn>
VpWmbpOnlineScheduleStats execute_vp_wmbp_online_schedule(
    cand_pos_t n,
    BeginColumnFn &&begin_column,
    VpParentStartsFn &&vp_parent_starts,
    BeParentsFn &&be_parents,
    UpdateItemFn &&update_item,
    CommitVpFn &&commit_vp,
    CommitBeFn &&commit_be,
    DiscoverWmbpParentsFn &&discover_wmbp_parents,
    CommitWmbpFn &&commit_wmbp,
    RemainingPhaseStartsFn &&remaining_phase_starts,
    CommitRemainingPhaseFn &&commit_remaining_phase) {
  VpWmbpOnlineScheduleStats stats;
  for (cand_pos_t j = 1; j <= n; ++j) {
    begin_column(j);
    for (const cand_pos_t i : vp_parent_starts(j)) {
      update_item(deductions::ItemKey{
          grammar::NonTerminal::VP, i, j, -1, -1});
      ++stats.vp_parent_visits;
    }
    commit_vp(j);

    // BE has no same-column dependency on WMBP. Every public BE rule is now
    // event-owned, so visit only the union of its pending anchored parents.
    for (const auto &parent : be_parents(j)) {
      update_item(parent);
      ++stats.be_parent_visits;
    }
    // WMB consumes complete BE items at this right boundary.  Commit that
    // consumer view only after every full BE item has been evaluated, and
    // before the later WMBP/remaining-family phases can observe the column.
    commit_be(j);
    const auto discovery = discover_wmbp_parents(j);
    stats.reverse_parent_vp_probes += discovery.parent_vp_probes;
    for (const auto &parent : discovery.parents) {
      update_item(parent);
      ++stats.wmbp_parent_visits;
    }
    commit_wmbp(j);

    for (const auto nonterminal : vp_wmbp_online_remaining_phases) {
      for (const cand_pos_t i : remaining_phase_starts(nonterminal, j)) {
        update_item(deductions::ItemKey{
            nonterminal, i, j, -1, -1});
        ++stats.family_parent_visits[static_cast<std::size_t>(nonterminal)];
      }
      commit_remaining_phase(nonterminal, j);
    }
  }
  return stats;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_WMBP_ONLINE_SCHEDULE_HH
