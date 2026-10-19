#ifndef CPARTY_SCFG2_BEAM_WMBP_VP_PHASE_BARRIER_PROFILE_HH
#define CPARTY_SCFG2_BEAM_WMBP_VP_PHASE_BARRIER_PROFILE_HH

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <array>
#include <cstddef>
#include <functional>
#include <map>
#include <queue>
#include <cstdint>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct WmbpVpPhaseBarrierProfile {
  std::size_t column_count = 0;
  std::size_t item_count = 0;
  std::size_t same_column_edge_count = 0;
  std::size_t vp_barrier_edge_count = 0;
  std::size_t barrier_wmbp_edge_count = 0;
  std::size_t cyclic_column_count = 0;
  std::size_t missing_child_count = 0;
  std::uint64_t nonterminal_before_barrier_mask = 0;
  std::uint64_t nonterminal_after_barrier_mask = 0;

  bool feasible() const {
    return cyclic_column_count == 0 && missing_child_count == 0;
  }

  std::uint64_t split_nonterminal_mask() const {
    return nonterminal_before_barrier_mask &
           nonterminal_after_barrier_mask;
  }
};

struct VpPrerequisiteBarrierSchedule {
  std::vector<deductions::ItemKey> schedule;
  std::size_t prerequisite_item_count = 0;
  std::size_t wmbp_prerequisite_count = 0;
  std::size_t missing_child_count = 0;
  std::uint64_t prerequisite_nonterminal_mask = 0;

  bool feasible(std::size_t expected_size) const {
    return missing_child_count == 0 && wmbp_prerequisite_count == 0 &&
           schedule.size() == expected_size;
  }
};

struct NonterminalPrerequisiteBarrierSchedule {
  std::vector<deductions::ItemKey> schedule;
  std::size_t prerequisite_item_count = 0;
  std::size_t target_item_count = 0;
  std::size_t wmbp_prerequisite_count = 0;
  std::size_t missing_child_count = 0;
  std::uint64_t prerequisite_nonterminal_mask = 0;
  std::array<std::size_t, 64> prerequisite_count_by_nonterminal{};

  bool feasible(std::size_t expected_size) const {
    return missing_child_count == 0 && schedule.size() == expected_size;
  }

  std::size_t prerequisite_count(grammar::NonTerminal nonterminal) const {
    return prerequisite_count_by_nonterminal[
        static_cast<unsigned>(nonterminal)];
  }
};

inline deductions::ItemKey wmbp_vp_dependency_child(
    const deductions::Deduction &deduction,
    deductions::ItemKey child) {
  return child;
}

template <class Context>
WmbpVpPhaseBarrierProfile profile_wmbp_vp_phase_barrier(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx,
    std::vector<deductions::ItemKey> *barrier_schedule = nullptr) {
  std::map<cand_pos_t, std::vector<deductions::ItemKey>> columns;
  for (const auto &item : schedule) columns[item.j].push_back(item);

  WmbpVpPhaseBarrierProfile out;
  out.column_count = columns.size();
  out.item_count = schedule.size();
  for (const auto &[j, items] : columns) {
    (void)j;
    std::unordered_map<deductions::ItemKey,
                       std::size_t,
                       deductions::ItemKeyHash,
                       deductions::ItemKeyEq>
        node_for;
    for (std::size_t idx = 0; idx < items.size(); ++idx) {
      node_for.emplace(items[idx], idx);
    }
    const std::size_t barrier = items.size();
    std::vector<std::vector<std::size_t>> edges(items.size() + 1);
    std::vector<std::size_t> indegree(items.size() + 1, 0);
    const auto add_edge = [&](std::size_t from, std::size_t to) {
      edges[from].push_back(to);
      ++indegree[to];
    };

    for (std::size_t parent_idx = 0; parent_idx < items.size(); ++parent_idx) {
      for (const auto &deduction :
           deductions::generate_exact_basic_deductions_for_item(
               items[parent_idx], ctx)) {
        engine::require_legacy_span_schedule_support(deduction);
        for (std::size_t child_idx = 0;
             child_idx < deduction.child_count;
             ++child_idx) {
          auto child = wmbp_vp_dependency_child(
              deduction, deduction.children[child_idx].item);
          if (child.j != deduction.parent.j) continue;
          const auto found = node_for.find(child);
          if (found == node_for.end()) {
            ++out.missing_child_count;
            continue;
          }
          add_edge(found->second, parent_idx);
          ++out.same_column_edge_count;
        }
      }
    }
    for (std::size_t idx = 0; idx < items.size(); ++idx) {
      if (items[idx].nonterminal == grammar::NonTerminal::VP) {
        add_edge(idx, barrier);
        ++out.vp_barrier_edge_count;
      }
      if (items[idx].nonterminal == grammar::NonTerminal::WMBP) {
        add_edge(barrier, idx);
        ++out.barrier_wmbp_edge_count;
      }
    }

    std::priority_queue<std::size_t,
                        std::vector<std::size_t>,
                        std::greater<std::size_t>>
        ready;
    for (std::size_t idx = 0; idx < indegree.size(); ++idx) {
      if (indegree[idx] == 0) ready.push(idx);
    }
    std::size_t visited = 0;
    bool barrier_seen = false;
    while (!ready.empty()) {
      const std::size_t from = ready.top();
      ready.pop();
      ++visited;
      if (from == barrier) {
        barrier_seen = true;
      } else {
        const auto bit = std::uint64_t{1}
            << static_cast<unsigned>(items[from].nonterminal);
        if (barrier_seen) {
          out.nonterminal_after_barrier_mask |= bit;
        } else {
          out.nonterminal_before_barrier_mask |= bit;
        }
        if (barrier_schedule) barrier_schedule->push_back(items[from]);
      }
      for (const auto to : edges[from]) {
        if (--indegree[to] == 0) ready.push(to);
      }
    }
    if (visited != items.size() + 1) ++out.cyclic_column_count;
  }
  return out;
}

template <class Context>
std::vector<deductions::ItemKey> build_wmbp_vp_phase_barrier_schedule(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx) {
  std::vector<deductions::ItemKey> out;
  out.reserve(schedule.size());
  const auto profile =
      profile_wmbp_vp_phase_barrier(schedule, ctx, &out);
  if (!profile.feasible() || out.size() != schedule.size()) return {};
  return out;
}

template <class Context>
NonterminalPrerequisiteBarrierSchedule
build_nonterminal_set_prerequisite_barrier_schedule(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx,
    std::uint64_t target_nonterminal_mask) {
  std::map<cand_pos_t, std::vector<deductions::ItemKey>> columns;
  for (const auto &item : schedule) columns[item.j].push_back(item);
  NonterminalPrerequisiteBarrierSchedule out;
  out.schedule.reserve(schedule.size());
  for (const auto &[j, items] : columns) {
    (void)j;
    std::unordered_map<deductions::ItemKey,
                       std::size_t,
                       deductions::ItemKeyHash,
                       deductions::ItemKeyEq>
        node_for;
    for (std::size_t idx = 0; idx < items.size(); ++idx) {
      node_for.emplace(items[idx], idx);
    }
    std::vector<std::vector<std::size_t>> dependencies(items.size());
    for (std::size_t parent_idx = 0; parent_idx < items.size(); ++parent_idx) {
      for (const auto &deduction :
           deductions::generate_exact_basic_deductions_for_item(
               items[parent_idx], ctx)) {
        engine::require_legacy_span_schedule_support(deduction);
        for (std::size_t child_idx = 0;
             child_idx < deduction.child_count;
             ++child_idx) {
          auto child = wmbp_vp_dependency_child(
              deduction, deduction.children[child_idx].item);
          if (child.j != deduction.parent.j) continue;
          const auto found = node_for.find(child);
          if (found == node_for.end()) {
            ++out.missing_child_count;
            continue;
          }
          dependencies[parent_idx].push_back(found->second);
        }
      }
    }
    std::vector<bool> prerequisite(items.size(), false);
    std::vector<std::size_t> pending;
    for (std::size_t idx = 0; idx < items.size(); ++idx) {
      const auto bit = std::uint64_t{1}
          << static_cast<unsigned>(items[idx].nonterminal);
      if ((target_nonterminal_mask & bit) != 0) {
        pending.push_back(idx);
        ++out.target_item_count;
      }
    }
    while (!pending.empty()) {
      const std::size_t current = pending.back();
      pending.pop_back();
      if (prerequisite[current]) continue;
      prerequisite[current] = true;
      for (const auto child : dependencies[current]) pending.push_back(child);
    }
    for (std::size_t idx = 0; idx < items.size(); ++idx) {
      if (!prerequisite[idx]) continue;
      out.schedule.push_back(items[idx]);
      ++out.prerequisite_item_count;
      out.prerequisite_nonterminal_mask |=
          std::uint64_t{1}
          << static_cast<unsigned>(items[idx].nonterminal);
      ++out.prerequisite_count_by_nonterminal[
          static_cast<unsigned>(items[idx].nonterminal)];
      if (items[idx].nonterminal == grammar::NonTerminal::WMBP) {
        ++out.wmbp_prerequisite_count;
      }
    }
    for (std::size_t idx = 0; idx < items.size(); ++idx) {
      if (!prerequisite[idx]) out.schedule.push_back(items[idx]);
    }
  }
  if (!out.feasible(schedule.size())) out.schedule.clear();
  return out;
}

template <class Context>
NonterminalPrerequisiteBarrierSchedule
build_nonterminal_prerequisite_barrier_schedule(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx,
    grammar::NonTerminal target) {
  return build_nonterminal_set_prerequisite_barrier_schedule(
      schedule,
      ctx,
      std::uint64_t{1} << static_cast<unsigned>(target));
}

template <class Context>
VpPrerequisiteBarrierSchedule build_vp_prerequisite_barrier_schedule(
    const std::vector<deductions::ItemKey> &schedule,
    Context &ctx) {
  const auto generic = build_nonterminal_prerequisite_barrier_schedule(
      schedule, ctx, grammar::NonTerminal::VP);
  VpPrerequisiteBarrierSchedule out;
  out.schedule = generic.schedule;
  out.prerequisite_item_count = generic.prerequisite_item_count;
  out.missing_child_count = generic.missing_child_count;
  out.prerequisite_nonterminal_mask =
      generic.prerequisite_nonterminal_mask;
  out.wmbp_prerequisite_count = generic.wmbp_prerequisite_count;
  if (!out.feasible(schedule.size())) out.schedule.clear();
  return out;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMBP_VP_PHASE_BARRIER_PROFILE_HH
