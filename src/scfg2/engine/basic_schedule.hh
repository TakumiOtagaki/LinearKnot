#ifndef CPARTY_SCFG2_ENGINE_BASIC_SCHEDULE_HH
#define CPARTY_SCFG2_ENGINE_BASIC_SCHEDULE_HH

#include "scfg2/deductions/item_key.hh"

#include <array>
#include <vector>

namespace cparty::scfg2::engine {

// Observation of the exact span schedule's fixed child-before-parent phases.
// Sparse derived consumers can use these phases without enumerating n² spans.
inline constexpr auto basic_exact_span_family_order() {
  using N = grammar::NonTerminal;
  return std::array{N::VM, N::V, N::WMv, N::VP, N::BE, N::VPL, N::VPR, N::WMBP, N::WMBW, N::WMB,
                    N::WMp, N::WM, N::WIP, N::WI, N::W};
}

inline std::vector<deductions::ItemKey> build_pkfree_basic_schedule(cand_pos_t n) {
  std::vector<deductions::ItemKey> schedule;
  if (n < 1) {
    return schedule;
  }

  const auto append_primary = [&](grammar::NonTerminal nt, cand_pos_t span, auto &out) {
    for (cand_pos_t i = 1; i + span - 1 <= n; ++i) {
      const cand_pos_t j = i + span - 1;
      out.push_back(deductions::ItemKey{nt, i, j, -1, -1});
    }
  };

  for (cand_pos_t span = 1; span <= n; ++span) {
    // Order same-span families after the families they read from so the
    // single-pass chart update sees all required dependencies.
    append_primary(grammar::NonTerminal::VM, span, schedule);
    append_primary(grammar::NonTerminal::V, span, schedule);
    append_primary(grammar::NonTerminal::WMv, span, schedule);
    append_primary(grammar::NonTerminal::VP, span, schedule);
    append_primary(grammar::NonTerminal::VPL, span, schedule);
    append_primary(grammar::NonTerminal::VPR, span, schedule);
    append_primary(grammar::NonTerminal::WMBP, span, schedule);
    append_primary(grammar::NonTerminal::WMBW, span, schedule);
    append_primary(grammar::NonTerminal::WMB, span, schedule);
    append_primary(grammar::NonTerminal::WMp, span, schedule);
    append_primary(grammar::NonTerminal::WM, span, schedule);
    append_primary(grammar::NonTerminal::WIP, span, schedule);
    append_primary(grammar::NonTerminal::WI, span, schedule);
    append_primary(grammar::NonTerminal::W, span, schedule);
  }

  return schedule;
}

template <class ExactContext>
inline std::vector<deductions::ItemKey> build_pkfree_basic_exact_schedule(
    cand_pos_t n,
    const ExactContext &ctx) {
  std::vector<deductions::ItemKey> schedule;
  if (n < 1) {
    return schedule;
  }

  const auto append_primary = [&](grammar::NonTerminal nt, cand_pos_t span) {
    for (cand_pos_t i = 1; i + span - 1 <= n; ++i) {
      const cand_pos_t j = i + span - 1;
      schedule.push_back(deductions::ItemKey{nt, i, j, -1, -1});
    }
  };

  for (cand_pos_t span = 1; span <= n; ++span) {
    // Keep the exact schedule on the same single-pass dependency order as the
    // basic chart, with BE inserted before the WMBP/WMB families that consume
    // same-span anchored supports.
    for (const auto family : basic_exact_span_family_order()) {
      if (family != grammar::NonTerminal::BE) {
        append_primary(family, span);
        continue;
      }
      for (cand_pos_t i = 1; i + span - 1 <= n; ++i) {
        const cand_pos_t j = i + span - 1;
        if (ctx.pair_partner(i) != j) continue;
        for (cand_pos_t ip = i; ip <= j; ++ip) {
          const cand_pos_t jp = ctx.pair_partner(ip);
          if (jp > ip && jp <= j)
            schedule.push_back(deductions::ItemKey{family, i, j, ip, jp});
        }
      }
    }
  }
  return schedule;
}

template <class ExactContext>
inline std::vector<deductions::ItemKey>
build_pkfree_basic_exact_right_boundary_schedule(cand_pos_t n,
                                                 const ExactContext &ctx) {
  std::vector<deductions::ItemKey> schedule;
  if (n < 1) {
    return schedule;
  }

  const auto append_primary = [&](grammar::NonTerminal nt,
                                  cand_pos_t i,
                                  cand_pos_t j) {
    schedule.push_back(deductions::ItemKey{nt, i, j, -1, -1});
  };

  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = j; i >= 1; --i) {
      // The right boundary is the future beam bucket. Inward same-column
      // children have larger i and are therefore complete before this cell.
      // Within one (i,j), retain the exact schedule's established
      // child-before-parent nonterminal phase order.
      append_primary(grammar::NonTerminal::VM, i, j);
      append_primary(grammar::NonTerminal::V, i, j);
      append_primary(grammar::NonTerminal::WMv, i, j);
      append_primary(grammar::NonTerminal::VP, i, j);
      if (ctx.pair_partner(i) == j) {
        for (cand_pos_t ip = i; ip <= j; ++ip) {
          const cand_pos_t jp = ctx.pair_partner(ip);
          if (jp > ip && jp <= j) {
            schedule.push_back(deductions::ItemKey{
                grammar::NonTerminal::BE, i, j, ip, jp});
          }
        }
      }
      append_primary(grammar::NonTerminal::VPL, i, j);
      append_primary(grammar::NonTerminal::VPR, i, j);
      append_primary(grammar::NonTerminal::WMBP, i, j);
      append_primary(grammar::NonTerminal::WMBW, i, j);
      append_primary(grammar::NonTerminal::WMB, i, j);
      append_primary(grammar::NonTerminal::WMp, i, j);
      append_primary(grammar::NonTerminal::WM, i, j);
      append_primary(grammar::NonTerminal::WIP, i, j);
      append_primary(grammar::NonTerminal::WI, i, j);
      append_primary(grammar::NonTerminal::W, i, j);
    }
  }
  return schedule;
}

template <class ExactContext>
inline std::vector<deductions::ItemKey>
build_pkfree_basic_exact_vp_barrier_schedule(cand_pos_t n,
                                             const ExactContext &ctx) {
  const auto right_boundary =
      build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  std::vector<deductions::ItemKey> schedule;
  schedule.reserve(right_boundary.size());
  std::size_t column_begin = 0;
  while (column_begin < right_boundary.size()) {
    const cand_pos_t j = right_boundary[column_begin].j;
    std::size_t column_end = column_begin;
    while (column_end < right_boundary.size() &&
           right_boundary[column_end].j == j) {
      ++column_end;
    }
    for (std::size_t idx = column_begin; idx < column_end; ++idx) {
      if (right_boundary[idx].nonterminal == grammar::NonTerminal::VP) {
        schedule.push_back(right_boundary[idx]);
      }
    }
    for (std::size_t idx = column_begin; idx < column_end; ++idx) {
      if (right_boundary[idx].nonterminal != grammar::NonTerminal::VP) {
        schedule.push_back(right_boundary[idx]);
      }
    }
    column_begin = column_end;
  }
  return schedule;
}

template <class ExactContext>
inline std::vector<deductions::ItemKey>
build_pkfree_basic_exact_wmbp_barrier_schedule(cand_pos_t n,
                                               const ExactContext &ctx) {
  const auto vp_barrier =
      build_pkfree_basic_exact_vp_barrier_schedule(n, ctx);
  std::vector<deductions::ItemKey> schedule;
  schedule.reserve(vp_barrier.size());
  std::size_t column_begin = 0;
  while (column_begin < vp_barrier.size()) {
    const cand_pos_t j = vp_barrier[column_begin].j;
    std::size_t column_end = column_begin;
    while (column_end < vp_barrier.size() &&
           vp_barrier[column_end].j == j) {
      ++column_end;
    }
    // VP is already the first phase. BE and WMBP form the completed WMBP
    // phase; preserve their established relative order within the phase.
    for (std::size_t idx = column_begin; idx < column_end; ++idx) {
      const auto nt = vp_barrier[idx].nonterminal;
      if (nt == grammar::NonTerminal::VP ||
          nt == grammar::NonTerminal::BE ||
          nt == grammar::NonTerminal::WMBP) {
        schedule.push_back(vp_barrier[idx]);
      }
    }
    for (std::size_t idx = column_begin; idx < column_end; ++idx) {
      const auto nt = vp_barrier[idx].nonterminal;
      if (nt != grammar::NonTerminal::VP &&
          nt != grammar::NonTerminal::BE &&
          nt != grammar::NonTerminal::WMBP) {
        schedule.push_back(vp_barrier[idx]);
      }
    }
    column_begin = column_end;
  }
  return schedule;
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_BASIC_SCHEDULE_HH
