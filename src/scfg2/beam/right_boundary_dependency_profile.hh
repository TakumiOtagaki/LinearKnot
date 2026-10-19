#ifndef CPARTY_SCFG2_BEAM_RIGHT_BOUNDARY_DEPENDENCY_PROFILE_HH
#define CPARTY_SCFG2_BEAM_RIGHT_BOUNDARY_DEPENDENCY_PROFILE_HH

#include "scfg2/deductions/deduction.hh"

#include <array>
#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

inline constexpr std::size_t nonterminal_count =
    static_cast<std::size_t>(grammar::NonTerminal::BE) + 1;

struct RightBoundaryDependencyProfile {
  std::size_t edge_count = 0;
  std::size_t prior_column_count = 0;
  std::size_t future_column_count = 0;
  std::size_t same_column_inward_count = 0;
  std::size_t same_column_outward_count = 0;
  std::size_t same_span_phase_count = 0;
  std::array<std::array<bool, nonterminal_count>, nonterminal_count>
      phase_edges{};
  // Every child -> parent edge whose two endpoints belong to the same right
  // boundary. Unlike phase_edges this also includes inward-span dependencies
  // (child.i > parent.i).  It is the dependency surface that a whole-column
  // family barrier must respect.
  std::array<std::array<bool, nonterminal_count>, nonterminal_count>
      same_column_family_edges{};

  void observe(const deductions::Deduction &deduction) {
    for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
      auto child = deduction.children[idx].item;
      ++edge_count;
      if (child.j < deduction.parent.j) {
        ++prior_column_count;
      } else if (child.j > deduction.parent.j) {
        ++future_column_count;
      } else if (child.i > deduction.parent.i) {
        ++same_column_inward_count;
      } else if (child.i < deduction.parent.i) {
        ++same_column_outward_count;
      } else {
        ++same_span_phase_count;
        phase_edges[static_cast<std::size_t>(child.nonterminal)]
                   [static_cast<std::size_t>(deduction.parent.nonterminal)] =
            true;
      }
      if (child.j == deduction.parent.j) {
        same_column_family_edges
            [static_cast<std::size_t>(child.nonterminal)]
            [static_cast<std::size_t>(deduction.parent.nonterminal)] = true;
      }
    }
  }

  bool phase_graph_is_acyclic() const {
    std::array<std::size_t, nonterminal_count> indegree{};
    for (std::size_t from = 0; from < nonterminal_count; ++from) {
      for (std::size_t to = 0; to < nonterminal_count; ++to) {
        if (phase_edges[from][to]) ++indegree[to];
      }
    }
    std::vector<std::size_t> ready;
    for (std::size_t node = 0; node < nonterminal_count; ++node) {
      if (indegree[node] == 0) ready.push_back(node);
    }
    std::size_t visited = 0;
    while (!ready.empty()) {
      const std::size_t from = ready.back();
      ready.pop_back();
      ++visited;
      for (std::size_t to = 0; to < nonterminal_count; ++to) {
        if (phase_edges[from][to] && --indegree[to] == 0) {
          ready.push_back(to);
        }
      }
    }
    return visited == nonterminal_count;
  }

  // Self edges are discharged by visiting starts from right to left inside a
  // family phase.  Cross-family cycles, however, prove that no single global
  // family ordering can replace the cell-interleaved exact schedule.
  bool same_column_family_graph_is_acyclic() const {
    std::array<std::size_t, nonterminal_count> indegree{};
    for (std::size_t from = 0; from < nonterminal_count; ++from) {
      for (std::size_t to = 0; to < nonterminal_count; ++to) {
        if (from != to && same_column_family_edges[from][to]) {
          ++indegree[to];
        }
      }
    }
    std::vector<std::size_t> ready;
    for (std::size_t node = 0; node < nonterminal_count; ++node) {
      if (indegree[node] == 0) ready.push_back(node);
    }
    std::size_t visited = 0;
    while (!ready.empty()) {
      const std::size_t from = ready.back();
      ready.pop_back();
      ++visited;
      for (std::size_t to = 0; to < nonterminal_count; ++to) {
        if (from != to && same_column_family_edges[from][to] &&
            --indegree[to] == 0) {
          ready.push_back(to);
        }
      }
    }
    return visited == nonterminal_count;
  }
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_RIGHT_BOUNDARY_DEPENDENCY_PROFILE_HH
