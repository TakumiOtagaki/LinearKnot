#ifndef CPARTY_SCFG2_BEAM_SCAFFOLD_BRANCHING_PATHS_HH
#define CPARTY_SCFG2_BEAM_SCAFFOLD_BRANCHING_PATHS_HH

#include "scfg2/types.hh"

#include <cstddef>
#include <vector>

namespace cparty::scfg2::beam {

struct ScaffoldBranchingPathEntry {
  cand_pos_t ancestor = -1;
  cand_pos_t child = -1;
};

struct ScaffoldBranchingPaths {
  std::vector<std::vector<ScaffoldBranchingPathEntry>> by_endpoint;
  // Endpoints that had more branching ancestors than max_entries. Zero means
  // the cap removed no CASE3 candidate.
  std::size_t truncated_endpoints = 0;
};

// For each unpaired endpoint inside a scaffold pair, the (ancestor, child)
// steps of its path to the virtual root (ancestor 0) whose ancestor has at
// least two pair children, innermost first and at most max_entries long.
//
// A VP_WI_CASE3 source needs B and b to be arms of two distinct pair children
// of LCA(i,j), so a single-child ancestor can never be that LCA. Every CASE3
// consumer combines an entry only with a different child of the same
// ancestor, so dropping single-child steps changes no query result, and they
// do not count toward max_entries.
template <class Context>
ScaffoldBranchingPaths scaffold_branching_paths(
    cand_pos_t n,
    const Context &ctx,
    std::size_t max_entries) {
  const auto size = static_cast<std::size_t>(n + 1);
  std::vector<std::size_t> pair_children(size, 0);
  for (cand_pos_t pos = 1; pos <= n; ++pos) {
    if (ctx.pair_partner(pos) <= pos) continue;
    const cand_pos_t parent = ctx.parent_index(pos);
    if (parent >= 0) ++pair_children[static_cast<std::size_t>(parent)];
  }
  const auto branching = [&](cand_pos_t ancestor) {
    return pair_children[static_cast<std::size_t>(ancestor)] >= 2;
  };

  // up[c]: the nearest proper ancestor of pair c that is branching or the
  // virtual root, with the child of that ancestor on c's path. A parent's
  // open position precedes its children, so one ascending pass suffices.
  std::vector<ScaffoldBranchingPathEntry> up(size);
  for (cand_pos_t c = 1; c <= n; ++c) {
    if (ctx.pair_partner(c) <= c) continue;
    const cand_pos_t parent = ctx.parent_index(c);
    if (parent <= 0 || branching(parent)) {
      up[static_cast<std::size_t>(c)] = ScaffoldBranchingPathEntry{parent, c};
    } else {
      up[static_cast<std::size_t>(c)] = up[static_cast<std::size_t>(parent)];
    }
  }

  ScaffoldBranchingPaths out;
  out.by_endpoint.resize(size);
  for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
    if (!ctx.is_unpaired(endpoint)) continue;
    cand_pos_t child = ctx.parent_index(endpoint);
    if (child <= 0 || ctx.pair_partner(child) <= endpoint) continue;
    auto &path = out.by_endpoint[static_cast<std::size_t>(endpoint)];
    bool finished = false;
    while (path.size() < max_entries) {
      const auto step = up[static_cast<std::size_t>(child)];
      if (step.ancestor < 0) {
        finished = true;
        break;
      }
      if (branching(step.ancestor)) path.push_back(step);
      if (step.ancestor == 0) {
        finished = true;
        break;
      }
      child = step.ancestor;
    }
    if (!finished) {
      // up[] lands on the next branching ancestor or the root, so one more
      // step decides whether the cap dropped a branching ancestor.
      const auto next = up[static_cast<std::size_t>(child)];
      if (next.ancestor >= 0 && branching(next.ancestor)) {
        ++out.truncated_endpoints;
      }
    }
  }
  return out;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_SCAFFOLD_BRANCHING_PATHS_HH
