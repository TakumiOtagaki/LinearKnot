#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_NONLEFT_LAST_CHILD_INDEX_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_NONLEFT_LAST_CHILD_INDEX_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/scaffold_branching_paths.hh"
#include "scfg2/beam/scaffold_endpoint_buckets.hh"
#include "scfg2/beam/virtual_candidate.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <cstddef>
#include <map>
#include <limits>
#include <optional>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase3NonleftIndexStats {
  std::size_t boundary_sources = 0;
  std::size_t other_sources = 0;
  std::size_t zero_child_sources = 0;
  std::size_t child_indexed_sources = 0;
  std::size_t endpoint_pairs = 0;
  std::size_t border_mismatches = 0;
  std::size_t duplicate_sources = 0;
};

// Exact characterization oracle. Construction deliberately scans public-VP
// parents; production may use the resulting readiness contract only after a
// scaffold/event-derived index proves the same source set.
template <class Context>
class ExactVpCase3NonleftLastChildIndex {
 public:
  static ExactVpCase3NonleftLastChildIndex build(cand_pos_t n,
                                                  Context &ctx) {
    ExactVpCase3NonleftLastChildIndex out;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    for (cand_pos_t j = 1; j <= n; ++j) {
      for (cand_pos_t i = 1; i <= j; ++i) {
        const deductions::ItemKey parent{
            grammar::NonTerminal::VP, i, j, -1, -1};
        for (auto source :
             deductions::generate_exact_vp_case3_deductions_for_item(
                 parent, ctx)) {
          engine::require_legacy_span_schedule_support(source);
          if (weights::is_zero(weight.local_weight_for(source))) {
            continue;
          }
          const auto partition = vp_case3_partition(source);
          if (!partition.has_value() ||
              *partition == VpCase3Partition::LeftNonempty) {
            continue;
          }
          if (*partition == VpCase3Partition::Boundary) {
            ++out.stats_.boundary_sources;
          } else {
            ++out.stats_.other_sources;
          }
          if (source.child_count == 0) {
            out.zero_child_by_column_[j].push_back(std::move(source));
            ++out.stats_.zero_child_sources;
          } else {
            const auto trigger = source.children[source.child_count - 1].item;
            out.by_last_child_[trigger].push_back(std::move(source));
            ++out.stats_.child_indexed_sources;
          }
        }
      }
    }
    return out;
  }

  const std::vector<deductions::Deduction> &sources_for_last_child(
      const deductions::ItemKey &child) const {
    const auto found = by_last_child_.find(child);
    if (found == by_last_child_.end()) {
      static const std::vector<deductions::Deduction> empty;
      return empty;
    }
    return found->second;
  }

  const std::vector<deductions::Deduction> &zero_child_sources(
      cand_pos_t j) const {
    const auto found = zero_child_by_column_.find(j);
    if (found == zero_child_by_column_.end()) {
      static const std::vector<deductions::Deduction> empty;
      return empty;
    }
    return found->second;
  }

  const VpCase3NonleftIndexStats &stats() const { return stats_; }

 private:
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      by_last_child_;
  std::unordered_map<cand_pos_t, std::vector<deductions::Deduction>>
      zero_child_by_column_;
  VpCase3NonleftIndexStats stats_{};
};

// Scaffold-tree characterization. Unlike the parent-scan oracle above, this
// discovers a non-left CASE3 parent from its left boundary endpoint and a
// later sibling-band endpoint. It still forms the complete endpoint product;
// retained-child online production is a later ownership step.
template <class Context>
class ScaffoldVpCase3NonleftLastChildIndex {
 public:
  static ScaffoldVpCase3NonleftLastChildIndex build(cand_pos_t n,
                                                     Context &ctx) {
    ScaffoldVpCase3NonleftLastChildIndex out;
    struct PathEntry {
      cand_pos_t ancestor = -1;
      cand_pos_t child = -1;
    };
    using EndpointsByChild = std::map<cand_pos_t, std::vector<cand_pos_t>>;
    std::map<cand_pos_t, EndpointsByChild> endpoints_by_ancestor;
    std::vector<std::vector<PathEntry>> paths_by_endpoint(
        static_cast<std::size_t>(n + 1));
    for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
      if (!ctx.is_unpaired(endpoint)) continue;
      cand_pos_t child = ctx.parent_index(endpoint);
      if (child <= 0 || ctx.pair_partner(child) <= endpoint) continue;
      for (cand_pos_t ancestor = ctx.parent_index(child);
           ancestor >= 0;
           ancestor = ctx.parent_index(child)) {
        endpoints_by_ancestor[ancestor][child].push_back(endpoint);
        paths_by_endpoint[static_cast<std::size_t>(endpoint)].push_back(
            PathEntry{ancestor, child});
        if (ancestor == 0) break;
        child = ancestor;
      }
    }

    engine::AllLocalWeightModel<Context> weight(&ctx);
    std::set<std::pair<cand_pos_t, cand_pos_t>> seen;
    for (cand_pos_t i = 1; i <= n; ++i) {
      if (!ctx.is_unpaired(i) || ctx.parent_index(i) <= 0 ||
          ctx.pair_partner(ctx.parent_index(i)) != i + 1) {
        continue;
      }
      for (const auto &path :
           paths_by_endpoint[static_cast<std::size_t>(i)]) {
        const auto ancestor = endpoints_by_ancestor.find(path.ancestor);
        if (ancestor == endpoints_by_ancestor.end()) continue;
        for (auto right = ancestor->second.upper_bound(path.child);
             right != ancestor->second.end();
             ++right) {
          for (const cand_pos_t j : right->second) {
            ++out.stats_.endpoint_pairs;
            if (i >= j) continue;
            const cand_pos_t p = ctx.pair_partner(ctx.parent_index(i));
            const cand_pos_t q = ctx.pair_partner(path.child);
            const cand_pos_t k = right->first;
            const cand_pos_t l = ctx.parent_index(j);
            const deductions::ItemKey parent{
                grammar::NonTerminal::VP, i, j, -1, -1};
            for (auto source :
                 deductions::generate_exact_vp_case3_deductions_for_item(
                     parent, ctx)) {
              engine::require_legacy_span_schedule_support(source);
              const auto partition = vp_case3_partition(source);
              if (!partition.has_value() ||
                  *partition == VpCase3Partition::LeftNonempty ||
                  weights::is_zero(weight.local_weight_for(source))) {
                continue;
              }
              if (source.split.p != p || source.split.q != q ||
                  source.split.k != k || source.split.l != l) {
                ++out.stats_.border_mismatches;
                continue;
              }
              if (!seen.emplace(i, j).second) {
                ++out.stats_.duplicate_sources;
                continue;
              }
              if (*partition == VpCase3Partition::Boundary) {
                ++out.stats_.boundary_sources;
              } else {
                ++out.stats_.other_sources;
              }
              if (source.child_count == 0) {
                out.zero_child_by_column_[j].push_back(std::move(source));
                ++out.stats_.zero_child_sources;
              } else {
                const auto trigger =
                    source.children[source.child_count - 1].item;
                out.by_last_child_[trigger].push_back(std::move(source));
                ++out.stats_.child_indexed_sources;
              }
            }
          }
        }
      }
    }
    return out;
  }

  const std::vector<deductions::Deduction> &sources_for_last_child(
      const deductions::ItemKey &child) const {
    const auto found = by_last_child_.find(child);
    if (found == by_last_child_.end()) {
      static const std::vector<deductions::Deduction> empty;
      return empty;
    }
    return found->second;
  }

  const std::vector<deductions::Deduction> &zero_child_sources(
      cand_pos_t j) const {
    const auto found = zero_child_by_column_.find(j);
    if (found == zero_child_by_column_.end()) {
      static const std::vector<deductions::Deduction> empty;
      return empty;
    }
    return found->second;
  }

  const VpCase3NonleftIndexStats &stats() const { return stats_; }

 private:
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      by_last_child_;
  std::unordered_map<cand_pos_t, std::vector<deductions::Deduction>>
      zero_child_by_column_;
  VpCase3NonleftIndexStats stats_{};
};

struct VpCase3NonleftChildEventQuery {
  std::vector<deductions::Deduction> sources;
  std::size_t border_mismatches = 0;
  std::size_t duplicate_sources = 0;
};

struct VpCase3NonleftPrefixGeometry {
  cand_pos_t i = -1;
  cand_pos_t p = -1;
  cand_pos_t q = -1;
  cand_pos_t k = -1;
  cand_pos_t target_close = -1;
};

inline bool operator<(const VpCase3NonleftPrefixGeometry &lhs,
                      const VpCase3NonleftPrefixGeometry &rhs) {
  return std::tie(lhs.i, lhs.p, lhs.q, lhs.k, lhs.target_close) <
         std::tie(rhs.i, rhs.p, rhs.q, rhs.k, rhs.target_close);
}

// Online inverse geometry. Construction stores only endpoint paths; it does
// not form left/right endpoint products. A retained last WI child, or an
// explicit childless target-column event, asks for the compatible sources.
// Paths keep only branching ancestors (see scaffold_branching_paths), and
// max_path_depth counts those.
class VpCase3NonleftChildEventIndex {
 public:
  template <class Context>
  static VpCase3NonleftChildEventIndex build(cand_pos_t n,
                                              const Context &ctx,
                                              std::size_t max_path_depth =
                                                  std::numeric_limits<std::size_t>::max()) {
    VpCase3NonleftChildEventIndex out;
    out.n_ = n;
    out.endpoints_by_ancestor_ = ScaffoldEndpointBuckets::build(n, ctx);
    out.prefix_endpoints_ = ScaffoldEndpointBuckets::build(n, ctx);
    out.empty_right_endpoints_ = ScaffoldEndpointBuckets::build(n, ctx);
    auto paths = scaffold_branching_paths(n, ctx, max_path_depth);
    out.paths_by_endpoint_ = std::move(paths.by_endpoint);
    out.truncated_endpoint_count_ = paths.truncated_endpoints;
    for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
      const bool empty_right = ctx.parent_index(endpoint) == endpoint - 1;
      const cand_pos_t parent = ctx.parent_index(endpoint);
      const bool prefix_endpoint =
          parent > 0 && ctx.pair_partner(parent) == endpoint + 1;
      for (const auto &path :
           out.paths_by_endpoint_[static_cast<std::size_t>(endpoint)]) {
        out.endpoints_by_ancestor_.append(path.child, endpoint);
        ++out.path_entry_count_;
        if (prefix_endpoint) {
          out.prefix_endpoints_.append(path.child, endpoint);
        }
        if (empty_right) {
          out.empty_right_endpoints_.append(path.child, endpoint);
        }
      }
    }
    return out;
  }

  std::size_t path_entry_count() const { return path_entry_count_; }
  std::size_t truncated_endpoint_count() const {
    return truncated_endpoint_count_;
  }

  template <class Context>
  VpCase3NonleftChildEventQuery query_completed_wi(
      const deductions::ItemKey &item,
      Context &ctx) const {
    VpCase3NonleftChildEventQuery out;
    if (item.nonterminal != grammar::NonTerminal::WI) return out;
    std::set<std::pair<cand_pos_t, cand_pos_t>> seen;

    // The completed WI is the nonempty right factor WI(l+1,j-1).
    const cand_pos_t l = item.i - 1;
    const cand_pos_t j = item.j + 1;
    if (j >= 1 && j <= n_ && ctx.is_unpaired(j) &&
        ctx.parent_index(j) == l) {
      for (const auto &path :
           paths_by_endpoint_[static_cast<std::size_t>(j)]) {
        for (const auto child : endpoints_by_ancestor_.children(path.ancestor)) {
          if (child >= path.child) break;
          const cand_pos_t q = ctx.pair_partner(child);
          for (const cand_pos_t i : endpoints_by_ancestor_.bucket(path.ancestor, child)) {
            const cand_pos_t p = ctx.pair_partner(ctx.parent_index(i));
            if (p != i + 1) continue;
            append_verified_source(i,
                                   j,
                                   p,
                                   q,
                                   path.child,
                                   l,
                                   item,
                                   ctx,
                                   seen,
                                   out);
          }
        }
      }
    }

    // The completed WI is the nonempty middle factor WI(q+1,k-1),
    // while the right factor is empty. The right endpoint may be inside a
    // descendant band below sibling k, so use its recorded ancestor path.
    const cand_pos_t q = item.i - 1;
    const cand_pos_t k = item.j + 1;
    const cand_pos_t left_child = ctx.pair_partner(q);
    if (left_child > 0 && ctx.pair_partner(left_child) == q) {
      const cand_pos_t ancestor = ctx.parent_index(left_child);
      const auto &left = endpoints_by_ancestor_.bucket(ancestor, left_child);
      const auto &right = empty_right_endpoints_.bucket(ancestor, k);
      for (const cand_pos_t i : left) {
        const cand_pos_t p = ctx.pair_partner(ctx.parent_index(i));
        if (p != i + 1) continue;
        for (const cand_pos_t right_j : right) {
          append_verified_source(i,
                                 right_j,
                                 p,
                                 q,
                                 k,
                                 right_j - 1,
                                 item,
                                 ctx,
                                 seen,
                                 out);
        }
      }
    }
    return out;
  }

  template <class Context>
  VpCase3NonleftChildEventQuery query_zero_child_column(
      cand_pos_t j,
      Context &ctx) const {
    VpCase3NonleftChildEventQuery out;
    if (j < 1 || j > n_ || !ctx.is_unpaired(j) ||
        ctx.parent_index(j) != j - 1) {
      return out;
    }
    std::set<std::pair<cand_pos_t, cand_pos_t>> seen;
    for (const auto &path :
         paths_by_endpoint_[static_cast<std::size_t>(j)]) {
      for (const auto child : endpoints_by_ancestor_.children(path.ancestor)) {
        if (child >= path.child) break;
        const cand_pos_t q = ctx.pair_partner(child);
        if (q + 1 != path.child) continue;
        for (const cand_pos_t i : endpoints_by_ancestor_.bucket(path.ancestor, child)) {
          const cand_pos_t p = ctx.pair_partner(ctx.parent_index(i));
          if (p != i + 1) continue;
          append_verified_source(i,
                                 j,
                                 p,
                                 q,
                                 path.child,
                                 j - 1,
                                 std::nullopt,
                                 ctx,
                                 seen,
                                 out);
        }
      }
    }
    return out;
  }

  template <class Context>
  std::vector<VpCase3NonleftPrefixGeometry> query_middle_prefixes(
      const deductions::ItemKey &middle,
      const Context &ctx,
      std::size_t limit = std::numeric_limits<std::size_t>::max()) const {
    if (middle.nonterminal != grammar::NonTerminal::WI) return {};
    return query_prefix_geometry(
        middle.i - 1, middle.j + 1, ctx, limit);
  }

  template <class Context>
  std::vector<VpCase3NonleftPrefixGeometry> query_empty_middle_prefixes(
      cand_pos_t k,
      const Context &ctx,
      std::size_t limit = std::numeric_limits<std::size_t>::max()) const {
    return query_prefix_geometry(k - 1, k, ctx, limit);
  }

 private:
  using PathEntry = ScaffoldBranchingPathEntry;
  template <class Context>
  std::vector<VpCase3NonleftPrefixGeometry> query_prefix_geometry(
      cand_pos_t q,
      cand_pos_t k,
      const Context &ctx,
      std::size_t limit) const {
    std::vector<VpCase3NonleftPrefixGeometry> out;
    if (q < 1 || k < 1 || k > n_) return out;
    const cand_pos_t left_child = ctx.pair_partner(q);
    const cand_pos_t target_close = ctx.pair_partner(k);
    if (left_child <= 0 || ctx.pair_partner(left_child) != q ||
        target_close <= k ||
        ctx.parent_index(left_child) != ctx.parent_index(k) ||
        left_child >= k) {
      return out;
    }
    const cand_pos_t ancestor = ctx.parent_index(left_child);
    // The endpoint-only predicate is evaluated once during construction.
    // Every entry in this bucket is eligible, so a width-b query examines
    // at most b endpoints, including on a long mostly-unpaired scaffold.
    for (const cand_pos_t i : prefix_endpoints_.bucket(ancestor, left_child)) {
      if (out.size() >= limit) break;
      out.push_back(VpCase3NonleftPrefixGeometry{
          i, i + 1, q, k, target_close});
    }
    return out;
  }

  template <class Context>
  static void append_verified_source(
      cand_pos_t i,
      cand_pos_t j,
      cand_pos_t p,
      cand_pos_t q,
      cand_pos_t k,
      cand_pos_t l,
      const std::optional<deductions::ItemKey> &trigger,
      Context &ctx,
      std::set<std::pair<cand_pos_t, cand_pos_t>> &seen,
      VpCase3NonleftChildEventQuery &out) {
    if (i < 1 || i >= j) return;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    const deductions::ItemKey parent{
        grammar::NonTerminal::VP, i, j, -1, -1};
    for (auto source :
         deductions::generate_exact_vp_case3_deductions_for_item(parent,
                                                                  ctx)) {
      engine::require_legacy_span_schedule_support(source);
      const auto partition = vp_case3_partition(source);
      if (!partition.has_value() ||
          *partition == VpCase3Partition::LeftNonempty ||
          weights::is_zero(weight.local_weight_for(source))) {
        continue;
      }
      if (source.split.p != p || source.split.q != q ||
          source.split.k != k || source.split.l != l) {
        ++out.border_mismatches;
        continue;
      }
      if (trigger.has_value()) {
        if (source.child_count == 0 ||
            !deductions::ItemKeyEq{}(
                source.children[source.child_count - 1].item, *trigger)) {
          continue;
        }
      } else if (source.child_count != 0) {
        continue;
      }
      if (!seen.emplace(i, j).second) {
        ++out.duplicate_sources;
        continue;
      }
      out.sources.push_back(std::move(source));
    }
  }

  cand_pos_t n_ = 0;
  ScaffoldEndpointBuckets endpoints_by_ancestor_, prefix_endpoints_, empty_right_endpoints_;
  std::vector<std::vector<PathEntry>> paths_by_endpoint_;
  std::size_t path_entry_count_ = 0;
  std::size_t truncated_endpoint_count_ = 0;
};

}  // namespace cparty::scfg2::beam

#endif
