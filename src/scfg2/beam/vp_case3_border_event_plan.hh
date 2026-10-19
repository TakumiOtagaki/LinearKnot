#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_BORDER_EVENT_PLAN_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_BORDER_EVENT_PLAN_HH

#include "scfg2/weights/local_weight.hh"
#include "scfg2/beam/scaffold_endpoint_buckets.hh"

#include "scfg2/beam/scaffold_branching_paths.hh"
#include "scfg2/beam/vp_case3_tail_column_plan.hh"

#include <array>
#include <cstddef>
#include <iterator>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

// Observation-only CASE3 discovery from the scaffold nesting tree.  For an
// unpaired endpoint, every entry records the immediate child band on its path
// below an ancestor.  Two endpoints can have all four CASE3 borders exactly
// when their path entries meet at one ancestor through distinct, ordered
// children.
class ExactVpCase3BorderEventPlan {
 public:
  template <class Context>
  static ExactVpCase3BorderEventPlan discover(cand_pos_t n,
                                               const Context &ctx) {
    ExactVpCase3BorderEventPlan plan;
    using EndpointsByChild = std::map<cand_pos_t, std::vector<cand_pos_t>>;
    std::map<cand_pos_t, EndpointsByChild> endpoints_by_ancestor;

    for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
      if (!ctx.is_unpaired(endpoint)) continue;
      cand_pos_t child = ctx.parent_index(endpoint);
      if (child <= 0 || ctx.pair_partner(child) <= endpoint) continue;
      for (cand_pos_t ancestor = ctx.parent_index(child);
           ancestor >= 0;
           ancestor = ctx.parent_index(child)) {
        endpoints_by_ancestor[ancestor][child].push_back(endpoint);
        if (ancestor == 0) break;
        child = ancestor;
      }
    }

    engine::AllLocalWeightModel<Context> weight(&ctx);
    using SourceIdentity =
        std::tuple<grammar::NonTerminal, cand_pos_t, cand_pos_t>;
    std::set<SourceIdentity> seen;
    for (const auto &[ancestor, by_child] : endpoints_by_ancestor) {
      (void)ancestor;
      for (auto left = by_child.begin(); left != by_child.end(); ++left) {
        for (auto right = std::next(left); right != by_child.end(); ++right) {
          for (const cand_pos_t i : left->second) {
            for (const cand_pos_t j : right->second) {
              ++plan.endpoint_pair_count_;
              if (i >= j) continue;
              const cand_pos_t p = ctx.pair_partner(ctx.parent_index(i));
              const cand_pos_t q = ctx.pair_partner(left->first);
              const cand_pos_t k = right->first;
              const cand_pos_t l = ctx.parent_index(j);
              for (const auto nt : {grammar::NonTerminal::VP}) {
                const deductions::ItemKey parent{nt, i, j, -1, -1};
                const auto generated =
                    deductions::generate_exact_vp_case3_deductions_for_item(
                        parent, ctx);
                for (const auto &deduction : generated) {
                  engine::require_legacy_span_schedule_support(deduction);
                  if (vp_case3_partition(deduction) !=
                          VpCase3Partition::LeftNonempty ||
                      weights::is_zero(weight.local_weight_for(deduction))) {
                    continue;
                  }
                  // Keep the border-event derivation independently checked
                  // against the semantic deduction generator while this is an
                  // observation surface rather than the runtime owner.
                  if (deduction.split.p != p || deduction.split.q != q ||
                      deduction.split.k != k || deduction.split.l != l) {
                    ++plan.border_mismatch_count_;
                    continue;
                  }
                  if (!seen.emplace(nt, i, j).second) {
                    ++plan.duplicate_source_count_;
                    continue;
                  }
                  const VpCase3TailKey tail{j, q, k, l};
                  plan.sources_.push_back(
                      VpCase3PlannedSource{tail, deduction});
                }
              }
            }
          }
        }
      }
    }
    return plan;
  }

  const auto &sources() const { return sources_; }
  std::size_t source_count() const { return sources_.size(); }
  std::size_t endpoint_pair_count() const { return endpoint_pair_count_; }
  std::size_t border_mismatch_count() const {
    return border_mismatch_count_;
  }
  std::size_t duplicate_source_count() const {
    return duplicate_source_count_;
  }

  ExactVpCase3TailColumnPlan into_tail_column_plan() && {
    return ExactVpCase3TailColumnPlan::from_sources(std::move(sources_));
  }

 private:
  std::vector<VpCase3PlannedSource> sources_;
  std::size_t endpoint_pair_count_ = 0;
  std::size_t border_mismatch_count_ = 0;
  std::size_t duplicate_source_count_ = 0;
};

struct VpCase3BorderEventChildQuery {
  std::vector<VpCase3PlannedSource> sources;
  std::size_t endpoint_pair_count = 0;
  std::size_t border_mismatch_count = 0;
  std::size_t duplicate_source_count = 0;
  std::size_t readiness_violation_count = 0;
};

using VpCase3EndpointPathEntry = ScaffoldBranchingPathEntry;

// Scaffold-only inverse index for online production. Unlike the complete
// event plan above, this surface does not form sibling-band cross products at
// construction time. A completed WI child identifies one left endpoint and
// queries only its compatible future right endpoints. Paths keep only
// branching ancestors (see scaffold_branching_paths), and max_path_depth
// counts those.
class ExactVpCase3BorderEventChildIndex {
 public:
  template <class Context>
  static ExactVpCase3BorderEventChildIndex build(cand_pos_t n,
                                                  const Context &ctx,
                                                  std::size_t max_path_depth =
                                                      std::numeric_limits<std::size_t>::max()) {
    ExactVpCase3BorderEventChildIndex index;
    index.endpoints_by_ancestor_ = ScaffoldEndpointBuckets::build(n, ctx);
    index.records_paths_ = max_path_depth > 0;
    auto paths = scaffold_branching_paths(n, ctx, max_path_depth);
    index.paths_by_endpoint_ = std::move(paths.by_endpoint);
    index.truncated_endpoint_count_ = paths.truncated_endpoints;
    for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
      for (const auto &path :
           index.paths_by_endpoint_[static_cast<std::size_t>(endpoint)]) {
        index.endpoints_by_ancestor_.append(path.child, endpoint);
        ++index.path_entry_count_;
      }
    }
    return index;
  }

  std::size_t path_entry_count() const { return path_entry_count_; }
  std::size_t truncated_endpoint_count() const {
    return truncated_endpoint_count_;
  }

  const std::vector<cand_pos_t> &child_bands(cand_pos_t ancestor) const {
    return endpoints_by_ancestor_.children(ancestor);
  }

  const std::vector<VpCase3EndpointPathEntry> &paths_for_endpoint(
      cand_pos_t endpoint) const {
    if (endpoint < 0 ||
        static_cast<std::size_t>(endpoint) >= paths_by_endpoint_.size()) {
      static const std::vector<VpCase3EndpointPathEntry> empty;
      return empty;
    }
    return paths_by_endpoint_[static_cast<std::size_t>(endpoint)];
  }

  // Recognition is structural only. An endpoint without a branching ancestor
  // can never match, but rejecting it here would change which WI items
  // compete for the finite CASE3 producer frontier.
  template <class Context>
  bool recognizes_completed_wi(const deductions::ItemKey &item,
                                const Context &ctx) const {
    if (item.nonterminal != grammar::NonTerminal::WI) return false;
    const cand_pos_t i = item.i - 1;
    const cand_pos_t p = item.j + 1;
    return records_paths_ && i >= 1 &&
           static_cast<std::size_t>(i) < paths_by_endpoint_.size() &&
           ctx.is_unpaired(i) && ctx.parent_index(i) > 0 &&
           ctx.pair_partner(ctx.parent_index(i)) == p;
  }

  template <class Context>
  VpCase3BorderEventChildQuery query_completed_wi(
      const deductions::ItemKey &item,
      const Context &ctx) const {
    VpCase3BorderEventChildQuery out;
    if (!recognizes_completed_wi(item, ctx)) return out;
    const cand_pos_t i = item.i - 1;
    const cand_pos_t p = item.j + 1;

    engine::AllLocalWeightModel<Context> weight(&ctx);
    using SourceIdentity =
        std::tuple<grammar::NonTerminal, cand_pos_t, cand_pos_t>;
    std::set<SourceIdentity> seen;
    for (const auto &path :
         paths_by_endpoint_[static_cast<std::size_t>(i)]) {
      const auto &children = endpoints_by_ancestor_.children(path.ancestor);
      for (auto right = std::upper_bound(children.begin(), children.end(), path.child);
           right != children.end();
           ++right) {
        for (const cand_pos_t j : endpoints_by_ancestor_.bucket(path.ancestor, *right)) {
          ++out.endpoint_pair_count;
          if (i >= j) continue;
          const cand_pos_t q = ctx.pair_partner(path.child);
          const cand_pos_t k = *right;
          const cand_pos_t l = ctx.parent_index(j);
          for (const auto nt : {grammar::NonTerminal::VP}) {
            const deductions::ItemKey parent{nt, i, j, -1, -1};
            const auto generated =
                deductions::generate_exact_vp_case3_deductions_for_item(
                    parent, ctx);
            for (const auto &deduction : generated) {
              engine::require_legacy_span_schedule_support(deduction);
              if (vp_case3_partition(deduction) !=
                      VpCase3Partition::LeftNonempty ||
                  weights::is_zero(weight.local_weight_for(deduction))) {
                continue;
              }
              if (deduction.split.p != p || deduction.split.q != q ||
                  deduction.split.k != k || deduction.split.l != l) {
                ++out.border_mismatch_count;
                continue;
              }
              if (!seen.emplace(nt, i, j).second) {
                ++out.duplicate_source_count;
                continue;
              }
              out.sources.push_back(VpCase3PlannedSource{
                  VpCase3TailKey{j, q, k, l}, deduction});
            }
          }
        }
      }
    }
    return out;
  }

  // Observation-only direction-reversal oracle. Rather than letting one left
  // WI scan all future endpoints, a target endpoint joins against the left WI
  // items that are already known to the caller. The initial implementation
  // deliberately scans that input inventory; a later retained-producer
  // frontier will replace the scan after source parity is established.
  template <class Context>
  VpCase3BorderEventChildQuery query_target_endpoint(
      cand_pos_t j,
      const std::vector<deductions::ItemKey> &left_wi_items,
      const Context &ctx) const {
    VpCase3BorderEventChildQuery out;
    if (j < 1 ||
        static_cast<std::size_t>(j) >= paths_by_endpoint_.size() ||
        !ctx.is_unpaired(j) || paths_by_endpoint_[static_cast<std::size_t>(j)].empty()) {
      return out;
    }

    engine::AllLocalWeightModel<Context> weight(&ctx);
    using SourceIdentity =
        std::tuple<grammar::NonTerminal, cand_pos_t, cand_pos_t>;
    std::set<SourceIdentity> seen;
    for (const auto &item : left_wi_items) {
      if (!recognizes_completed_wi(item, ctx)) continue;
      const cand_pos_t i = item.i - 1;
      if (i >= j) continue;
      const cand_pos_t p = item.j + 1;
      for (const auto &left_path :
           paths_by_endpoint_[static_cast<std::size_t>(i)]) {
        for (const auto &right_path :
             paths_by_endpoint_[static_cast<std::size_t>(j)]) {
          if (left_path.ancestor != right_path.ancestor ||
              left_path.child >= right_path.child) {
            continue;
          }
          ++out.endpoint_pair_count;
          if (item.j >= j) ++out.readiness_violation_count;
          const cand_pos_t q = ctx.pair_partner(left_path.child);
          const cand_pos_t k = right_path.child;
          const cand_pos_t l = ctx.parent_index(j);
          for (const auto nt : {grammar::NonTerminal::VP}) {
            const deductions::ItemKey parent{nt, i, j, -1, -1};
            const auto generated =
                deductions::generate_exact_vp_case3_deductions_for_item(
                    parent, ctx);
            for (const auto &deduction : generated) {
              engine::require_legacy_span_schedule_support(deduction);
              if (vp_case3_partition(deduction) !=
                      VpCase3Partition::LeftNonempty ||
                  weights::is_zero(weight.local_weight_for(deduction))) {
                continue;
              }
              if (deduction.split.p != p || deduction.split.q != q ||
                  deduction.split.k != k || deduction.split.l != l) {
                ++out.border_mismatch_count;
                continue;
              }
              if (!seen.emplace(nt, i, j).second) {
                ++out.duplicate_source_count;
                continue;
              }
              out.sources.push_back(VpCase3PlannedSource{
                  VpCase3TailKey{j, q, k, l}, deduction});
            }
          }
        }
      }
    }
    return out;
  }

 private:
  ScaffoldEndpointBuckets endpoints_by_ancestor_;
  std::vector<std::vector<VpCase3EndpointPathEntry>> paths_by_endpoint_;
  std::size_t path_entry_count_ = 0;
  std::size_t truncated_endpoint_count_ = 0;
  bool records_paths_ = true;
};

template <class Context>
ExactVpCase3TailColumnPlan build_vp_case3_tail_column_plan_from_border_events(
    cand_pos_t n,
    const Context &ctx) {
  return ExactVpCase3BorderEventPlan::discover(n, ctx)
      .into_tail_column_plan();
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_BORDER_EVENT_PLAN_HH
