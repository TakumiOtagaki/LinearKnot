#ifndef CPARTY_SCFG2_BEAM_VP_CASE12_LAST_CHILD_INDEX_HH
#define CPARTY_SCFG2_BEAM_VP_CASE12_LAST_CHILD_INDEX_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <limits>
#include <set>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

inline constexpr std::size_t vp_case12_rule_index(grammar::RuleId rule) {
  return rule == grammar::RuleId::VP_WI_CASE1 ? 0 : 1;
}

struct VpCase12LastChildStats {
  std::array<std::size_t, 2> sources{};
  std::array<std::array<std::size_t, 3>, 2> child_counts{};
  std::array<std::size_t, 2> left_empty{};
  std::array<std::size_t, 2> right_empty{};
  std::array<std::size_t, 2> zero_child{};
  std::array<std::size_t, 2> child_indexed{};
  std::size_t endpoint_candidates = 0;
  std::size_t border_mismatches = 0;
  std::size_t duplicate_sources = 0;
};

// Parent-scan characterization oracle. Production ownership must not use this
// index until an independent scaffold/retained-child inverse surface proves
// the same source set.
template <class Context>
class ExactVpCase12LastChildIndex {
 public:
  static ExactVpCase12LastChildIndex build(cand_pos_t n, Context &ctx) {
    ExactVpCase12LastChildIndex out;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    for (cand_pos_t j = 1; j <= n; ++j) {
      for (cand_pos_t i = 1; i <= j; ++i) {
        const deductions::ItemKey parent{
            grammar::NonTerminal::VP, i, j, -1, -1};
        for (const auto rule : {grammar::RuleId::VP_WI_CASE1,
                                grammar::RuleId::VP_WI_CASE2}) {
          auto generated = deductions::generate_exact_basic_deductions(
              parent, rule, ctx);
          for (auto &source : generated) {
            engine::require_legacy_span_schedule_support(source);
            if (weights::is_zero(weight.local_weight_for(source))) {
              continue;
            }
            const auto idx = vp_case12_rule_index(rule);
            ++out.stats_.sources[idx];
            ++out.stats_.child_counts[idx][source.child_count];
            if (parent.i + 1 > source.split.p - 1) {
              ++out.stats_.left_empty[idx];
            }
            if (source.split.q + 1 > parent.j - 1) {
              ++out.stats_.right_empty[idx];
            }
            if (source.child_count == 0) {
              out.zero_child_by_column_[j].push_back(std::move(source));
              ++out.stats_.zero_child[idx];
            } else {
              const auto trigger =
                  source.children[source.child_count - 1].item;
              out.by_last_child_[trigger].push_back(std::move(source));
              ++out.stats_.child_indexed[idx];
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

  const VpCase12LastChildStats &stats() const { return stats_; }

 private:
  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      by_last_child_;
  std::unordered_map<cand_pos_t, std::vector<deductions::Deduction>>
      zero_child_by_column_;
  VpCase12LastChildStats stats_{};
};

struct VpCase12PrefixGeometry {
  grammar::RuleId rule = grammar::RuleId::VP_WI_CASE1;
  cand_pos_t i = -1;
  cand_pos_t p = -1;
  cand_pos_t q = -1;
  cand_pos_t target_begin = -1;
  cand_pos_t target_end = -1;
};

inline bool operator==(const VpCase12PrefixGeometry &lhs,
                       const VpCase12PrefixGeometry &rhs) {
  return lhs.rule == rhs.rule && lhs.i == rhs.i && lhs.p == rhs.p &&
      lhs.q == rhs.q && lhs.target_begin == rhs.target_begin &&
      lhs.target_end == rhs.target_end;
}

inline bool operator<(const VpCase12PrefixGeometry &lhs,
                      const VpCase12PrefixGeometry &rhs) {
  return std::tie(lhs.rule,
                  lhs.i,
                  lhs.p,
                  lhs.q,
                  lhs.target_begin,
                  lhs.target_end) <
         std::tie(rhs.rule,
                  rhs.i,
                  rhs.p,
                  rhs.q,
                  rhs.target_begin,
                  rhs.target_end);
}

struct VpCase12PrefixPlanStats {
  std::array<std::size_t, 2> sources{};
  std::array<std::size_t, 2> prefixes{};
};

// Parent-scan factorization oracle.  This proves the interval prefix shape and
// its algebraic reconstruction, but must not become production ownership: a
// retained-child event surface has to create the same prefixes independently.
template <class Context>
class ExactVpCase12PrefixPlan {
 public:
  static ExactVpCase12PrefixPlan build(cand_pos_t n, Context &ctx) {
    ExactVpCase12PrefixPlan out;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    std::set<VpCase12PrefixGeometry> unique;
    for (cand_pos_t j = 1; j <= n; ++j) {
      for (cand_pos_t i = 1; i <= j; ++i) {
        const deductions::ItemKey parent{
            grammar::NonTerminal::VP, i, j, -1, -1};
        for (const auto rule : {grammar::RuleId::VP_WI_CASE1,
                                grammar::RuleId::VP_WI_CASE2}) {
          for (const auto &source :
               deductions::generate_exact_basic_deductions(
                   parent, rule, ctx)) {
            engine::require_legacy_span_schedule_support(source);
            if (weights::is_zero(weight.local_weight_for(source))) {
              continue;
            }
            const auto idx = vp_case12_rule_index(rule);
            ++out.stats_.sources[idx];
            VpCase12PrefixGeometry geometry{};
            geometry.rule = rule;
            geometry.i = i;
            geometry.p = source.split.p;
            if (rule == grammar::RuleId::VP_WI_CASE1) {
              geometry.q = source.split.q;
              geometry.target_begin = source.split.q + 1;
              const cand_pos_t child = ctx.pair_partner(source.split.q);
              const cand_pos_t ancestor = ctx.parent_index(child);
              geometry.target_end = ancestor == 0
                  ? n
                  : ctx.pair_partner(ancestor) - 1;
            } else {
              geometry.q = -1;
              geometry.target_begin = source.split.p;
              geometry.target_end = ctx.pair_partner(source.split.p);
            }
            unique.insert(geometry);
          }
        }
      }
    }
    out.prefixes_.assign(unique.begin(), unique.end());
    for (const auto &prefix : out.prefixes_) {
      ++out.stats_.prefixes[vp_case12_rule_index(prefix.rule)];
    }
    return out;
  }

  const std::vector<VpCase12PrefixGeometry> &prefixes() const {
    return prefixes_;
  }

  const VpCase12PrefixPlanStats &stats() const { return stats_; }

 private:
  std::vector<VpCase12PrefixGeometry> prefixes_;
  VpCase12PrefixPlanStats stats_{};
};

// Independent scaffold-tree oracle for the same prefix geometries.  Unlike
// ExactVpCase12PrefixPlan, this does not invoke a deduction generator or a
// local-weight model.
template <class Context>
class ScaffoldVpCase12PrefixPlan {
 public:
  static ScaffoldVpCase12PrefixPlan build(cand_pos_t n, Context &ctx) {
    ScaffoldVpCase12PrefixPlan out;
    struct PathEntry {
      cand_pos_t ancestor = -1;
      cand_pos_t child = -1;
    };
    std::vector<std::vector<PathEntry>> paths(
        static_cast<std::size_t>(n + 1));
    using EndpointsByBranch =
        std::map<cand_pos_t, std::vector<cand_pos_t>>;
    std::map<cand_pos_t, EndpointsByBranch> endpoints_by_ancestor;
    for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
      if (!ctx.is_unpaired(endpoint)) continue;
      const cand_pos_t immediate = ctx.parent_index(endpoint);
      if (immediate >= 0) {
        endpoints_by_ancestor[immediate][endpoint].push_back(endpoint);
      }
      if (immediate <= 0 || ctx.pair_partner(immediate) <= endpoint) continue;
      cand_pos_t child = immediate;
      for (cand_pos_t ancestor = ctx.parent_index(child);
           ancestor >= 0;
           ancestor = ctx.parent_index(child)) {
        paths[static_cast<std::size_t>(endpoint)].push_back(
            PathEntry{ancestor, child});
        endpoints_by_ancestor[ancestor][child].push_back(endpoint);
        if (ancestor == 0) break;
        child = ancestor;
      }
    }

    std::set<VpCase12PrefixGeometry> unique;
    for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
      if (!ctx.is_unpaired(endpoint)) continue;
      const cand_pos_t immediate = ctx.parent_index(endpoint);
      for (const auto &entry : paths[static_cast<std::size_t>(endpoint)]) {
        const auto branches = endpoints_by_ancestor.find(entry.ancestor);
        if (branches == endpoints_by_ancestor.end()) continue;

        // CASE1: endpoint is inside the left child; the existence of a later
        // nonempty branch makes its ancestor suffix a useful target interval.
        bool has_later_direct_endpoint = false;
        for (auto branch = branches->second.upper_bound(entry.child);
             branch != branches->second.end() &&
                 !has_later_direct_endpoint;
             ++branch) {
          has_later_direct_endpoint = std::any_of(
              branch->second.begin(),
              branch->second.end(),
              [&](cand_pos_t right_endpoint) {
                return ctx.parent_index(right_endpoint) == entry.ancestor;
              });
        }
        if (has_later_direct_endpoint) {
          unique.insert(VpCase12PrefixGeometry{
              grammar::RuleId::VP_WI_CASE1,
              endpoint,
              ctx.pair_partner(immediate),
              ctx.pair_partner(entry.child),
              ctx.pair_partner(entry.child) + 1,
              entry.ancestor == 0
                  ? n
                  : ctx.pair_partner(entry.ancestor) - 1});
        }

        // CASE2: endpoint supplies evidence that this right child has a
        // possible final endpoint.  Every endpoint in an earlier ordered
        // branch supplies one reusable (i,p) prefix.
        for (auto branch = branches->second.begin();
             branch != branches->second.lower_bound(entry.child);
             ++branch) {
          for (const cand_pos_t left_endpoint : branch->second) {
            if (ctx.parent_index(left_endpoint) != entry.ancestor) continue;
            unique.insert(VpCase12PrefixGeometry{
                grammar::RuleId::VP_WI_CASE2,
                left_endpoint,
                entry.child,
                -1,
                entry.child,
                ctx.pair_partner(entry.child)});
          }
        }
      }
    }
    out.prefixes_.assign(unique.begin(), unique.end());
    for (const auto &prefix : out.prefixes_) {
      ++out.stats_.prefixes[vp_case12_rule_index(prefix.rule)];
    }
    return out;
  }

  const std::vector<VpCase12PrefixGeometry> &prefixes() const {
    return prefixes_;
  }

  const VpCase12PrefixPlanStats &stats() const { return stats_; }

 private:
  std::vector<VpCase12PrefixGeometry> prefixes_;
  VpCase12PrefixPlanStats stats_{};
};

struct VpCase12LeftEventStats {
  std::array<std::size_t, 2> nonempty{};
  std::array<std::size_t, 2> empty{};
};

// Trigger-shape oracle for prefix production.  It deliberately indexes the
// independently derived scaffold plan rather than deductions.  The next
// production slice may replace its eager build with online geometry queries
// without changing the completed-WI / empty-column event contract.
template <class Context>
class VpCase12LeftEventIndex {
 public:
  static VpCase12LeftEventIndex build(cand_pos_t n, Context &ctx) {
    VpCase12LeftEventIndex out;
    const auto plan = ScaffoldVpCase12PrefixPlan<Context>::build(n, ctx);
    for (const auto &prefix : plan.prefixes()) {
      const auto idx = vp_case12_rule_index(prefix.rule);
      if (prefix.i + 1 <= prefix.p - 1) {
        const deductions::ItemKey child{
            grammar::NonTerminal::WI,
            prefix.i + 1,
            prefix.p - 1,
            -1,
            -1};
        out.by_completed_left_[child].push_back(prefix);
        ++out.stats_.nonempty[idx];
      } else {
        out.empty_by_activation_[prefix.target_begin].push_back(prefix);
        ++out.stats_.empty[idx];
      }
    }
    return out;
  }

  const std::vector<VpCase12PrefixGeometry> &query_completed_left(
      const deductions::ItemKey &child) const {
    const auto found = by_completed_left_.find(child);
    if (found == by_completed_left_.end()) return empty();
    return found->second;
  }

  const std::vector<VpCase12PrefixGeometry> &query_empty_activation(
      cand_pos_t column) const {
    const auto found = empty_by_activation_.find(column);
    if (found == empty_by_activation_.end()) return empty();
    return found->second;
  }

  const VpCase12LeftEventStats &stats() const { return stats_; }

 private:
  static const std::vector<VpCase12PrefixGeometry> &empty() {
    static const std::vector<VpCase12PrefixGeometry> value;
    return value;
  }

  std::unordered_map<deductions::ItemKey,
                     std::vector<VpCase12PrefixGeometry>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      by_completed_left_;
  std::unordered_map<cand_pos_t, std::vector<VpCase12PrefixGeometry>>
      empty_by_activation_;
  VpCase12LeftEventStats stats_{};
};

// Direct nonempty-left geometry query.  No prefix product or deduction source
// is materialized at construction.  CASE1 walks only the completed child's
// ancestor chain and may stop inside that walk at finite width; CASE2 has at
// most one geometry for a fixed (i,p) left span.
class VpCase12DirectNonemptyLeftQuery {
 public:
  template <class Context>
  static VpCase12DirectNonemptyLeftQuery build(cand_pos_t n,
                                                Context &ctx) {
    VpCase12DirectNonemptyLeftQuery out;
    out.n_ = n;
    out.last_direct_endpoint_.assign(
        static_cast<std::size_t>(n + 1), -1);
    out.unpaired_count_.assign(static_cast<std::size_t>(n + 1), 0);
    out.next_case1_child_.assign(static_cast<std::size_t>(n + 1), -1);
    for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
      out.unpaired_count_[static_cast<std::size_t>(endpoint)] =
          out.unpaired_count_[static_cast<std::size_t>(endpoint - 1)] +
          (ctx.is_unpaired(endpoint) ? 1 : 0);
      if (!ctx.is_unpaired(endpoint)) continue;
      const cand_pos_t ancestor = ctx.parent_index(endpoint);
      if (ancestor >= 0) {
        out.last_direct_endpoint_[static_cast<std::size_t>(ancestor)] =
            endpoint;
      }
    }
    // Parent band openings precede child band openings, so ascending order
    // makes the nearest valid ancestor jump available in one pass.
    for (cand_pos_t child = 1; child <= n; ++child) {
      if (ctx.pair_partner(child) <= child) continue;
      const cand_pos_t ancestor = ctx.parent_index(child);
      if (ancestor >= 0 &&
          out.last_direct_endpoint_[static_cast<std::size_t>(ancestor)] >
              child) {
        out.next_case1_child_[static_cast<std::size_t>(child)] = child;
      } else if (ancestor > 0) {
        out.next_case1_child_[static_cast<std::size_t>(child)] =
            out.next_case1_child_[static_cast<std::size_t>(ancestor)];
      }
    }
    return out;
  }

  template <class Context>
  std::vector<VpCase12PrefixGeometry> query_completed_left(
      const deductions::ItemKey &left,
      Context &ctx,
      std::size_t width = std::numeric_limits<std::size_t>::max()) const {
    std::vector<VpCase12PrefixGeometry> out;
    if (width == 0 || left.nonterminal != grammar::NonTerminal::WI) {
      return out;
    }
    const cand_pos_t i = left.i - 1;
    const cand_pos_t p = left.j + 1;
    if (i < 1 || p > n_ || !ctx.is_unpaired(i)) return out;

    const cand_pos_t immediate = ctx.parent_index(i);
    if (immediate > 0 && ctx.pair_partner(immediate) == p) {
      cand_pos_t child =
          next_case1_child_[static_cast<std::size_t>(immediate)];
      while (child > 0 && out.size() < width) {
        const cand_pos_t ancestor = ctx.parent_index(child);
        out.push_back(VpCase12PrefixGeometry{
            grammar::RuleId::VP_WI_CASE1,
            i,
            p,
            ctx.pair_partner(child),
            ctx.pair_partner(child) + 1,
            ancestor == 0 ? n_ : ctx.pair_partner(ancestor) - 1});
        child = ancestor > 0
            ? next_case1_child_[static_cast<std::size_t>(ancestor)]
            : -1;
      }
    }

    const cand_pos_t p_close = ctx.pair_partner(p);
    if (p_close > p && immediate >= 0 &&
        ctx.parent_index(p) == immediate && i < p &&
        unpaired_in_open_band(p, p_close)) {
      out.push_back(VpCase12PrefixGeometry{
          grammar::RuleId::VP_WI_CASE2,
          i,
          p,
          -1,
          p,
          p_close});
    }
    std::sort(out.begin(), out.end());
    return out;
  }

 private:
  bool unpaired_in_open_band(cand_pos_t open, cand_pos_t close) const {
    if (open + 1 > close - 1) return false;
    return unpaired_count_[static_cast<std::size_t>(close - 1)] >
        unpaired_count_[static_cast<std::size_t>(open)];
  }

  cand_pos_t n_ = 0;
  std::vector<cand_pos_t> last_direct_endpoint_;
  std::vector<cand_pos_t> unpaired_count_;
  std::vector<cand_pos_t> next_case1_child_;
};

// Direct empty-left geometry query.  CASE1 summarizes qualifying endpoints
// bottom-up in each scaffold-band subtree; a finite instance stores at most
// width endpoints per band.  CASE2 is a single activation-column check.
class VpCase12DirectEmptyLeftQuery {
 public:
  template <class Context>
  static VpCase12DirectEmptyLeftQuery build(
      cand_pos_t n,
      Context &ctx,
      std::size_t width = std::numeric_limits<std::size_t>::max()) {
    VpCase12DirectEmptyLeftQuery out;
    out.n_ = n;
    out.width_ = width;
    out.last_direct_endpoint_.assign(
        static_cast<std::size_t>(n + 1), -1);
    out.unpaired_count_.assign(static_cast<std::size_t>(n + 1), 0);
    out.case1_subtree_endpoints_.resize(static_cast<std::size_t>(n + 1));
    std::vector<std::vector<cand_pos_t>> child_bands(
        static_cast<std::size_t>(n + 1));
    for (cand_pos_t pos = 1; pos <= n; ++pos) {
      out.unpaired_count_[static_cast<std::size_t>(pos)] =
          out.unpaired_count_[static_cast<std::size_t>(pos - 1)] +
          (ctx.is_unpaired(pos) ? 1 : 0);
      if (ctx.is_unpaired(pos)) {
        const cand_pos_t ancestor = ctx.parent_index(pos);
        if (ancestor >= 0) {
          out.last_direct_endpoint_[static_cast<std::size_t>(ancestor)] = pos;
        }
      }
      if (ctx.pair_partner(pos) > pos) {
        const cand_pos_t ancestor = ctx.parent_index(pos);
        if (ancestor >= 0) {
          child_bands[static_cast<std::size_t>(ancestor)].push_back(pos);
        }
      }
    }

    // Child band openings are greater than their parents, hence descending
    // coordinate order is a bottom-up scaffold-tree traversal.
    for (cand_pos_t band = n; band >= 1; --band) {
      const cand_pos_t close = ctx.pair_partner(band);
      if (close <= band) continue;
      auto &summary =
          out.case1_subtree_endpoints_[static_cast<std::size_t>(band)];
      // Child subtrees are disjoint and listed by increasing open coordinate.
      // Their sorted summaries can be concatenated up to width, with the
      // band's own terminal endpoint last. No global candidate sort needed.
      for (const cand_pos_t child :
           child_bands[static_cast<std::size_t>(band)]) {
        const auto &child_summary =
            out.case1_subtree_endpoints_[static_cast<std::size_t>(child)];
        const auto count = std::min(width - summary.size(), child_summary.size());
        summary.insert(summary.end(),
                       child_summary.begin(),
                       child_summary.begin() + static_cast<std::ptrdiff_t>(count));
        if (summary.size() == width) break;
      }
      const cand_pos_t own_endpoint = close - 1;
      if (summary.size() < width && own_endpoint > band && ctx.is_unpaired(own_endpoint) &&
          ctx.parent_index(own_endpoint) == band) summary.push_back(own_endpoint);
    }
    return out;
  }

  template <class Context>
  std::vector<VpCase12PrefixGeometry> query_activation(
      cand_pos_t column,
      Context &ctx) const {
    std::vector<VpCase12PrefixGeometry> out;
    if (width_ == 0 || column < 2 || column > n_) return out;

    const cand_pos_t q = column - 1;
    const cand_pos_t child = ctx.pair_partner(q);
    if (child > 0 && child < q && ctx.pair_partner(child) == q) {
      const cand_pos_t ancestor = ctx.parent_index(child);
      if (ancestor >= 0 &&
          last_direct_endpoint_[static_cast<std::size_t>(ancestor)] > child) {
        for (const cand_pos_t i :
             case1_subtree_endpoints_[static_cast<std::size_t>(child)]) {
          out.push_back(VpCase12PrefixGeometry{
              grammar::RuleId::VP_WI_CASE1,
              i,
              i + 1,
              q,
              column,
              ancestor == 0 ? n_ : ctx.pair_partner(ancestor) - 1});
        }
      }
    }

    const cand_pos_t p = column;
    const cand_pos_t i = p - 1;
    const cand_pos_t p_close = ctx.pair_partner(p);
    const cand_pos_t ancestor = ctx.parent_index(i);
    if (p_close > p && ctx.is_unpaired(i) && ancestor >= 0 &&
        ctx.parent_index(p) == ancestor && unpaired_in_open_band(p, p_close)) {
      out.push_back(VpCase12PrefixGeometry{
          grammar::RuleId::VP_WI_CASE2, i, p, -1, p, p_close});
    }
    std::sort(out.begin(), out.end());
    return out;
  }

  std::size_t case1_summary_entry_count() const {
    std::size_t count = 0;
    for (const auto &summary : case1_subtree_endpoints_) {
      count += summary.size();
    }
    return count;
  }

 private:
  bool unpaired_in_open_band(cand_pos_t open, cand_pos_t close) const {
    if (open + 1 > close - 1) return false;
    return unpaired_count_[static_cast<std::size_t>(close - 1)] >
        unpaired_count_[static_cast<std::size_t>(open)];
  }

  cand_pos_t n_ = 0;
  std::size_t width_ = std::numeric_limits<std::size_t>::max();
  std::vector<cand_pos_t> last_direct_endpoint_;
  std::vector<cand_pos_t> unpaired_count_;
  std::vector<std::vector<cand_pos_t>> case1_subtree_endpoints_;
};

template <class Context>
class ScaffoldVpCase12LastChildIndex {
 public:
  static ScaffoldVpCase12LastChildIndex build(cand_pos_t n, Context &ctx) {
    ScaffoldVpCase12LastChildIndex out;
    struct PathEntry {
      cand_pos_t ancestor = -1;
      cand_pos_t child = -1;
    };
    std::vector<std::vector<PathEntry>> paths(
        static_cast<std::size_t>(n + 1));
    using EndpointsByBranch =
        std::map<cand_pos_t, std::vector<cand_pos_t>>;
    std::map<cand_pos_t, EndpointsByBranch> endpoints_by_ancestor;
    for (cand_pos_t endpoint = 1; endpoint <= n; ++endpoint) {
      if (!ctx.is_unpaired(endpoint)) continue;
      const cand_pos_t immediate = ctx.parent_index(endpoint);
      if (immediate >= 0) {
        endpoints_by_ancestor[immediate][endpoint].push_back(endpoint);
      }
      if (immediate <= 0 || ctx.pair_partner(immediate) <= endpoint) continue;
      cand_pos_t child = immediate;
      for (cand_pos_t ancestor = ctx.parent_index(child);
           ancestor >= 0;
           ancestor = ctx.parent_index(child)) {
        paths[static_cast<std::size_t>(endpoint)].push_back(
            PathEntry{ancestor, child});
        endpoints_by_ancestor[ancestor][child].push_back(endpoint);
        if (ancestor == 0) break;
        child = ancestor;
      }
    }

    engine::AllLocalWeightModel<Context> weight(&ctx);
    std::set<std::tuple<grammar::RuleId, cand_pos_t, cand_pos_t>> seen;
    auto append = [&](grammar::RuleId rule,
                      cand_pos_t i,
                      cand_pos_t j,
                      cand_pos_t p,
                      cand_pos_t q) {
      ++out.stats_.endpoint_candidates;
      if (i < 1 || i >= j || j > n) return;
      const deductions::ItemKey parent{
          grammar::NonTerminal::VP, i, j, -1, -1};
      for (auto source :
           deductions::generate_exact_basic_deductions(parent, rule, ctx)) {
        engine::require_legacy_span_schedule_support(source);
        if (weights::is_zero(weight.local_weight_for(source))) {
          continue;
        }
        if (source.split.p != p || source.split.q != q) {
          ++out.stats_.border_mismatches;
          continue;
        }
        if (!seen.emplace(rule, i, j).second) {
          ++out.stats_.duplicate_sources;
          continue;
        }
        out.insert(std::move(source));
      }
    };

    // CASE1: i lies below child, while j is a later unpaired endpoint
    // directly owned by ancestor.
    for (cand_pos_t i = 1; i <= n; ++i) {
      if (paths[static_cast<std::size_t>(i)].empty()) continue;
      const cand_pos_t p = ctx.pair_partner(ctx.parent_index(i));
      for (const auto &path : paths[static_cast<std::size_t>(i)]) {
        const cand_pos_t q = ctx.pair_partner(path.child);
        const auto branches = endpoints_by_ancestor.find(path.ancestor);
        if (branches == endpoints_by_ancestor.end()) continue;
        for (auto right = branches->second.upper_bound(path.child);
             right != branches->second.end();
             ++right) {
          for (const cand_pos_t j : right->second) {
            if (j <= q) continue;
            append(grammar::RuleId::VP_WI_CASE1, i, j, p, q);
          }
        }
      }
    }

    // CASE2 is the mirror: i is directly owned by ancestor before child,
    // while j lies below child.
    for (cand_pos_t j = 1; j <= n; ++j) {
      if (paths[static_cast<std::size_t>(j)].empty()) continue;
      const cand_pos_t q = ctx.parent_index(j);
      for (const auto &path : paths[static_cast<std::size_t>(j)]) {
        const cand_pos_t p = path.child;
        const auto branches = endpoints_by_ancestor.find(path.ancestor);
        if (branches == endpoints_by_ancestor.end()) continue;
        for (auto left = branches->second.begin();
             left != branches->second.lower_bound(path.child);
             ++left) {
          for (const cand_pos_t i : left->second) {
            if (i >= p) continue;
            append(grammar::RuleId::VP_WI_CASE2, i, j, p, q);
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

  const VpCase12LastChildStats &stats() const { return stats_; }

 private:
  void insert(deductions::Deduction source) {
    const auto idx = vp_case12_rule_index(source.rule);
    ++stats_.sources[idx];
    ++stats_.child_counts[idx][source.child_count];
    if (source.parent.i + 1 > source.split.p - 1) {
      ++stats_.left_empty[idx];
    }
    if (source.split.q + 1 > source.parent.j - 1) {
      ++stats_.right_empty[idx];
    }
    if (source.child_count == 0) {
      zero_child_by_column_[source.parent.j].push_back(std::move(source));
      ++stats_.zero_child[idx];
    } else {
      const auto trigger = source.children[source.child_count - 1].item;
      by_last_child_[trigger].push_back(std::move(source));
      ++stats_.child_indexed[idx];
    }
  }

  std::unordered_map<deductions::ItemKey,
                     std::vector<deductions::Deduction>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      by_last_child_;
  std::unordered_map<cand_pos_t, std::vector<deductions::Deduction>>
      zero_child_by_column_;
  VpCase12LastChildStats stats_{};
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE12_LAST_CHILD_INDEX_HH
