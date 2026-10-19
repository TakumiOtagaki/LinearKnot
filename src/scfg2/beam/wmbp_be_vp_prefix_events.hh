#ifndef CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_EVENTS_HH
#define CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_EVENTS_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/wmbp_be_vp_prefix_selection.hh"

#include <algorithm>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::beam {

template <class Value>
struct BasicWmbpBeVpPrefixParentDiscovery {
  std::vector<deductions::Deduction> deductions;
  std::vector<BasicScoredWmbpBeVpPrefix<Value>> scored;
};

struct WmbpSourceParentKeyDiscovery {
  std::vector<deductions::ItemKey> parents;
  std::size_t parent_vp_probes = 0;
};

template <class Semiring, class Chart, class Context>
std::vector<BasicScoredWmbpBeVpPrefix<typename Semiring::value_type>>
score_wmbp_be_vp_prefixes_from_deductions(
    const std::vector<deductions::Deduction> &deductions,
    const Chart &chart,
    Context &ctx) {
  engine::AllLocalWeightModel<Context> weight(&ctx);
  std::vector<BasicScoredWmbpBeVpPrefix<typename Semiring::value_type>> out;
  for (const auto &deduction : deductions) {
    if (!is_wmbp_be_vp_rule(deduction.rule)) continue;
    engine::require_legacy_span_schedule_support(deduction);
    const auto local = weight.local_weight_for(deduction);
    if (weights::is_zero(local)) continue;
    const auto factor = evaluate_wmbp_be_vp_prefix<Semiring>(deduction, chart, local);
    out.push_back({wmbp_be_vp_prefix_frame_key(deduction), *factor});
  }
  return out;
}

// Cache-taking overloads below remain standalone historical test references.
template <class Semiring, class Chart, class Context>
std::vector<BasicScoredWmbpBeVpPrefix<typename Semiring::value_type>>
score_wmbp_be_vp_prefixes_from_deductions(
    const std::vector<deductions::Deduction> &deductions,
    const Chart &chart,
    Context &ctx,
    BasicExactWmbpBeVpPrefixCache<typename Semiring::value_type> &cache) {
  engine::AllLocalWeightModel<Context> weight(&ctx);
  using Scored =
      BasicScoredWmbpBeVpPrefix<typename Semiring::value_type>;
  std::vector<Scored> out;
  for (const auto &deduction : deductions) {
    if (!is_wmbp_be_vp_rule(deduction.rule)) {
      continue;
    }
    engine::require_legacy_span_schedule_support(deduction);
    const auto local = weight.local_weight_for(deduction);
    if (weights::is_zero(local)) continue;
    const auto factor = cache.template materialize_prefix<Semiring>(
        deduction, chart, local);
    if (!factor.has_value()) continue;
    out.push_back(
        Scored{wmbp_be_vp_prefix_frame_key(deduction), *factor});
  }
  return out;
}

template <class Context>
WmbpSourceParentKeyDiscovery
discover_wmbp_source_parents_from_retained_vp_starts(
    cand_pos_t j,
    const std::vector<cand_pos_t> &vp_starts,
    Context &ctx) {
  WmbpSourceParentKeyDiscovery out;
  std::unordered_set<deductions::ItemKey,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      seen;
  for (const cand_pos_t vp_i : vp_starts) {
    const std::vector<cand_pos_t> one_start{vp_i};
    const deductions::ItemKey direct_parent{
        grammar::NonTerminal::WMBP, vp_i, j, -1, -1};
    ++out.parent_vp_probes;
    if (deductions::exact_wmbp_item_active(direct_parent, ctx) &&
        seen.insert(direct_parent).second) {
      out.parents.push_back(direct_parent);
    }

    const cand_pos_t split_parent_max =
        std::min(j, ctx.parent_index(vp_i));
    for (cand_pos_t parent_i = 1; parent_i <= split_parent_max; ++parent_i) {
      if (parent_i == vp_i) continue;
      ++out.parent_vp_probes;
      const deductions::ItemKey parent{
          grammar::NonTerminal::WMBP, parent_i, j, -1, -1};
      if (!deductions::exact_wmbp_item_active(parent, ctx)) continue;
      const auto family = deductions::
          generate_exact_wmbp_be_vp_family_deductions_for_vp_starts(
              parent, ctx, one_start, true, true, true);
      const bool has_source = !family.wmbp.empty() || !family.wmbw.empty() ||
                              !family.wi.empty();
      if (has_source && seen.insert(parent).second) {
        out.parents.push_back(parent);
      }
    }
  }
  std::sort(out.parents.begin(), out.parents.end(), [](const auto &lhs,
                                                        const auto &rhs) {
    return lhs.i > rhs.i;
  });
  return out;
}

template <class Context>
std::vector<deductions::Deduction>
generate_wmbp_parent_deductions_for_retained_vp_starts(
    const deductions::ItemKey &parent,
    const std::vector<cand_pos_t> &vp_starts,
    Context &ctx) {
  if (!deductions::exact_wmbp_item_active(parent, ctx)) return {};
  auto family = deductions::
      generate_exact_wmbp_be_vp_family_deductions_for_vp_starts(
          parent, ctx, vp_starts, true, true, true);
  std::vector<deductions::Deduction> out;
  out.insert(out.end(),
             std::make_move_iterator(family.wmbp.begin()),
             std::make_move_iterator(family.wmbp.end()));
  out.insert(out.end(),
             std::make_move_iterator(family.wmbw.begin()),
             std::make_move_iterator(family.wmbw.end()));
  if (std::binary_search(vp_starts.begin(), vp_starts.end(), parent.i)) {
    auto direct = deductions::generate_exact_basic_deductions(
        parent, grammar::RuleId::WMBP_DIRECT_VP, ctx);
    out.insert(out.end(),
               std::make_move_iterator(direct.begin()),
               std::make_move_iterator(direct.end()));
  }
  out.insert(out.end(),
             std::make_move_iterator(family.wi.begin()),
             std::make_move_iterator(family.wi.end()));
  return out;
}

// Discover one parent's admissible WMBP BE-prefix candidates at the point
// where that parent is ready. The shared exact deduction generator remains the
// sole grammar/admissibility source; this helper only derives a scored finite
// selection surface from its output.
template <class Semiring, class Chart, class Context>
BasicWmbpBeVpPrefixParentDiscovery<typename Semiring::value_type>
discover_wmbp_be_vp_prefixes_for_parent(
    const deductions::ItemKey &parent,
    const Chart &chart,
    Context &ctx,
    BasicExactWmbpBeVpPrefixCache<typename Semiring::value_type> &cache) {
  BasicWmbpBeVpPrefixParentDiscovery<typename Semiring::value_type> out;
  if (parent.nonterminal != grammar::NonTerminal::WMBP) return out;
  out.deductions =
      deductions::generate_exact_basic_deductions_for_item(parent, ctx);
  out.scored = score_wmbp_be_vp_prefixes_from_deductions<Semiring>(
      out.deductions, chart, ctx, cache);
  return out;
}

inline bool same_scored_wmbp_be_vp_prefixes(
    std::vector<ScoredWmbpBeVpPrefix> lhs,
    std::vector<ScoredWmbpBeVpPrefix> rhs) {
  std::sort(lhs.begin(), lhs.end(), wmbp_be_vp_prefix_rank_before);
  std::sort(rhs.begin(), rhs.end(), wmbp_be_vp_prefix_rank_before);
  if (lhs.size() != rhs.size()) return false;
  for (std::size_t idx = 0; idx < lhs.size(); ++idx) {
    if (!(lhs[idx].key == rhs[idx].key) ||
        lhs[idx].prefix_factor != rhs[idx].prefix_factor) {
      return false;
    }
  }
  return true;
}

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_EVENTS_HH
