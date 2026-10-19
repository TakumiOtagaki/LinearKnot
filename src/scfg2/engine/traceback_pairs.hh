#ifndef CPARTY_SCFG2_ENGINE_TRACEBACK_PAIRS_HH
#define CPARTY_SCFG2_ENGINE_TRACEBACK_PAIRS_HH

#include "scfg2/deductions/exact_context.hh"
#include "scfg2/engine/traceback.hh"

#include <algorithm>
#include <set>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace cparty::scfg2::engine {

enum class TracebackPairLane : unsigned char { Round, Square };

struct TracebackPair {
  cand_pos_t i = -1;
  cand_pos_t j = -1;
  TracebackPairLane lane = TracebackPairLane::Round;
};

inline bool traceback_pair_equal(const TracebackPair &lhs, const TracebackPair &rhs) {
  return lhs.i == rhs.i && lhs.j == rhs.j && lhs.lane == rhs.lane;
}

inline auto traceback_pair_identity(const TracebackPair &pair) {
  return std::make_tuple(pair.i, pair.j, pair.lane);
}

struct TracebackPairHash {
  std::size_t operator()(const TracebackPair &pair) const {
    std::size_t h = std::hash<cand_pos_t>{}(pair.i);
    h ^= std::hash<cand_pos_t>{}(pair.j) + 0x9e3779b9U + (h << 6) + (h >> 2);
    h ^= static_cast<unsigned char>(pair.lane) + 0x9e3779b9U + (h << 6) + (h >> 2);
    return h;
  }
};

struct TracebackPairEqual {
  bool operator()(const TracebackPair &a, const TracebackPair &b) const {
    return traceback_pair_equal(a, b);
  }
};

// Preserve first-occurrence order and the full lane identity, with expected
// O(A) work for A occurrences and O(P) storage for P unique pairs.
// Bulk consumers must not repeatedly scan an ever-growing output vector.
inline void stable_unique_traceback_pairs(std::vector<TracebackPair> *pairs) {
  if (!pairs) return;
  std::unordered_set<TracebackPair, TracebackPairHash, TracebackPairEqual> seen;
  std::size_t kept = 0;
  for (const auto &pair : *pairs) {
    if (seen.insert(pair).second) {
      (*pairs)[kept++] = pair;
    }
  }
  pairs->resize(kept);
}

inline bool contains_traceback_pair(const std::vector<TracebackPair> &pairs,
                                    const TracebackPair &needle) {
  for (const auto &pair : pairs) {
    if (traceback_pair_equal(pair, needle)) return true;
  }
  return false;
}

inline bool same_traceback_pair_set(const std::vector<TracebackPair> &lhs,
                                    const std::vector<TracebackPair> &rhs) {
  if (lhs.size() != rhs.size()) return false;
  std::set<std::tuple<cand_pos_t, cand_pos_t, TracebackPairLane>> right;
  for (const auto &pair : rhs) right.insert(traceback_pair_identity(pair));
  for (const auto &pair : lhs) {
    if (right.find(traceback_pair_identity(pair)) == right.end()) return false;
  }
  return true;
}

inline void append_unique_traceback_pair(std::vector<TracebackPair> *out,
                                         const TracebackPair &pair) {
  if (out && !contains_traceback_pair(*out, pair)) out->push_back(pair);
}

inline void append_traceback_pair(std::vector<TracebackPair> *out,
                                  cand_pos_t i,
                                  cand_pos_t j,
                                  TracebackPairLane lane) {
  if (!out || i <= 0 || j <= 0 || i >= j) return;
  const TracebackPair pair{std::min(i, j), std::max(i, j), lane};
  for (const auto &existing : *out) {
    if (existing.i == pair.i && existing.j == pair.j && existing.lane == pair.lane) return;
  }
  out->push_back(pair);
}

template <class PairabilityContext>
TracebackPairLane classify_traceback_pair_lane(const PairabilityContext &ctx,
                                               cand_pos_t i,
                                               cand_pos_t j,
                                               TracebackPairLane fallback) {
  if (ctx.pair_round(i) == j || ctx.pair_round(j) == i) return TracebackPairLane::Round;
  if (ctx.pair_square(i) == j || ctx.pair_square(j) == i) return TracebackPairLane::Square;
  return fallback;
}

template <class PairabilityContext>
void append_origin_pair(std::vector<TracebackPair> *out,
                        const PairabilityContext &ctx,
                        cand_pos_t i,
                        cand_pos_t j,
                        TracebackPairLane fallback = TracebackPairLane::Square) {
  append_traceback_pair(out, i, j, classify_traceback_pair_lane(ctx, i, j, fallback));
}

template <class PairabilityContext>
void append_partner_pair(std::vector<TracebackPair> *out,
                         const PairabilityContext &ctx,
                         cand_pos_t pos,
                         TracebackPairLane fallback = TracebackPairLane::Square) {
  if (pos <= 0) return;
  const cand_pos_t mate = ctx.pair_partner(pos);
  if (mate > 0 && mate != pos) append_origin_pair(out, ctx, pos, mate, fallback);
}

// Observation-only output surface for one deduction. It must not decide
// admissibility or enumerate rules; constrained consumers may inspect this
// result without changing the underlying green deduction stream.
template <class PairabilityContext>
void append_deduction_output_pairs(const deductions::Deduction &deduction,
                                   const PairabilityContext &ctx,
                                   std::vector<TracebackPair> *out) {
  const auto &parent = deduction.parent;
  switch (parent.nonterminal) {
    case grammar::NonTerminal::V:
      append_traceback_pair(out, parent.i, parent.j, TracebackPairLane::Round);
      break;
    case grammar::NonTerminal::BE:
      append_traceback_pair(out, parent.i, parent.j, TracebackPairLane::Round);
      append_traceback_pair(out, parent.ip, parent.jp, TracebackPairLane::Round);
      break;
    case grammar::NonTerminal::VP:
    append_traceback_pair(out, parent.i, parent.j, TracebackPairLane::Square);
      if (deduction.rule == grammar::RuleId::VP_WI_CASE1 ||
          deduction.rule == grammar::RuleId::VP_WI_CASE3) {
        append_partner_pair(out, ctx, ctx.parent_index(parent.i), TracebackPairLane::Round);
      }
      if (deduction.rule == grammar::RuleId::VP_WI_CASE2 ||
          deduction.rule == grammar::RuleId::VP_WI_CASE3) {
        append_partner_pair(out, ctx, ctx.parent_index(parent.j), TracebackPairLane::Round);
      }
      break;
    case grammar::NonTerminal::WMBP:
      if (deduction.rule == grammar::RuleId::WMBP_SPLIT_BE_WI_VP) {
        append_partner_pair(out, ctx, parent.i);
        append_partner_pair(out, ctx, deduction.split.p);
      } else if (deduction.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP ||
                 deduction.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP) {
        append_partner_pair(out, ctx, deduction.split.p);
        append_partner_pair(out, ctx, deduction.split.q);
      }
      break;
    case grammar::NonTerminal::WMB:
      if (deduction.rule == grammar::RuleId::WMB_SPLIT_BE_WMBP_WI) {
        append_origin_pair(out, ctx, deduction.split.p, parent.j);
        append_partner_pair(out, ctx, deduction.split.q);
      }
      break;
    default:
      break;
  }
}

template <class PairabilityContext>
std::vector<TracebackPair> deduction_output_pairs(
    const deductions::Deduction &deduction,
    const PairabilityContext &ctx) {
  std::vector<TracebackPair> out;
  append_deduction_output_pairs(deduction, ctx, &out);
  return out;
}

template <class PairabilityContext>
bool is_fixed_scaffold_pair(const TracebackPair &pair,
                            const PairabilityContext &ctx);

template <class PairabilityContext>
void append_deduction_generated_pairs(const deductions::Deduction &deduction,
                                      const PairabilityContext &ctx,
                                      std::vector<TracebackPair> *out);

template <class PairabilityContext>
bool deduction_output_is_subset_of(const deductions::Deduction &deduction,
                                   const PairabilityContext &ctx,
                                   const std::vector<TracebackPair> &target) {
  std::vector<TracebackPair> emitted;
  append_deduction_generated_pairs(deduction, ctx, &emitted);
  for (const auto &pair : emitted) {
    if (!contains_traceback_pair(target, pair)) return false;
  }
  return true;
}

template <class PairabilityContext>
bool is_fixed_scaffold_pair(const TracebackPair &pair,
                            const PairabilityContext &ctx) {
  return ctx.pair_partner(pair.i) == pair.j || ctx.pair_partner(pair.j) == pair.i;
}

template <class PairabilityContext>
std::vector<TracebackPair> generated_traceback_pairs(
    const std::vector<TracebackPair> &pairs,
    const PairabilityContext &ctx) {
  std::vector<TracebackPair> out;
  for (const auto &pair : pairs) {
    if (!is_fixed_scaffold_pair(pair, ctx)) {
      out.push_back(TracebackPair{pair.i, pair.j, TracebackPairLane::Square});
    }
  }
  stable_unique_traceback_pairs(&out);
  return out;
}

template <class PairabilityContext>
void append_deduction_generated_pairs(const deductions::Deduction &deduction,
                                      const PairabilityContext &ctx,
                                      std::vector<TracebackPair> *out) {
  const auto nt = deduction.parent.nonterminal;
  if (nt != grammar::NonTerminal::V &&
      nt != grammar::NonTerminal::VP) {
    return;
  }
  const TracebackPair pair{
      deduction.parent.i, deduction.parent.j, TracebackPairLane::Square};
  if (!is_fixed_scaffold_pair(pair, ctx)) append_unique_traceback_pair(out, pair);
}

inline bool target_has_endpoint(const std::vector<TracebackPair> &target,
                                cand_pos_t pos) {
  for (const auto &pair : target) {
    if (pair.i == pos || pair.j == pos) return true;
  }
  return false;
}

inline bool target_has_endpoint_in_range(
    const std::vector<TracebackPair> &target,
    cand_pos_t first,
    cand_pos_t last) {
  if (last < first) return false;
  for (cand_pos_t pos = first; pos <= last; ++pos) {
    if (target_has_endpoint(target, pos)) return true;
  }
  return false;
}

template <class PairabilityContext>
bool deduction_locally_matches_generated_target(
    const deductions::Deduction &deduction,
    const PairabilityContext &ctx,
    const std::vector<TracebackPair> &target) {
  std::vector<TracebackPair> emitted;
  append_deduction_generated_pairs(deduction, ctx, &emitted);
  for (const auto &pair : emitted) {
    if (!contains_traceback_pair(target, pair)) return false;
  }

  const auto range_is_clear = [&](cand_pos_t first, cand_pos_t last) {
    return !target_has_endpoint_in_range(target, first, last);
  };
  const auto &parent = deduction.parent;
  const auto &split = deduction.split;
  switch (deduction.rule) {
    case grammar::RuleId::V_HAIRPIN:
      return range_is_clear(parent.i + 1, parent.j - 1);
    case grammar::RuleId::V_INTERNAL:
      return range_is_clear(parent.i + 1, split.k - 1) &&
             range_is_clear(split.l + 1, parent.j - 1);
    case grammar::RuleId::WI_BASE_SINGLE:
      return range_is_clear(parent.i, parent.i);
    case grammar::RuleId::WI_EXTEND_UNPAIRED:
    case grammar::RuleId::W_EXTEND_UNPAIRED:
    case grammar::RuleId::WMv_EXTEND_UNPAIRED:
    case grammar::RuleId::WMp_EXTEND_UNPAIRED:
    case grammar::RuleId::WM_EXTEND_UNPAIRED:
    case grammar::RuleId::WIP_EXTEND_UNPAIRED:
      return range_is_clear(parent.j, parent.j);
    case grammar::RuleId::VM_SPLIT_WMp_BASE:
      return range_is_clear(parent.i + 1, split.k - 1);
    case grammar::RuleId::WM_START_V:
    case grammar::RuleId::WM_START_WMB:
    case grammar::RuleId::WIP_BASEPAIR_V:
    case grammar::RuleId::WIP_BASEPAIR_WMB:
    case grammar::RuleId::VPL_SPLIT_VP:
      return range_is_clear(parent.i, split.k - 1);
    case grammar::RuleId::VPR_SPLIT_VP_BASEPAIR:
      return range_is_clear(split.k + 1, parent.j);
    case grammar::RuleId::VP_INTERNAL_LOOP:
      return range_is_clear(parent.i + 1, split.k - 1) &&
             range_is_clear(split.l + 1, parent.j - 1);
    case grammar::RuleId::BE_INTERNAL_LOOP:
      return range_is_clear(parent.i + 1, split.k - 1) &&
             range_is_clear(split.l + 1, parent.j - 1);
    case grammar::RuleId::BE_WIP_BASEPAIR:
      return range_is_clear(split.l + 1, parent.j - 1);
    case grammar::RuleId::BE_BASEPAIR_WIP:
      return range_is_clear(parent.i + 1, split.k - 1);
    default:
      return true;
  }
}

// Each deduction emits a fixed-size local list. Uniqueness is applied once
// after the full traversal, retaining the historical first-occurrence order.
template <class Value, class PairabilityContext>
void collect_traceback_pair_occurrences_into(
    const DeductionTraceNode<Value> &node,
    const PairabilityContext &ctx,
    std::vector<TracebackPair> *out,
    bool generated_only) {
  std::vector<const DeductionTraceNode<Value> *> pending{&node};
  while (!pending.empty()) {
    const auto *current = pending.back();
    pending.pop_back();
    if (current->chosen_deduction) {
      std::vector<TracebackPair> local;
      if (generated_only) {
        append_deduction_generated_pairs(*current->chosen_deduction, ctx, &local);
      } else {
        append_deduction_output_pairs(*current->chosen_deduction, ctx, &local);
      }
      out->insert(out->end(), local.begin(), local.end());
    }
    for (auto child = current->children.rbegin();
         child != current->children.rend(); ++child)
      pending.push_back(&*child);
  }
}

template <class Value, class PairabilityContext>
void collect_generated_traceback_pairs_into(
    const DeductionTraceNode<Value> &node,
    const PairabilityContext &ctx,
    std::vector<TracebackPair> *out) {
  collect_traceback_pair_occurrences_into(node, ctx, out, true);
  stable_unique_traceback_pairs(out);
}

template <class Value, class PairabilityContext>
std::vector<TracebackPair> collect_generated_traceback_pairs(
    const DeductionTraceNode<Value> &root,
    const PairabilityContext &ctx) {
  std::vector<TracebackPair> out;
  collect_generated_traceback_pairs_into(root, ctx, &out);
  return out;
}

template <class PairabilityContext>
std::vector<TracebackPair> generated_target_pairs_for_item(
    const deductions::ItemKey &item,
    const PairabilityContext &ctx,
    const std::vector<TracebackPair> &target) {
  std::vector<TracebackPair> out;
  for (const auto &pair : target) {
    if (is_fixed_scaffold_pair(pair, ctx)) continue;
    bool owned = item.i <= pair.i && pair.j <= item.j;
    if (item.nonterminal == grammar::NonTerminal::VM) {
      owned = item.i < pair.i && pair.j < item.j;
    } else if (item.nonterminal == grammar::NonTerminal::BE) {
      const bool in_left_wing = item.i < pair.i && pair.j < item.ip;
      const bool in_right_wing = item.jp < pair.i && pair.j < item.j;
      owned = in_left_wing || in_right_wing;
    }
    if (owned) {
      append_unique_traceback_pair(&out, pair);
    }
  }
  return out;
}

// Intersects one deduction with a fixed output pair map. The target generated
// inside a primary item span must be partitioned exactly between this
// deduction's own emissions and its children. BE is a fixed-scaffold carrier,
// so its round pairs are input geometry rather than generated output.
template <class PairabilityContext>
bool deduction_exactly_partitions_generated_target(
    const deductions::Deduction &deduction,
    const PairabilityContext &ctx,
    const std::vector<TracebackPair> &target) {
  const auto required = generated_target_pairs_for_item(deduction.parent, ctx, target);
  std::vector<TracebackPair> covered;

  std::vector<TracebackPair> emitted;
  append_deduction_generated_pairs(deduction, ctx, &emitted);
  for (const auto &pair : emitted) {
    if (!contains_traceback_pair(target, pair)) return false;
    append_unique_traceback_pair(&covered, pair);
  }
  for (std::size_t idx = 0; idx < deduction.child_count; ++idx) {
    for (const auto &pair : generated_target_pairs_for_item(
             deduction.children[idx].item, ctx, target)) {
      append_unique_traceback_pair(&covered, pair);
    }
  }
  return same_traceback_pair_set(required, covered);
}

template <class Value, class PairabilityContext>
void collect_traceback_pairs_into(const DeductionTraceNode<Value> &node,
                                  const PairabilityContext &ctx,
                                  std::vector<TracebackPair> *out) {
  collect_traceback_pair_occurrences_into(node, ctx, out, false);
  stable_unique_traceback_pairs(out);
}

template <class Value, class PairabilityContext>
std::vector<TracebackPair> collect_traceback_pairs(const DeductionTraceNode<Value> &root,
                                                   const PairabilityContext &ctx) {
  std::vector<TracebackPair> out;
  collect_traceback_pairs_into(root, ctx, &out);
  return out;
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_TRACEBACK_PAIRS_HH
