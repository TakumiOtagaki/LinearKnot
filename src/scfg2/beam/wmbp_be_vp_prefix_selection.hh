#ifndef CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_SELECTION_HH
#define CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_SELECTION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/wmbp_be_vp_prefix_cache.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct WmbpBeVpPrefixFrameKey {
  deductions::ItemKey parent{};
  WmbpBeVpPrefixKey prefix{};
};

inline bool operator==(const WmbpBeVpPrefixFrameKey &lhs,
                       const WmbpBeVpPrefixFrameKey &rhs) {
  return deductions::ItemKeyEq{}(lhs.parent, rhs.parent) &&
         lhs.prefix == rhs.prefix;
}

struct WmbpBeVpPrefixFrameKeyHash {
  std::size_t operator()(const WmbpBeVpPrefixFrameKey &key) const noexcept {
    return (WmbpBeVpPrefixKeyHash{}(key.prefix) * 1315423911u) ^
           deductions::ItemKeyHash{}(key.parent);
  }
};

inline auto wmbp_be_vp_prefix_frame_tie_key(
    const WmbpBeVpPrefixFrameKey &key) {
  return std::make_tuple(item_tie_key(key.parent),
                         wmbp_be_vp_prefix_tie_key(key.prefix));
}

inline WmbpBeVpPrefixFrameKey wmbp_be_vp_prefix_frame_key(
    const deductions::Deduction &deduction) {
  return {deduction.parent, wmbp_be_vp_prefix_key(deduction)};
}

struct WmbpBeVpPrefixPlanRecord {
  deductions::Deduction representative{};
  weights::AnyLocalWeight local_weight = 0.0;
  std::size_t source_count = 0;
};

template <class Value>
struct BasicScoredWmbpBeVpPrefix {
  WmbpBeVpPrefixFrameKey key{};
  Value prefix_factor{};
};

using ScoredWmbpBeVpPrefix = BasicScoredWmbpBeVpPrefix<double>;

inline bool wmbp_be_vp_prefix_rank_before(
    const ScoredWmbpBeVpPrefix &lhs,
    const ScoredWmbpBeVpPrefix &rhs) {
  if (lhs.prefix_factor != rhs.prefix_factor) {
    return lhs.prefix_factor > rhs.prefix_factor;
  }
  return wmbp_be_vp_prefix_frame_tie_key(lhs.key) <
         wmbp_be_vp_prefix_frame_tie_key(rhs.key);
}

class ExactWmbpBeVpPrefixParentPlan {
 private:
  struct FrameLess {
    bool operator()(const WmbpBeVpPrefixFrameKey &lhs,
                    const WmbpBeVpPrefixFrameKey &rhs) const {
      return wmbp_be_vp_prefix_frame_tie_key(lhs) <
             wmbp_be_vp_prefix_frame_tie_key(rhs);
    }
  };

 public:
  template <class Context>
  static ExactWmbpBeVpPrefixParentPlan build(
      const std::vector<deductions::ItemKey> &schedule,
      Context &ctx) {
    ExactWmbpBeVpPrefixParentPlan out;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    for (const auto &item : schedule) {
      if (item.nonterminal != grammar::NonTerminal::WMBP) continue;
      for (const auto &deduction :
           deductions::generate_exact_basic_deductions_for_item(item, ctx)) {
        if (!is_wmbp_be_vp_rule(deduction.rule)) {
          continue;
        }
        engine::require_legacy_span_schedule_support(deduction);
        const auto local = weight.local_weight_for(deduction);
        if (weights::is_zero(local)) continue;
        const auto key = wmbp_be_vp_prefix_frame_key(deduction);
        auto &record = out.parents_[item][key];
        if (record.source_count == 0) {
          record = WmbpBeVpPrefixPlanRecord{deduction, local, 1};
          ++out.prefix_frame_count_;
        } else {
          ++record.source_count;
        }
        ++out.source_count_;
      }
    }
    return out;
  }

  template <class Semiring, class Chart>
  std::vector<BasicScoredWmbpBeVpPrefix<typename Semiring::value_type>>
  materialize_parent(
      const deductions::ItemKey &parent,
      const Chart &chart,
      BasicExactWmbpBeVpPrefixCache<typename Semiring::value_type> &cache)
      const {
    using Scored =
        BasicScoredWmbpBeVpPrefix<typename Semiring::value_type>;
    std::vector<Scored> scored;
    const auto found = parents_.find(parent);
    if (found == parents_.end()) return scored;
    scored.reserve(found->second.size());
    for (const auto &[key, record] : found->second) {
      const auto factor = cache.template materialize_prefix<Semiring>(
          record.representative, chart, record.local_weight);
      if (factor.has_value()) scored.push_back({key, *factor});
    }
    return scored;
  }

  std::size_t prefix_frame_count() const { return prefix_frame_count_; }
  std::size_t source_count() const { return source_count_; }

 private:
  struct ItemLess {
    bool operator()(const deductions::ItemKey &lhs,
                    const deductions::ItemKey &rhs) const {
      return item_tie_key(lhs) < item_tie_key(rhs);
    }
  };

  std::map<deductions::ItemKey,
           std::map<WmbpBeVpPrefixFrameKey,
                    WmbpBeVpPrefixPlanRecord,
                    FrameLess>,
           ItemLess>
      parents_;
  std::size_t prefix_frame_count_ = 0;
  std::size_t source_count_ = 0;
};

class TopWidthWmbpBeVpPrefixSelection {
 public:
  explicit TopWidthWmbpBeVpPrefixSelection(std::size_t width)
      : width_(width) {}

  template <class Semiring, class Value>
  void commit_parent(const deductions::ItemKey &parent,
                     std::vector<BasicScoredWmbpBeVpPrefix<Value>> scored) {
    std::sort(scored.begin(), scored.end(), [](const auto &lhs,
                                               const auto &rhs) {
      const double lhs_priority = Semiring::priority(lhs.prefix_factor);
      const double rhs_priority = Semiring::priority(rhs.prefix_factor);
      if (lhs_priority != rhs_priority) return lhs_priority > rhs_priority;
      return wmbp_be_vp_prefix_frame_tie_key(lhs.key) <
             wmbp_be_vp_prefix_frame_tie_key(rhs.key);
    });
    const std::size_t kept = std::min(width_, scored.size());
    for (std::size_t idx = 0; idx < kept; ++idx) {
      selected_.insert(scored[idx].key);
    }
    candidate_per_parent_[parent] = scored.size();
    selected_per_parent_[parent] = kept;
  }

  void commit_parent(const deductions::ItemKey &parent,
                     std::vector<ScoredWmbpBeVpPrefix> scored) {
    commit_parent<semirings::SumProductSemiring>(parent, std::move(scored));
  }

  bool allows(const WmbpBeVpPrefixFrameKey &key) const {
    return selected_.find(key) != selected_.end();
  }

  bool allows(const deductions::Deduction &deduction) const {
    if (!is_wmbp_be_vp_rule(deduction.rule)) return true;
    return allows(wmbp_be_vp_prefix_frame_key(deduction));
  }

  std::size_t candidate_count() const {
    std::size_t count = completed_candidates_;
    for (const auto &[parent, value] : candidate_per_parent_) {
      (void)parent;
      count += value;
    }
    return count;
  }
  std::size_t selected_count() const {
    return completed_selected_ + selected_.size();
  }
  std::size_t pruned_count() const {
    return candidate_count() - selected_count();
  }
  bool respects_width() const {
    if (!completed_respects_width_) return false;
    for (const auto &[parent, value] : selected_per_parent_) {
      (void)parent;
      if (value > width_) return false;
    }
    return true;
  }

  // Forward consumers only ask about a parent's branches while evaluating
  // that parent. The graph-based outside does not read this selection history.
  // Release its keys after consumption, retaining aggregate diagnostics only.
  // Callers that need historical lookups (standalone references) omit this.
  void release_completed_parents() {
    completed_candidates_ = candidate_count();
    completed_selected_ = selected_count();
    completed_respects_width_ = respects_width();
    selected_.clear();
    candidate_per_parent_.clear();
    selected_per_parent_.clear();
  }

  std::size_t live_selection_count() const { return selected_.size(); }
  std::size_t live_parent_count() const { return candidate_per_parent_.size(); }

 private:
  std::size_t width_ = 0;
  std::size_t completed_candidates_ = 0;
  std::size_t completed_selected_ = 0;
  bool completed_respects_width_ = true;
  std::unordered_set<WmbpBeVpPrefixFrameKey,
                     WmbpBeVpPrefixFrameKeyHash>
      selected_;
  std::unordered_map<deductions::ItemKey,
                     std::size_t,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      candidate_per_parent_;
  std::unordered_map<deductions::ItemKey,
                     std::size_t,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      selected_per_parent_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_SELECTION_HH
