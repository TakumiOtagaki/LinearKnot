#ifndef CPARTY_SCFG2_BEAM_BE_LEFT_RIGHT_FRAME_SELECTION_HH
#define CPARTY_SCFG2_BEAM_BE_LEFT_RIGHT_FRAME_SELECTION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/bounded_ranked_bucket.hh"

#include "scfg2/beam/be_left_right_frame_profile.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/all_local_weight.hh"
#include "scfg2/engine/legacy_span_support.hh"

#include <algorithm>
#include <cstddef>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct BeRightShellFramePlanRecord {
  BeRightShellFrameKey key{};
  deductions::Deduction source{};
  weights::AnyLocalWeight local_weight = 0.0;
};

template <class Value>
struct BasicScoredBeRightShellFrame {
  BeRightShellFrameKey key{};
  Value complete_contribution{};
};

using ScoredBeRightShellFrame = BasicScoredBeRightShellFrame<double>;

inline bool be_right_shell_frame_rank_before(
    const ScoredBeRightShellFrame &lhs,
    const ScoredBeRightShellFrame &rhs) {
  if (lhs.complete_contribution != rhs.complete_contribution) {
    return lhs.complete_contribution > rhs.complete_contribution;
  }
  return be_right_shell_frame_tie_key(lhs.key) <
         be_right_shell_frame_tie_key(rhs.key);
}

class ExactBeRightShellFrameColumnPlan {
 public:
  template <class Context>
  static ExactBeRightShellFrameColumnPlan build(
      const std::vector<deductions::ItemKey> &schedule,
      Context &ctx) {
    ExactBeRightShellFrameColumnPlan out;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    for (const auto &item : schedule) {
      if (item.nonterminal != grammar::NonTerminal::BE) continue;
      for (const auto &deduction :
           deductions::generate_exact_basic_deductions_for_item(item, ctx)) {
        if (!recursive_be_rule_has_right_wip(deduction.rule)) {
          continue;
        }
        engine::require_legacy_span_schedule_support(deduction);
        const auto local = weight.local_weight_for(deduction);
        if (weights::is_zero(local)) continue;
        out.columns_[deduction.parent.j].push_back(
            BeRightShellFramePlanRecord{
                be_right_shell_frame_key(deduction), deduction, local});
        ++out.source_count_;
      }
    }
    return out;
  }

  template <class Semiring, class Chart, class TailAllowsFn>
  std::vector<BasicScoredBeRightShellFrame<typename Semiring::value_type>>
  score_column(
      cand_pos_t j,
      const Chart &chart,
      BasicExactBeFrameCache<typename Semiring::value_type> &cache,
      TailAllowsFn &&tail_allows) const {
    using Scored =
        BasicScoredBeRightShellFrame<typename Semiring::value_type>;
    std::vector<Scored> scored;
    const auto found = columns_.find(j);
    if (found == columns_.end()) return scored;
    for (const auto &record : found->second) {
      if (!tail_allows(record.key.tail)) continue;
      auto factor = cache.template materialize_right_tail<Semiring>(
          record.key.tail, chart);
      if (record.key.mode == BeRightShellFrameMode::LeftRight) {
        factor = Semiring::times(
            factor, chart.get(record.source.children[0].item));
      }
      scored.push_back(Scored{
          record.key,
          Semiring::times(
              weights::inject_local<Semiring>(record.local_weight), factor)});
    }
    return scored;
  }

  std::size_t source_count() const { return source_count_; }

 private:
  std::map<cand_pos_t, std::vector<BeRightShellFramePlanRecord>> columns_;
  std::size_t source_count_ = 0;
};

class TopWidthBeRightShellFrameSelection {
 public:
  explicit TopWidthBeRightShellFrameSelection(std::size_t width)
      : width_(width) {}

  template <class Semiring, class Value>
  void commit_column(cand_pos_t j,
                     std::vector<BasicScoredBeRightShellFrame<Value>> scored) {
    using Cell = BasicScoredBeRightShellFrame<Value>;
    struct Rank {
      bool operator()(const Cell &lhs, const Cell &rhs) const {
        const double x = Semiring::priority(lhs.complete_contribution);
        const double y = Semiring::priority(rhs.complete_contribution);
        return x != y ? x > y : be_right_shell_frame_tie_key(lhs.key) < be_right_shell_frame_tie_key(rhs.key);
      }
    };
    BoundedRankedBucket<Cell, Rank> bucket;
    for (auto &cell : scored) bucket.insert(std::move(cell), width_);
    for (const auto &cell : bucket.ranked()) selected_.insert(cell.key);
    candidate_per_column_[j] = scored.size();
    selected_per_column_[j] = bucket.size();
  }

  bool allows(const BeRightShellFrameKey &key) const {
    return selected_.find(key) != selected_.end();
  }

  bool allows(const deductions::Deduction &deduction) const {
    if (!recursive_be_rule_has_right_wip(deduction.rule)) return true;
    return allows(be_right_shell_frame_key(deduction));
  }

  std::size_t candidate_count() const {
    std::size_t count = 0;
    for (const auto &[j, value] : candidate_per_column_) {
      (void)j;
      count += value;
    }
    return count;
  }
  std::size_t selected_count() const { return selected_.size(); }
  std::size_t pruned_count() const {
    return candidate_count() - selected_count();
  }
  bool respects_width() const {
    for (const auto &[j, value] : selected_per_column_) {
      (void)j;
      if (value > width_) return false;
    }
    return true;
  }

 private:
  std::size_t width_ = 0;
  std::unordered_set<BeRightShellFrameKey,
                     BeRightShellFrameKeyHash>
      selected_;
  std::unordered_map<cand_pos_t, std::size_t> candidate_per_column_;
  std::unordered_map<cand_pos_t, std::size_t> selected_per_column_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_LEFT_RIGHT_FRAME_SELECTION_HH
