#ifndef CPARTY_SCFG2_BEAM_BE_RIGHT_TAIL_SELECTION_HH
#define CPARTY_SCFG2_BEAM_BE_RIGHT_TAIL_SELECTION_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/bounded_ranked_bucket.hh"

#include "scfg2/beam/be_frame_cache.hh"
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

struct BeRightTailPlanRecord {
  BeRightTailKey key{};
  deductions::Deduction representative{};
};

template <class Value>
struct BasicScoredBeRightTail {
  BeRightTailKey key{};
  Value factor{};
};

using ScoredBeRightTail = BasicScoredBeRightTail<double>;

inline auto be_right_tail_tie_key(const BeRightTailKey &key) {
  return std::make_tuple(
      key.core.k, key.core.l, key.core.ip, key.core.jp, key.j);
}

class ExactBeRightTailColumnPlan {
 public:
  template <class Context>
  static ExactBeRightTailColumnPlan build(
      const std::vector<deductions::ItemKey> &schedule,
      Context &ctx) {
    ExactBeRightTailColumnPlan out;
    engine::AllLocalWeightModel<Context> weight(&ctx);
    for (const auto &item : schedule) {
      if (item.nonterminal != grammar::NonTerminal::BE) continue;
      for (const auto &deduction :
           deductions::generate_exact_basic_deductions_for_item(item, ctx)) {
        if (!recursive_be_rule_has_right_wip(deduction.rule)) continue;
        engine::require_legacy_span_schedule_support(deduction);
        if (weights::is_zero(weight.local_weight_for(deduction))) {
          continue;
        }
        const BeRightTailKey key{deduction.parent.j,
                                 be_core_key(deduction)};
        auto &column = out.columns_[deduction.parent.j];
        if (column.find(key) == column.end()) {
          column.emplace(key, BeRightTailPlanRecord{key, deduction});
          ++out.tail_count_;
        }
        ++out.source_count_;
      }
    }
    return out;
  }

  template <class Semiring, class Chart>
  std::vector<BasicScoredBeRightTail<typename Semiring::value_type>>
  materialize_column(
      cand_pos_t j,
      const Chart &chart,
      BasicExactBeFrameCache<typename Semiring::value_type> &cache) const {
    using Scored = BasicScoredBeRightTail<typename Semiring::value_type>;
    std::vector<Scored> scored;
    const auto column = columns_.find(j);
    if (column == columns_.end()) return scored;
    scored.reserve(column->second.size());
    for (const auto &[key, record] : column->second) {
      const auto factor = cache.template materialize_right_tail_for_source<
          Semiring>(
          record.representative, chart);
      if (factor.has_value()) scored.push_back({key, *factor});
    }
    return scored;
  }

  template <class Chart>
  std::vector<ScoredBeRightTail> materialize_column(
      cand_pos_t j,
      const Chart &chart,
      ExactBeFrameCache &cache) const {
    return materialize_column<semirings::SumProductSemiring>(j, chart, cache);
  }

  std::size_t tail_count() const { return tail_count_; }
  std::size_t source_count() const { return source_count_; }
  std::size_t column_count() const { return columns_.size(); }

  template <class Fn>
  void for_each_record(Fn &&fn) const {
    for (const auto &[j, column] : columns_) {
      (void)j;
      for (const auto &[key, record] : column) {
        (void)key;
        fn(record);
      }
    }
  }

 private:
  std::map<cand_pos_t,
           std::unordered_map<BeRightTailKey,
                              BeRightTailPlanRecord,
                              BeRightTailKeyHash>>
      columns_;
  std::size_t tail_count_ = 0;
  std::size_t source_count_ = 0;
};

class TopWidthBeRightTailSelection {
 public:
  explicit TopWidthBeRightTailSelection(std::size_t width) : width_(width) {}

  template <class Semiring, class Value>
  void commit_column(cand_pos_t j,
                     std::vector<BasicScoredBeRightTail<Value>> scored) {
    using Cell = BasicScoredBeRightTail<Value>;
    struct Rank {
      bool operator()(const Cell &lhs, const Cell &rhs) const {
        const double x = Semiring::priority(lhs.factor), y = Semiring::priority(rhs.factor);
        return x != y ? x > y : be_right_tail_tie_key(lhs.key) < be_right_tail_tie_key(rhs.key);
      }
    };
    BoundedRankedBucket<Cell, Rank> bucket;
    for (auto &cell : scored) bucket.insert(std::move(cell), width_);
    for (const auto &cell : bucket.ranked()) selected_.insert(cell.key);
    candidate_per_column_[j] = scored.size();
    selected_per_column_[j] = bucket.size();
  }

  bool allows(const BeRightTailKey &key) const {
    return selected_.find(key) != selected_.end();
  }

  bool allows(const deductions::Deduction &deduction) const {
    if (!recursive_be_rule_has_right_wip(deduction.rule)) return true;
    return allows(BeRightTailKey{deduction.parent.j,
                                 be_core_key(deduction)});
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
  std::unordered_set<BeRightTailKey,
                     BeRightTailKeyHash>
      selected_;
  std::unordered_map<cand_pos_t, std::size_t> candidate_per_column_;
  std::unordered_map<cand_pos_t, std::size_t> selected_per_column_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_RIGHT_TAIL_SELECTION_HH
