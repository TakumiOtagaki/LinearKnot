#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_LEFT_CHILD_INDEX_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_LEFT_CHILD_INDEX_HH

#include "scfg2/beam/vp_case3_border_event_plan.hh"
#include "scfg2/beam/vp_case3_left_cache.hh"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase3LeftChildAttachment {
  VpCase3LeftKey left{};
  unsigned carrier_mask = 0;
};

struct VpCase3LeftFutureSegment {
  deductions::ItemKey child{};
  cand_pos_t i = -1;
  cand_pos_t p = -1;
  cand_pos_t j_begin = -1;
  cand_pos_t j_end = -1;
  cand_pos_t q = -1;
  cand_pos_t k = -1;
  cand_pos_t l = -1;
  unsigned carrier_mask = 0;
};

struct VpCase3LeftChildIndexProfile {
  std::size_t child_count = 0;
  std::size_t attachment_count = 0;
  std::size_t child_target_group_count = 0;
  std::size_t target_run_count = 0;
  std::size_t signature_segment_count = 0;
  std::size_t maximum_same_target_count = 0;
  std::size_t maximum_target_runs_per_child = 0;
  std::size_t maximum_signature_segments_per_child = 0;
};

inline std::vector<VpCase3LeftFutureSegment>
build_vp_case3_left_future_segments_for_child(
    const deductions::ItemKey &child,
    const std::vector<VpCase3PlannedSource> &sources) {
  struct Geometry {
    cand_pos_t i = -1;
    cand_pos_t p = -1;
    cand_pos_t q = -1;
    cand_pos_t k = -1;
    cand_pos_t l = -1;
    unsigned carrier_mask = 0;
  };
  std::map<cand_pos_t, Geometry> by_target;
  for (const auto &source : sources) {
    auto &geometry = by_target[source.tail.j];
    geometry.i = source.deduction.parent.i;
    geometry.p = source.deduction.split.p;
    geometry.q = source.tail.q;
    geometry.k = source.tail.k;
    geometry.l = source.tail.l;
    geometry.carrier_mask |=
        vp_case3_carrier_bit(source.deduction.parent.nonterminal);
  }

  std::vector<VpCase3LeftFutureSegment> out;
  VpCase3LeftFutureSegment current{};
  bool active = false;
  for (const auto &[j, geometry] : by_target) {
    const bool continues =
        active && j == current.j_end + 1 && geometry.i == current.i &&
        geometry.p == current.p && geometry.q == current.q &&
        geometry.k == current.k && geometry.l == current.l &&
        geometry.carrier_mask == current.carrier_mask;
    if (!continues) {
      if (active) out.push_back(current);
      current = VpCase3LeftFutureSegment{child,
                                         geometry.i,
                                         geometry.p,
                                         j,
                                         j,
                                         geometry.q,
                                         geometry.k,
                                         geometry.l,
                                         geometry.carrier_mask};
      active = true;
    } else {
      current.j_end = j;
    }
  }
  if (active) out.push_back(current);
  return out;
}

class ExactVpCase3LeftChildIndex {
 public:
  explicit ExactVpCase3LeftChildIndex(
      const std::vector<VpCase3TailAttachmentSegment> &segments) {
    for (const auto &segment : segments) {
      for (cand_pos_t i = segment.i_begin; i <= segment.i_end; ++i) {
        const deductions::ItemKey child{
            grammar::NonTerminal::WI, i + 1, segment.p - 1, -1, -1};
        by_child_[child].push_back(VpCase3LeftChildAttachment{
            VpCase3LeftKey{segment.tail, i, segment.p},
            segment.carrier_mask});
        ++attachment_count_;
      }
    }
    build_future_segments();
  }

  const std::vector<VpCase3LeftFutureSegment> *future_segments_for(
      const deductions::ItemKey &child) const {
    const auto it = future_segments_by_child_.find(child);
    return it == future_segments_by_child_.end() ? nullptr : &it->second;
  }

  std::size_t future_segment_count() const { return future_segment_count_; }

  std::size_t child_count() const { return by_child_.size(); }
  std::size_t attachment_count() const { return attachment_count_; }
  std::size_t max_attachments_per_child() const {
    std::size_t maximum = 0;
    for (const auto &[child, attachments] : by_child_) {
      (void)child;
      maximum = std::max(maximum, attachments.size());
    }
    return maximum;
  }

  VpCase3LeftChildIndexProfile profile() const {
    VpCase3LeftChildIndexProfile out;
    out.child_count = child_count();
    out.attachment_count = attachment_count();
    for (const auto &[child, attachments] : by_child_) {
      (void)child;
      auto ordered = attachments;
      std::sort(ordered.begin(), ordered.end(), [](const auto &lhs,
                                                   const auto &rhs) {
        return std::make_tuple(lhs.left.tail.j,
                               lhs.left.tail.q,
                               lhs.left.tail.k,
                               lhs.left.tail.l) <
               std::make_tuple(rhs.left.tail.j,
                               rhs.left.tail.q,
                               rhs.left.tail.k,
                               rhs.left.tail.l);
      });
      std::size_t child_runs = 0;
      std::size_t child_segments = 0;
      cand_pos_t previous_j = -2;
      cand_pos_t previous_q = -2;
      cand_pos_t previous_k = -2;
      cand_pos_t previous_l = -2;
      unsigned previous_mask = 0;
      std::unordered_map<cand_pos_t, std::size_t> per_target;
      for (const auto &attachment : ordered) {
        const auto &tail = attachment.left.tail;
        const bool new_run = tail.j != previous_j + 1;
        if (new_run) ++child_runs;
        if (new_run || tail.q != previous_q || tail.k != previous_k ||
            tail.l != previous_l ||
            attachment.carrier_mask != previous_mask) {
          ++child_segments;
        }
        ++per_target[tail.j];
        previous_j = tail.j;
        previous_q = tail.q;
        previous_k = tail.k;
        previous_l = tail.l;
        previous_mask = attachment.carrier_mask;
      }
      out.child_target_group_count += per_target.size();
      out.target_run_count += child_runs;
      out.signature_segment_count += child_segments;
      out.maximum_target_runs_per_child =
          std::max(out.maximum_target_runs_per_child, child_runs);
      out.maximum_signature_segments_per_child =
          std::max(out.maximum_signature_segments_per_child, child_segments);
      for (const auto &[j, count] : per_target) {
        (void)j;
        out.maximum_same_target_count =
            std::max(out.maximum_same_target_count, count);
      }
    }
    return out;
  }

 private:
  void build_future_segments() {
    for (const auto &[child, attachments] : by_child_) {
      auto ordered = attachments;
      std::sort(ordered.begin(), ordered.end(), [](const auto &lhs,
                                                   const auto &rhs) {
        return lhs.left.tail.j < rhs.left.tail.j;
      });
      VpCase3LeftFutureSegment current{};
      bool active = false;
      for (const auto &attachment : ordered) {
        const auto &tail = attachment.left.tail;
        const bool continues =
            active && tail.j == current.j_end + 1 &&
            tail.q == current.q && tail.k == current.k &&
            tail.l == current.l &&
            attachment.carrier_mask == current.carrier_mask;
        if (!continues) {
          if (active) {
            future_segments_by_child_[child].push_back(current);
            ++future_segment_count_;
          }
          current = VpCase3LeftFutureSegment{child,
                                             attachment.left.i,
                                             attachment.left.p,
                                             tail.j,
                                             tail.j,
                                             tail.q,
                                             tail.k,
                                             tail.l,
                                             attachment.carrier_mask};
          active = true;
        } else {
          current.j_end = tail.j;
        }
      }
      if (active) {
        future_segments_by_child_[child].push_back(current);
        ++future_segment_count_;
      }
    }
  }

  std::unordered_map<deductions::ItemKey,
                     std::vector<VpCase3LeftChildAttachment>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      by_child_;
  std::unordered_map<deductions::ItemKey,
                     std::vector<VpCase3LeftFutureSegment>,
                     deductions::ItemKeyHash,
                     deductions::ItemKeyEq>
      future_segments_by_child_;
  std::size_t attachment_count_ = 0;
  std::size_t future_segment_count_ = 0;
};

template <class Value>
class BasicPendingVpCase3LeftCandidates {
 public:
  explicit BasicPendingVpCase3LeftCandidates(
      std::size_t active_width = std::numeric_limits<std::size_t>::max())
      : active_width_(active_width) {}

  template <class Semiring, class Chart>
  void observe_completed_child(const deductions::ItemKey &item,
                               const Chart &chart,
                               const ExactVpCase3LeftChildIndex &index) {
    if (item.nonterminal != grammar::NonTerminal::WI) return;
    const auto *segments = index.future_segments_for(item);
    if (!segments) return;
    const Value left_factor = chart.get(item);
    const double priority = Semiring::priority(left_factor);
    for (const auto &segment : *segments) {
      starts_[segment.j_begin].push_back(
          PendingSegment{segment, left_factor, priority, next_segment_id_++});
      observed_attachment_count_ +=
          static_cast<std::size_t>(segment.j_end - segment.j_begin + 1);
      ++registered_segment_count_;
    }
  }

  template <class Chart>
  void observe_completed_child(const deductions::ItemKey &item,
                               const Chart &chart,
                               const ExactVpCase3LeftChildIndex &index) {
    observe_completed_child<semirings::SumProductSemiring>(
        item, chart, index);
  }

  template <class Semiring, class Chart, class Context>
  void observe_completed_child(
      const deductions::ItemKey &item,
      const Chart &chart,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    if (item.nonterminal != grammar::NonTerminal::WI) return;
    const auto query = index.query_completed_wi(item, ctx);
    queried_endpoint_pair_count_ += query.endpoint_pair_count;
    border_mismatch_count_ += query.border_mismatch_count;
    duplicate_source_count_ += query.duplicate_source_count;
    queried_source_count_ += query.sources.size();
    const auto segments =
        build_vp_case3_left_future_segments_for_child(item, query.sources);
    const Value left_factor = chart.get(item);
    const double priority = Semiring::priority(left_factor);
    for (const auto &segment : segments) {
      starts_[segment.j_begin].push_back(
          PendingSegment{segment, left_factor, priority, next_segment_id_++});
      observed_attachment_count_ +=
          static_cast<std::size_t>(segment.j_end - segment.j_begin + 1);
      ++registered_segment_count_;
    }
  }

  template <class Chart, class Context>
  void observe_completed_child(
      const deductions::ItemKey &item,
      const Chart &chart,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    observe_completed_child<semirings::SumProductSemiring>(
        item, chart, index, ctx);
  }

  template <class Semiring, class Chart>
  void observe_target_sources(
      const std::vector<VpCase3PlannedSource> &sources,
      const Chart &chart) {
    using SourcesByChild =
        std::unordered_map<deductions::ItemKey,
                           std::vector<VpCase3PlannedSource>,
                           deductions::ItemKeyHash,
                           deductions::ItemKeyEq>;
    SourcesByChild by_child;
    for (const auto &source : sources) {
      if (source.deduction.child_count == 0) continue;
      by_child[source.deduction.children[0].item].push_back(source);
    }
    for (const auto &[child, child_sources] : by_child) {
      const auto segments =
          build_vp_case3_left_future_segments_for_child(child, child_sources);
      const Value left_factor = chart.get(child);
      const double priority = Semiring::priority(left_factor);
      for (const auto &segment : segments) {
        starts_[segment.j_begin].push_back(
            PendingSegment{
                segment, left_factor, priority, next_segment_id_++});
        observed_attachment_count_ +=
            static_cast<std::size_t>(segment.j_end - segment.j_begin + 1);
        ++registered_segment_count_;
      }
    }
  }

  template <class Chart>
  void observe_target_sources(
      const std::vector<VpCase3PlannedSource> &sources,
      const Chart &chart) {
    observe_target_sources<semirings::SumProductSemiring>(sources, chart);
  }

  template <class Semiring, class Chart>
  void observe_target_query(
      const VpCase3BorderEventChildQuery &query,
      const Chart &chart) {
    queried_endpoint_pair_count_ += query.endpoint_pair_count;
    border_mismatch_count_ += query.border_mismatch_count;
    duplicate_source_count_ += query.duplicate_source_count;
    queried_source_count_ += query.sources.size();
    observe_target_sources<Semiring>(query.sources, chart);
  }

  template <class Chart>
  void observe_target_query(
      const VpCase3BorderEventChildQuery &query,
      const Chart &chart) {
    observe_target_query<semirings::SumProductSemiring>(query, chart);
  }

  template <class Semiring>
  std::vector<BasicVpCase3LeftRecord<Value>> finalize_column(
      cand_pos_t j,
      const BasicExactVpCase3TailCache<Value> &tails) {
    activate_column(j);
    return score_active<Semiring>(j, tails);
  }

  std::vector<BasicVpCase3LeftRecord<Value>> finalize_column(
      cand_pos_t j,
      const BasicExactVpCase3TailCache<Value> &tails) {
    return finalize_column<semirings::SumProductSemiring>(j, tails);
  }

  template <class Semiring, class Chart>
  std::vector<BasicVpCase3LeftRecord<Value>> finalize_column(
      cand_pos_t j,
      const Chart &chart,
      BasicExactVpCase3TailCache<Value> &tails) {
    prepare_column<Semiring>(j, chart, tails);
    return score_prepared_column<Semiring>(
        j, tails, [](const VpCase3TailKey &) { return true; });
  }

  template <class Chart>
  std::vector<BasicVpCase3LeftRecord<Value>> finalize_column(
      cand_pos_t j,
      const Chart &chart,
      BasicExactVpCase3TailCache<Value> &tails) {
    return finalize_column<semirings::SumProductSemiring>(j, chart, tails);
  }

  template <class Semiring, class Chart>
  void prepare_column(cand_pos_t j,
                      const Chart &chart,
                      BasicExactVpCase3TailCache<Value> &tails) {
    activate_column(j);
    for (const auto &pending : active_) {
      ++tail_materialization_probe_count_;
      (void)tails.template materialize_tail<Semiring>(
          VpCase3TailKey{j,
                         pending.segment.q,
                         pending.segment.k,
                         pending.segment.l},
          chart);
    }
  }

  template <class Chart>
  void prepare_column(cand_pos_t j,
                      const Chart &chart,
                      BasicExactVpCase3TailCache<Value> &tails) {
    prepare_column<semirings::SumProductSemiring>(j, chart, tails);
  }

  template <class Semiring, class TailAllowedFn>
  std::vector<BasicVpCase3LeftRecord<Value>> score_prepared_column(
      cand_pos_t j,
      const BasicExactVpCase3TailCache<Value> &tails,
      TailAllowedFn &&tail_allowed) {
    std::vector<BasicVpCase3LeftRecord<Value>> finalized;
    finalized.reserve(active_.size());
    for (const auto &pending : active_) {
      ++prepared_score_probe_count_;
      const VpCase3TailKey tail{j,
                                pending.segment.q,
                                pending.segment.k,
                                pending.segment.l};
      if (!tail_allowed(tail)) continue;
      const auto tail_factor = tails.tail_factor(tail);
      if (!tail_factor.has_value()) continue;
      finalized.push_back(BasicVpCase3LeftRecord<Value>{
          VpCase3LeftKey{tail, pending.segment.i, pending.segment.p},
          Semiring::times(pending.left_factor, *tail_factor),
          pending.segment.carrier_mask});
    }
    return finalized;
  }

  template <class TailAllowedFn>
  std::vector<BasicVpCase3LeftRecord<Value>> score_prepared_column(
      cand_pos_t j,
      const BasicExactVpCase3TailCache<Value> &tails,
      TailAllowedFn &&tail_allowed) {
    return score_prepared_column<semirings::SumProductSemiring>(
        j, tails, std::forward<TailAllowedFn>(tail_allowed));
  }

  std::size_t queried_endpoint_pair_count() const {
    return queried_endpoint_pair_count_;
  }
  std::size_t queried_source_count() const { return queried_source_count_; }
  std::size_t border_mismatch_count() const { return border_mismatch_count_; }
  std::size_t duplicate_source_count() const {
    return duplicate_source_count_;
  }

  std::size_t observed_attachment_count() const {
    return observed_attachment_count_;
  }
  std::size_t registered_segment_count() const {
    return registered_segment_count_;
  }
  std::size_t active_expiry_probe_count() const {
    return active_expiry_probe_count_;
  }
  std::size_t tail_materialization_probe_count() const {
    return tail_materialization_probe_count_;
  }
  std::size_t prepared_score_probe_count() const {
    return prepared_score_probe_count_;
  }
  std::size_t max_active_segment_count() const {
    return max_active_segment_count_;
  }
  std::size_t finite_registered_segment_count() const {
    return finite_registered_segment_count_;
  }
  std::size_t finite_expired_segment_count() const {
    return finite_expired_segment_count_;
  }
  std::size_t finite_max_heap_size() const {
    return finite_max_heap_size_;
  }

  std::size_t live_segment_records() const {
    std::size_t total = active_.size() + finite_active_.size();
    for (const auto &[column, entries] : starts_) {
      (void)column;
      total += entries.size();
    }
    return total;
  }

  // The production endpoint-query path registers only the current target
  // column. Release its expired reservations and lazy heap entries together,
  // preserving aggregate work counters. Future-segment reference callers
  // must not discard a live later target through this operation.
  void release_target_column(cand_pos_t j) {
    for (const auto &[column, entries] : starts_) {
      (void)column;
      for (const auto &entry : entries)
        if (entry.segment.j_end > j)
          throw std::logic_error("CASE3 target release would discard a future segment");
    }
    starts_.clear();
    active_.clear();
    finite_active_ = {};
  }

 private:
  void activate_column(cand_pos_t j) {
    if (active_width_ != std::numeric_limits<std::size_t>::max()) {
      const auto starts = starts_.find(j);
      if (starts != starts_.end()) {
        for (const auto &entry : starts->second) finite_active_.push(entry);
        finite_registered_segment_count_ += starts->second.size();
      }
      active_.clear();
      while (active_.size() < active_width_ && !finite_active_.empty()) {
        ++active_expiry_probe_count_;
        const auto entry = finite_active_.top();
        finite_active_.pop();
        if (entry.segment.j_end < j) {
          ++finite_expired_segment_count_;
          continue;
        }
        active_.push_back(entry);
      }
      for (const auto &entry : active_) finite_active_.push(entry);
      finite_max_heap_size_ =
          std::max(finite_max_heap_size_, finite_active_.size());
      max_active_segment_count_ =
          std::max(max_active_segment_count_, active_.size());
      return;
    }
    active_expiry_probe_count_ += active_.size();
    active_.erase(
        std::remove_if(active_.begin(), active_.end(), [&](const auto &entry) {
          return entry.segment.j_end < j;
        }),
        active_.end());
    const auto starts = starts_.find(j);
    if (starts != starts_.end()) {
      active_.insert(active_.end(), starts->second.begin(), starts->second.end());
    }
    max_active_segment_count_ =
        std::max(max_active_segment_count_, active_.size());
  }

  template <class Semiring>
  std::vector<BasicVpCase3LeftRecord<Value>> score_active(
      cand_pos_t j,
      const BasicExactVpCase3TailCache<Value> &tails) const {
    std::vector<BasicVpCase3LeftRecord<Value>> finalized;
    finalized.reserve(active_.size());
    for (const auto &pending : active_) {
      const VpCase3TailKey tail{j,
                                pending.segment.q,
                                pending.segment.k,
                                pending.segment.l};
      const auto tail_factor = tails.tail_factor(tail);
      if (!tail_factor.has_value()) continue;
      finalized.push_back(BasicVpCase3LeftRecord<Value>{
          VpCase3LeftKey{
              tail, pending.segment.i, pending.segment.p},
          Semiring::times(pending.left_factor, *tail_factor),
          pending.segment.carrier_mask});
    }
    return finalized;
  }
  struct PendingSegment {
    VpCase3LeftFutureSegment segment{};
    Value left_factor{};
    double priority = 0.0;
    std::size_t id = 0;
  };
  struct PendingSegmentLowerPriority {
    bool operator()(const PendingSegment &lhs,
                    const PendingSegment &rhs) const {
      if (lhs.priority != rhs.priority) {
        return lhs.priority < rhs.priority;
      }
      const auto tie = [](const PendingSegment &entry) {
        return std::make_tuple(entry.segment.i,
                               entry.segment.p,
                               entry.segment.q,
                               entry.segment.k,
                               entry.segment.l,
                               entry.segment.j_begin,
                               entry.segment.j_end,
                               entry.id);
      };
      return tie(lhs) > tie(rhs);
    }
  };
  std::unordered_map<cand_pos_t, std::vector<PendingSegment>> starts_;
  std::vector<PendingSegment> active_;
  std::priority_queue<PendingSegment,
                      std::vector<PendingSegment>,
                      PendingSegmentLowerPriority>
      finite_active_;
  std::size_t active_width_ = std::numeric_limits<std::size_t>::max();
  std::size_t next_segment_id_ = 0;
  std::size_t observed_attachment_count_ = 0;
  std::size_t registered_segment_count_ = 0;
  std::size_t queried_endpoint_pair_count_ = 0;
  std::size_t queried_source_count_ = 0;
  std::size_t border_mismatch_count_ = 0;
  std::size_t duplicate_source_count_ = 0;
  std::size_t active_expiry_probe_count_ = 0;
  std::size_t tail_materialization_probe_count_ = 0;
  std::size_t prepared_score_probe_count_ = 0;
  std::size_t max_active_segment_count_ = 0;
  std::size_t finite_registered_segment_count_ = 0;
  std::size_t finite_expired_segment_count_ = 0;
  std::size_t finite_max_heap_size_ = 0;
};

using PendingVpCase3LeftCandidates =
    BasicPendingVpCase3LeftCandidates<double>;

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_LEFT_CHILD_INDEX_HH
