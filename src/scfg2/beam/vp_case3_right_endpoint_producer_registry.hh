#ifndef CPARTY_SCFG2_BEAM_VP_CASE3_RIGHT_ENDPOINT_PRODUCER_REGISTRY_HH
#define CPARTY_SCFG2_BEAM_VP_CASE3_RIGHT_ENDPOINT_PRODUCER_REGISTRY_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/beam/vp_case3_border_event_plan.hh"

#include <array>
#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace cparty::scfg2::beam {

struct VpCase3RegisteredLeftProducer {
  deductions::ItemKey item{};
  double left_factor = 0.0;
  cand_pos_t endpoint = -1;
  cand_pos_t p = -1;
  cand_pos_t ancestor = -1;
  cand_pos_t left_child = -1;
};

inline auto vp_case3_registered_producer_tie_key(
    const VpCase3RegisteredLeftProducer &producer) {
  return std::make_tuple(producer.item.i,
                         producer.item.j,
                         producer.ancestor,
                         producer.left_child);
}

inline bool vp_case3_registered_producer_rank_before(
    const VpCase3RegisteredLeftProducer &lhs,
    const VpCase3RegisteredLeftProducer &rhs) {
  if (lhs.left_factor != rhs.left_factor) {
    return lhs.left_factor > rhs.left_factor;
  }
  return vp_case3_registered_producer_tie_key(lhs) <
      vp_case3_registered_producer_tie_key(rhs);
}

struct VpCase3RightEndpointProducerMatch {
  VpCase3RegisteredLeftProducer producer{};
  VpCase3EndpointPathEntry right_path{};
};

template <class Context>
VpCase3BorderEventChildQuery build_vp_case3_sources_from_endpoint_matches(
    cand_pos_t j,
    const std::vector<VpCase3RightEndpointProducerMatch> &matches,
    const Context &ctx) {
  VpCase3BorderEventChildQuery out;
  engine::AllLocalWeightModel<Context> weight(&ctx);
  using SourceIdentity =
      std::tuple<grammar::NonTerminal, cand_pos_t, cand_pos_t>;
  std::set<SourceIdentity> seen;
  for (const auto &match : matches) {
    ++out.endpoint_pair_count;
    const auto &producer = match.producer;
    if (producer.item.j >= j) ++out.readiness_violation_count;
    const cand_pos_t i = producer.endpoint;
    const cand_pos_t p = producer.p;
    const cand_pos_t q = ctx.pair_partner(producer.left_child);
    const cand_pos_t k = match.right_path.child;
    const cand_pos_t l = ctx.parent_index(j);
    for (const auto nt : {grammar::NonTerminal::VP}) {
      const deductions::ItemKey parent{nt, i, j, -1, -1};
      const auto generated =
          deductions::generate_exact_vp_case3_deductions_for_item(parent, ctx);
      for (const auto &deduction : generated) {
        engine::require_legacy_span_schedule_support(deduction);
        if (vp_case3_partition(deduction) != VpCase3Partition::LeftNonempty ||
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
  return out;
}

// Observation-first online inverse of query_completed_wi(). Left WI items are
// registered when ready; a target endpoint then visits compatible producers
// from strictly earlier sibling bands. This exact registry intentionally has
// no finite selection yet.
class ExactVpCase3RightEndpointProducerRegistry {
 public:
  template <class Context>
  void observe_completed_wi(
      const deductions::ItemKey &item,
      double left_factor,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    if (!index.recognizes_completed_wi(item, ctx)) return;
    const cand_pos_t endpoint = item.i - 1;
    const cand_pos_t p = item.j + 1;
    for (const auto &path : index.paths_for_endpoint(endpoint)) {
      by_ancestor_child_[path.ancestor][path.child].push_back(
          VpCase3RegisteredLeftProducer{
              item,
              left_factor,
              endpoint,
              p,
              path.ancestor,
              path.child});
    }
  }

  template <class Context>
  VpCase3BorderEventChildQuery query_target_endpoint(
      cand_pos_t j,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) const {
    VpCase3BorderEventChildQuery out;
    if (!ctx.is_unpaired(j)) return out;

    engine::AllLocalWeightModel<Context> weight(&ctx);
    using SourceIdentity =
        std::tuple<grammar::NonTerminal, cand_pos_t, cand_pos_t>;
    std::set<SourceIdentity> seen;
    for (const auto &right_path : index.paths_for_endpoint(j)) {
      const auto ancestor = by_ancestor_child_.find(right_path.ancestor);
      if (ancestor == by_ancestor_child_.end()) continue;
      const auto &by_child = ancestor->second;
      for (auto left = by_child.begin();
           left != by_child.lower_bound(right_path.child);
           ++left) {
        for (const auto &producer : left->second) {
          ++out.endpoint_pair_count;
          if (producer.item.j >= j) ++out.readiness_violation_count;
          const cand_pos_t i = producer.endpoint;
          const cand_pos_t p = producer.p;
          const cand_pos_t q = ctx.pair_partner(producer.left_child);
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
  using ProducersByChild =
      std::map<cand_pos_t, std::vector<VpCase3RegisteredLeftProducer>>;
  std::map<cand_pos_t, ProducersByChild> by_ancestor_child_;
};

// Finite preceding-sibling frontier. Each child bucket and each accumulated
// prefix retains at most b producers. Queries for one target endpoint return
// at most b producers per endpoint ancestor path; later semantic filtering
// and column-global selection remain separate boundaries.
// Production queries sibling children monotonically after all producers of
// earlier children are ready. A static-child cursor visits each child once.
class TopWidthVpCase3PrecedingSiblingProducerFrontier {
 public:
  explicit TopWidthVpCase3PrecedingSiblingProducerFrontier(std::size_t width)
      : width_(width) {}

  template <class Context>
  void observe_completed_wi(
      const deductions::ItemKey &item,
      double left_factor,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    if (!index.recognizes_completed_wi(item, ctx)) return;
    const cand_pos_t endpoint = item.i - 1;
    const cand_pos_t p = item.j + 1;
    for (const auto &path : index.paths_for_endpoint(endpoint)) {
      if (static_cast<std::size_t>(path.ancestor) >= states_.size()) states_.resize(path.ancestor + 1);
      auto &state = states_[path.ancestor];
      state.registered = true;
      if (state.has_merged_child && path.child <= state.merged_through_child) {
        ++late_registration_count_;
      }
      if (static_cast<std::size_t>(path.child) >= by_child_.size()) {
        by_child_.resize(path.child + 1);
        registered_child_.resize(path.child + 1);
      }
      registered_child_[path.child] = 1;
      auto &bucket = by_child_[path.child];
      retain(bucket,
             VpCase3RegisteredLeftProducer{item,
                                           left_factor,
                                           endpoint,
                                           p,
                                           path.ancestor,
                                           path.child});
      ++candidate_registration_count_;
    }
  }

  template <class Context>
  std::vector<VpCase3RightEndpointProducerMatch> query_target_endpoint(
      cand_pos_t j,
      const ExactVpCase3BorderEventChildIndex &index,
      const Context &ctx) {
    std::vector<VpCase3RightEndpointProducerMatch> out;
    if (!ctx.is_unpaired(j)) return out;
    for (const auto &right_path : index.paths_for_endpoint(j)) {
      if (static_cast<std::size_t>(right_path.ancestor) >= states_.size()) continue;
      auto &state = states_[right_path.ancestor];
      if (!state.registered) continue;
      if (state.has_query && right_path.child < state.last_query_child) {
        ++nonmonotone_query_count_;
      }
      state.has_query = true;
      state.last_query_child = right_path.child;
      const auto &children = index.child_bands(right_path.ancestor);
      while (state.next_child < children.size() && children[state.next_child] < right_path.child) {
        const auto child = children[state.next_child++];
        ++visited_child_count_;
        if (static_cast<std::size_t>(child) >= by_child_.size() || !registered_child_[child]) continue;
        for (const auto &producer : by_child_[child]) retain(state.prefix, producer);
        state.has_merged_child = true;
        state.merged_through_child = child;
      }
      maximum_prefix_size_ = std::max(maximum_prefix_size_, state.prefix.size());
      maximum_returned_per_path_ =
          std::max(maximum_returned_per_path_, state.prefix.size());
      for (const auto &producer : state.prefix) {
        out.push_back(VpCase3RightEndpointProducerMatch{
            producer, right_path});
      }
    }
    return out;
  }

  std::size_t candidate_registration_count() const {
    return candidate_registration_count_;
  }
  std::size_t late_registration_count() const {
    return late_registration_count_;
  }
  std::size_t nonmonotone_query_count() const {
    return nonmonotone_query_count_;
  }
  std::size_t visited_child_count() const { return visited_child_count_; }
  std::size_t maximum_prefix_size() const { return maximum_prefix_size_; }
  std::size_t maximum_returned_per_path() const {
    return maximum_returned_per_path_;
  }

 private:
  void retain(std::vector<VpCase3RegisteredLeftProducer> &ranked,
              const VpCase3RegisteredLeftProducer &candidate) {
    ranked.push_back(candidate);
    std::sort(ranked.begin(),
              ranked.end(),
              vp_case3_registered_producer_rank_before);
    if (ranked.size() > width_) ranked.resize(width_);
  }

  struct AncestorState {
    std::vector<VpCase3RegisteredLeftProducer> prefix;
    std::size_t next_child = 0;
    cand_pos_t merged_through_child = -1;
    cand_pos_t last_query_child = -1;
    bool has_merged_child = false;
    bool has_query = false;
    bool registered = false;
  };

  std::size_t width_ = 0;
  std::vector<AncestorState> states_;
  std::vector<std::vector<VpCase3RegisteredLeftProducer>> by_child_;
  std::vector<unsigned char> registered_child_;
  std::size_t candidate_registration_count_ = 0;
  std::size_t late_registration_count_ = 0;
  std::size_t nonmonotone_query_count_ = 0;
  std::size_t visited_child_count_ = 0;
  std::size_t maximum_prefix_size_ = 0;
  std::size_t maximum_returned_per_path_ = 0;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VP_CASE3_RIGHT_ENDPOINT_PRODUCER_REGISTRY_HH
