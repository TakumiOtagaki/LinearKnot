#ifndef CPARTY_SCFG2_BEAM_VIRTUAL_PROFILE_HH
#define CPARTY_SCFG2_BEAM_VIRTUAL_PROFILE_HH

#include "scfg2/beam/virtual_candidate.hh"

#include <algorithm>
#include <cstddef>
#include <optional>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::beam {

struct VirtualBeamBucketKeyHash {
  std::size_t operator()(const VirtualBeamBucketKey &key) const noexcept {
    return (static_cast<std::size_t>(key.family) * 1315423911u) ^
           static_cast<std::size_t>(key.right_boundary + 3);
  }
};

struct VirtualBeamSourceObservation {
  VirtualBeamFullIdentity identity{};
  double source_contribution = 0.0;
  double ranking_value = 0.0;
  bool selected_source = false;
};

// Observation hypothesis only: the CASE3 private-tail product depends on the
// middle/right spans (q,k,l,j), not on the left parent coordinates (i,p).
// Compatibility provenance remains in full identity even when this key is
// used to measure potential algebraic sharing.
inline bool same_proposed_algebraic_identity(
    const VirtualBeamSourceObservation &lhs,
    const VirtualBeamSourceObservation &rhs) {
  if (lhs.identity.family != rhs.identity.family) return false;
  if (lhs.identity.family != VirtualBeamFamily::VpCase3LeftTail) {
    return lhs.identity == rhs.identity;
  }
  return lhs.identity.source_parent.j == rhs.identity.source_parent.j &&
         lhs.identity.source_split.q == rhs.identity.source_split.q &&
         lhs.identity.source_split.k == rhs.identity.source_split.k &&
         lhs.identity.source_split.l == rhs.identity.source_split.l;
}

struct VirtualBeamPriorityOverrides {
  std::optional<double> vp_case3_left_tail;
  std::optional<double> vp_case3_left;
};

inline bool virtual_beam_observation_precedes(
    const VirtualBeamSourceObservation &lhs,
    const VirtualBeamSourceObservation &rhs) {
  if (lhs.ranking_value != rhs.ranking_value) {
    return lhs.ranking_value > rhs.ranking_value;
  }
  const auto key = [](const VirtualBeamSourceObservation &entry) {
    const auto &parent = entry.identity.source_parent;
    const auto &split = entry.identity.source_split;
    return std::tuple{static_cast<unsigned char>(entry.identity.family),
                      static_cast<unsigned char>(parent.nonterminal),
                      parent.i,
                      parent.j,
                      parent.ip,
                      parent.jp,
                      split.k,
                      split.l,
                      split.p,
                      split.q};
  };
  return key(lhs) < key(rhs);
}

struct VirtualBeamBucketSummary {
  std::size_t candidate_count = 0;
  std::size_t parent_facing_count = 0;
  std::size_t selected_source_count = 0;
  double source_contribution_sum = 0.0;
  double max_source_contribution = 0.0;
  double max_ranking_value = 0.0;
  bool uses_factor_priority = false;
  std::vector<VirtualBeamSourceObservation> source_observations;

  std::vector<VirtualBeamSourceObservation> ranked_source_observations() const {
    auto ranked = source_observations;
    std::sort(ranked.begin(), ranked.end(), virtual_beam_observation_precedes);
    return ranked;
  }

  double top_source_contribution_sum(std::size_t width) const {
    if (width == 0 || source_observations.empty()) return 0.0;
    const auto ranked = ranked_source_observations();
    const std::size_t kept = std::min(width, ranked.size());
    double total = 0.0;
    for (std::size_t idx = 0; idx < kept; ++idx) {
      total += ranked[idx].source_contribution;
    }
    return total;
  }

  std::size_t selected_source_count_in_top(std::size_t width) const {
    if (width == 0 || source_observations.empty()) return 0;
    const auto ranked = ranked_source_observations();
    const std::size_t kept = std::min(width, ranked.size());
    std::size_t selected = 0;
    for (std::size_t idx = 0; idx < kept; ++idx) {
      if (ranked[idx].selected_source) ++selected;
    }
    return selected;
  }

  std::size_t proposed_unique_identity_count() const {
    return proposed_unique_representatives().size();
  }

  std::vector<std::size_t> proposed_unique_representatives() const {
    std::vector<std::size_t> representatives;
    for (std::size_t idx = 0; idx < source_observations.size(); ++idx) {
      bool seen = false;
      for (const std::size_t representative : representatives) {
        if (same_proposed_algebraic_identity(
                source_observations[idx], source_observations[representative])) {
          seen = true;
          break;
        }
      }
      if (!seen) representatives.push_back(idx);
    }
    return representatives;
  }

  std::vector<std::size_t> ranked_proposed_unique_representatives() const {
    auto representatives = proposed_unique_representatives();
    std::sort(representatives.begin(), representatives.end(),
              [&](std::size_t lhs, std::size_t rhs) {
                return virtual_beam_observation_precedes(
                    source_observations[lhs], source_observations[rhs]);
              });
    return representatives;
  }

  bool belongs_to_top_proposed_identity(
      const VirtualBeamSourceObservation &observation,
      const std::vector<std::size_t> &ranked_representatives,
      std::size_t kept) const {
    for (std::size_t idx = 0; idx < kept; ++idx) {
      if (same_proposed_algebraic_identity(
              observation, source_observations[ranked_representatives[idx]])) {
        return true;
      }
    }
    return false;
  }

  std::size_t selected_source_count_in_top_unique(std::size_t width) const {
    if (width == 0) return 0;
    const auto ranked = ranked_proposed_unique_representatives();
    const std::size_t kept = std::min(width, ranked.size());
    std::size_t selected = 0;
    for (const auto &observation : source_observations) {
      if (observation.selected_source &&
          belongs_to_top_proposed_identity(observation, ranked, kept)) {
        ++selected;
      }
    }
    return selected;
  }

  double top_unique_source_contribution_sum(std::size_t width) const {
    if (width == 0) return 0.0;
    const auto ranked = ranked_proposed_unique_representatives();
    const std::size_t kept = std::min(width, ranked.size());
    double total = 0.0;
    for (const auto &observation : source_observations) {
      if (belongs_to_top_proposed_identity(observation, ranked, kept)) {
        total += observation.source_contribution;
      }
    }
    return total;
  }

  std::size_t max_sources_per_proposed_identity() const {
    std::size_t maximum = 0;
    for (std::size_t idx = 0; idx < source_observations.size(); ++idx) {
      std::size_t matches = 0;
      for (const auto &candidate : source_observations) {
        if (same_proposed_algebraic_identity(
                source_observations[idx], candidate)) {
          ++matches;
        }
      }
      maximum = std::max(maximum, matches);
    }
    return maximum;
  }

  bool proposed_identities_have_consistent_ranking_values() const {
    for (std::size_t lhs = 0; lhs < source_observations.size(); ++lhs) {
      for (std::size_t rhs = lhs + 1; rhs < source_observations.size(); ++rhs) {
        if (same_proposed_algebraic_identity(source_observations[lhs],
                                             source_observations[rhs]) &&
            source_observations[lhs].ranking_value !=
                source_observations[rhs].ranking_value) {
          return false;
        }
      }
    }
    return true;
  }

  bool full_source_identities_are_unique() const {
    for (std::size_t lhs = 0; lhs < source_observations.size(); ++lhs) {
      for (std::size_t rhs = lhs + 1; rhs < source_observations.size(); ++rhs) {
        if (source_observations[lhs].identity ==
            source_observations[rhs].identity) {
          return false;
        }
      }
    }
    return true;
  }
};

// Observation-only profile. source_contribution is the contribution of the
// original direct deduction. For private factors such as CASE3_LTAIL it is a
// coverage/provenance metric, not yet the factor's algebraic partial product.
class VirtualBeamProfile {
 public:
  void observe(const deductions::Deduction &deduction,
               double source_contribution,
               bool selected_source = false,
               const VirtualBeamPriorityOverrides &priority_overrides = {}) {
    // Exact generation visits structurally admissible deductions whose child
    // product is semiring zero. They are not live frontier candidates and
    // must not consume a beam slot.
    if (!(source_contribution > 0.0)) return;
    const auto candidates = virtual_beam_candidates_for_deduction(deduction);
    for (std::size_t idx = 0; idx < candidates.size; ++idx) {
      const auto &candidate = candidates.candidates[idx];
      auto &summary = buckets_[candidate.bucket];
      double ranking_value = source_contribution;
      if (candidate.full_identity.family ==
              VirtualBeamFamily::VpCase3LeftTail &&
          priority_overrides.vp_case3_left_tail.has_value()) {
        ranking_value = *priority_overrides.vp_case3_left_tail;
        summary.uses_factor_priority = true;
      } else if (candidate.full_identity.family ==
                     VirtualBeamFamily::VpCase3Left &&
                 priority_overrides.vp_case3_left.has_value()) {
        ranking_value = *priority_overrides.vp_case3_left;
        summary.uses_factor_priority = true;
      }
      ++summary.candidate_count;
      if (candidate.parent_facing) ++summary.parent_facing_count;
      if (selected_source) ++summary.selected_source_count;
      summary.source_contribution_sum += source_contribution;
      summary.max_source_contribution =
          std::max(summary.max_source_contribution, source_contribution);
      summary.max_ranking_value =
          std::max(summary.max_ranking_value, ranking_value);
      summary.source_observations.push_back(VirtualBeamSourceObservation{
          candidate.full_identity,
          source_contribution,
          ranking_value,
          selected_source});
      ++total_candidate_observations_;
    }
  }

  const VirtualBeamBucketSummary *find(const VirtualBeamBucketKey &key) const {
    const auto it = buckets_.find(key);
    return it == buckets_.end() ? nullptr : &it->second;
  }

  std::size_t bucket_count() const { return buckets_.size(); }
  std::size_t total_candidate_observations() const {
    return total_candidate_observations_;
  }
  const auto &buckets() const { return buckets_; }

 private:
  std::unordered_map<VirtualBeamBucketKey,
                     VirtualBeamBucketSummary,
                     VirtualBeamBucketKeyHash>
      buckets_;
  std::size_t total_candidate_observations_ = 0;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_VIRTUAL_PROFILE_HH
