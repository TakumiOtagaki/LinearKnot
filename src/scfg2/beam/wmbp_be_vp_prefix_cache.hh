#ifndef CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_CACHE_HH
#define CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_CACHE_HH

#include "scfg2/beam/wmbp_be_vp_prefix_value.hh"
#include "scfg2/semirings/semiring.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <cmath>
#include <cstddef>
#include <optional>
#include <type_traits>
#include <unordered_map>

namespace cparty::scfg2::beam {

// Historical cache reference for standalone differential/profile tests only.
// Production recomputes these simple products and never instantiates this map.
template <class Value>
class BasicExactWmbpBeVpPrefixCache {
 public:
  template <class Semiring, class Chart, class LocalWeight>
  std::optional<Value> materialize_prefix(
      const deductions::Deduction &deduction,
      const Chart &chart,
      LocalWeight local_weight) {
    const auto evaluated = [&] {
      if constexpr (std::is_same_v<LocalWeight, weights::AnyLocalWeight>) {
        return std::visit([&](auto local) {
          return evaluate_wmbp_be_vp_prefix<Semiring>(deduction, chart, local);
        }, local_weight);
      } else {
        return evaluate_wmbp_be_vp_prefix<Semiring>(deduction, chart, local_weight);
      }
    }();
    if (!evaluated.has_value()) return std::nullopt;
    const auto key = wmbp_be_vp_prefix_key(deduction);
    const Value &candidate = *evaluated;

    const auto [found, inserted] = prefix_factors_.emplace(
        key, semirings::checkpoint_value(candidate));
    if (!inserted) {
      ++reuse_count_;
      const double found_priority = Semiring::priority(found->second);
      const double candidate_priority = Semiring::priority(candidate);
      const double tolerance =
          std::max(1e-14,
                   std::max(std::fabs(found_priority),
                            std::fabs(candidate_priority)) *
                       1e-12);
      if (std::fabs(found_priority - candidate_priority) > tolerance) {
        ++inconsistent_reuse_count_;
      }
    }
    return found->second;
  }

  template <class Semiring, class Chart, class LocalWeight>
  std::optional<Value> reconstruct_source(
      const deductions::Deduction &deduction,
      const Chart &chart,
      LocalWeight local_weight) {
    const auto prefix =
        materialize_prefix<Semiring>(deduction, chart, local_weight);
    if (!prefix.has_value()) return std::nullopt;
    ++reconstructed_source_count_;
    return Semiring::times(
        *prefix, chart.get(wmbp_be_vp_right_child(deduction)));
  }

  std::size_t prefix_count() const { return prefix_factors_.size(); }
  std::size_t materialization_count() const {
    return prefix_factors_.size();
  }
  std::size_t reuse_count() const { return reuse_count_; }
  std::size_t inconsistent_reuse_count() const {
    return inconsistent_reuse_count_;
  }
  std::size_t reconstructed_source_count() const {
    return reconstructed_source_count_;
  }

 private:
  std::unordered_map<WmbpBeVpPrefixKey,
                     Value,
                     WmbpBeVpPrefixKeyHash>
      prefix_factors_;
  std::size_t reuse_count_ = 0;
  std::size_t inconsistent_reuse_count_ = 0;
  std::size_t reconstructed_source_count_ = 0;
};

using ExactWmbpBeVpPrefixCache = BasicExactWmbpBeVpPrefixCache<double>;

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_CACHE_HH
