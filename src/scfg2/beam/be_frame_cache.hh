#ifndef CPARTY_SCFG2_BEAM_BE_FRAME_CACHE_HH
#define CPARTY_SCFG2_BEAM_BE_FRAME_CACHE_HH

#include "scfg2/deductions/deduction.hh"
#include "scfg2/semirings/semiring.hh"
#include "scfg2/semirings/value_boundary.hh"

#include <cstddef>
#include <optional>
#include <unordered_map>

namespace cparty::scfg2::beam {

struct BeCoreKey {
  cand_pos_t k = -1;
  cand_pos_t l = -1;
  cand_pos_t ip = -1;
  cand_pos_t jp = -1;
};

inline constexpr bool operator==(const BeCoreKey &lhs, const BeCoreKey &rhs) {
  return lhs.k == rhs.k && lhs.l == rhs.l && lhs.ip == rhs.ip &&
         lhs.jp == rhs.jp;
}

struct BeCoreKeyHash {
  std::size_t operator()(const BeCoreKey &key) const noexcept {
    std::size_t out = static_cast<std::size_t>(key.k + 11);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.l + 13);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.ip + 17);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.jp + 19);
    return out;
  }
};

struct BeRightTailKey {
  cand_pos_t j = -1;
  BeCoreKey core{};
};

inline constexpr bool operator==(const BeRightTailKey &lhs,
                                 const BeRightTailKey &rhs) {
  return lhs.j == rhs.j && lhs.core == rhs.core;
}

struct BeRightTailKeyHash {
  std::size_t operator()(const BeRightTailKey &key) const noexcept {
    return (BeCoreKeyHash{}(key.core) * 1315423911u) ^
           static_cast<std::size_t>(key.j + 23);
  }
};

inline bool is_recursive_be_frame_rule(grammar::RuleId rule) {
  return rule == grammar::RuleId::BE_INTERNAL_LOOP ||
         rule == grammar::RuleId::BE_WIP_WIP ||
         rule == grammar::RuleId::BE_WIP_BASEPAIR ||
         rule == grammar::RuleId::BE_BASEPAIR_WIP;
}

inline bool recursive_be_rule_has_left_wip(grammar::RuleId rule) {
  return rule == grammar::RuleId::BE_WIP_WIP ||
         rule == grammar::RuleId::BE_WIP_BASEPAIR;
}

inline bool recursive_be_rule_has_right_wip(grammar::RuleId rule) {
  return rule == grammar::RuleId::BE_WIP_WIP ||
         rule == grammar::RuleId::BE_BASEPAIR_WIP;
}

inline BeCoreKey be_core_key(const deductions::Deduction &deduction) {
  return BeCoreKey{deduction.split.k,
                   deduction.split.l,
                   deduction.parent.ip,
                   deduction.parent.jp};
}

// Exact/private partial-product cache shared by every recursive BE side-shell
// mode. It changes only factor ownership: the original rule-local weight stays
// on the final source contribution.
template <class Value>
class BasicExactBeFrameCache {
 public:
  template <class Semiring, class Chart>
  Value materialize_right_tail(const BeRightTailKey &key,
                               const Chart &chart) {
    const Value core = materialize_core<Semiring>(key.core, chart);
    const auto found = right_tail_factors_.find(key);
    if (found != right_tail_factors_.end()) {
      ++right_tail_reuse_count_;
      return found->second;
    }
    const Value factor = semirings::checkpoint_value(Semiring::times(
        core,
        chart.get(deductions::ItemKey{grammar::NonTerminal::WIP,
                                      key.core.l + 1,
                                      key.j - 1,
                                      -1,
                                      -1})));
    right_tail_factors_.emplace(key, factor);
    return factor;
  }

  template <class Semiring, class Chart>
  std::optional<Value> materialize_right_tail_for_source(
      const deductions::Deduction &deduction,
      const Chart &chart) {
    if (!recursive_be_rule_has_right_wip(deduction.rule)) {
      return std::nullopt;
    }
    return materialize_right_tail<Semiring>(
        BeRightTailKey{deduction.parent.j, be_core_key(deduction)}, chart);
  }

  template <class Semiring, class Chart, class LocalWeight>
  std::optional<Value> reconstruct_source(
      const deductions::Deduction &deduction,
      const Chart &chart,
      LocalWeight local_weight) {
    if (!is_recursive_be_frame_rule(deduction.rule)) return std::nullopt;
    Value factor = Semiring::zero();
    if (recursive_be_rule_has_right_wip(deduction.rule)) {
      factor = *materialize_right_tail_for_source<Semiring>(deduction, chart);
    } else {
      factor = materialize_core<Semiring>(deduction, chart);
    }
    if (recursive_be_rule_has_left_wip(deduction.rule)) {
      factor = Semiring::times(
          factor, chart.get(deduction.children[0].item));
    }
    ++reconstructed_source_count_;
    return Semiring::times(Semiring::from_weight(local_weight), factor);
  }

  std::size_t core_factor_count() const { return core_factors_.size(); }
  std::size_t core_reuse_count() const { return core_reuse_count_; }
  std::size_t right_tail_factor_count() const {
    return right_tail_factors_.size();
  }
  std::size_t right_tail_reuse_count() const {
    return right_tail_reuse_count_;
  }
  std::size_t reconstructed_source_count() const {
    return reconstructed_source_count_;
  }

 private:
  template <class Semiring, class Chart>
  Value materialize_core(const deductions::Deduction &deduction,
                         const Chart &chart) {
    const BeCoreKey key = be_core_key(deduction);
    return materialize_core<Semiring>(key, chart);
  }

  template <class Semiring, class Chart>
  Value materialize_core(const BeCoreKey &key, const Chart &chart) {
    const auto found = core_factors_.find(key);
    if (found != core_factors_.end()) {
      ++core_reuse_count_;
      return found->second;
    }
    const Value factor =
        chart.get(deductions::ItemKey{grammar::NonTerminal::BE,
                                      key.k,
                                      key.l,
                                      key.ip,
                                      key.jp});
    core_factors_.emplace(key, factor);
    return factor;
  }

  std::unordered_map<BeCoreKey, Value, BeCoreKeyHash> core_factors_;
  std::unordered_map<BeRightTailKey, Value, BeRightTailKeyHash>
      right_tail_factors_;
  std::size_t core_reuse_count_ = 0;
  std::size_t right_tail_reuse_count_ = 0;
  std::size_t reconstructed_source_count_ = 0;
};

using ExactBeFrameCache = BasicExactBeFrameCache<double>;

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_BE_FRAME_CACHE_HH
