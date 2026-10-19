#ifndef CPARTY_SCFG2_SEMIRINGS_SEMIRING_HH
#define CPARTY_SCFG2_SEMIRINGS_SEMIRING_HH

#include "scfg2/weights/local_weight.hh"

#include <algorithm>
#include <cmath>
#include <limits>

namespace cparty::scfg2::semirings {

struct SumProductSemiring {
  using value_type = double;

  static constexpr value_type zero() { return 0.0; }
  static constexpr value_type one() { return 1.0; }
  static constexpr value_type from_weight(double weight) { return weight; }
  static value_type from_log_weight(double value) { return std::exp(value); }
  static value_type from_weight(weights::LogWeight weight) { return from_log_weight(weight.value); }
  static constexpr bool is_zero(value_type value) { return value == zero(); }
  static constexpr double priority(value_type value) { return value; }
  static double traceback_priority(value_type value) {
    return is_zero(value) ? -std::numeric_limits<double>::infinity()
                          : std::log(value);
  }
  static constexpr value_type plus(value_type lhs, value_type rhs) { return lhs + rhs; }
  static constexpr value_type times(value_type lhs, value_type rhs) { return lhs * rhs; }
  static constexpr bool better(value_type lhs, value_type rhs) { return lhs > rhs; }
};

struct MaxProductSemiring {
  using value_type = double;

  static constexpr value_type zero() { return 0.0; }
  static constexpr value_type one() { return 1.0; }
  static constexpr value_type from_weight(double weight) { return weight; }
  static value_type from_log_weight(double value) { return std::exp(value); }
  static value_type from_weight(weights::LogWeight weight) { return from_log_weight(weight.value); }
  static constexpr bool is_zero(value_type value) { return value == zero(); }
  static constexpr double priority(value_type value) { return value; }
  static double traceback_priority(value_type value) {
    return is_zero(value) ? -std::numeric_limits<double>::infinity()
                          : std::log(value);
  }
  static constexpr value_type plus(value_type lhs, value_type rhs) { return std::max(lhs, rhs); }
  static constexpr value_type times(value_type lhs, value_type rhs) { return lhs * rhs; }
  static constexpr bool better(value_type lhs, value_type rhs) { return lhs > rhs; }
};

struct LogMaxSemiring {
  using value_type = double;

  static value_type zero() {
    return -std::numeric_limits<value_type>::infinity();
  }
  static constexpr value_type one() { return 0.0; }
  static value_type from_weight(double weight) { return std::log(weight); }
  static constexpr value_type from_log_weight(double value) { return value; }
  static value_type from_weight(weights::LogWeight weight) { return from_log_weight(weight.value); }
  static bool is_zero(value_type value) { return value == zero(); }
  static constexpr double priority(value_type value) { return value; }
  static constexpr double traceback_priority(value_type value) {
    return value;
  }
  static constexpr value_type plus(value_type lhs, value_type rhs) {
    return std::max(lhs, rhs);
  }
  static constexpr value_type times(value_type lhs, value_type rhs) {
    return lhs + rhs;
  }
  static constexpr bool better(value_type lhs, value_type rhs) {
    return lhs > rhs;
  }
};

struct LogSumExpSemiring {
  using value_type = double;

  static value_type zero() {
    return -std::numeric_limits<value_type>::infinity();
  }
  static constexpr value_type one() { return 0.0; }
  static value_type from_weight(double weight) { return std::log(weight); }
  static constexpr value_type from_log_weight(double value) { return value; }
  static value_type from_weight(weights::LogWeight weight) { return from_log_weight(weight.value); }
  static bool is_zero(value_type value) { return value == zero(); }
  static constexpr double priority(value_type value) { return value; }
  static constexpr double traceback_priority(value_type value) {
    return value;
  }
  static value_type plus(value_type lhs, value_type rhs) {
    if (is_zero(lhs)) return rhs;
    if (is_zero(rhs)) return lhs;
    const value_type high = std::max(lhs, rhs);
    const value_type low = std::min(lhs, rhs);
    return high + std::log1p(std::exp(low - high));
  }
  static constexpr value_type times(value_type lhs, value_type rhs) {
    return lhs + rhs;
  }
  static constexpr bool better(value_type lhs, value_type rhs) {
    return lhs > rhs;
  }
};

struct CoupledLogSumMaxValue {
  double log_sum = -std::numeric_limits<double>::infinity();
  double log_max = -std::numeric_limits<double>::infinity();
};

// Direct-product semiring for a single retained deduction DAG. Beam ranking
// follows log_sum; traceback follows log_max. Consequently log_max <= log_sum
// is preserved on that shared DAG.
struct CoupledLogSumMaxSemiring {
  using value_type = CoupledLogSumMaxValue;

  static value_type zero() { return {}; }
  static constexpr value_type one() { return {0.0, 0.0}; }
  static value_type from_weight(double weight) {
    const double log_weight = std::log(weight);
    return {log_weight, log_weight};
  }
  static value_type from_log_weight(double value) { return {value, value}; }
  static value_type from_weight(weights::LogWeight weight) { return from_log_weight(weight.value); }
  static bool is_zero(value_type value) {
    return value.log_sum == -std::numeric_limits<double>::infinity();
  }
  static constexpr double priority(value_type value) {
    return value.log_sum;
  }
  static constexpr double traceback_priority(value_type value) {
    return value.log_max;
  }
  static value_type plus(value_type lhs, value_type rhs) {
    return {LogSumExpSemiring::plus(lhs.log_sum, rhs.log_sum),
            std::max(lhs.log_max, rhs.log_max)};
  }
  static constexpr value_type times(value_type lhs, value_type rhs) {
    return {lhs.log_sum + rhs.log_sum, lhs.log_max + rhs.log_max};
  }
  static constexpr bool better(value_type lhs, value_type rhs) {
    return lhs.log_max > rhs.log_max;
  }
};

struct TargetedCoupledLogSumMaxValue {
  double log_sum = -std::numeric_limits<double>::infinity();
  double log_max = -std::numeric_limits<double>::infinity();
  double log_target = -std::numeric_limits<double>::infinity();
};

// Extends the shared inside/Viterbi retained-DAG value with the mass of one
// externally constrained output. A deduction-aware transform may clear only
// log_target; log_sum remains the beam priority and log_max remains the
// unconstrained traceback priority.
struct TargetedCoupledLogSumMaxSemiring {
  using value_type = TargetedCoupledLogSumMaxValue;

  static value_type zero() { return {}; }
  static constexpr value_type one() { return {0.0, 0.0, 0.0}; }
  static value_type from_weight(double weight) {
    const double log_weight = std::log(weight);
    return {log_weight, log_weight, log_weight};
  }
  static value_type from_log_weight(double value) { return {value, value, value}; }
  static value_type from_weight(weights::LogWeight weight) { return from_log_weight(weight.value); }
  static bool is_zero(value_type value) {
    return value.log_sum == -std::numeric_limits<double>::infinity();
  }
  static constexpr double priority(value_type value) {
    return value.log_sum;
  }
  static constexpr double traceback_priority(value_type value) {
    return value.log_max;
  }
  static value_type plus(value_type lhs, value_type rhs) {
    return {LogSumExpSemiring::plus(lhs.log_sum, rhs.log_sum),
            std::max(lhs.log_max, rhs.log_max),
            LogSumExpSemiring::plus(lhs.log_target, rhs.log_target)};
  }
  static constexpr value_type times(value_type lhs, value_type rhs) {
    return {lhs.log_sum + rhs.log_sum,
            lhs.log_max + rhs.log_max,
            lhs.log_target + rhs.log_target};
  }
  static constexpr bool better(value_type lhs, value_type rhs) {
    return lhs.log_max > rhs.log_max;
  }
  static value_type reject_target(value_type value) {
    value.log_target = -std::numeric_limits<double>::infinity();
    return value;
  }
};

}  // namespace cparty::scfg2::semirings

#endif  // CPARTY_SCFG2_SEMIRINGS_SEMIRING_HH
