#ifndef CPARTY_SCFG2_WEIGHTS_LOCAL_WEIGHT_HH
#define CPARTY_SCFG2_WEIGHTS_LOCAL_WEIGHT_HH

#include <cmath>
#include <limits>
#include <type_traits>
#include <variant>

namespace cparty::scfg2::weights {

// No implicit conversion: a log weight must never enter a raw-weight API.
struct LogWeight { double value; };

using AnyLocalWeight = std::variant<double, LogWeight>;

template <class Semiring, class = void>
struct HasLogEntry : std::false_type {};
template <class Semiring>
struct HasLogEntry<Semiring, std::void_t<decltype(Semiring::from_log_weight(0.))>> : std::true_type {};
template <class Semiring>
auto inject_local(const AnyLocalWeight &weight) {
  return std::visit([](auto local) {
    if constexpr (std::is_same_v<decltype(local), LogWeight>) {
      if constexpr (HasLogEntry<Semiring>::value) return Semiring::from_log_weight(local.value);
      else return Semiring::from_weight(std::exp(local.value)); // Historical test semirings only.
    } else return Semiring::from_weight(local);
  }, weight);
}

inline constexpr double log_zero = -std::numeric_limits<double>::infinity();
inline bool is_zero(double weight) { return weight == 0.0; }
inline bool is_zero(LogWeight weight) { return weight.value == log_zero; }
inline double log_weight(double weight) { return std::log(weight); }
inline double log_weight(LogWeight weight) { return weight.value; }
inline double raw_weight(double weight) { return weight; }
inline double raw_weight(LogWeight weight) { return std::exp(weight.value); }

template <class Context, class = void>
struct LogContext : std::false_type {};
template <class Context>
struct LogContext<Context, std::void_t<decltype(Context::log_local_weights)>>
    : std::bool_constant<Context::log_local_weights> {};

template <class Context>
using LocalWeight = std::conditional_t<LogContext<Context>::value, LogWeight, double>;

template <class Context>
LocalWeight<Context> local_weight(const Context &, double value) {
  if constexpr (LogContext<Context>::value) return LogWeight{value};
  else return value;
}
template <class Context>
constexpr double impossible_weight() {
  return LogContext<Context>::value ? log_zero : 0.0;
}

}  // namespace cparty::scfg2::weights
#endif
