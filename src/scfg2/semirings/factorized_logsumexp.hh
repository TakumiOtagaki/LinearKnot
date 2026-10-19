#ifndef CPARTY_SCFG2_SEMIRINGS_FACTORIZED_LOGSUMEXP_HH
#define CPARTY_SCFG2_SEMIRINGS_FACTORIZED_LOGSUMEXP_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/semirings/semiring.hh"

#include <cstdint>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cparty::scfg2::semirings {

inline constexpr std::uint32_t invalid_outside_state =
    std::numeric_limits<std::uint32_t>::max();

class FactorizedInsideGraph;
struct FactorizedExpression;

struct FactorizedLogSumExpValue {
  double log_value = -std::numeric_limits<double>::infinity();
  std::uint32_t node = invalid_outside_state;
  std::shared_ptr<FactorizedExpression> expression;

  FactorizedLogSumExpValue checkpointed() const;
  bool discard_chart_entry() const {
    return log_value == -std::numeric_limits<double>::infinity();
  }
};

inline bool operator<(const FactorizedLogSumExpValue &lhs,
                      const FactorizedLogSumExpValue &rhs) {
  return lhs.log_value < rhs.log_value;
}

// Each state is an item or a consumed/observed private aggregate:
//   alpha(P) = logsumexp_t c(t)
//   c(t) = local(t) + sum_C alpha(C).
// Unobserved local constants, multiplication chains, and logadd chains are
// flattened; a semantic leaf is a term with no child factors.
struct FactorizedInsideState {
  double log_value;
  std::uint32_t first_term;
  std::uint32_t term_count;
};

struct FactorizedInsideTerm {
  double log_value;
  std::uint32_t first_factor;
  std::uint32_t factor_count;
  // Observed local factor, independent of child inside sums. Thermodynamic
  // consumers must not use log_value as an edge weight or recover this by
  // subtracting large inside values.
  double local_log_weight;
};

class FactorizedInsideGraph {
 public:
  const std::vector<FactorizedInsideState> &states() const { return states_; }
  const std::vector<FactorizedInsideTerm> &terms() const { return terms_; }
  const std::vector<std::uint32_t> &factors() const { return factors_; }

  std::size_t storage_bytes() const {
    return states_.capacity() * sizeof(FactorizedInsideState) +
           terms_.capacity() * sizeof(FactorizedInsideTerm) +
           factors_.capacity() * sizeof(std::uint32_t);
  }

  // Analytic beta recursion on the frozen forward sums/products. Repeated
  // factors and shared aggregates receive every incoming contribution once.
  std::vector<double> outside_log_beta(std::uint32_t root) const {
    if (root >= states_.size()) {
      throw std::invalid_argument("factorized outside root is invalid");
    }
    std::vector<double> beta(states_.size(), LogSumExpSemiring::zero());
    beta[root] = 0.0;
    for (std::size_t cursor = states_.size(); cursor > 0; --cursor) {
      const auto index = cursor - 1;
      if (LogSumExpSemiring::is_zero(beta[index])) continue;
      const auto &parent = states_[index];
      for (std::size_t t = parent.first_term;
           t < std::size_t(parent.first_term) + parent.term_count; ++t) {
        const auto &term = terms_[t];
        for (std::size_t f = term.first_factor;
             f < std::size_t(term.first_factor) + term.factor_count; ++f) {
          const auto child = factors_[f];
          beta[child] = LogSumExpSemiring::plus(
              beta[child], beta[index] + term.log_value - states_[child].log_value);
        }
      }
    }
    return beta;
  }

  std::vector<double> posteriors(std::uint32_t root) const {
    if (root >= states_.size()) {
      throw std::invalid_argument("factorized outside root is invalid");
    }
    // Normalized beta: post(P) * exp(c(t)-alpha(P)) is the term's
    // posterior, passed to each factor. This avoids subtracting a large
    // root log Z at every state and computes one exponential per term.
    std::vector<double> out(states_.size(), 0.0);
    out[root] = 1.0;
    for (std::size_t cursor = states_.size(); cursor > 0; --cursor) {
      const auto index = cursor - 1;
      if (out[index] == 0.0) continue;
      const auto &state = states_[index];
      for (std::size_t t = state.first_term;
           t < std::size_t(state.first_term) + state.term_count; ++t) {
        const auto &term = terms_[t];
        const double probability = out[index] * std::exp(term.log_value - state.log_value);
        for (std::size_t f = term.first_factor;
             f < std::size_t(term.first_factor) + term.factor_count; ++f)
          out[factors_[f]] += probability;
      }
    }
    return out;
  }

  // Frozen-graph finite differences: selection is never re-run.
  double evaluate_with_log_shift(std::uint32_t root, std::uint32_t shifted,
                                 double delta) const {
    if (root >= states_.size() || shifted >= states_.size()) {
      throw std::invalid_argument("factorized outside shift state is invalid");
    }
    std::vector<double> alpha(states_.size());
    for (std::size_t i = 0; i < states_.size(); ++i) {
      if (i == shifted) {
        alpha[i] = states_[i].log_value + delta;
        continue;
      }
      double total = LogSumExpSemiring::zero();
      const auto &state = states_[i];
      for (std::size_t t = state.first_term;
           t < std::size_t(state.first_term) + state.term_count; ++t) {
        const auto &term = terms_[t];
        double value = term.log_value;
        for (std::size_t f = term.first_factor;
             f < std::size_t(term.first_factor) + term.factor_count; ++f) {
          const auto child = factors_[f];
          value += alpha[child] - states_[child].log_value;
        }
        total = LogSumExpSemiring::plus(total, value);
      }
      alpha[i] = total;
    }
    return alpha[root];
  }

 private:
  friend struct FactorizedLogSumExpValue;
  std::uint32_t materialize(FactorizedExpression &expression);

  static std::uint32_t checked_index(std::size_t size) {
    if (size >= invalid_outside_state) {
      throw std::length_error("factorized outside exceeds 32-bit record capacity");
    }
    return static_cast<std::uint32_t>(size);
  }

  std::vector<FactorizedInsideState> states_;
  std::vector<FactorizedInsideTerm> terms_;
  std::vector<std::uint32_t> factors_;
};

// Only unsealed candidates use expressions. Products contain state IDs;
// a running sum shares a vector of candidate terms with a prefix length,
// preserving old snapshots without copying its entire prefix at every +.
struct FactorizedExpression {
  FactorizedInsideGraph *graph = nullptr;
  double log_value = LogSumExpSemiring::zero();
  double local_log_weight = 0.0;
  std::uint32_t state = invalid_outside_state;
  std::vector<std::uint32_t> factors;
  std::shared_ptr<std::vector<FactorizedLogSumExpValue>> sum;
  std::size_t term_count = 0;
};

inline FactorizedLogSumExpValue FactorizedLogSumExpValue::checkpointed() const {
  if (node != invalid_outside_state || !expression) return *this;
  return {log_value, expression->graph->materialize(*expression), {}};
}

inline std::uint32_t FactorizedInsideGraph::materialize(
    FactorizedExpression &expression) {
  if (expression.state != invalid_outside_state) return expression.state;
  std::size_t term_count = expression.sum ? 0 : 1;
  if (expression.sum) {
    // Seal every nested/shared sum in the old left-to-right order before
    // recording this parent's first term. Recursive child terms must not
    // split the parent's contiguous term range.
    for (std::size_t i = 0; i < expression.term_count; ++i) {
      const auto &value = (*expression.sum)[i];
      if (value.discard_chart_entry()) continue;
      ++term_count;
      if (value.node == invalid_outside_state &&
          (value.expression->state != invalid_outside_state || value.expression->sum))
        (void)value.checkpointed();
    }
  }
  const auto first_term = checked_index(terms_.size());
  checked_index(terms_.size() + term_count);
  const auto append_product = [&](double log_value,
      const std::vector<std::uint32_t> &factors, double local_log_weight) {
    const auto first_factor = checked_index(factors_.size());
    checked_index(factors_.size() + factors.size());
    factors_.insert(factors_.end(), factors.begin(), factors.end());
    terms_.push_back({log_value, first_factor, checked_index(factors.size()), local_log_weight});
  };
  const auto append_state = [&](double log_value, std::uint32_t node) {
    const auto first_factor = checked_index(factors_.size());
    checked_index(factors_.size() + 1);
    factors_.push_back(node);
    terms_.push_back({log_value, first_factor, 1, 0.0});
  };
  if (expression.sum) {
    for (std::size_t i = 0; i < expression.term_count; ++i) {
      const auto &value = (*expression.sum)[i];
      if (value.discard_chart_entry()) continue;
      if (value.node != invalid_outside_state) {
        append_state(value.log_value, value.node);
      } else if (value.expression->state != invalid_outside_state) {
        append_state(value.log_value, value.expression->state);
      } else {
        // Plain products were never sealed by the first pass; their factors
        // stay owned by the sum's values until this parent is complete.
        append_product(value.log_value, value.expression->factors,
                       value.expression->local_log_weight);
      }
    }
  } else {
    append_product(expression.log_value, expression.factors, expression.local_log_weight);
  }
  const auto state = checked_index(states_.size());
  states_.push_back({expression.log_value, first_term, checked_index(term_count)});
  expression.state = state;
  expression.sum.reset();
  std::vector<std::uint32_t>().swap(expression.factors);
  return state;
}

struct FactorizedLogSumExpSemiring {
  using value_type = FactorizedLogSumExpValue;

  class Scope {
   public:
    explicit Scope(FactorizedInsideGraph &graph) : previous_(active_graph_) {
      if (previous_) throw std::logic_error("nested factorized inside scopes are unsupported");
      active_graph_ = &graph;
    }
    ~Scope() { active_graph_ = previous_; }
    Scope(const Scope &) = delete;
    Scope &operator=(const Scope &) = delete;
   private:
    FactorizedInsideGraph *previous_;
  };

  static value_type zero() { return {}; }
  static value_type one() { return constant(0.0); }
  static value_type from_weight(double weight) { return constant(std::log(weight)); }
  static value_type from_log_weight(double value) { return constant(value); }
  static value_type from_weight(weights::LogWeight weight) { return from_log_weight(weight.value); }
  static bool is_zero(const value_type &value) { return value.discard_chart_entry(); }
  static double priority(const value_type &value) { return value.log_value; }
  static double traceback_priority(const value_type &value) { return value.log_value; }
  static bool better(const value_type &lhs, const value_type &rhs) {
    return lhs.log_value > rhs.log_value;
  }

  static value_type times(const value_type &lhs, const value_type &rhs) {
    if (is_zero(lhs) || is_zero(rhs)) return zero();
    auto out = constant(lhs.log_value + rhs.log_value);
    out.expression->local_log_weight = 0.0;
    const auto append = [&](const value_type &value) {
      auto &factors = out.expression->factors;
      if (value.node != invalid_outside_state) {
        factors.push_back(value.node);
      } else if (value.expression->state != invalid_outside_state ||
                 value.expression->sum) {
        factors.push_back(value.checkpointed().node);
      } else {
        factors.insert(factors.end(), value.expression->factors.begin(),
                       value.expression->factors.end());
        out.expression->local_log_weight += value.expression->local_log_weight;
      }
    };
    append(lhs);
    append(rhs);
    return out;
  }

  static value_type plus(const value_type &lhs, const value_type &rhs) {
    if (is_zero(lhs)) return rhs;
    if (is_zero(rhs)) return lhs;
    auto out = constant(LogSumExpSemiring::plus(lhs.log_value, rhs.log_value));
    auto &expression = *out.expression;
    if (lhs.expression && lhs.expression->state == invalid_outside_state &&
        lhs.expression->sum &&
        lhs.expression->term_count == lhs.expression->sum->size() &&
        lhs.expression.use_count() == 1) {
      expression.sum = lhs.expression->sum;
    } else {
      expression.sum = std::make_shared<std::vector<value_type>>();
      // A reused sum is a private state, not a second copy of its seeds.
      // Plain product terms stay transient until the whole sum is sealed.
      expression.sum->push_back(lhs.expression && lhs.expression->sum
                                    ? lhs.checkpointed() : lhs);
    }
    expression.sum->push_back(rhs);
    expression.term_count = expression.sum->size();
    return out;
  }

 private:
  static value_type constant(double log_value) {
    if (LogSumExpSemiring::is_zero(log_value)) return zero();
    if (!active_graph_) throw std::logic_error("factorized semiring has no active scope");
    auto expression = std::make_shared<FactorizedExpression>();
    expression->graph = active_graph_;
    expression->log_value = log_value;
    expression->local_log_weight = log_value;
    return {log_value, invalid_outside_state, std::move(expression)};
  }
  inline static thread_local FactorizedInsideGraph *active_graph_ = nullptr;
};

}  // namespace cparty::scfg2::semirings

#endif
