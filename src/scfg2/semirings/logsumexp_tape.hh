#ifndef CPARTY_SCFG2_SEMIRINGS_LOGSUMEXP_TAPE_HH
#define CPARTY_SCFG2_SEMIRINGS_LOGSUMEXP_TAPE_HH

#include "scfg2/weights/local_weight.hh"

#include "scfg2/semirings/semiring.hh"

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cparty::scfg2::semirings {

inline constexpr std::size_t invalid_tape_node =
    std::numeric_limits<std::size_t>::max();

struct LogSumExpTapeValue {
  double log_value = -std::numeric_limits<double>::infinity();
  std::size_t node = invalid_tape_node;
};

inline constexpr bool operator<(LogSumExpTapeValue lhs,
                                LogSumExpTapeValue rhs) {
  return lhs.log_value < rhs.log_value;
}

struct LogSumExpTapeNode {
  enum class Operation : unsigned char { Constant, Product, LogAdd };

  Operation operation = Operation::Constant;
  double log_value = -std::numeric_limits<double>::infinity();
  std::size_t lhs = invalid_tape_node;
  std::size_t rhs = invalid_tape_node;
};

class LogSumExpTape {
 public:
  LogSumExpTapeValue constant(double log_value) {
    return append(LogSumExpTapeNode::Operation::Constant,
                  log_value,
                  invalid_tape_node,
                  invalid_tape_node);
  }

  LogSumExpTapeValue product(LogSumExpTapeValue lhs,
                             LogSumExpTapeValue rhs) {
    require_node(lhs);
    require_node(rhs);
    return append(LogSumExpTapeNode::Operation::Product,
                  lhs.log_value + rhs.log_value,
                  lhs.node,
                  rhs.node);
  }

  LogSumExpTapeValue log_add(LogSumExpTapeValue lhs,
                             LogSumExpTapeValue rhs) {
    require_node(lhs);
    require_node(rhs);
    return append(LogSumExpTapeNode::Operation::LogAdd,
                  LogSumExpSemiring::plus(lhs.log_value, rhs.log_value),
                  lhs.node,
                  rhs.node);
  }

  std::vector<double> reverse(std::size_t root) const {
    if (root >= nodes_.size()) {
      throw std::invalid_argument("outside tape root is invalid");
    }
    std::vector<double> adjoints(nodes_.size(), 0.0);
    adjoints[root] = 1.0;
    for (std::size_t cursor = nodes_.size(); cursor > 0; --cursor) {
      const std::size_t index = cursor - 1;
      const auto &node = nodes_[index];
      const double adjoint = adjoints[index];
      if (adjoint == 0.0 ||
          node.operation == LogSumExpTapeNode::Operation::Constant) {
        continue;
      }
      if (node.operation == LogSumExpTapeNode::Operation::Product) {
        adjoints[node.lhs] += adjoint;
        adjoints[node.rhs] += adjoint;
        continue;
      }
      adjoints[node.lhs] +=
          adjoint * std::exp(nodes_[node.lhs].log_value - node.log_value);
      adjoints[node.rhs] +=
          adjoint * std::exp(nodes_[node.rhs].log_value - node.log_value);
    }
    return adjoints;
  }

  double evaluate_with_log_shift(std::size_t root,
                                 std::size_t shifted_node,
                                 double delta) const {
    if (root >= nodes_.size() || shifted_node >= nodes_.size()) {
      throw std::invalid_argument("outside tape shift node is invalid");
    }
    std::vector<double> values;
    values.reserve(nodes_.size());
    for (std::size_t index = 0; index < nodes_.size(); ++index) {
      const auto &node = nodes_[index];
      if (index == shifted_node) {
        values.push_back(node.log_value + delta);
      } else if (node.operation == LogSumExpTapeNode::Operation::Constant) {
        values.push_back(node.log_value);
      } else if (node.operation == LogSumExpTapeNode::Operation::Product) {
        values.push_back(values[node.lhs] + values[node.rhs]);
      } else {
        values.push_back(
            LogSumExpSemiring::plus(values[node.lhs], values[node.rhs]));
      }
    }
    return values[root];
  }

  const std::vector<LogSumExpTapeNode> &nodes() const { return nodes_; }

 private:
  static void require_node(LogSumExpTapeValue value) {
    if (value.node == invalid_tape_node) {
      throw std::logic_error("finite outside tape value has no node");
    }
  }

  LogSumExpTapeValue append(LogSumExpTapeNode::Operation operation,
                            double log_value,
                            std::size_t lhs,
                            std::size_t rhs) {
    const std::size_t node = nodes_.size();
    nodes_.push_back(LogSumExpTapeNode{operation, log_value, lhs, rhs});
    return {log_value, node};
  }

  std::vector<LogSumExpTapeNode> nodes_;
};

struct TapeLogSumExpSemiring {
  using value_type = LogSumExpTapeValue;

  class Scope {
   public:
    explicit Scope(LogSumExpTape &tape) : previous_(active_tape_) {
      if (previous_) {
        throw std::logic_error("nested outside tape scopes are unsupported");
      }
      active_tape_ = &tape;
    }
    ~Scope() { active_tape_ = previous_; }
    Scope(const Scope &) = delete;
    Scope &operator=(const Scope &) = delete;

   private:
    LogSumExpTape *previous_ = nullptr;
  };

  static value_type zero() { return {}; }
  static value_type one() { return tape().constant(0.0); }
  static value_type from_weight(double weight) {
    return tape().constant(std::log(weight));
  }
  static value_type from_log_weight(double value) { return tape().constant(value); }
  static value_type from_weight(weights::LogWeight weight) { return from_log_weight(weight.value); }
  static bool is_zero(value_type value) {
    return value.log_value ==
           -std::numeric_limits<double>::infinity();
  }
  static double priority(value_type value) { return value.log_value; }
  static double traceback_priority(value_type value) {
    return value.log_value;
  }
  static value_type plus(value_type lhs, value_type rhs) {
    if (is_zero(lhs)) return rhs;
    if (is_zero(rhs)) return lhs;
    return tape().log_add(lhs, rhs);
  }
  static value_type times(value_type lhs, value_type rhs) {
    if (is_zero(lhs) || is_zero(rhs)) return zero();
    return tape().product(lhs, rhs);
  }
  static bool better(value_type lhs, value_type rhs) {
    return lhs.log_value > rhs.log_value;
  }

 private:
  static LogSumExpTape &tape() {
    if (!active_tape_) {
      throw std::logic_error("outside tape semiring has no active scope");
    }
    return *active_tape_;
  }

  inline static thread_local LogSumExpTape *active_tape_ = nullptr;
};

}  // namespace cparty::scfg2::semirings

#endif  // CPARTY_SCFG2_SEMIRINGS_LOGSUMEXP_TAPE_HH
