#include "scfg2/semirings/logsumexp_tape.hh"

#include <cassert>
#include <cmath>
#include <iostream>

namespace {

bool close(double lhs, double rhs, double tolerance = 1e-12) {
  return std::fabs(lhs - rhs) <= tolerance;
}

double expression(double x, double y, double z) {
  return std::log(std::exp(x + y) + std::exp(x + z));
}

}  // namespace

int main() {
  using cparty::scfg2::semirings::LogSumExpTape;
  using cparty::scfg2::semirings::TapeLogSumExpSemiring;

  LogSumExpTape tape;
  TapeLogSumExpSemiring::value_type root;
  TapeLogSumExpSemiring::value_type x;
  TapeLogSumExpSemiring::value_type y;
  TapeLogSumExpSemiring::value_type z;
  {
    TapeLogSumExpSemiring::Scope scope(tape);
    x = TapeLogSumExpSemiring::from_weight(std::exp(0.3));
    y = TapeLogSumExpSemiring::from_weight(std::exp(-0.7));
    z = TapeLogSumExpSemiring::from_weight(std::exp(1.1));
    root = TapeLogSumExpSemiring::plus(
        TapeLogSumExpSemiring::times(x, y),
        TapeLogSumExpSemiring::times(x, z));
  }

  assert(close(root.log_value, expression(0.3, -0.7, 1.1)));
  const auto adjoints = tape.reverse(root.node);
  assert(close(adjoints[x.node], 1.0));

  constexpr double epsilon = 1e-6;
  const auto finite_difference = [&](double x0, double y0, double z0,
                                     int variable) {
    double *selected = variable == 0 ? &x0 : variable == 1 ? &y0 : &z0;
    *selected += epsilon;
    const double high = expression(x0, y0, z0);
    *selected -= 2.0 * epsilon;
    const double low = expression(x0, y0, z0);
    return (high - low) / (2.0 * epsilon);
  };
  assert(close(adjoints[x.node],
               finite_difference(0.3, -0.7, 1.1, 0),
               1e-9));
  assert(close(adjoints[y.node],
               finite_difference(0.3, -0.7, 1.1, 1),
               1e-9));
  assert(close(adjoints[z.node],
               finite_difference(0.3, -0.7, 1.1, 2),
               1e-9));
  assert(close(adjoints[y.node] + adjoints[z.node], 1.0));
  const auto shifted_difference = [&](std::size_t node) {
    constexpr double shift = 1e-5;
    return (tape.evaluate_with_log_shift(root.node, node, shift) -
            tape.evaluate_with_log_shift(root.node, node, -shift)) /
           (2.0 * shift);
  };
  assert(close(adjoints[x.node], shifted_difference(x.node), 1e-9));
  assert(close(adjoints[y.node], shifted_difference(y.node), 1e-9));
  assert(close(adjoints[z.node], shifted_difference(z.node), 1e-9));

  std::cout << "scfg2_logsumexp_tape_test\t1\n";
  return 0;
}
