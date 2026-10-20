#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <optional>
#include <random>
#include <string>
#include <utility>
#include <vector>

// One-shot faults isolate the optional support allocations from authoritative
// value/presence allocations, without relying on library allocation counts.
namespace allocation_fault {
thread_local long remaining = -1;
thread_local bool fired = false;
void arm(long successful_allocations) {
  remaining = successful_allocations;
  fired = false;
}
void disable() { remaining = -1; }
void *allocate(std::size_t bytes) {
  if (remaining >= 0 && remaining-- == 0) {
    remaining = -1;
    fired = true;
    throw std::bad_alloc{};
  }
  if (void *memory = std::malloc(bytes ? bytes : 1)) return memory;
  throw std::bad_alloc{};
}
}  // namespace allocation_fault
void *operator new(std::size_t bytes) { return allocation_fault::allocate(bytes); }
void *operator new[](std::size_t bytes) { return allocation_fault::allocate(bytes); }
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete(void *memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory, std::size_t) noexcept { std::free(memory); }

namespace {
namespace sc = cparty::scfg2;
using Chart = sc::engine::Chart<double>;
using Key = sc::deductions::ItemKey;
using NT = sc::grammar::NonTerminal;
using R = sc::grammar::RuleId;
using S = sc::semirings::LogSumExpSemiring;
using M = sc::semirings::LogMaxSemiring;
using Scratch = sc::replay::ExactNumericSplitScratch;
constexpr double zero = -std::numeric_limits<double>::infinity();
constexpr int families = static_cast<int>(NT::BE);
std::size_t compared_values = 0, compared_rows = 0, compared_selections = 0;

void require(bool ok, const char *message) {
  if (!ok) {
    allocation_fault::disable();
    std::cerr << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(a)) == 0; }
bool discard_marker(const double &value) { return value == 99.0; }

const std::array<double, 8> values{{zero, -1e300, 0.0, -0.0,
    std::numeric_limits<double>::quiet_NaN(),
    std::numeric_limits<double>::infinity(), -7.125, 99.0}};

void compare_key(const Chart &expected, const Chart &actual, const Key &key) {
  require(expected.contains(key) == actual.contains(key), "support index changed presence");
  require(same(expected.get(key), actual.get(key)), "support index changed value bits");
  ++compared_values;
}

void compare_entries(const Chart &expected, const Chart &actual) {
  require(expected.size() == actual.size(), "support index changed chart size");
  std::size_t seen = 0;
  actual.for_each_entry([&](const Key &key, double value) {
    require(expected.contains(key) && same(expected.get(key), value),
            "support index changed full entry identity/value");
    ++seen;
  });
  require(seen == expected.size(), "support index changed entry count");
  const auto &view = actual.values();
  require(view.size() == expected.size(), "support index changed materialized size");
  expected.for_each_entry([&](const Key &key, double value) {
    const auto found = view.find(key);
    require(found != view.end() && same(found->second, value),
            "support index changed materialized full keys/value bits");
  });
}

void compare_start_row(const Chart &expected, const Chart &actual, cand_pos_t n,
                       NT first_family, cand_pos_t first_end,
                       NT second_family, cand_pos_t second_end,
                       cand_pos_t begin, cand_pos_t end) {
  std::vector<cand_pos_t> want, got;
  for (cand_pos_t i = begin; i <= end; ++i)
    if ((i <= first_end && !S::is_zero(expected.get({first_family, i, first_end, -1, -1}))) ||
        (i <= second_end && !S::is_zero(expected.get({second_family, i, second_end, -1, -1}))))
      want.push_back(i);
  require(actual.for_each_live_primary_start_union(first_family, first_end,
              second_family, second_end, begin, end,
              [&](cand_pos_t i) { got.push_back(i); }), "canonical start row requested fallback");
  require(got == want, "support start union omitted/duplicated/reordered a live coordinate");
  require(std::is_sorted(got.begin(), got.end()) &&
              std::adjacent_find(got.begin(), got.end()) == got.end(),
          "support start union is not strictly ascending");
  require(actual.live_primary_index_length() == n, "support query changed index ownership");
  ++compared_rows;
}

void compare_end_row(const Chart &expected, const Chart &actual, NT family,
                     cand_pos_t start, cand_pos_t begin, cand_pos_t end) {
  std::vector<cand_pos_t> want, got;
  for (cand_pos_t j = begin; j <= end; ++j)
    if (j >= start && !S::is_zero(expected.get({family, start, j, -1, -1}))) want.push_back(j);
  require(actual.for_each_live_primary_end(family, start, begin, end,
              [&](cand_pos_t j) { got.push_back(j); }), "canonical end row requested fallback");
  require(got == want, "support end row omitted/duplicated/reordered a live coordinate");
  ++compared_rows;
}

void compare_rows(const Chart &expected, const Chart &actual, cand_pos_t n) {
  for (int f = 0; f < families; ++f) {
    const auto nt = static_cast<NT>(f), next = static_cast<NT>((f + 1) % families);
    for (cand_pos_t fixed = 1; fixed <= n; ++fixed) {
      compare_start_row(expected, actual, n, nt, fixed, next, std::min(n, cand_pos_t(fixed + 1)), 1, n);
      compare_end_row(expected, actual, nt, fixed, 1, n);
    }
    // Independent partial masks on both sides of byte and machine-word edges.
    for (cand_pos_t edge : {1, 7, 8, 9, 62, 63, 64, 65, 66, 126, 127, 128, 129}) {
      if (edge > n) continue;
      compare_start_row(expected, actual, n, nt, n, next, n, edge, n);
      compare_start_row(expected, actual, n, nt, n, nt, n, 1, edge);
      compare_start_row(expected, actual, n, nt, n, next, n, edge, edge);
      compare_end_row(expected, actual, nt, 1, edge, n);
      compare_end_row(expected, actual, nt, 1, 1, edge);
      compare_end_row(expected, actual, nt, 1, edge, edge);
    }
  }
  std::size_t called = 0;
  auto count = [&](cand_pos_t) { ++called; };
  require(actual.for_each_live_primary_start_union(NT::V, n, NT::WMB, n, 2, 1, count) &&
              actual.for_each_live_primary_end(NT::VP, 1, 2, 1, count) && !called,
          "empty support range emitted coordinates");
  for (const auto &invalid : std::vector<std::pair<NT, cand_pos_t>>{
           {NT::BE, 1}, {static_cast<NT>(255), 1}, {NT::V, 0}, {NT::V, n + 1}}) {
    require(!actual.for_each_live_primary_start_union(invalid.first, invalid.second,
                NT::V, n, 1, n, count) &&
                !actual.for_each_live_primary_end(invalid.first, invalid.second, 1, n, count),
            "invalid support row did not request fallback");
  }
  require(!actual.for_each_live_primary_start_union(NT::V, n, NT::WMB, n, 0, n, count) &&
              !actual.for_each_live_primary_end(NT::VP, 1, 1, n + 1, count) && !called,
          "unsupported support range emitted a partial result");
}

void storage_boundaries(cand_pos_t n) {
  Chart expected(zero, discard_marker), actual(zero, discard_marker, n);
  require(actual.live_primary_index_length() == n && expected.live_primary_index_bytes() == 0,
          "dense/sparse support opt-in configuration changed");
  require(actual.live_primary_index_bytes() ==
              2 * static_cast<std::size_t>(families) * n * ((n + 63) / 64) * sizeof(std::uint64_t),
          "support index is not bounded by the documented two bit matrices");
  for (int f = 0; f < families; ++f)
    for (cand_pos_t j = 1; j <= n; ++j)
      for (cand_pos_t i = 1; i <= j; ++i) {
        const Key key{static_cast<NT>(f), i, j, -1, -1};
        const auto value = values[(i * 7 + j * 3 + f) % values.size()];
        expected.set(key, value); actual.set(key, value);
        compare_key(expected, actual, key);
      }
  compare_rows(expected, actual, n);
  // Full keys outside the primary route must never contaminate canonical bits.
  const std::vector<Key> fallback{{NT::BE, 1, n, 2, n - 1}, {NT::BE, 1, n, -1, -1},
      {NT::V, 1, n, 0, -1}, {NT::V, 1, n, -1, 0}, {NT::V, 1, n, 2, n - 1},
      {NT::WM, 0, n, -1, -1}, {NT::WI, -2, n, -1, -1}, {NT::WIP, 2, 1, -1, -1},
      {NT::VP, 1, n + 1, -1, -1}, {static_cast<NT>(255), 1, n, -1, -1}};
  for (const auto &key : fallback) {
    expected.set(key, -1e300); actual.set(key, -1e300);
    compare_key(expected, actual, key);
  }
  compare_rows(expected, actual, n);
  compare_entries(expected, actual);
  const auto &materialized = actual.values();
  // Clear and revive an endpoint in every family, preserving explicit-zero
  // presence and the distinct signed-zero bit patterns in the map view.
  for (int f = 0; f < families; ++f) {
    const Key key{static_cast<NT>(f), 1, n, -1, -1};
    for (const auto value : values) {
      expected.set(key, value); actual.set(key, value);
      compare_key(expected, actual, key);
      require((materialized.find(key) != materialized.end()) == expected.contains(key),
              "support mutation stopped updating the compatibility view");
      if (expected.contains(key)) require(same(materialized.at(key), expected.get(key)),
          "support mutation changed compatibility-view bits");
      compare_start_row(expected, actual, n, static_cast<NT>(f), n, static_cast<NT>(f), n, 1, n);
      compare_end_row(expected, actual, static_cast<NT>(f), 1, 1, n);
    }
  }
  compare_entries(expected, actual);
}

void copy_move_contract() {
  constexpr cand_pos_t n = 65;
  Chart expected(zero), source(zero, nullptr, n);
  for (int f = 0; f < families; ++f)
    for (cand_pos_t j = 1; j <= n; ++j) {
      const Key key{static_cast<NT>(f), 1, j, -1, -1};
      expected.set(key, values[(f + j) % (values.size() - 1)]);
      source.set(key, expected.get(key));
    }
  source.set({NT::BE, 1, n, 2, 63}, -1e300);
  expected.set({NT::BE, 1, n, 2, 63}, -1e300);
  Chart copied(source), assigned(17.0, nullptr, 3);
  assigned = source;
  source.set({NT::V, 1, n, -1, -1}, zero);
  source.set({NT::VP, 2, n, -1, -1}, 42.0);
  compare_entries(expected, copied); compare_rows(expected, copied, n);
  compare_entries(expected, assigned); compare_rows(expected, assigned, n);
  Chart moved(std::move(copied));
  require(copied.live_primary_index_length() == 0 && copied.live_primary_index_bytes() == 0 &&
              copied.dense_length() == 0, "move construction retained support ownership");
  compare_entries(expected, moved); compare_rows(expected, moved, n);
  Chart destination(-0.0, nullptr, 7);
  destination = std::move(assigned);
  require(assigned.live_primary_index_length() == 0 && assigned.live_primary_index_bytes() == 0,
          "move assignment retained support ownership");
  compare_entries(expected, destination); compare_rows(expected, destination, n);
  destination = destination;
  destination = std::move(destination);
  compare_rows(expected, destination, n);
  assigned.set({NT::V, 1, 1, -1, -1}, -0.0);
  require(assigned.get({NT::V, 1, 1, -1, -1}) == 0.0 &&
              assigned.live_primary_index_bytes() == 0, "moved-from sparse chart cannot be reused");
}

template <class NumericChart>
void require_no_index(const NumericChart &chart) {
  require(chart.live_primary_index_length() == 0 && chart.live_primary_index_bytes() == 0,
          "unsupported chart allocated dense support rows");
  std::size_t calls = 0;
  require(!chart.for_each_live_primary_start_union(NT::V, 3, NT::WMB, 3, 1, 3,
              [&](cand_pos_t) { ++calls; }) &&
              !chart.for_each_live_primary_end(NT::VP, 1, 1, 3,
              [&](cand_pos_t) { ++calls; }) && !calls,
          "unsupported chart did not request a complete fallback");
}

void default_contract() {
  for (const double default_value : {0.0, -0.0, 17.0,
           std::numeric_limits<double>::quiet_NaN()}) {
    Chart chart(default_value, nullptr, 65);
    require_no_index(chart);
    require(same(chart.get({NT::V, 1, 65, -1, -1}), default_value),
            "support policy changed a custom default");
  }
  Chart sparse(zero), beam(zero, nullptr, 0, 65);
  require_no_index(sparse); require_no_index(beam);
  sc::engine::Chart<long double> non_double(-std::numeric_limits<long double>::infinity(), nullptr, 65);
  require_no_index(non_double);
}

void allocation_contract() {
  constexpr cand_pos_t n = 65;
  std::size_t constructor_fallbacks = 0, copy_fallbacks = 0, assignment_fallbacks = 0;
  for (long fail = 0; fail < 10; ++fail) {
    std::optional<Chart> chart;
    allocation_fault::arm(fail);
    try { chart.emplace(zero, nullptr, n); } catch (const std::bad_alloc &) {}
    const bool fired = allocation_fault::fired;
    allocation_fault::disable();
    if (chart && fired) {
      ++constructor_fallbacks;
      require_no_index(*chart);
      require(chart->dense_length() == n && chart->size() == 0,
              "support allocation failure damaged dense storage");
      chart->set({NT::V, 1, n, -1, -1}, -1e300);
      require(same(chart->get({NT::V, 1, n, -1, -1}), -1e300),
              "support allocation fallback changed subsequent dense updates");
    }
  }
  Chart original(zero, nullptr, n);
  original.set({NT::V, 1, n, -1, -1}, -1e300);
  original.set({NT::VP, 64, n, -1, -1}, -0.0);
  for (long fail = 0; fail < 10; ++fail) {
    std::optional<Chart> copied;
    allocation_fault::arm(fail);
    try { copied.emplace(original); } catch (const std::bad_alloc &) {}
    const bool fired = allocation_fault::fired;
    allocation_fault::disable();
    if (copied && fired) {
      ++copy_fallbacks;
      require_no_index(*copied);
      require(copied->dense_length() == n, "optional support copy failure lost dense storage");
      compare_entries(original, *copied);
    }
    Chart assigned(zero, nullptr, 3);
    allocation_fault::arm(fail);
    bool completed = false;
    try { assigned = original; completed = true; } catch (const std::bad_alloc &) {}
    const bool assignment_fired = allocation_fault::fired;
    allocation_fault::disable();
    if (completed && assignment_fired) {
      ++assignment_fallbacks;
      require_no_index(assigned);
      require(assigned.dense_length() == n, "optional support assignment failure lost dense storage");
      compare_entries(original, assigned);
    }
  }
  require(constructor_fallbacks == 2 && copy_fallbacks == 2 && assignment_fallbacks == 2,
          "allocation faults did not isolate both optional bit vectors in construct/copy/assignment");
  require(original.live_primary_index_length() == n,
          "failed optional deepcopy mutated the source support index");
}

struct ChartView {
  const Chart &chart;
  double get(const Key &key) const { return chart.get(key); }
};
struct CustomLogSemiring : S {};

template <class Semiring, class Context>
void compare_selections(const Key &parent, const Chart &sparse, const Chart &dense,
                         const Context &ctx) {
  Scratch expected, actual;
  expected.prepare<Semiring>(parent, sparse, ctx);
  actual.prepare<Semiring>(parent, dense, ctx);
  require(expected.selected == actual.selected,
          "live support changed one of the 16 selected split arrays/order");
  for (const auto &selection : actual.selected)
    require(std::is_sorted(selection.begin(), selection.end()) &&
                std::adjacent_find(selection.begin(), selection.end()) == selection.end(),
            "live support selected a duplicate or out-of-order split");
  ++compared_selections;
}

void integration(cand_pos_t n, sc::WmbClosurePolicy policy) {
  std::string sequence(n, 'G'), scaffold(n, '.');
  for (cand_pos_t i = 0; i < n; ++i) sequence[i] = "GCAU"[i % 4];
  const auto add_pair = [&](cand_pos_t i, cand_pos_t j) {
    scaffold[i] = '('; scaffold[j] = ')';
    sequence[i] = 'G'; sequence[j] = 'C';
  };
  add_pair(0, n - 1); add_pair(1, n - 2);
  add_pair(7, 20); add_pair(8, 19);
  if (n > 50) { add_pair(32, 47); add_pair(33, 46); }
  sc::weights::RuntimeSnapshot owner(sequence, scaffold);
  sc::replay::WFinalLogWeightContext ctx(owner, policy);
  Chart sparse(zero), dense(zero, nullptr, n);
  std::mt19937 rng(20261020 + n);
  for (int f = 0; f < families; ++f)
    for (cand_pos_t j = 1; j <= n; ++j)
      for (cand_pos_t i = 1; i <= j; ++i) {
        if (rng() % 3) continue;
        const auto value = values[rng() % (values.size() - 1)];
        const Key key{static_cast<NT>(f), i, j, -1, -1};
        sparse.set(key, value); dense.set(key, value);
      }
  for (const auto nt : {NT::WM, NT::VM, NT::WI, NT::WIP, NT::VPR, NT::WMBP})
    for (cand_pos_t j = 1; j <= n; ++j)
      for (cand_pos_t i = 1; i <= j; ++i) {
        const Key parent{nt, i, j, -1, -1};
        require(Scratch::standard_live_index<S>(parent, dense, ctx) &&
                    Scratch::standard_live_index<M>(parent, dense, ctx),
                "standard canonical integration failed to opt into the support index");
        compare_selections<S>(parent, sparse, dense, ctx);
        compare_selections<M>(parent, sparse, dense, ctx);
      }
  // Static interface identity matters: a custom chart/context/semiring must
  // keep all old get/predicate and diagnostic behavior, even with the same data.
  const Key parent{NT::WM, 3, n - 3, -1, -1};
  const sc::deductions::PairabilityContext &custom_context = ctx;
  ChartView view{dense};
  require(!Scratch::standard_live_index<S>(parent, dense, custom_context) &&
              !Scratch::standard_live_index<CustomLogSemiring>(parent, dense, ctx) &&
              !Scratch::standard_live_index<sc::semirings::SumProductSemiring>(parent, dense, ctx) &&
              !Scratch::standard_live_index<S>(parent, view, ctx),
          "custom interface incorrectly enabled live support");
  compare_selections<S>(parent, sparse, dense, custom_context);
  compare_selections<CustomLogSemiring>(parent, sparse, dense, ctx);
  compare_selections<sc::semirings::SumProductSemiring>(parent, sparse, dense, ctx);
  Scratch custom_chart, ordinary;
  custom_chart.prepare<S>(parent, view, ctx);
  ordinary.prepare<S>(parent, sparse, ctx);
  require(custom_chart.selected == ordinary.selected, "custom chart full-range fallback changed selections");
  for (const Key malformed : std::vector<Key>{{NT::WM, 3, n - 3, 0, -1},
           {NT::WM, 3, n - 3, -1, 0}, {NT::VM, 1, n + 1, -1, -1}}) {
    require(!Scratch::standard_live_index<S>(malformed, dense, ctx),
            "noncanonical parent incorrectly enabled live support");
    compare_selections<S>(malformed, sparse, dense, ctx);
    compare_selections<M>(malformed, sparse, dense, ctx);
  }
  Chart custom_default(0.0, nullptr, n), sparse_default(0.0);
  custom_default.set({NT::V, 5, n - 3, -1, -1}, zero);
  sparse_default.set({NT::V, 5, n - 3, -1, -1}, zero);
  compare_selections<S>(parent, sparse_default, custom_default, ctx);
  require(!Scratch::standard_live_index<S>(parent, custom_default, ctx),
          "non-log default incorrectly enabled live support");

  // Out-of-range split diagnostics must survive even when every child is zero.
  const Key diagnostic_parent{NT::WIP, 1, 9, -1, -1};
  Chart empty(zero, nullptr, n);
  bool threw = false;
  try {
    Scratch::for_each_right_live_start<S, R::WIP_SPLIT_V, R::WIP_SPLIT_WMB>(
        diagnostic_parent, empty, ctx, 10, 10, [&](cand_pos_t k) {
          Scratch::children_nonzero<S, R::WIP_SPLIT_V>(diagnostic_parent, k, empty);
        });
  } catch (const std::logic_error &) { threw = true; }
  require(threw, "live support concealed an invalid dependency outside the canonical split range");
}
}  // namespace

int main() {
  for (cand_pos_t n : {1, 7, 8, 9, 63, 64, 65, 127, 128, 129}) storage_boundaries(n);
  copy_move_contract();
  default_contract();
  allocation_contract();
  for (const auto policy : {sc::WmbClosurePolicy::LegacyCompatible,
                            sc::WmbClosurePolicy::RequireWeakClosure}) {
    integration(24, policy);
    integration(70, policy);
  }
  std::cout << "scfg2_exact_live_support=ok values=" << compared_values
            << " rows=" << compared_rows << " selections=" << compared_selections << '\n';
}
