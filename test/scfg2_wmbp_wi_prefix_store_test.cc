#include "scfg2/beam/wmbp_wi_prefix_store.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <type_traits>

namespace {
using namespace cparty::scfg2;
using deductions::ItemKey;

void require(bool ok, const char *message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}

template <class Semiring>
struct Chart {
  using Value = typename Semiring::value_type;
  std::unordered_map<ItemKey, Value, deductions::ItemKeyHash, deductions::ItemKeyEq> values;
  Value get(const ItemKey &item) const {
    const auto found = values.find(item);
    return found == values.end() ? Semiring::zero() : found->second;
  }
};

// Deliberately expensive, independent full-sort oracle for the original
// primary ranks, with FIFO resolving the formerly unspecified complete ties.
template <class Semiring>
struct SortReference {
  using Value = typename Semiring::value_type;
  using Prefix = beam::WmbpWiPrefixCell<Value>;
  struct Wi { ItemKey item; Value value; };
  std::map<cand_pos_t, std::vector<Wi>> histories;
  std::map<cand_pos_t, std::vector<Prefix>> prefixes;
  beam::WmbpWiPrefixStoreStats stats;

  void wi(cand_pos_t j, const std::vector<cand_pos_t> &starts,
          const Chart<Semiring> &chart, std::size_t width) {
    for (cand_pos_t i : starts) {
      const ItemKey item{grammar::NonTerminal::WI, i, j, -1, -1};
      const Value value = chart.get(item);
      if (Semiring::is_zero(value)) continue;
      auto &v = histories[i];
      v.push_back({item, value});
      ++stats.wi_history_insertions;
      std::stable_sort(v.begin(), v.end(), [](const Wi &a, const Wi &b) {
        const double x = Semiring::priority(a.value), y = Semiring::priority(b.value);
        return x != y ? x > y : a.item.j < b.item.j;
      });
      if (v.size() > width) v.resize(width);
    }
  }
  void add(Prefix cell, std::size_t width) {
    auto &v = prefixes[cell.vp_i];
    v.push_back(std::move(cell));
    ++stats.prefix_candidates;
    std::stable_sort(v.begin(), v.end(), [](const Prefix &a, const Prefix &b) {
      const double x = Semiring::priority(a.factor), y = Semiring::priority(b.factor);
      if (x != y) return x > y;
      return a.be.i != b.be.i ? a.be.i < b.be.i : a.be.j < b.be.j;
    });
    if (v.size() > width) v.resize(width);
  }
  void be(const std::vector<ItemKey> &items, const Chart<Semiring> &chart, std::size_t width) {
    for (const auto &item : items) {
      const Value value = chart.get(item);
      if (Semiring::is_zero(value)) continue;
      ++stats.empty_probes;
      add(Prefix{item, {}, false, item.ip + 1, value}, width);
      const auto found = histories.find(item.ip + 1);
      if (found == histories.end()) continue;
      for (const auto &entry : found->second) {
        ++stats.be_wi_probes;
        if (entry.item.j + 1 >= item.j) continue;
        add(Prefix{item, entry.item, true, entry.item.j + 1,
                   Semiring::times(value, entry.value)}, width);
      }
    }
  }
};

template <class Semiring>
void compare(const beam::WmbpWiPrefixStore<Semiring> &store,
             const SortReference<Semiring> &reference, cand_pos_t upto) {
  for (cand_pos_t vp = 1; vp <= upto; ++vp) {
    const auto &got = store.retained(vp);
    const auto found = reference.prefixes.find(vp);
    const std::size_t wanted = found == reference.prefixes.end() ? 0 : found->second.size();
    require(got.size() == wanted, "prefix count differs from full-sort oracle");
    for (std::size_t x = 0; x < wanted; ++x) {
      const auto &a = got[x], &b = found->second[x];
      require(deductions::ItemKeyEq{}(a.be, b.be) &&
              deductions::ItemKeyEq{}(a.wi, b.wi) && a.has_wi == b.has_wi && a.vp_i == b.vp_i &&
              Semiring::priority(a.factor) == Semiring::priority(b.factor) &&
              Semiring::traceback_priority(a.factor) == Semiring::traceback_priority(b.factor),
              "prefix identity/value/order differs from full-sort oracle");
    }
  }
  const auto &a = store.stats(), &b = reference.stats;
  require(a.wi_history_insertions == b.wi_history_insertions &&
          a.be_wi_probes == b.be_wi_probes && a.empty_probes == b.empty_probes &&
          a.prefix_candidates == b.prefix_candidates, "candidate/timing accounting changed");
}

std::size_t differential_columns = 0;
template <class Semiring>
void differential() {
  const std::vector<std::size_t> widths{0, 1, 3, 17, std::numeric_limits<std::size_t>::max()};
  std::mt19937 rng(20261009);
  for (std::size_t history_width : widths) for (std::size_t prefix_width : widths) {
    beam::WmbpWiPrefixStore<Semiring> store;
    SortReference<Semiring> reference;
    Chart<Semiring> chart;
    for (cand_pos_t j = 12; j < 44; ++j) {
      for (cand_pos_t i : {3, 5, 8}) {
        auto value = Semiring::from_weight(static_cast<double>(rng() % 5));
        if constexpr (std::is_same_v<Semiring, semirings::CoupledLogSumMaxSemiring>)
          value.log_max -= i % 3;
        chart.values[ItemKey{grammar::NonTerminal::WI, i, j, -1, -1}] = value;
      }
      store.commit_wi_column(j, {3, 5, 8}, chart, history_width);
      reference.wi(j, {3, 5, 8}, chart, history_width);
      const std::vector<ItemKey> be_items{
          {grammar::NonTerminal::BE, 1, j, 4, j - 2},
          {grammar::NonTerminal::BE, 1, j, 7, j - 2},
          {grammar::NonTerminal::BE, 2, j, 2, j - 2},
          {grammar::NonTerminal::BE, 3, j, 3, j - 2}};
      const double shared = static_cast<double>(rng() % 5);
      chart.values[be_items[0]] = Semiring::from_weight(shared);
      chart.values[be_items[1]] = Semiring::from_weight(shared);  // complete primary ties
      chart.values[be_items[2]] = Semiring::from_weight(static_cast<double>(rng() % 5));
      // be_items[3] is absent/zero; no WI history exists for its left endpoint.
      store.commit_be_column(be_items, chart, prefix_width);
      reference.be(be_items, chart, prefix_width);
      compare(store, reference, j + 1);
      const auto before = store.stats();
      compare(store, reference, j + 1);
      const auto &after = store.stats();
      require(after.bucket_view_builds == before.bucket_view_builds &&
              after.bucket_view_cells == before.bucket_view_cells &&
              after.bucket_view_comparisons == before.bucket_view_comparisons,
              "unchanged prefix reread rebuilt ranked view");
      ++differential_columns;
    }
  }
}

void scaling(std::size_t width, std::size_t columns) {
  using Semiring = semirings::SumProductSemiring;
  beam::WmbpWiPrefixStore<Semiring> store;
  Chart<Semiring> chart;
  const cand_pos_t left = static_cast<cand_pos_t>(width + 2);
  for (std::size_t x = 0; x < width; ++x) {
    const cand_pos_t j = static_cast<cand_pos_t>(2 * width + 2 + x);
    chart.values[ItemKey{grammar::NonTerminal::WI, left, j, -1, -1}] = 1.0;
    store.commit_wi_column(j, {left}, chart, width);
  }
  for (std::size_t col = 0; col < columns; ++col) {
    const cand_pos_t q = static_cast<cand_pos_t>(4 * width + 3 + col);
    std::vector<ItemKey> items;
    for (std::size_t x = 1; x <= width; ++x) {
      items.push_back(ItemKey{grammar::NonTerminal::BE, static_cast<cand_pos_t>(x), q,
                             static_cast<cand_pos_t>(width + 1), q - 1});
      chart.values[items.back()] = static_cast<double>((col + 1) * (width + 1) + x);
    }
    store.commit_be_column(items, chart, width);
    for (std::size_t x = 0; x <= width; ++x) {
      const cand_pos_t vp = x == width ? left : static_cast<cand_pos_t>(2 * width + 3 + x);
      const auto &v = store.retained(vp);
      require(v.size() == width, "scaling prefix width violated");
      for (std::size_t r = 0; r < width; ++r)
        require(v[r].be.j == q && v[r].be.i == static_cast<cand_pos_t>(width - r) &&
                v[r].has_wi == (x != width), "scaling selected incorrect latest prefix");
    }
    const auto before = store.stats();
    for (std::size_t x = 0; x <= width; ++x)
      (void)store.retained(x == width ? left : static_cast<cand_pos_t>(2 * width + 3 + x));
    require(store.stats().bucket_view_builds == before.bucket_view_builds,
            "scaling reread materialized output twice");
  }
  std::size_t lg = 0;
  for (std::size_t x = width + 2; x > 1; x >>= 1) ++lg;
  const auto &s = store.stats();
  require(s.be_wi_probes == columns * width * width &&
          s.prefix_candidates == columns * (width * width + width), "b^2 candidate budget incorrect");
  require(s.bucket_update_comparisons <=
          (s.wi_history_insertions + s.prefix_candidates) * (6 * lg + 6),
          "prefix updates exceeded logarithmic budget");
  require(s.bucket_view_builds == 1 + columns * (width + 1) &&
          s.bucket_view_cells == width * s.bucket_view_builds,
          "ranked-view count exceeded production output phases");
  require(s.bucket_view_comparisons <= s.bucket_view_cells * (8 * lg + 8),
          "prefix views exceeded b^2 log b budget");
  std::cout << "b=" << width << " columns=" << columns << " candidates=" << s.prefix_candidates
            << " update_comparisons=" << s.bucket_update_comparisons
            << " view_comparisons=" << s.bucket_view_comparisons
            << " view_builds=" << s.bucket_view_builds << '\n';
}
}  // namespace

int main() {
  differential<semirings::SumProductSemiring>();
  differential<semirings::MaxProductSemiring>();
  differential<semirings::LogSumExpSemiring>();
  differential<semirings::LogMaxSemiring>();
  differential<semirings::CoupledLogSumMaxSemiring>();
  {
    semirings::FactorizedInsideGraph graph;
    semirings::FactorizedLogSumExpSemiring::Scope scope(graph);
    differential<semirings::FactorizedLogSumExpSemiring>();
  }
  for (std::size_t width : {1u, 4u, 16u, 64u, 256u}) scaling(width, 16);
  for (std::size_t columns : {8u, 32u, 128u}) scaling(16, columns);
  std::cout << "differential_columns=" << differential_columns << "\nscfg2_wmbp_wi_prefix_store=ok\n";
}
