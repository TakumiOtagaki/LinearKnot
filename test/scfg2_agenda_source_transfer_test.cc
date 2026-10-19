#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <algorithm>
#include <array>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace sc = cparty::scfg2;
using Item = sc::deductions::ItemKey;
using Deduction = sc::deductions::Deduction;
using NT = sc::grammar::NonTerminal;

void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
bool same(const Deduction &a, const Deduction &b) {
  const sc::deductions::ItemKeyEq equal;
  if (!equal(a.parent, b.parent) || a.rule != b.rule ||
      a.split.k != b.split.k || a.split.l != b.split.l ||
      a.split.p != b.split.p || a.split.q != b.split.q ||
      a.child_count != b.child_count) return false;
  for (std::size_t c = 0; c < a.children.size(); ++c)
    if (!equal(a.children[c].item, b.children[c].item)) return false;
  return true;
}
bool same(const std::vector<Deduction> &a, const std::vector<Deduction> &b) {
  return a.size() == b.size() &&
      std::equal(a.begin(), a.end(), b.begin(), [](const auto &x, const auto &y) {
        return same(x, y);
      });
}

struct Counts { std::size_t parents = 0, sources = 0, aggregates = 0; };

// The reference keeps its source views until the next column. The consumed
// agenda must preserve both the source order and all future column results.
template <class S, class Agenda, class Build, class Aggregate>
Counts exercise(NT family, int n, Build build, Aggregate aggregate) {
  Agenda reference, consumed;
  std::vector<std::pair<std::vector<Deduction>, std::vector<Deduction>>> held;
  Counts counts;
  for (int j = 1; j <= n; ++j) {
    build(reference, j);
    build(consumed, j);
    for (const auto &[owned, expected] : held)
      require(same(owned, expected), "transferred storage did not survive next-column mutation");
    held.clear();
    const auto starts = reference.parent_starts();
    require(starts == consumed.parent_starts(), "consumption changed future parent order");
    require(consumed.take_for({family, n + 1, j, -1, -1}).empty(),
            "missing source lookup was not empty");
    for (int i : starts) {
      const Item parent{family, i, j, -1, -1};
      const auto expected = reference.deductions_for(parent);
      const auto &view = consumed.deductions_for(parent);
      require(same(expected, view), "consumption changed future deductions");
      const auto *storage = view.data();
      auto owned = consumed.take_for(parent);
      require(same(owned, expected), "transfer changed deduction fields/order");
      if (!owned.empty()) require(owned.data() == storage, "source storage was copied");
      require(consumed.deductions_for(parent).empty() && consumed.take_for(parent).empty(),
              "completed source storage was consumed twice");
      require(consumed.parent_starts() == starts, "transfer invalidated the completed parent list");
      for (const auto &source : expected) {
        const auto before = aggregate(reference, source);
        const auto after = aggregate(consumed, source);
        require(before.has_value() == after.has_value(), "transfer retired an aggregate owner");
        if (before) {
          require(S::priority(*before) == S::priority(*after) &&
                  S::traceback_priority(*before) == S::traceback_priority(*after),
                  "transfer changed a retained aggregate value");
          ++counts.aggregates;
        }
      }
      ++counts.parents;
      counts.sources += owned.size();
      held.emplace_back(std::move(owned), expected);
    }
  }
  for (const auto &[owned, expected] : held)
    require(same(owned, expected), "transferred storage changed after later transfers");
  return counts;
}

template <class S>
std::array<Counts, 9> check(const std::string &input_sequence,
                          const std::string &input_scaffold, std::size_t width) {
  auto sequence = input_sequence, scaffold = input_scaffold;
  sc::weights::RuntimeSnapshot owner(sequence, scaffold);
  sc::replay::WFinalLogWeightContext ctx(owner);
  const int n = sequence.size();
  using Value = typename S::value_type;
  sc::engine::Chart<Value> chart(S::zero());
  std::vector<std::vector<int>> starts(n + 1);
  std::vector<std::vector<Item>> be(n + 1);
  // An immutable synthetic child view exercises every agenda independently
  // of numerical DP reachability. Geometry/weights use the real scaffold.
  for (int j = 1; j <= n; ++j) for (int i = 1; i <= j; ++i) {
    starts[j].push_back(i);
    for (auto family : sc::engine::basic_exact_span_family_order()) {
      if (family == NT::BE) continue;
      chart.set({family, i, j, -1, -1}, S::from_log_weight(0.01 * (i + j)));
    }
    if (ctx.pair_partner(i) != j) continue;
    for (int ip = i; ip < j; ++ip) {
      const int jp = ctx.pair_partner(ip);
      if (jp <= ip || jp > j) continue;
      const Item item{NT::BE, i, j, ip, jp};
      be[j].push_back(item);
      chart.set(item, S::one());
    }
  }
  const auto ends = [&](int j) -> const std::vector<int> & {
    static const std::vector<int> empty;
    return j >= 1 && j <= n ? starts[j] : empty;
  };
  const auto retained = [&](NT, int j) -> const auto & { return ends(j); };
  sc::beam::WmbBtailStore<S> btails;
  for (int t = 1; t <= n; ++t) btails.build(t, chart, ctx, ends(t - 1), ends, width);
  const auto none = [](const auto &, const auto &) -> std::optional<Value> { return {}; };
  std::array<Counts, 9> out;
  out[0] = exercise<S, sc::beam::VmRetainedJoinAgenda<S>>(NT::VM, n,
      [&](auto &a, int j) { a.build_column(j, chart, ctx, retained, width); },
      [](const auto &a, const auto &d) -> std::optional<Value> {
        return a.reconstructs_base(d) ? std::optional<Value>{a.base_contribution(d.parent)} : std::nullopt;
      });
  out[1] = exercise<S, sc::beam::WmRetainedJoinAgenda<S>>(NT::WM, n,
      [&](auto &a, int j) { a.build_column(j, chart, ctx, retained, width); },
      [](const auto &a, const auto &d) -> std::optional<Value> {
        return a.reconstructs_start(d) ? std::optional<Value>{a.start_contribution(d.parent)} : std::nullopt;
      });
  out[2] = exercise<S, sc::beam::WiRetainedJoinAgenda>(NT::WI, n,
      [&](auto &a, int j) { a.build_column(j, ctx, retained); }, none);
  out[3] = exercise<S, sc::beam::WipRetainedJoinAgenda<S>>(NT::WIP, n,
      [&](auto &a, int j) { a.build_column(j, chart, ctx, retained, width); },
      [](const auto &a, const auto &d) -> std::optional<Value> {
        return a.reconstructs_start(d) ? std::optional<Value>{a.start_contribution(d.parent)} : std::nullopt;
      });
  out[4] = exercise<S, sc::beam::WmbwRetainedJoinAgenda>(NT::WMBW, n,
      [&](auto &a, int j) { a.build_column(j, ctx, ends, ends); }, none);
  out[5] = exercise<S, sc::beam::VplRetainedAgenda<S>>(NT::VPL, n,
      [&](auto &a, int j) { a.build_column(j, chart, ctx, ends(j), width); },
      [](const auto &a, const auto &d) -> std::optional<Value> {
        return a.reconstructs(d) ? std::optional<Value>{a.contribution(d.parent)} : std::nullopt;
      });
  out[6] = exercise<S, sc::beam::VprRetainedAgenda<S>>(NT::VPR, n,
      [&](auto &a, int j) { a.build_column(j, chart, ctx, ends(j), ends, width); },
      [](const auto &a, const auto &d) -> std::optional<Value> {
        return a.reconstructs_base(d) ? std::optional<Value>{a.base_contribution(d.parent)} : std::nullopt;
      });
  out[7] = exercise<S, sc::beam::WRetainedJoinAgenda>(NT::W, n,
      [&](auto &a, int j) { a.build_column(j, ctx, retained); }, none);
  out[8] = exercise<S, sc::beam::WmbRetainedAgenda<S>>(NT::WMB, n,
      [&](auto &a, int j) { a.build_column(j, chart, ctx, ends(j), be[j], btails); },
      [](const auto &a, const auto &d) -> std::optional<Value> {
        return a.reconstructs_split(d) ? std::optional<Value>{a.split_contribution(d.parent)} : std::nullopt;
      });
  return out;
}

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", SCFG2_TEST_PARAM_FILE, 1);
    std::array<Counts, 9> total{};
    for (std::size_t width : {std::size_t{3}, std::numeric_limits<std::size_t>::max()}) {
      for (auto [seq, g] : std::vector<std::pair<std::string, std::string>>{
          {"GGGGGAACCAACCC", "((.....))....."},
          {"GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
           "....(...)((....)..)((.(.((...))))(.(....)((...))))"}}) {
        for (bool maximum : {false, true}) {
          const auto got = maximum ? check<sc::semirings::LogMaxSemiring>(seq, g, width)
                                   : check<sc::semirings::LogSumExpSemiring>(seq, g, width);
          for (std::size_t k = 0; k < total.size(); ++k) {
            total[k].parents += got[k].parents;
            total[k].sources += got[k].sources;
            total[k].aggregates += got[k].aggregates;
          }
        }
      }
    }
    for (std::size_t k = 0; k < total.size(); ++k) {
      require(total[k].parents && total[k].sources, "agenda transfer coverage was vacuous");
      if (k == 0 || k == 1 || k == 3 || k == 5 || k == 6 || k == 8)
        require(total[k].aggregates, "aggregate lifetime coverage was vacuous");
    }
    std::cout << "scfg2_agenda_source_transfer=ok families=9 semirings=2 widths=2\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
