#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/wmb_regenerated_btail_cells.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

namespace {
using namespace cparty::scfg2;
using S = semirings::LogSumExpSemiring;
using T = semirings::TapeLogSumExpSemiring;
using Key = deductions::ItemKey;
void require(bool ok, const char *message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
bool close(double a, double b) {
  return a == b || std::fabs(a - b) <= 1e-10 * std::max({1.0, std::fabs(a), std::fabs(b)});
}
std::size_t cases = 0, cells = 0, splits = 0, empty_seeds = 0, wi_seeds = 0;

void representative_tie_witness() {
  struct Context {
    cand_pos_t pair_partner(cand_pos_t) const { return 1; }
    cand_pos_t border_Bp(cand_pos_t, cand_pos_t t) const { return t; }
  } ctx;
  struct Chart {
    double score7, score4;
    double get(const Key &key) const {
      if (key.nonterminal == grammar::NonTerminal::WI) return 0.0;
      if (key.nonterminal != grammar::NonTerminal::WMBP || key.i != 1)
        return S::zero();
      if (key.j == 9) return 0.0;
      if (key.j == 7) return score7;
      if (key.j == 4) return score4;
      return S::zero();
    }
  };
  const std::vector<cand_pos_t> one{1}, empty;
  const auto starts = [&](cand_pos_t k) -> const auto & {
    return k == 4 || k == 7 || k == 9 ? one : empty;
  };
  struct Fixture { double score7, score4; cand_pos_t expected_k; };
  for (const auto fixture : {Fixture{0.0, 0.0, 4},
                             Fixture{0.5e-12, 0.25e-12, 4},
                             Fixture{2.0e-12, 0.0, 7}}) {
    const Chart chart{fixture.score7, fixture.score4};
    const double total = std::log(1.0 + std::exp(fixture.score7) +
                                 std::exp(fixture.score4));
    for (const auto &wi : {std::vector<cand_pos_t>{8, 5},
                           std::vector<cand_pos_t>{5, 8}}) {
      for (const auto width : {std::size_t{0}, std::size_t{1},
                               std::numeric_limits<std::size_t>::max()}) {
        beam::WmbBtailStore<S> store;
        store.build(10, chart, ctx, wi, starts, width);
        const auto &kept = store.retained(10);
        require(kept.size() == (width == 0 ? 0u : 1u),
                "BTAIL tie witness changed selected cardinality");
        if (kept.empty()) continue;
        require(kept[0].i == 1 && kept[0].best_k == fixture.expected_k,
                "BTAIL exact/near-tie representative changed");
        require(std::fabs(kept[0].total - total) <= 1e-14,
                "BTAIL tie witness changed the sum of its three terms");
      }
    }
  }
}

// Storage witness with b completed WI children and b distinct WMBP starts
// per child: b^2 candidate starts, but only b persistent cells per close.
// The small context supplies the BTAIL seed contract, not a folding oracle.
void retained_capacity_witness() {
  constexpr cand_pos_t width = 20, candidate_count = width * width;
  struct Context {
    cand_pos_t pair_partner(cand_pos_t) const { return 1; }
    cand_pos_t border_Bp(cand_pos_t, cand_pos_t t) const { return t; }
  } ctx;
  struct Chart {
    double get(const Key &key) const {
      return key.nonterminal == grammar::NonTerminal::WMBP ? -double(key.i) : 0.;
    }
  } chart;
  std::vector<std::vector<cand_pos_t>> starts(width);
  std::vector<cand_pos_t> wi;
  for (cand_pos_t k = 0; k < width; ++k) {
    wi.push_back(candidate_count + k + 1);
    for (cand_pos_t a = 1; a <= width; ++a) starts[k].push_back(k * width + a);
  }
  const std::vector<cand_pos_t> empty;
  const auto wmbp = [&](cand_pos_t k) -> const auto & {
    return k >= candidate_count && k < candidate_count + width
        ? starts[k - candidate_count] : empty;
  };
  for (const std::size_t retained_width : {0u, 1u, 20u}) {
    beam::WmbBtailStore<S> store;
    for (cand_pos_t t = candidate_count + 2 * width; t < candidate_count + 4 * width; ++t) {
      store.build(t, chart, ctx, wi, wmbp, retained_width);
      const auto &kept = store.retained(t);
      require(kept.size() == retained_width, "BTAIL capacity witness did not prune candidates");
      require(kept.capacity() <= retained_width,
              "BTAIL history retained its quadratic candidate capacity");
      for (std::size_t x = 0; x < kept.size(); ++x)
        require(kept[x].i == cand_pos_t(x + 1) && kept[x].total == -double(x + 1) &&
                    kept[x].best_k == candidate_count + cand_pos_t(x) / width,
                "BTAIL compact storage changed rank, value or representative");
    }
  }
  beam::WmbBtailStore<S> unlimited;
  unlimited.build(candidate_count + 2 * width, chart, ctx, wi, wmbp,
                  std::numeric_limits<std::size_t>::max());
  require(unlimited.retained(candidate_count + 2 * width).size() == candidate_count,
          "BTAIL compact storage pruned KeepAll");
}
void check(std::string seq, std::string g, std::size_t width, std::size_t tail_width) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(seq.size());
  auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(width);
  widths.wmb_btail = tail_width;
  const auto run = replay::run_w_final_vp_wmbp_online<S>(ctx, n, widths);
  std::array<std::vector<std::vector<cand_pos_t>>, replay::w_final_online_nonterminal_count> starts;
  std::vector<std::vector<Key>> be(n + 1);
  for (auto &family : starts) family.resize(n + 1);
  for (const auto &[item, value] : run.chart.values()) {
    (void)value;
    starts[static_cast<std::size_t>(item.nonterminal)][item.j].push_back(item.i);
    if (item.nonterminal == grammar::NonTerminal::BE) be[item.j].push_back(item);
  }
  for (auto &family : starts) for (auto &column : family) std::sort(column.begin(), column.end());
  for (auto &column : be) std::sort(column.begin(), column.end(), [&](const auto &a, const auto &b) {
    const auto x = run.chart.get(a), y = run.chart.get(b);
    return x != y ? x > y : beam::item_tie_key(a) < beam::item_tie_key(b);
  });
  const auto kept = [&](grammar::NonTerminal nt, cand_pos_t j) -> const auto & {
    return starts[static_cast<std::size_t>(nt)][j];
  };
  const auto regenerated = beam::WmbRegeneratedBtailCells<S>::from_retained_chart(
      n, run.chart, ctx, kept, tail_width);
  if (tail_width != std::numeric_limits<std::size_t>::max())
    require(regenerated.size() <= std::size_t(n) * tail_width, "BTAIL history exceeds nb cells");

  // Original forward private aggregates, with independent tape leaves for
  // retained public children. No reverse code or regenerated term list is
  // used to form the reference objective.
  semirings::LogSumExpTape tape;
  T::Scope scope(tape);
  engine::Chart<T::value_type> leaves(T::zero());
  for (const auto &[item, value] : run.chart.values()) leaves.set(item, tape.constant(value));
  beam::WmbBtailStore<T> reference;
  for (cand_pos_t t = 2; t <= n; ++t) {
    reference.build(t, leaves, ctx, kept(grammar::NonTerminal::WI, t - 1),
        [&](cand_pos_t k) -> const auto & { return kept(grammar::NonTerminal::WMBP, k); }, tail_width);
    const auto &actual = regenerated.store().retained(t);
    const auto &wanted = reference.retained(t);
    require(actual.size() == wanted.size(), "BTAIL retained count changed");
    for (std::size_t x = 0; x < actual.size(); ++x)
      require(actual[x].i == wanted[x].i && actual[x].best_k == wanted[x].best_k &&
              close(actual[x].total, wanted[x].total.log_value), "BTAIL alpha/rank/representative changed");
    cells += actual.size();
  }
  auto objective = T::zero();
  beam::WmbRetainedAgenda<T> agenda;
  beam::WmbRetainedAgenda<S> scalar;
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);
  engine::Chart<double> public_beta(S::zero());
  for (cand_pos_t j = 1; j <= n; ++j) {
    agenda.build_column(j, leaves, ctx, kept(grammar::NonTerminal::WMBP, j), be[j], reference);
    scalar.build_column(j, run.chart, ctx, kept(grammar::NonTerminal::WMBP, j), be[j], regenerated.store());
    for (const auto i : kept(grammar::NonTerminal::WMB, j)) {
      const Key parent{grammar::NonTerminal::WMB, i, j, -1, -1};
      auto alpha = S::zero();
      for (const auto &source : scalar.deductions_for(parent)) {
        if (!engine::legacy_span_schedule_supports(source) || weight.weight_for(source) == 0.0) continue;
        alpha = S::plus(alpha, scalar.reconstructs_split(source)
            ? scalar.split_contribution(parent)
            : replay::runtime_deduction_contribution<S>(run.chart, source, weight.weight_for(source)));
      }
      require(close(alpha, run.chart.get(parent)), "WMB forward aggregate does not rebuild from final chart");
      const auto term = agenda.split_contribution(parent);
      if (T::is_zero(term)) continue;
      const double w = 1 + (i + j) % 5;
      objective = T::plus(objective, T::times(term, T::from_weight(w)));
      public_beta.set(parent, S::from_weight(w));
      ++splits;
    }
  }
  auto private_beta = regenerated.empty_outside();
  for (cand_pos_t j = n; j >= 1; --j) {
    regenerated.propagate_wmb_splits(j, run.chart, ctx, be[j], &public_beta, &private_beta);
    regenerated.propagate_column(j, run.chart, ctx, kept, &private_beta, &public_beta);
  }
  if (!T::is_zero(objective)) {
    const auto posterior = tape.reverse(objective.node);
    for (const auto &[item, value] : leaves.values()) {
      if (item.nonterminal != grammar::NonTerminal::BE &&
          item.nonterminal != grammar::NonTerminal::WI &&
          item.nonterminal != grammar::NonTerminal::WMBP) continue;
      require(close(std::exp(public_beta.get(item) + run.chart.get(item) - objective.log_value),
                    posterior[value.node]), "WMB/BTAIL public derivative differs from tape");
    }
    for (cand_pos_t t = 2; t <= n; ++t) {
      for (const auto &cell : reference.retained(t)) {
        const auto found = private_beta[t].find(cell.i);
        const double b = found == private_beta[t].end() ? S::zero() : found->second;
        require(close(std::exp(b + cell.total.log_value - objective.log_value), posterior[cell.total.node]),
                "BTAIL private derivative differs from tape");
      }
      beam::for_each_wmb_btail_seed(t, ctx, kept(grammar::NonTerminal::WI, t - 1),
          [&](cand_pos_t k) -> const auto & { return kept(grammar::NonTerminal::WMBP, k); },
          [&](cand_pos_t i, cand_pos_t, bool has_wi) {
            if (private_beta[t].count(i)) { if (has_wi) ++wi_seeds; else ++empty_seeds; }
          });
    }
  }
  ++cases;
}
}  // namespace

int main() {
  representative_tie_witness();
  retained_capacity_witness();
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  std::vector<std::pair<std::string, std::string>> inputs{
      {"GCGCGCGCGCGC", "..(....)...."},
      {"AGGUACCCUACAAUUACGGCAUCAGUCUAAAAACUCUGGGAUGGGCA",
       "(.(....).....)..(..(....)........(....)).(....)"},
      {"GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
       "....(...)((....)..)((.(.((...))))(.(....)((...))))"}};
  std::mt19937 rng(20261009);
  for (int x = 0; x < 20; ++x) {
    auto [seq, g] = inputs[1 + x % 2];
    for (auto &c : seq) c = "ACGU"[rng() % 4];
    std::vector<std::size_t> open;
    for (std::size_t i = 0; i < g.size(); ++i) {
      if (g[i] == '(') open.push_back(i);
      if (g[i] == ')') { const auto k = open.back(); open.pop_back(); seq[i] = "UGCA"[std::string("ACGU").find(seq[k])]; }
    }
    inputs.push_back({seq, g});
  }
  for (const auto &[seq, g] : inputs) {
    for (const std::size_t b : {1, 3, 20, 100}) check(seq, g, b, b);
    check(seq, g, std::numeric_limits<std::size_t>::max(), std::numeric_limits<std::size_t>::max());
    for (const std::size_t b : {0, 1, 3}) check(seq, g, 20, b);
  }
  require(splits > 0 && empty_seeds > 0 && wi_seeds > 0, "missing BTAIL empty/nonempty WI witnesses");
  std::cout << "scfg2_wmb_regenerated_btail_cells=ok cases=" << cases << " cells=" << cells
            << " splits=" << splits << " empty=" << empty_seeds << " wi=" << wi_seeds << '\n';
}
