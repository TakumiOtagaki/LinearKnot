#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/retained_rule_admission.hh"
#include "scfg2/beam/vpr_regenerated_base_cells.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <set>

namespace {
using namespace cparty::scfg2;
using S = semirings::LogSumExpSemiring;
using Key = deductions::ItemKey;
using Cells = beam::VprRegeneratedBaseCells<S>;
void require(bool ok, const char *message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
bool close(double a, double b) {
  return a == b || std::fabs(a - b) <= 1e-10 * std::max({1.0, std::fabs(a), std::fabs(b)});
}
struct Capture : beam::RunTermObserver<double> {
  std::vector<Cells::Column> cells;
  beam::RetainedRuleAdmission<S> admission;
  std::unordered_map<Key, std::set<grammar::RuleId>, deductions::ItemKeyHash,
                     deductions::ItemKeyEq> rules;
  explicit Capture(std::size_t n) : cells(n + 1) {}
  void on_run_term(grammar::RuleId, cand_pos_t, cand_pos_t, const double &) override {}
  void on_left_run_unit(grammar::RuleId, cand_pos_t, cand_pos_t,
                        const double &, const double *, bool) override {}
  void on_vpr_run_cell(cand_pos_t i, cand_pos_t j, const double &v, bool carried) override {
    cells[j].emplace(i, Cells::Cell{v, carried, false});
  }
  template <class Chart>
  void operator()(const Key &parent, const deductions::Deduction &d,
                   double candidate, const Chart &chart, double local) {
    admission(parent, d, candidate, chart, local);
    if (!S::is_zero(candidate)) rules[parent].insert(d.rule);
  }
};
std::size_t checked = 0, cells_checked = 0, private_without_public = 0;
void check(std::string seq, std::string g, std::size_t width) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(seq.size());
  auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(width);
  Capture observed(n);
  const auto run = [&] {
    beam::ScopedRunTermObserver<double> scope(&observed);
    return replay::run_w_final_vp_wmbp_online_observed<S>(ctx, n, widths, observed);
  }();
  observed.admission.finish(run.chart);
  require(observed.admission.pending_size() == 0 &&
          observed.admission.size() == run.chart.values().size(), "rule admission size/flush");
  for (const auto &[parent, rules] : observed.rules)
    for (const auto rule : grammar::all_rule_ids())
      require(observed.admission.has(parent, rule) ==
              (!S::is_zero(run.chart.get(parent)) && rules.count(rule)),
              "rule admission retained a pruned state or lost a contributing rule");
  std::array<std::vector<std::vector<cand_pos_t>>, replay::w_final_online_nonterminal_count> starts;
  for (auto &family : starts) family.resize(n + 1);
  for (const auto &[item, value] : run.chart.values()) {
    (void)value;
    starts[static_cast<std::size_t>(item.nonterminal)][item.j].push_back(item.i);
  }
  for (auto &family : starts) for (auto &column : family) std::sort(column.begin(), column.end());
  const auto kept = [&](grammar::NonTerminal nt, cand_pos_t j) -> const auto & {
    return starts[static_cast<std::size_t>(nt)][j];
  };
  const auto rebuilt = Cells::from_retained_chart(n, run.chart, ctx, kept, width);
  using T = semirings::TapeLogSumExpSemiring;
  semirings::LogSumExpTape tape;
  T::Scope scope(tape);
  engine::Chart<T::value_type> seeds(T::zero());
  std::vector<std::unordered_map<cand_pos_t, T::value_type>> private_nodes(n + 1);
  auto objective = T::zero();
  const auto cp = T::from_weight(ctx.vpr_split_vp_basepair_weight(1));
  for (cand_pos_t j = 1; j <= n; ++j) {
    require(rebuilt.cells(j).size() == observed.cells[j].size(), "private cell count changed");
    for (auto &[i, recorded] : observed.cells[j]) {
      const auto &cell = rebuilt.cells(j).at(i);
      const Key parent{grammar::NonTerminal::VPR, i, j, -1, -1};
      recorded.public_base = observed.admission.has(parent, grammar::RuleId::VPR_SPLIT_VP_BASEPAIR);
      require(close(cell.alpha, recorded.alpha) && cell.carried == recorded.carried &&
              cell.public_base == recorded.public_base, "private alpha/topology/admission changed");
      auto total = T::zero();
      if (recorded.carried) total = T::times(private_nodes[j - 1].at(i), cp);
      const Key vp{grammar::NonTerminal::VP, i, j - 1, -1, -1};
      if (i < j - 1 && !S::is_zero(run.chart.get(vp))) {
        if (!seeds.contains(vp)) seeds.set(vp, tape.constant(run.chart.get(vp)));
        total = T::plus(total, T::times(seeds.get(vp), cp));
      }
      require(close(total.log_value, recorded.alpha), "observed private recurrence does not replay");
      private_nodes[j][i] = total;
      if (recorded.public_base)
        objective = T::plus(objective, T::times(total, T::from_weight(1 + (i + j) % 5)));
      ++cells_checked;
    }
  }
  auto beta = rebuilt.empty_outside();
  engine::Chart<double> public_beta(S::zero());
  for (cand_pos_t j = n; j >= 1; --j) {
    for (const auto &[i, cell] : observed.cells[j]) {
      (void)cell;
      rebuilt.add_public_outside({grammar::NonTerminal::VPR, i, j, -1, -1},
                                  S::from_weight(1 + (i + j) % 5), &beta);
    }
    rebuilt.propagate_column(j, run.chart, ctx, &beta, &public_beta);
  }
  if (!T::is_zero(objective)) {
    const auto posterior = tape.reverse(objective.node);
    for (const auto &[vp, value] : seeds.values())
      require(close(std::exp(public_beta.get(vp) + run.chart.get(vp) - objective.log_value),
                    posterior[value.node]), "private outside VP seed derivative differs from tape");
    for (cand_pos_t j = 1; j <= n; ++j) for (const auto &[i, cell] : observed.cells[j]) {
      const auto found = beta[j].find(i);
      const double b = found == beta[j].end() ? S::zero() : found->second;
      const double p = std::exp(b + cell.alpha - objective.log_value);
      require(close(p, posterior[private_nodes[j].at(i).node]), "private carry derivative differs from tape");
      if (!cell.public_base && p > 0) ++private_without_public;
    }
  }
  if (width != std::numeric_limits<std::size_t>::max())
    require(rebuilt.size() <= std::size_t(n) * width, "private history exceeds nb cells");
  ++checked;
}
}  // namespace

int main() {
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
    for (const std::size_t width : {0, 1, 3, 20, 100}) check(seq, g, width);
    check(seq, g, std::numeric_limits<std::size_t>::max());
  }
  require(private_without_public > 0, "missing carry through a private cell without a public base");
  std::cout << "scfg2_vpr_regenerated_base_cells=ok cases=" << checked << " cells="
            << cells_checked << " private_without_public=" << private_without_public << '\n';
}
