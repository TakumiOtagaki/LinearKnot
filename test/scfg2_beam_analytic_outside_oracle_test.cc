#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/engine/outside_marginals.hh"
#include "scfg2/profiles/beam_profiles.hh"

#include <cstdlib>
#include <functional>
#include <iostream>
#include <random>
#include <tuple>

namespace {
using namespace cparty::scfg2;
using Key = deductions::ItemKey;

bool close(double a, double b, double tolerance = 1e-9) {
  return a == b || std::fabs(a - b) <= tolerance * std::max({1.0, std::fabs(a), std::fabs(b)});
}

auto key(const deductions::Deduction &d) {
  return std::make_tuple(d.parent.nonterminal, d.parent.i, d.parent.j,
                        d.parent.ip, d.parent.jp, d.rule,
                        d.split.k, d.split.l, d.split.p, d.split.q);
}

template <class Value>
struct Capture {
  struct Candidate { deductions::Deduction deduction; Value value; };
  struct Item { Key key; Value value; };
  std::vector<Candidate> candidates;
  std::vector<Item> items;
  template <class Chart>
  void operator()(const Key &, const deductions::Deduction &d, Value v,
                  const Chart &, double) {
    candidates.push_back({d, semirings::checkpoint_value(v)});
  }
  template <class Chart>
  void after_item(const Key &k, const Chart &chart) {
    items.push_back({k, chart.get(k)});
  }
};

std::size_t checked = 0, captured_positive = 0, shared_positive = 0;

bool check_case(std::string sequence, std::string scaffold,
                const replay::WFinalVpWmbpOnlineWidths &widths) {
  using S = semirings::FactorizedLogSumExpSemiring;
  using T = semirings::TapeLogSumExpSemiring;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const Key root{grammar::NonTerminal::W, 1, n, -1, -1};
  Capture<S::value_type> a;
  Capture<T::value_type> t;
  profiles::detail::RunTermRecorder<S::value_type> ar;
  profiles::detail::RunTermRecorder<T::value_type> tr;
  semirings::FactorizedInsideGraph graph;
  semirings::LogSumExpTape tape;
  const auto analytic = [&] {
    S::Scope scope(graph);
    beam::ScopedRunTermObserver<S::value_type> observe(&ar);
    return replay::run_w_final_vp_wmbp_online_observed<S>(ctx, n, widths, a);
  }();
  const auto taped = [&] {
    T::Scope scope(tape);
    beam::ScopedRunTermObserver<T::value_type> observe(&tr);
    return replay::run_w_final_vp_wmbp_online_observed<T>(ctx, n, widths, t);
  }();
  const auto av = analytic.chart.get(root);
  const auto tv = taped.chart.get(root);
  const auto fail = [&](const std::string &what) {
    std::cerr << what << " n=" << n << " b=" << widths.v << " seq=" << sequence
              << " G=" << scaffold << '\n';
    return false;
  };
  if (!close(av.log_value, tv.log_value, 1e-11) ||
      analytic.item_visits != taped.item_visits ||
      analytic.positive_items != taped.positive_items ||
      analytic.frontier_selected != taped.frontier_selected ||
      analytic.frontier_candidates != taped.frontier_candidates ||
      analytic.vp_selected != taped.vp_selected ||
      analytic.wmbp_selected != taped.wmbp_selected ||
      analytic.case3_left_selected != taped.case3_left_selected ||
      analytic.be_right_shell_selected != taped.be_right_shell_selected ||
      analytic.wmbp_prefix_selected != taped.wmbp_prefix_selected) return fail("root/selection");
  if (!std::isfinite(av.log_value)) {
    ++checked;
    return true;  // Identical disconnected frozen beams, with forcing off.
  }
  const auto ap = graph.posteriors(av.node);
  const auto tp = tape.reverse(tv.node);
  const auto pa = [&](const S::value_type &v) { return v.node < ap.size() ? ap[v.node] : 0.0; };
  const auto pt = [&](const T::value_type &v) { return v.node < tp.size() ? tp[v.node] : 0.0; };
  if (a.items.size() != t.items.size() || a.candidates.size() != t.candidates.size() ||
      ar.terms.size() != tr.terms.size() || ar.vpr_cells.size() != tr.vpr_cells.size() ||
      ar.left_units.size() != tr.left_units.size()) return fail("observer topology");
  for (std::size_t i = 0; i < a.items.size(); ++i) {
    const auto &x = a.items[i]; const auto &y = t.items[i];
    if (!deductions::ItemKeyEq{}(x.key, y.key) || !close(x.value.log_value, y.value.log_value) ||
        !close(pa(x.value), pt(y.value))) {
      std::cerr << " item=" << int(x.key.nonterminal) << ',' << x.key.i << ',' << x.key.j
                << " p=" << pa(x.value) << '/' << pt(y.value) << '\n';
      return fail("pre-commit item posterior");
    }
    if (!close(analytic.chart.get(x.key).log_value, taped.chart.get(x.key).log_value))
      return fail("final chart value");
    if (!std::isfinite(taped.chart.get(x.key).log_value) && pt(y.value) > 1e-12) {
      if (captured_positive == 0) {
        std::cout << "capture_witness b=" << widths.v << " seq=" << sequence
                  << " G=" << scaffold << " nt=" << int(x.key.nonterminal)
                  << " i=" << x.key.i << " j=" << x.key.j
                  << " posterior=" << pt(y.value) << '\n';
      }
      ++captured_positive;
      return fail("pruned public item still contributes to root");
    }
  }
  for (std::size_t i = 0; i < a.candidates.size(); ++i) {
    const auto &x = a.candidates[i]; const auto &y = t.candidates[i];
    if (key(x.deduction) != key(y.deduction) || !close(x.value.log_value, y.value.log_value) ||
        !close(pa(x.value), pt(y.value))) {
      std::cerr << " rule=" << int(x.deduction.rule) << " p=" << pa(x.value) << '/' << pt(y.value) << '\n';
      return fail("candidate node posterior");
    }
  }
  for (std::size_t i = 0; i < ar.terms.size(); ++i) {
    const auto &x = ar.terms[i]; const auto &y = tr.terms[i];
    if (std::tie(x.rule, x.lo, x.hi) != std::tie(y.rule, y.lo, y.hi) ||
        !close(x.value.log_value, y.value.log_value) || !close(pa(x.value), pt(y.value)))
      return fail("run evaluation instance");
  }
  for (std::size_t i = 0; i < ar.vpr_cells.size(); ++i) {
    const auto &x = ar.vpr_cells[i]; const auto &y = tr.vpr_cells[i];
    if (std::tie(x.i, x.t, x.carried) != std::tie(y.i, y.t, y.carried) ||
        !close(x.total.log_value, y.total.log_value) || !close(pa(x.total), pt(y.total)))
      return fail("VPR cell");
    if (pt(y.total) > 0) ++shared_positive;
  }
  for (std::size_t i = 0; i < ar.left_units.size(); ++i) {
    const auto &x = ar.left_units[i]; const auto &y = tr.left_units[i];
    if (std::tie(x.rule, x.j, x.base, x.has_alias, x.linked) !=
        std::tie(y.rule, y.j, y.base, y.has_alias, y.linked) ||
        !close(pa(x.unit), pt(y.unit)) || !close(pa(x.alias), pt(y.alias)))
      return fail("KeepAll shared left unit");
  }
  const auto am = engine::collect_outside_marginals(analytic.chart, ap, n, ctx);
  const auto tm = engine::collect_outside_marginals(taped.chart, tp, n, ctx);
  if (!am.ok || !tm.ok || am.pairs.size() != tm.pairs.size()) return fail("BPP setup");
  for (std::size_t i = 0; i < am.pairs.size(); ++i) {
    const auto &x = am.pairs[i]; const auto &y = tm.pairs[i];
    if (x.i != y.i || x.j != y.j || !close(x.round_probability, y.round_probability) ||
        !close(x.square_probability, y.square_probability) ||
        !close(x.union_probability, y.union_probability)) return fail("BPP");
  }
  for (std::size_t i = 1; i < am.unpaired_probability.size(); ++i)
    if (!close(am.unpaired_probability[i], tm.unpaired_probability[i])) return fail("unpaired");
  std::size_t derivatives = 0;
  for (const auto &[item, v] : analytic.chart.values()) {
    (void)item;
    if (pa(v) <= 1e-12) continue;
    const auto taped_v = taped.chart.get(item);
    for (const double eps : {1e-3, 1e-4, 1e-5}) {
      const double fd = (tape.evaluate_with_log_shift(tv.node, taped_v.node, eps) -
                         tape.evaluate_with_log_shift(tv.node, taped_v.node, -eps)) / (2 * eps);
      if (!close(pa(v), fd, 1e-6)) return fail("frozen tape derivative");
    }
    if (++derivatives == 8) break;
  }
  const std::vector<cand_pos_t> windows{1, 2, 4, 7};
  const auto aprofile = profiles::beam_profiles(sequence, scaffold, widths, windows);
  const auto tprofile = profiles::beam_profiles_taped(sequence, scaffold, widths, windows);
  if (aprofile.ok != tprofile.ok || !close(aprofile.log_partition_weight, tprofile.log_partition_weight))
    return fail("profiles setup");
  for (std::size_t i = 1; i < aprofile.context.size(); ++i)
    for (std::size_t c = 0; c < profiles::loop_context_count; ++c)
      if (!close(aprofile.context[i][c], tprofile.context[i][c])) return fail("contexts");
  for (std::size_t w = 0; w < aprofile.accessibility.size(); ++w)
    for (std::size_t i = 1; i < aprofile.accessibility[w].size(); ++i)
      if (!close(aprofile.accessibility[w][i], tprofile.accessibility[w][i])) return fail("accessibility");
  ++checked;
  return true;
}

std::pair<std::string, std::string> random_case(std::mt19937 &rng, int n) {
  std::string g(n, '.');
  std::function<void(int, int)> fill = [&](int l, int r) {
    if (r - l < 4) return;
    if (std::uniform_int_distribution<int>(0, 2)(rng) == 0) {
      const int k = std::uniform_int_distribution<int>(l + 4, r)(rng);
      g[l] = '('; g[k] = ')'; fill(l + 1, k - 1); fill(k + 1, r);
    } else fill(l + 1, r);
  };
  fill(0, n - 1);
  std::string seq(n, 'A');
  for (auto &c : seq) c = "ACGU"[std::uniform_int_distribution<int>(0, 3)(rng)];
  std::vector<int> stack;
  for (int i = 0; i < n; ++i) {
    if (g[i] == '(') stack.push_back(i);
    if (g[i] == ')') {
      const int j = stack.back(); stack.pop_back();
      seq[i] = "UGCA"[std::string("ACGU").find(seq[j])];
    }
  }
  return {seq, g};
}
}  // namespace

int main(int argc, char **argv) {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  using cparty::scfg2::replay::WFinalVpWmbpOnlineWidths;
  std::vector<std::pair<std::string, std::string>> cases{
      {"GAGAAACAAAAAAGAAACC", "..................."},
      {"GCGCGCGCGCGC", "..(....)...."},
      // Historical pre-prune WIP witness: WIP(26,39) used to be zero in the
      // final b=3 chart while its old node still contributed to the root.
      {"AGGUACCCUACAAUUACGGCAUCAGUCUAAAAACUCUGGGAUGGGCA",
       "(.(....).....)..(..(....)........(....)).(....)"},
      {"GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
       "....(...)((....)..)((.(.((...))))(.(....)((...))))"}};
  std::mt19937 rng(20261008);
  const int random_count = argc > 1 ? std::atoi(argv[1]) : 20;
  for (int i = 0; i < random_count; ++i) cases.push_back(random_case(rng, 20 + rng() % 41));
  for (const auto &[seq, g] : cases) {
    for (std::size_t width : {1, 2, 3, 5, 10, 20})
      if (!check_case(seq, g, WFinalVpWmbpOnlineWidths::uniform(width))) return EXIT_FAILURE;
    if (!check_case(seq, g, WFinalVpWmbpOnlineWidths{})) return EXIT_FAILURE;
    auto unforced = WFinalVpWmbpOnlineWidths::uniform(3);
    unforced.force_scaffold_derivation = false;
    if (!check_case(seq, g, unforced)) return EXIT_FAILURE;
  }
  if (captured_positive != 0 || shared_positive == 0) return EXIT_FAILURE;
  std::cout << "scfg2_beam_analytic_outside_oracle=ok cases=" << checked
            << " captured_positive=" << captured_positive
            << " shared_vpr_positive=" << shared_positive << '\n';
}
