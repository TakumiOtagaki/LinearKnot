#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/wmbp_regenerated_sources.hh"
#include "scfg2/beam/be_regenerated_sources.hh"
#include "scfg2/engine/outside.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <array>
#include <cmath>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <type_traits>

namespace {
using namespace cparty::scfg2;
using Key = deductions::ItemKey;
using D = deductions::Deduction;

void require(bool ok, const char *message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
bool close(double a, double b) {
  return a == b || std::fabs(a - b) <= 1e-10 * std::max({1.0, std::fabs(a), std::fabs(b)});
}
auto id(const D &d) {
  return std::make_tuple(d.rule, d.parent.i, d.parent.j,
                        d.split.k, d.split.l, d.split.p, d.split.q,
                        beam::item_tie_key(d.children[0].item),
                        beam::item_tie_key(d.children[1].item),
                        beam::item_tie_key(d.children[2].item), d.child_count);
}

template <class S>
struct Capture {
  beam::WmbpWiPrefixAdmission admission;
  std::unordered_map<Key, std::vector<D>, deductions::ItemKeyHash,
                     deductions::ItemKeyEq> sources;
  template <class Chart>
  void after_be_commit(const std::vector<Key> &items, const Chart &chart,
                        const beam::WmbpWiPrefixStore<S> &store) {
    admission.template observe_be_column<S>(items, chart, store);
  }
  void before_wmbp_column(cand_pos_t j, const std::vector<cand_pos_t> &starts,
                          const beam::WmbpWiPrefixStore<S> &store) {
    admission.template observe_vp_column<S>(j, starts, store);
  }
  template <class Chart>
  void operator()(const Key &parent, const D &d, typename S::value_type v,
                   const Chart &, double) {
    if ((parent.nonterminal == grammar::NonTerminal::WMBP ||
         parent.nonterminal == grammar::NonTerminal::BE) && !S::is_zero(v))
      sources[parent].push_back(d);
  }
};

std::size_t cases = 0, parents = 0;
std::array<std::size_t, 5> modes{};
std::array<std::size_t, 6> be_modes{};

// Save queries, then check them only after later insertions have evicted the
// original selections. Constant values deliberately exercise complete FIFO
// ties between BE anchors with the same outer coordinates.
void check_historical_cutoffs() {
  using S = semirings::LogSumExpSemiring;
  using Cell = beam::WmbpWiPrefixCell<double>;
  struct Query { cand_pos_t j, l; std::vector<Cell> selected; };
  for (const std::size_t history_width : {0, 1, 3})
    for (const std::size_t prefix_width : {0, 1, 3}) {
      engine::Chart<double> chart(S::zero());
      beam::WmbpWiPrefixStore<S> store;
      beam::WmbpWiPrefixAdmission admission;
      std::vector<Query> queries;
      std::vector<Key> all_be, all_wi;
      for (cand_pos_t j = 10; j < 24; ++j) {
        const std::vector<Key> cores{
            {grammar::NonTerminal::BE, 1, j, 4, j - 2},
            {grammar::NonTerminal::BE, 1, j, 7, j - 2}};
        for (const auto &be : cores) {
          chart.set(be, j % 3 == 0 ? 1.0 : 0.0);
          all_be.push_back(be);
        }
        store.commit_be_column(cores, chart, prefix_width);
        admission.observe_be_column<S>(cores, chart, store);
        std::vector<cand_pos_t> starts{5, 8};
        for (cand_pos_t l = 11; l < j; ++l) starts.push_back(l);
        admission.observe_vp_column<S>(j, starts, store);
        for (const auto l : starts) queries.push_back({j, l, store.retained(l)});
        for (const cand_pos_t i : {5, 8}) {
          const Key wi{grammar::NonTerminal::WI, i, j, -1, -1};
          chart.set(wi, j % 4 == 0 ? 2.0 : 0.0);
          all_wi.push_back(wi);
        }
        store.commit_wi_column(j, {5, 8}, chart, history_width);
      }
      for (const auto &query : queries) {
        const auto test = [&](const Key &be, const Key *wi) {
          const auto l = wi ? wi->j + 1 : be.ip + 1;
          if (l != query.l) return;
          D source;
          source.rule = grammar::RuleId::WMBP_SPLIT_BE_WI_VP;
          source.parent = {grammar::NonTerminal::WMBP, be.i, query.j, -1, -1};
          source.child_count = wi ? 3 : 2;
          source.children[0].item = be;
          if (wi) source.children[1].item = *wi;
          source.children[source.child_count - 1].item =
              {grammar::NonTerminal::VP, l, query.j, -1, -1};
          const bool wanted = std::any_of(query.selected.begin(), query.selected.end(),
              [&](const auto &cell) {
                return deductions::ItemKeyEq{}(cell.be, be) &&
                       cell.has_wi == bool(wi) &&
                       (!wi || deductions::ItemKeyEq{}(cell.wi, *wi));
              });
          require(admission.allows<S>(source, chart) == wanted,
                  "historical cutoff changed a prefix admission/FIFO tie");
        };
        for (const auto &be : all_be) {
          test(be, nullptr);
          for (const auto &wi : all_wi)
            if (wi.i == be.ip + 1) test(be, &wi);
        }
      }
      require(admission.be_boundary_count() == all_be.size() &&
              admission.vp_boundary_count() == queries.size(),
              "cutoff history grew with the number of selected prefixes");
    }
}

template <class S>
void check(std::string sequence, std::string scaffold,
            const replay::WFinalVpWmbpOnlineWidths &widths) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  Capture<S> capture;
  const auto run = replay::run_w_final_vp_wmbp_online_observed<S>(
      ctx, n, widths, capture);
  const auto rebuilt = beam::WmbpWiPrefixAdmission::from_retained_chart<S>(
      n, run.chart, widths.wmbp_wi_history, widths.wmbp_wi_prefix);
  std::vector<std::vector<cand_pos_t>> vp(n + 1);
  std::size_t be_count = 0, vp_count = 0;
  for (const auto &[item, value] : run.chart.values()) {
    (void)value;
    if (item.nonterminal == grammar::NonTerminal::VP) {
      vp[item.j].push_back(item.i);
      ++vp_count;
    }
    if (item.nonterminal == grammar::NonTerminal::BE) ++be_count;
  }
  for (auto &starts : vp) std::sort(starts.begin(), starts.end());
  require(capture.admission.be_boundary_count() == be_count &&
          capture.admission.vp_boundary_count() == vp_count &&
          rebuilt.be_boundary_count() == be_count &&
          rebuilt.vp_boundary_count() == vp_count,
          "admission storage is not one boundary per retained BE/VP");

  using L = semirings::LogSumExpSemiring;
  using T = semirings::TapeLogSumExpSemiring;
  semirings::LogSumExpTape tape;
  T::Scope scope(tape);
  engine::Chart<double> alpha(L::zero()), beta(L::zero());
  engine::Chart<T::value_type> leaves(T::zero());
  for (const auto &[item, value] : run.chart.values()) {
    const double log_alpha = std::is_same_v<S, semirings::SumProductSemiring> ||
                             std::is_same_v<S, semirings::MaxProductSemiring>
        ? std::log(S::priority(value)) : S::priority(value);
    alpha.set(item, log_alpha);
    leaves.set(item, tape.constant(log_alpha));
  }
  auto objective = T::zero();
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);
  auto be_provider = beam::BeRegeneratedSources<S>::from_retained_chart(n, run.chart, ctx);
  decltype(capture.sources) regenerated_be;
  for (cand_pos_t j = n; j >= 1; --j) {
    be_provider.prepare_column(j, run.chart, ctx, widths.be_right_tail, widths.be_right_shell);
    for (const auto &[item, value] : run.chart.values()) {
      (void)value;
      if (item.nonterminal == grammar::NonTerminal::BE && item.j == j)
        regenerated_be[item] = be_provider.deductions_for(item, run.chart, ctx);
    }
  }
  for (const auto &[parent, value] : run.chart.values()) {
    const bool is_be = parent.nonterminal == grammar::NonTerminal::BE;
    if (parent.nonterminal != grammar::NonTerminal::WMBP && !is_be) continue;
    const auto actual = is_be ? regenerated_be.at(parent) : beam::regenerate_wmbp_sources<S>(
        parent, vp[parent.j], run.chart, ctx, rebuilt, widths.wmbp_prefix);
    const auto observed = is_be ? actual : beam::regenerate_wmbp_sources<S>(
        parent, vp[parent.j], run.chart, ctx, capture.admission, widths.wmbp_prefix);
    std::multiset<decltype(id(D{}))> wanted, got, seen;
    auto total = S::zero();
    for (const auto &d : capture.sources[parent]) wanted.insert(id(d));
    for (const auto &d : observed) seen.insert(id(d));
    for (const auto &d : actual) {
      got.insert(id(d));
      total = S::plus(total, replay::runtime_deduction_contribution<S>(
          run.chart, d, weight.weight_for(d)));
      engine::propagate_outside_deduction<L>(d, alpha, L::one(), &beta,
          [&](const auto &source) { return weight.weight_for(source); },
          engine::RuntimeDeductionChildItem{});
      if (is_be) ++be_modes[static_cast<std::size_t>(d.rule) -
                            static_cast<std::size_t>(grammar::RuleId::BE_BASE_SAMEPAIR)];
      else if (d.rule == grammar::RuleId::WMBP_DIRECT_VP) ++modes[4];
      else ++modes[static_cast<std::size_t>(beam::wmbp_be_vp_prefix_key(d).mode)];
    }
    if (got != wanted || seen != wanted || !close(S::priority(total), S::priority(value))) {
      std::cerr << "n=" << n << " b=" << widths.vp << " parent=" << parent.i
                << ',' << parent.j << " nt=" << int(parent.nonterminal) << " regenerated=" << got.size()
                << " forward=" << wanted.size() << '\n';
      require(false, "BE/WMBP regeneration changed admitted branches/alpha");
    }
    // The independent arithmetic tape reads observed forward branches, not
    // the regenerated provider whose derivatives are being checked.
    for (const auto &d : capture.sources[parent])
      objective = T::plus(objective, replay::runtime_deduction_contribution<T>(
          leaves, d, weight.weight_for(d)));
    ++parents;
  }
  if (!T::is_zero(objective)) {
    const auto expected = tape.reverse(objective.node);
    for (const auto &[item, value] : leaves.values()) {
      const double actual = std::exp(beta.get(item) + alpha.get(item) - objective.log_value);
      require(close(actual, expected[value.node]), "BE/WMBP local outside differs from tape");
    }
  }
  ++cases;
}
std::pair<std::string, std::string> random_case(std::mt19937 &rng, int n) {
  std::string g(n, '.');
  std::function<void(int, int)> fill = [&](int l, int r) {
    if (r - l < 4) return;
    if (rng() % 3 == 0) {
      const int k = l + 4 + rng() % (r - l - 3);
      g[l] = '('; g[k] = ')'; fill(l + 1, k - 1); fill(k + 1, r);
    } else fill(l + 1, r);
  };
  fill(0, n - 1);
  std::string seq(n, 'A');
  for (auto &c : seq) c = "ACGU"[rng() % 4];
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
  check_historical_cutoffs();
  using replay::WFinalVpWmbpOnlineWidths;
  std::vector<std::pair<std::string, std::string>> inputs{
      {"GCGCGCGCGCGC", "..(....)...."},
      {"AGGUACCCUACAAUUACGGCAUCAGUCUAAAAACUCUGGGAUGGGCA",
       "(.(....).....)..(..(....)........(....)).(....)"},
      {"GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
       "....(...)((....)..)((.(.((...))))(.(....)((...))))"}};
  std::mt19937 rng(20261009);
  const int random_count = argc > 1 ? std::atoi(argv[1]) : 10;
  for (int x = 0; x < random_count; ++x) inputs.push_back(random_case(rng, 20 + rng() % 41));
  for (const auto &[seq, g] : inputs) {
    std::vector<WFinalVpWmbpOnlineWidths> settings;
    for (const std::size_t b : {1, 3, 20, 100}) settings.push_back(WFinalVpWmbpOnlineWidths::uniform(b));
    settings.push_back({});
    auto asymmetric = WFinalVpWmbpOnlineWidths::uniform(20);
    asymmetric.wmbp_wi_history = 1;
    asymmetric.wmbp_wi_prefix = 3;
    asymmetric.wmbp_prefix = 2;
    settings.push_back(asymmetric);
    asymmetric.wmbp_wi_prefix = 0;
    settings.push_back(asymmetric);
    for (const auto &widths : settings) {
      check<semirings::LogSumExpSemiring>(seq, g, widths);
      check<semirings::LogMaxSemiring>(seq, g, widths);
    }
  }
  for (const auto count : modes) require(count > 0, "a WMBP mode has no witness");
  for (const auto count : be_modes) require(count > 0, "a BE rule has no witness");
  std::cout << "scfg2_be_wmbp_regenerated_sources=ok cases=" << cases
            << " parents=" << parents << " modes=";
  for (const auto count : modes) std::cout << count << ',';
  std::cout << " be_rules=";
  for (const auto count : be_modes) std::cout << count << ',';
  std::cout << '\n';
}
