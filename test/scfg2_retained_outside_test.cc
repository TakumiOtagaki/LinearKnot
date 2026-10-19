#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/semirings/logsumexp_tape.hh"
#include "scfg2/semirings/value_boundary.hh"

// Oracle-only identity nodes distinguish a public chart assignment from a
// private prefix/carry value that happens to have the identical alpha. Plain
// tape aliases otherwise report the private value's later reuse as public
// item posterior (notably KeepAll WIP and finite VPR).
namespace cparty::scfg2::semirings {
template <> struct ValueBoundary<LogSumExpTapeValue> {
  static LogSumExpTapeValue checkpoint(LogSumExpTapeValue value) {
    return TapeLogSumExpSemiring::is_zero(value) ? value :
        TapeLogSumExpSemiring::times(value, TapeLogSumExpSemiring::one());
  }
  static bool discard(const LogSumExpTapeValue &) { return false; }
};
}

#include "scfg2/beam/retained_outside.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>

namespace {
using namespace cparty::scfg2;
using S = semirings::LogSumExpSemiring;
using T = semirings::TapeLogSumExpSemiring;
using NT = grammar::NonTerminal;
using Key = deductions::ItemKey;
std::size_t cases = 0, parents = 0;
std::size_t retired_segments = 0;
bool close(double a, double b) {
  return a == b || std::fabs(a - b) <= 2e-9 * std::max({1.0, std::fabs(a), std::fabs(b)});
}
void require(bool ok, const char *message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}

void check_case3_target_retirement() {
  using Pending = beam::BasicPendingVpCase3LeftCandidates<double>;
  const auto source = [](cand_pos_t i, cand_pos_t j) {
    deductions::Deduction d;
    d.rule = grammar::RuleId::VP_WI_CASE3;
    d.parent = {NT::VP, i, j, -1, -1};
    d.split.p = i + 3;
    d.split.q = i + 4;
    d.split.k = i + 6;
    d.split.l = i + 7;
    d.child_count = 1;
    d.children[0].item = {NT::WI, i + 1, i + 2, -1, -1};
    return beam::VpCase3PlannedSource{{j, d.split.q, d.split.k, d.split.l}, d};
  };
  for (const auto width : std::vector<std::size_t>{0, 1, 3, 20, std::numeric_limits<std::size_t>::max()}) {
    Pending pending(width);
    engine::Chart<double> chart(S::zero());
    for (cand_pos_t j = 100; j < 164; ++j) {
      std::vector<beam::VpCase3PlannedSource> query;
      for (cand_pos_t i = 1; i <= 16; ++i) {
        query.push_back(source(i, j));
        chart.set(query.back().deduction.children[0].item, S::one());
      }
      pending.observe_target_sources<S>(query, chart);
      beam::BasicExactVpCase3TailCache<double> tails;
      pending.prepare_column<S>(j, chart, tails);
      require(pending.live_segment_records() <= 48, "CASE3 reservations grow across retired target columns");
      const auto registered = pending.registered_segment_count();
      pending.release_target_column(j);
      require(pending.live_segment_records() == 0 && pending.registered_segment_count() == registered,
              "CASE3 retirement retained entries or lost aggregate counters");
      retired_segments += query.size();
    }
  }
  Pending pending(3);
  engine::Chart<double> chart(S::zero());
  const auto first = source(1, 100);
  chart.set(first.deduction.children[0].item, S::one());
  pending.observe_target_sources<S>({first, source(1, 101)}, chart);
  const auto live = pending.live_segment_records();
  bool rejected = false;
  try { pending.release_target_column(100); } catch (const std::logic_error &) { rejected = true; }
  require(rejected && pending.live_segment_records() == live,
          "CASE3 retirement discarded a live future target");
  pending.release_target_column(101);
  require(pending.live_segment_records() == 0, "completed future target did not retire");
}
void check(std::string seq, std::string scaffold,
           const replay::WFinalVpWmbpOnlineWidths &widths) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(seq.size());
  beam::RetainedRuleAdmission<S> rules;
  const auto run = replay::run_w_final_vp_wmbp_online_observed<S>(ctx, n, widths, rules);
  rules.finish(run.chart);
  const auto viterbi = replay::run_w_final_vp_wmbp_online_viterbi_with_semiring<
      semirings::LogMaxSemiring>(ctx, n, widths);
  require(viterbi.choices.values().size() == viterbi.dp.chart.values().size(),
          "Viterbi choices retain pruned parents or lose retained parents");
  for (const auto &[key, choice] : viterbi.choices.values()) {
    (void)choice;
    require(std::isfinite(viterbi.dp.chart.get(key)), "pruned Viterbi choice survived its barrier");
  }
  beam::RetainedOutsideProgram<S, replay::WFinalExactInsideContext> program(n, run.chart, ctx, rules, widths);
  if (widths.vpr != std::numeric_limits<std::size_t>::max() &&
      widths.wmb_btail != std::numeric_limits<std::size_t>::max())
    require(program.private_cell_count() <= std::size_t(n) * (widths.vpr + widths.wmb_btail),
            "private reverse cells exceed the two nb budgets");
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);
  for (cand_pos_t j = n; j >= 1; --j) {
    program.prepare_column(j);
    for (auto nt : beam::vp_wmbp_online_column_phases) {
      for (const auto &parent : program.views().items(nt, j)) {
        auto alpha = program.private_contribution(parent);
        const auto sources = program.deductions_for(parent);
        for (const auto &d : sources) {
          auto term = S::from_weight(weight.weight_for(d));
          for (std::size_t c = 0; c < d.child_count; ++c)
            term = S::times(term, run.chart.get(engine::RuntimeDeductionChildItem{}(d, c)));
          alpha = S::plus(alpha, term);
        }
        if (!close(alpha, run.chart.get(parent))) {
          std::cerr << "alpha mismatch case=" << cases << " nt=" << grammar::nonterminal_name(nt)
                    << " i=" << parent.i << " j=" << parent.j << " ip=" << parent.ip
                    << " jp=" << parent.jp << " forward=" << run.chart.get(parent)
                    << " rebuilt=" << alpha << " terms=" << sources.size() << '\n';
          for (const auto &d : sources) std::cerr << grammar::rule_name(d.rule) << " k=" << d.split.k << '\n';
          std::exit(EXIT_FAILURE);
        }
        ++parents;
      }
    }
  }
  const Key root{NT::W, 1, n, -1, -1};
  const auto reverse = program.outside(root);
  semirings::LogSumExpTape tape;
  const auto reference = [&] {
    T::Scope scope(tape);
    return replay::run_w_final_vp_wmbp_online<T>(ctx, n, widths);
  }();
  require(close(reference.chart.get(root).log_value, run.chart.get(root)), "tape/numeric root mismatch");
  if (!S::is_zero(run.chart.get(root))) {
    const auto posterior = tape.reverse(reference.chart.get(root).node);
    for (const auto &[key, value] : reference.chart.values()) {
      const double expected = posterior[value.node];
      const double actual = std::exp(reverse.beta.get(key) + run.chart.get(key) - run.chart.get(root));
      if (!close(expected, actual)) {
        std::cerr << "outside mismatch case=" << cases << " nt=" << grammar::nonterminal_name(key.nonterminal)
                  << " i=" << key.i << " j=" << key.j << " ip=" << key.ip << " jp=" << key.jp
                  << " expected=" << expected << " actual=" << actual << '\n';
        std::exit(EXIT_FAILURE);
      }
    }
    for (const auto &[key, value] : reverse.beta.values()) {
      (void)value;
      require(!S::is_zero(run.chart.get(key)), "outside reached a removed public state");
    }
  }
  ++cases;
}
}  // namespace

int main(int argc, char **argv) {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  check_case3_target_retirement();
  {
    // A root-only WMB winner must not become a new WM/WIP start seed when
    // those families are KeepAll and WMB remains finite. The dense prefixes
    // previously read the entire WMB union and disagreed with regeneration.
    const std::string seq =
        "CUUACUGUGUGAGCCUCACCACGGUACCACUGUGCAACCAAUGGCUUUCCAACGUUCAGGAAAUAAGUAGAGA";
    const std::string g =
        ".(x((...(.(...x..x....(.x.)).)....((....x)))..).(....).(....)..(....)...)";
    for (const auto family : {grammar::NonTerminal::WM, grammar::NonTerminal::WIP}) {
      auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(20);
      widths.force_scaffold_derivation = false;
      widths.wmb = 1;
      if (family == grammar::NonTerminal::WM)
        widths.wm = std::numeric_limits<std::size_t>::max();
      else
        widths.wip = std::numeric_limits<std::size_t>::max();
      check(seq, g, widths);
    }
  }
  std::vector<std::pair<std::string, std::string>> inputs{
      {"GGGAAACCCAGGGAAACCC", "(((...))).........."},
      {"GCGCGCGCGCGC", "..(....)...."},
      {"AGGUACCCUACAAUUACGGCAUCAGUCUAAAAACUCUGGGAUGGGCA",
       "(.(....).....)..(..(....)........(....)).(....)"},
      {"GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
       "....(...)((....)..)((.(.((...))))(.(....)((...))))"},
      {"GGGAAAAACGGAAAAACAAAACCAC", "(.........(.....)....)..."}};
  {
    std::ifstream fixture("test/fixtures/scfg2_wmbp_leftward_be_20261008/stem20.txt");
    std::string seq, g;
    require(bool(std::getline(fixture, seq) && std::getline(fixture, g)), "missing large BE fixture");
    inputs.emplace_back(seq, g);
  }
  std::mt19937 rng(20261009);
  const int random_count = argc > 1 ? std::atoi(argv[1]) : 4;
  for (int x = 0; x < random_count; ++x) {
    auto [seq, g] = inputs[2 + x % 2];
    for (auto &c : seq) c = "ACGU"[rng() % 4];
    std::vector<std::size_t> open;
    for (std::size_t i = 0; i < g.size(); ++i) {
      if (g[i] == '(') open.push_back(i);
      if (g[i] == ')') { const auto k = open.back(); open.pop_back(); seq[i] = "UGCA"[std::string("ACGU").find(seq[k])]; }
    }
    inputs.emplace_back(seq, g);
  }
  for (const auto &[seq, g] : inputs) {
    for (const std::size_t b : {0, 1, 3, 20, 100}) check(seq, g, replay::WFinalVpWmbpOnlineWidths::uniform(b));
    check(seq, g, replay::WFinalVpWmbpOnlineWidths{});
    for (const std::size_t private_width : {0, 1, 3}) {
      auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(20);
      widths.wmbp_wi_history = private_width;
      widths.wmbp_wi_prefix = private_width;
      widths.wmbp_prefix = private_width;
      widths.be_right_tail = private_width;
      widths.be_right_shell = private_width;
      widths.wmb_btail = private_width;
      widths.vpr = private_width;
      check(seq, g, widths);
    }
    for (const std::size_t b : {1, 20}) {
      auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(b);
      widths.force_scaffold_derivation = false;
      check(seq, g, widths);
    }
  }
  std::cout << "scfg2_retained_outside=ok cases=" << cases << " parents=" << parents
            << " retired_segments=" << retired_segments << '\n';
}
