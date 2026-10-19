#include "scfg2/weights/runtime_snapshot.hh"
// Empty WIP spans are impossible; an empty unpaired BE wing uses EMPTY.
#include "scfg2/beam/retained_outside.hh"
#include "scfg2/geometry/static_scaffold_geometry.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <tuple>
#include <vector>

namespace {
using namespace cparty::scfg2;
using R = grammar::RuleId;
using NT = grammar::NonTerminal;
using Key = deductions::ItemKey;
using S = semirings::LogSumExpSemiring;
using M = semirings::LogMaxSemiring;
std::size_t scaffolds = 0, parents = 0, sources = 0, false_empty_queries = 0;
std::size_t numeric_checks = 0, beam_checks = 0;

void require(bool ok, const std::string &message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
bool close(double a, double b) { return std::fabs(a - b) < 1e-10; }

class GeometryContext : public deductions::PairabilityContext {
 public:
  GeometryContext(const geometry::StaticScaffoldGeometry &g, bool empty_is_true)
      : g_(g), empty_is_true_(empty_is_true) {}
  cand_pos_t turn() const override { return 3; }
  cand_pos_t sequence_length() const override { return g_.size(); }
  // Geometry-only comparison; real sequence/energy checks are separate below.
  bool can_pair(cand_pos_t, cand_pos_t) const override { return true; }
  cand_pos_t pair_partner(cand_pos_t p) const override { return g_.partner(p); }
  cand_pos_t unpaired_prefix(cand_pos_t p) const override { return g_.unpaired_prefix(p); }
  bool is_unpaired(cand_pos_t p) const override {
    return p >= 1 && p <= g_.size() && g_.partner(p) < 0;
  }
  bool weakly_closed(cand_pos_t i, cand_pos_t j) const override {
    return empty_is_true_ && i == j + 1 ? true : g_.weakly_closed(i, j);
  }
 private:
  const geometry::StaticScaffoldGeometry &g_;
  bool empty_is_true_;
};

auto signatures(const std::vector<deductions::Deduction> &ds) {
  std::vector<std::tuple<R, cand_pos_t, cand_pos_t>> out;
  for (const auto &d : ds) out.emplace_back(d.rule, d.split.k, d.split.l);
  return out;
}

void check_geometry(const std::string &scaffold) {
  const geometry::StaticScaffoldGeometry g(scaffold, scaffold.size());
  require(g.usable(), "invalid BE geometry witness");
  for (int p = 1; p <= g.size(); ++p)
    if (g.partner(p) > p && g.partner(p) - p < 4) return;
  ++scaffolds;
  GeometryContext current(g, false), normalized(g, true);
  // The mathematical WC equivalence is for nonempty intervals, including
  // singletons. Check it against partner containment independently.
  for (int i = 1; i <= g.size(); ++i) for (int j = i; j <= g.size(); ++j) {
    bool expected = true;
    for (int p = i; p <= j; ++p)
      if (g.partner(p) > 0 && (g.partner(p) < i || g.partner(p) > j)) expected = false;
    require(g.weakly_closed(i, j) == expected, "nonempty WC containment mismatch");
  }
  for (int p = 1; p <= g.size() + 1; ++p)
    false_empty_queries += !g.weakly_closed(p, p - 1);
  for (int i = 1; i <= g.size(); ++i) {
    const int j = g.partner(i);
    if (j <= i) continue;
    for (int ip = i; ip < j; ++ip) {
      const int jp = g.partner(ip);
      if (jp <= ip || jp > j) continue;
      const Key parent{NT::BE, i, j, ip, jp};
      ++parents;
      for (const auto rule : {R::BE_INTERNAL_LOOP, R::BE_WIP_WIP,
                              R::BE_WIP_BASEPAIR, R::BE_BASEPAIR_WIP}) {
        const auto actual = deductions::generate_exact_basic_deductions(parent, rule, current);
        const auto changed = deductions::generate_exact_basic_deductions(parent, rule, normalized);
        require(signatures(actual) == signatures(changed),
                "normalizing adjacent-empty WC changed valid BE sources: " + scaffold);
        sources += actual.size();
      }
    }
  }
}

void enumerate(std::string &g, int balance, int remaining) {
  if (!remaining) { if (!balance) check_geometry(g); return; }
  for (char c : {'.', '(', ')'}) {
    const int next = balance + (c == '(') - (c == ')');
    if (next < 0 || next > remaining - 1) continue;
    g.push_back(c); enumerate(g, next, remaining - 1); g.pop_back();
  }
}

void check_wing(std::string seq, std::string g, int ip, int jp,
                R expected_rule, double energy, WmbClosurePolicy policy) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
  replay::WFinalLogWeightContext ctx(owner, policy);
  const auto n = static_cast<cand_pos_t>(seq.size());
  const Key parent{NT::BE, 1, n, ip, jp};
  const bool left_empty = ip == 2;
  require(!ctx.weakly_closed(left_empty ? 2 : n, left_empty ? 1 : n - 1),
          "BE empty-wing WC should be false at paired endpoints");
  require(ctx.is_empty_region(left_empty ? 1 : jp, left_empty ? ip : n),
          "BE empty wing failed EMPTY");
  const grammar::RuleSplit split{ip, jp};
  for (const auto rule : {R::BE_WIP_WIP,
                          left_empty ? R::BE_WIP_BASEPAIR : R::BE_BASEPAIR_WIP})
    require(!deductions::is_valid(deductions::expand_basic_rule(parent, rule, split)),
            "empty WIP child was accepted as a unit");
  require(deductions::exact_be_split_allows(parent, expected_rule, ip, jp, ctx),
          "valid opposite-wing rule was rejected");

  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  const auto allow = [](const auto &) { return true; };
  const auto ignore = [](const auto &, const auto &, const auto &, const auto &, auto) {};
  const auto inside = replay::run_w_final_exact_dp_schedule<S>(schedule, ctx, allow, ignore);
  const auto best = replay::run_w_final_exact_dp_schedule<M>(schedule, ctx, allow, ignore);
  require(close(owner.partition_free_energy_from_log_weight(inside.get(parent)), energy) &&
              close(owner.partition_free_energy_from_log_weight(best.get(parent)), energy),
          "BE empty-wing independent loop-sum energy mismatch: " + g);
  engine::AllLocalWeightModel<replay::WFinalLogWeightContext> weights(&ctx);
  std::size_t positive = 0;
  for (const auto &d : deductions::generate_exact_basic_deductions_for_item(parent, ctx)) {
    const double term = replay::runtime_deduction_contribution<S>(inside, d, weights.weight_for(d));
    if (S::is_zero(term)) continue;
    ++positive;
    require(d.rule == expected_rule && d.split.k == ip && d.split.l == jp,
            "unexpected positive BE empty-wing source");
  }
  require(positive == 1, "BE empty-wing witness has missing or duplicate mass");
  ++numeric_checks;

  for (const auto &widths : {replay::WFinalVpWmbpOnlineWidths::uniform(20),
                            replay::WFinalVpWmbpOnlineWidths{}}) {
    beam::RetainedRuleAdmission<S> admission;
    const auto run = replay::run_w_final_vp_wmbp_online_observed<S>(ctx, n, widths, admission);
    admission.finish(run.chart);
    require(close(run.chart.get(parent), inside.get(parent)), "beam lost BE empty-wing mass");
    beam::RetainedOutsideProgram<S, replay::WFinalLogWeightContext> program(n, run.chart, ctx, admission, widths);
    program.prepare_column(n);
    const auto regenerated = program.deductions_for(parent);
    require(regenerated.size() == 1 && regenerated.front().rule == expected_rule,
            "outside did not regenerate the BE empty-wing rule");
    const auto reverse = program.outside(parent);
    const Key wing{NT::WIP, left_empty ? jp + 1 : 2, left_empty ? n - 1 : ip - 1, -1, -1};
    const double occurrence = std::exp(reverse.beta.get(wing) + run.chart.get(wing) - run.chart.get(parent));
    require(close(occurrence, 1.0), "outside lost the nonempty BE wing");
    ++beam_checks;
  }
}

void check_stack() {
  std::string seq = "GGAAACC", g = "((...))";
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
  replay::WFinalLogWeightContext ctx(owner);
  const Key parent{NT::BE, 1, 7, 2, 6};
  require(!ctx.weakly_closed(2, 1) && !ctx.weakly_closed(7, 6),
          "stack's adjacent-empty WC convention changed");
  engine::AllLocalWeightModel<replay::WFinalLogWeightContext> weights(&ctx);
  std::size_t stack = 0, internal = 0;
  for (const auto &d : deductions::generate_exact_basic_deductions_for_item(parent, ctx)) {
    if (d.rule == R::BE_STACK) { ++stack; require(weights.weight_for(d) > 0., "stack lost local weight"); }
    if (d.rule == R::BE_INTERNAL_LOOP) { ++internal; require(weights.weight_for(d) == 0., "stack double counted as internal loop"); }
  }
  require(stack == 1 && internal == 1, "both-empty stack control was not exercised");
}
}  // namespace

int main() {
  for (int n = 1; n <= 13; ++n) { std::string g; enumerate(g, 0, n); }
  require(false_empty_queries > 0, "empty WC convention was not exercised");
  for (const auto *params : {"params/rna_Turner04.par", "params/rna_DirksPierce09.par"}) {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", params, 1);
    // H3 + b' + c' + a' + 2b', with a'=3.41, b'=0.56, c'=0.12.
    const double energy = std::string(params).find("Turner04") != std::string::npos ? 10.61 : 8.90;
    for (const auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible}) {
      check_wing("GGAAACGAAACAC", "((...)(...).)", 2, 6, R::BE_BASEPAIR_WIP, energy, policy);
      check_wing("GAGAAACGAAACC", "(.(...)(...))", 8, 12, R::BE_WIP_BASEPAIR, energy, policy);
    }
    check_stack();
  }
  std::cout << "scfg2_be_empty_guard=ok scaffolds=" << scaffolds << " BE_parents=" << parents
            << " split_sources=" << sources << " false_empty_queries=" << false_empty_queries
            << " numeric_checks=" << numeric_checks << " beam_checks=" << beam_checks << '\n';
}
