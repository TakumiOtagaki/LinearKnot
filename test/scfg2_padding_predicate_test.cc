#include "scfg2/weights/runtime_snapshot.hh"
// Compare LP/RP with independent scaffold endpoint counts, then check that
// first-column WM/WIP/VPL auxiliaries cannot enter the public W(1,n) root.
#include "scfg2/deductions/exact_context.hh"
#include "scfg2/geometry/static_scaffold_geometry.hh"
#include "scfg2/beam/retained_outside.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {
using namespace cparty::scfg2;
using NT = grammar::NonTerminal;
using Key = deductions::ItemKey;
using S = semirings::LogSumExpSemiring;
std::size_t scaffolds = 0, queries = 0, leading_differences = 0;
std::size_t reachable_items = 0, unused_positive_items = 0, beam_checks = 0;

void require(bool ok, const std::string &message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}

class GeometryContext : public deductions::PairabilityContext {
 public:
  explicit GeometryContext(const geometry::StaticScaffoldGeometry &g) : g_(g) {}
  cand_pos_t turn() const override { return 3; }
  cand_pos_t sequence_length() const override { return g_.size(); }
  bool can_pair(cand_pos_t, cand_pos_t) const override { return false; }
  cand_pos_t unpaired_prefix(cand_pos_t p) const override { return g_.unpaired_prefix(p); }
 private:
  const geometry::StaticScaffoldGeometry &g_;
};

void check_geometry(const std::string &g) {
  const int n = static_cast<int>(g.size());
  geometry::StaticScaffoldGeometry index(g, n);
  GeometryContext ctx(index);
  require(index.usable(), "invalid geometry in padding witness");
  // This oracle reads only the input characters, independently of up[],
  // sparse_tree and the production geometry implementation.
  std::vector<int> endpoints(n + 1);
  int last_endpoint = 0;
  for (int p = 1; p <= n; ++p) {
    const bool paired = g[p - 1] == '(' || g[p - 1] == ')';
    endpoints[p] = endpoints[p - 1] + paired;
    if (paired) last_endpoint = p;
    const int expected_up = paired ? 0 : p - (last_endpoint ? last_endpoint : 1);
    require(index.unpaired_prefix(p) == expected_up, "up[] differs from endpoint-count oracle: " + g);
  }
  const auto free = [&](int lo, int hi) {
    return lo > hi || endpoints[hi] == endpoints[lo - 1];
  };
  for (int a = 1; a <= n; ++a) for (int b = 1; b <= n; ++b) {
    const bool mathematical_left = b >= a && free(a, b - 1);
    const bool legacy_leading = a == 1 && b > a && mathematical_left;
    const bool mathematical_right = b >= a && free(a + 1, b);
    const bool left = ctx.can_pair_left_span(a, b);
    require(left == (mathematical_left && !legacy_leading),
            "unexpected LP mismatch G=" + g + " args=" + std::to_string(a) + "," + std::to_string(b));
    require(ctx.can_pair_right_span(a, b) == mathematical_right,
            "RP differs from G-free scan G=" + g + " args=" + std::to_string(a) + "," + std::to_string(b));
    leading_differences += legacy_leading;
    ++queries;
  }
  ++scaffolds;
}

void enumerate(std::string &g, int balance, int remaining) {
  if (!remaining) {
    if (balance) return;
    check_geometry(g);
    auto marked = g;
    bool changed = false;
    for (std::size_t p = 0; p < marked.size(); ++p)
      if (marked[p] == '.') { marked[p] = p % 2 ? 'X' : 'x'; changed = true; }
    if (changed) check_geometry(marked);
    return;
  }
  for (char c : {'.', '(', ')'}) {
    const int next = balance + (c == '(') - (c == ')');
    if (next < 0 || next > remaining - 1) continue;
    g.push_back(c); enumerate(g, next, remaining - 1); g.pop_back();
  }
}

bool unused_left_auxiliary(const Key &key) {
  return key.i == 1 && (key.nonterminal == NT::WM ||
      key.nonterminal == NT::WIP || key.nonterminal == NT::VPL);
}

template <class Semiring>
void check_exact(replay::WFinalLogWeightContext &ctx, cand_pos_t n) {
  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  const auto ignore = [](const auto &, const auto &, const auto &, const auto &, auto) {};
  const auto full = replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule, ctx, [](const auto &) { return true; }, ignore);
  // Removing entire unused auxiliary states is a stronger intervention than
  // changing only their padded LP(1,k) sources. Do not apply this filter to
  // the production beam: changing the competing pool needs separate study.
  const auto removed = replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule, ctx, [](const auto &d) { return !unused_left_auxiliary(d.parent); }, ignore);
  const Key root{NT::W, 1, n, -1, -1};
  require(std::isfinite(full.get(root)) && full.get(root) == removed.get(root),
          "unused leading auxiliaries affected the exact root");
  for (const auto &[key, value] : full.values())
    unused_positive_items += unused_left_auxiliary(key) && !Semiring::is_zero(value);
}

void check_root(std::string seq, std::string g, WmbClosurePolicy policy) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
  replay::WFinalLogWeightContext ctx(owner, policy);
  const auto n = static_cast<cand_pos_t>(seq.size());
  const Key root{NT::W, 1, n, -1, -1};
  std::vector<Key> pending{root};
  std::unordered_set<Key, deductions::ItemKeyHash, deductions::ItemKeyEq> seen;
  while (!pending.empty()) {
    const auto parent = pending.back(); pending.pop_back();
    if (!seen.insert(parent).second) continue;
    require(!unused_left_auxiliary(parent), "root reaches a leading LP auxiliary");
    for (const auto &d : deductions::generate_exact_basic_deductions_for_item(parent, ctx)) {
      if (!engine::legacy_span_schedule_supports(d)) continue;
      for (std::size_t c = 0; c < d.child_count; ++c)
        pending.push_back(engine::RuntimeDeductionChildItem{}(d, c));
    }
  }
  reachable_items += seen.size();
  check_exact<S>(ctx, n);
  check_exact<semirings::LogMaxSemiring>(ctx, n);
  for (const std::size_t width : {1, 3, 20}) {
    const auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(width);
    beam::RetainedRuleAdmission<S> admission;
    const auto run = replay::run_w_final_vp_wmbp_online_observed<S>(ctx, n, widths, admission);
    admission.finish(run.chart);
    require(std::isfinite(run.chart.get(root)), "beam lost padding audit root");
    beam::RetainedOutsideProgram<S, replay::WFinalLogWeightContext> program(n, run.chart, ctx, admission, widths);
    const auto reverse = program.outside(root);
    for (const auto &[key, value] : reverse.beta.values())
      require(!unused_left_auxiliary(key) || S::is_zero(value),
              "beam outside reached an unused leading auxiliary");
    ++beam_checks;
  }
}
}  // namespace

int main() {
  for (int n = 1; n <= 9; ++n) { std::string g; enumerate(g, 0, n); }
  require(leading_differences > 0, "leading LP exception was not exercised");
  for (const auto *params : {"params/rna_Turner04.par", "params/rna_DirksPierce09.par"}) {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", params, 1);
    for (const auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible}) {
      check_root("GCAAAGC", ".......", policy);
      check_root("GAGAAACAAGAAGAAAAAAAAACCAAAAC", "(.(...)..(............)).....", policy);
      check_root("GGGAAAAACGGAAAAACAAAACCAC", "(.........(.....)....)...", policy);
      check_root("UUAUAUGGGUAUUAUUCUCA", "((x...))((...))(...)", policy);
    }
  }
  require(unused_positive_items > 0, "root-removal check was vacuous");
  std::cout << "scfg2_padding_predicate=ok scaffolds=" << scaffolds << " queries=" << queries
            << " leading_lp_differences=" << leading_differences
            << " reachable_items=" << reachable_items << " unused_positive_items=" << unused_positive_items
            << " beam_checks=" << beam_checks << '\n';
}
