#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/engine/legacy_span_support.hh"
#include "scfg2/geometry/static_scaffold_geometry.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cstdlib>
#include <iostream>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace cparty::scfg2;
using D = deductions::Deduction;
using Key = deductions::ItemKey;
using NT = grammar::NonTerminal;
using R = grammar::RuleId;
std::size_t scaffolds = 0, items = 0, candidates = 0;
std::set<R> covered_rules;

void require(bool ok, const std::string &message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}

// Geometry-only superset: every ordered pair is sequence-compatible. This
// exercises split bounds independently of nucleotide/weight zero filters.
class GeometryContext : public deductions::PairabilityContext {
 public:
  GeometryContext(const geometry::StaticScaffoldGeometry &g, WmbClosurePolicy p)
      : g_(g), policy_(p) {}
  cand_pos_t turn() const override { return 3; }
  cand_pos_t sequence_length() const override { return g_.size(); }
  bool can_pair(cand_pos_t i, cand_pos_t j) const override {
    return i >= 1 && i < j && j <= g_.size();
  }
  cand_pos_t pair_partner(cand_pos_t p) const override { return g_.partner(p); }
  cand_pos_t parent_index(cand_pos_t p) const override { return g_.parent(p); }
  cand_pos_t unpaired_prefix(cand_pos_t p) const override { return g_.unpaired_prefix(p); }
  cand_pos_t border_b(cand_pos_t i, cand_pos_t j) const override { return g_.b(i, j); }
  cand_pos_t border_bp(cand_pos_t i, cand_pos_t j) const override { return g_.bp(i, j); }
  cand_pos_t border_B(cand_pos_t i, cand_pos_t j) const override { return g_.B(i, j); }
  cand_pos_t border_Bp(cand_pos_t i, cand_pos_t j) const override { return g_.Bp(i, j); }
  bool is_unpaired(cand_pos_t p) const override {
    return p >= 1 && p <= g_.size() && g_.partner(p) < 0;
  }
  bool weakly_closed(cand_pos_t i, cand_pos_t j) const override { return g_.weakly_closed(i, j); }
  WmbClosurePolicy wmb_closure_policy() const override { return policy_; }
 private:
  const geometry::StaticScaffoldGeometry &g_;
  WmbClosurePolicy policy_;
};

template <class Context>
void check_generated(cand_pos_t n, const Context &ctx) {
  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  for (const auto &parent : schedule) {
    ++items;
    for (const auto &d : deductions::generate_exact_basic_deductions_for_item(parent, ctx)) {
      require(engine::legacy_span_schedule_supports(d), "generated candidate violates canonical rank");
      engine::require_legacy_span_schedule_support(d);
      ++candidates;
      covered_rules.insert(d.rule);
    }
  }
}

void check_scaffold(const std::string &scaffold, bool include_short_pairs = false) {
  const geometry::StaticScaffoldGeometry g(scaffold, scaffold.size());
  if (!include_short_pairs)
    for (int p = 1; p <= g.size(); ++p)
      if (g.partner(p) > p && g.partner(p) - p < 4) return;
  for (const auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible}) {
    GeometryContext ctx(g, policy);
    check_generated(g.size(), ctx);
  }
  ++scaffolds;
}

void enumerate(std::string &g, int balance, int remaining) {
  if (!remaining) {
    if (balance) return;
    // Input syntax permits short fixed pairs even when the thermal weight
    // vanishes. They must not become internal dependency-order errors.
    check_scaffold(g, g.size() <= 9);
    auto marked = g;
    bool changed = false;
    for (std::size_t p = 0; p < marked.size(); ++p)
      if (marked[p] == '.') { marked[p] = p % 2 ? 'X' : 'x'; changed = true; }
    if (changed) check_scaffold(marked, marked.size() <= 9);
    return;
  }
  for (char c : {'.', '(', ')'}) {
    const int next = balance + (c == '(') - (c == ')');
    if (next < 0 || next > remaining - 1) continue;
    g.push_back(c); enumerate(g, next, remaining - 1); g.pop_back();
  }
}

void expect_violation(const D &d, const std::string &child_text) {
  require(!engine::legacy_span_schedule_supports(d), "invalid dependency passed observation predicate");
  try { engine::require_legacy_span_schedule_support(d); }
  catch (const std::logic_error &e) {
    const std::string message(e.what());
    for (const auto &part : {std::string("internal scfg2 schedule invariant"),
                             std::string(grammar::rule_name(d.rule)),
                             std::string("parent="), std::string("runtime_child[0]="),
                             child_text, std::string("split=(k=2,l=6,p=-1,q=-1)"),
                             std::string("expected runtime child rank < parent rank")})
      require(message.find(part) != std::string::npos, "incomplete dependency diagnostic: " + message);
    return;
  }
  require(false, "invalid dependency was silently accepted");
}

void check_diagnostics() {
  const Key parent{NT::V, 1, 7, -1, -1};
  auto d = deductions::expand_basic_rule(parent, R::V_INTERNAL, grammar::RuleSplit{2, 6});
  engine::require_legacy_span_schedule_support(d);  // Strictly smaller child.
  d.children[0].item = parent;
  expect_violation(d, "V(1,7;-1,-1)");  // Self-dependency.
  d.children[0].item = Key{NT::WI, 1, 7, -1, -1};
  expect_violation(d, "WI(1,7;-1,-1)");  // Same span, later phase.
  d.children[0].item = Key{NT::V, 1, 8, -1, -1};
  expect_violation(d, "V(1,8;-1,-1)");  // Wider child.
  engine::require_legacy_span_schedule_support(
      deductions::expand_basic_rule(parent, R::V_HAIRPIN, {}));
  engine::require_legacy_span_schedule_support(
      deductions::expand_basic_rule(parent, R::V_VM, {}));  // Earlier same-span phase.
  const auto direct = deductions::expand_basic_rule(Key{NT::WMBP, 1, 7, -1, -1}, R::WMBP_DIRECT_VP, {});
  require(engine::runtime_child_item_for_support(direct, 0).nonterminal == NT::VP,
          "dependency check ignored canonical VP child");
  engine::require_legacy_span_schedule_support(direct);
  d.child_count = d.children.size() + 1;
  require(!engine::legacy_span_schedule_supports(d), "invalid child count passed observation predicate");
  bool invalid_count = false;
  try { engine::require_legacy_span_schedule_support(d); }
  catch (const std::logic_error &e) {
    invalid_count = std::string(e.what()).find("invalid child count") != std::string::npos;
  }
  require(invalid_count, "invalid child count diagnostic changed");
}

struct ZeroWeight {
  bool called = false;
  weights::LogWeight local_weight_for(const D &) {
    called = true;
    return {weights::log_zero};
  }
};

template <class Semiring>
void check_executor_violation(bool caller_allows) {
  const Key parent{NT::V, 1, 7, -1, -1};
  auto corrupt = deductions::expand_basic_rule(parent, R::V_INTERNAL, grammar::RuleSplit{2, 6});
  corrupt.children[0].item = parent;
  auto provider = [&](const Key &) { return std::vector<D>{corrupt}; };
  bool allows_called = false, contribution_called = false, committed = false;
  auto allows = [&](const D &) { allows_called = true; return caller_allows; };
  auto observe = [](const auto &, const auto &, const auto &, const auto &, auto) {};
  auto before = [](const auto &, const auto &) {};
  auto after = [&](const auto &, const auto &) { committed = true; };
  auto contribution = [&](const auto &, const auto &, auto) {
    contribution_called = true;
    return Semiring::one();
  };
  ZeroWeight weight;
  engine::Chart<typename Semiring::value_type> chart(Semiring::zero());
  bool threw = false;
  try {
    replay::update_w_final_exact_dp_item_with_provider_and_contribution<Semiring>(
        parent, chart, weight, provider, allows, observe, before, after, contribution);
  } catch (const std::logic_error &e) {
    threw = std::string(e.what()).find("V_INTERNAL") != std::string::npos;
  }
  require(threw && !allows_called && !weight.called && !contribution_called &&
              !committed && !chart.contains(parent),
          "executor hid a dependency bug behind caller admission or zero weight");
}

void check_real(std::string seq, std::string g) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
  for (const auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible}) {
    replay::WFinalLogWeightContext ctx(owner, policy);
    check_generated(static_cast<cand_pos_t>(seq.size()), ctx);
  }
}
}  // namespace

int main() {
  check_diagnostics();
  for (bool caller_allows : {false, true}) {
    check_executor_violation<semirings::LogSumExpSemiring>(caller_allows);
    check_executor_violation<semirings::LogMaxSemiring>(caller_allows);
  }
  for (int n = 1; n <= 13; ++n) { std::string g; enumerate(g, 0, n); }
  std::mt19937 rng(20261013);
  for (int n : {20, 32, 50}) for (int repeat = 0; repeat < 20; ++repeat) {
    std::string g(n, '.');
    std::vector<int> stack;
    for (int p = 0; p < n; ++p) {
      if (!stack.empty() && p - stack.back() >= 4 && rng() % 3 == 0) {
        g[p] = ')'; stack.pop_back();
      } else if (rng() % 3 == 0) { g[p] = '('; stack.push_back(p); }
      else if (rng() % 5 == 0) g[p] = rng() % 2 ? 'x' : 'X';
    }
    for (int p : stack) g[p] = '.';
    check_scaffold(g);
  }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_DirksPierce09.par", 1);
  check_real("GAGAAACAAGAAGAAAAAAAAACCAAAAC", "(.(...)..(............)).....");
  check_real("GGGAAAAACGGAAAAACAAAACCAC", "(.........(.....)....)...");
  check_real("UUAUAUGGGUAUUAUUCUCA", "((x...))((...))(...)");
  check_real("GGAAACGAAACAC", "((...)(...).)");
  check_real("GAGAAACGAAACC", "(.(...)(...))");
  require(!covered_rules.count(R::VM_SCALE2) && !covered_rules.count(R::W_ROOT_WMB_DECOMPOSITION),
          "inactive exact rule unexpectedly generated");
  require(covered_rules.size() == 55, "dependency audit missed an effective exact rule");
  std::cout << "scfg2_schedule_invariant=ok scaffolds=" << scaffolds << " items=" << items
            << " candidates=" << candidates << " covered_rules=" << covered_rules.size() << '\n';
}
