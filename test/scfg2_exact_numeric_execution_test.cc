#include "scfg2/replay/w_final_exact_outside.hh"

#include <cstring>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>

namespace {
using namespace cparty::scfg2;
using Key = deductions::ItemKey;
using D = deductions::Deduction;
using NT = grammar::NonTerminal;
using S = semirings::LogSumExpSemiring;
using Options = replay::ExactExecutionOptions;
std::size_t unused_live_vm = 0, compared = 0;

void require(bool ok, const char *message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(a)) == 0; }
bool equal(const Key &a, const Key &b) { return deductions::ItemKeyEq{}(a, b); }
bool equal(const D &a, const D &b) {
  if (!equal(a.parent, b.parent) || a.rule != b.rule || a.child_count != b.child_count ||
      a.split.k != b.split.k || a.split.l != b.split.l ||
      a.split.p != b.split.p || a.split.q != b.split.q) return false;
  for (std::size_t c = 0; c < a.children.size(); ++c)
    if (!equal(a.children[c].item, b.children[c].item)) return false;
  return true;
}

struct CheckpointValue {
  int value = 0;
  inline static int calls = 0;
  CheckpointValue checkpointed() const { ++calls; return {value + 1}; }
  bool discard_chart_entry() const { return value < 0; }
};

void storage_contract() {
  const auto discard = +[](const double &v) { return S::is_zero(v); };
  engine::Chart<double> sparse(S::zero(), discard), dense(S::zero(), discard, 5);
  std::vector<Key> keys;
  for (int f = 0; f < static_cast<int>(NT::BE); ++f)
    for (int j = 1; j <= 5; ++j)
      for (int i = 1; i <= j; ++i) keys.push_back({NT(f), i, j, -1, -1});
  for (const auto &key : std::vector<Key>{{NT::BE, 1, 5, 2, 4},
           {NT::BE, 1, 5, 1, 4}, {NT::W, 1, 6, -1, -1},
           {NT::V, 1, 5, 2, 4}, {NT::WI, 0, 3, -1, -1}}) keys.push_back(key);
  int value = 0;
  for (const auto &key : keys) {
    require(!dense.contains(key) && same(dense.get(key), S::zero()), "dense missing value");
    sparse.set(key, ++value); dense.set(key, value);
  }
  const auto &view = dense.values();
  auto copy = dense;
  auto moved = std::move(copy);
  require(copy.size() == 0 && !copy.contains(keys.front()), "moved chart is not empty");
  copy.set(keys.front(), 13.0);
  require(copy.get(keys.front()) == 13.0, "moved chart cannot be reused");
  auto assigned = dense;
  assigned = std::move(copy);
  require(assigned.size() == 1 && copy.size() == 0, "chart move assignment failed");
  for (const auto &key : keys) {
    require(same(sparse.get(key), moved.get(key)), "dense key alias or copy failure");
    dense.set(key, 0.0);
    require(same(view.at(key), 0.0) && dense.contains(key), "dense live map view");
    dense.set(key, S::zero());
    require(!dense.contains(key) && !view.count(key), "dense discard failed");
  }
  require(dense.size() == 0 && moved.size() == keys.size(), "dense count or copy ownership");
  engine::Chart<double> zero(S::zero(), nullptr, 3);
  const Key k{NT::W, 1, 3, -1, -1};
  zero.set(k, S::zero());
  require(zero.contains(k) && zero.size() == 1, "explicit zero lost presence");
  std::size_t entries = 0;
  moved.for_each_entry([&](const Key &key, double v) {
    ++entries; require(same(sparse.get(key), v), "dense iteration changed values");
  });
  require(entries == keys.size(), "dense iteration count");
  engine::Chart<CheckpointValue> boundary({}, nullptr, 3);
  boundary.values();
  CheckpointValue::calls = 0;
  boundary.set(k, {7});
  require(CheckpointValue::calls == 1 && boundary.get(k).value == 8 &&
      boundary.values().at(k).value == 8, "checkpoint executed twice");
  boundary.set(k, {-1});
  require(!boundary.contains(k) && boundary.values().empty(), "boundary discard failed");
}

void product_contract() {
  const Key parent{NT::W, 1, 9, -1, -1};
  D d; d.parent = parent; d.child_count = 3;
  d.children[0].item = {NT::W, 1, 3, -1, -1};
  d.children[1].item = {NT::W, 4, 6, -1, -1};
  d.children[2].item = {NT::W, 7, 8, -1, -1};
  struct Weight {
    int calls = 0;
    weights::LogWeight local_weight_for(const D &) { ++calls; return weights::LogWeight{1.25}; }
  } weight;
  auto provider = [&](const Key &) { return std::vector<D>{d}; };
  auto allows = [](const D &) { return true; };
  int observations = 0, before = 0, after = 0;
  auto observe = [&](const auto &...) { ++observations; };
  auto before_item = [&](const auto &...) { ++before; };
  auto after_item = [&](const auto &...) { ++after; };
  replay::RuntimeProductContribution<S, true> fast;
  replay::RuntimeProductContribution<S> reference;
  engine::Chart<double> chart(S::zero(), nullptr, 9);
  auto update = [&](auto &contribution) {
    replay::update_w_final_exact_dp_item_with_provider_and_contribution<S>(
        parent, chart, weight, provider, allows, observe, before_item, after_item, contribution);
  };
  update(fast);
  require(weight.calls == 0 && observations == 0 && before == 1 && after == 1 &&
      chart.contains(parent) && S::is_zero(chart.get(parent)), "dead product was evaluated");
  update(reference);
  require(weight.calls == 1 && observations == 1, "reference zero observations changed");
  auto custom = [](const auto &, const auto &, auto) { return 17.0; };
  update(custom);
  require(chart.get(parent) == 17.0 && observations == 2, "custom product short-circuited");
  chart.set(d.children[0].item, 1e16);
  chart.set(d.children[1].item, -1e16);
  chart.set(d.children[2].item, 1e-16);
  update(reference); const double old = chart.get(parent);
  update(fast);
  require(same(old, chart.get(parent)), "live arithmetic order changed");
  chart.set(d.children[0].item, S::zero());
  // A dead first child must not conceal a later dependency-order violation.
  d.children[2].item = parent;
  bool threw = false;
  try { update(fast); } catch (const std::logic_error &) { threw = true; }
  require(threw, "dead product concealed schedule invariant");

  // Two distinct live alternatives with equal scores must keep the first
  // winner in both paths, rather than relying on incidental ties in fixtures.
  using M = semirings::LogMaxSemiring;
  D first = d; first.child_count = 1; first.split.k = 11;
  D second = first; second.split.k = 12;
  chart.set(first.children[0].item, 0.0);
  const auto check_tie = [&](auto product) {
    engine::ChoiceMap<double> choices;
    int events = 0;
    auto tied_provider = [&](const Key &) { return std::vector<D>{first, second}; };
    auto choose = [&](const Key &item, const D &candidate, double contribution, const auto &, auto) {
      ++events;
      choices.observe(item, candidate, contribution, [](double a, double b) { return a > b; });
    };
    replay::update_w_final_exact_dp_item_with_provider_and_contribution<M>(parent,
        chart, weight, tied_provider, allows, choose, before_item, after_item, product);
    require(events == 2 && choices.find(parent)->deduction.split.k == 11,
        "equal-score first-winner tie changed");
  };
  check_tie(replay::RuntimeProductContribution<M>{});
  check_tie(replay::RuntimeProductContribution<M, true>{});

  // A zero inside product can have a nonzero derivative with respect to its
  // zero child. The forward optimization must never suppress this adjoint.
  d.child_count = 2;
  chart.set(d.children[0].item, S::zero());
  chart.set(d.children[1].item, 4.0);
  engine::Chart<double> outside(S::zero(), nullptr, 9);
  engine::propagate_outside_deduction<S>(d, chart, 0.0, &outside,
      [](const D &) { return weights::LogWeight{2.0}; });
  require(outside.get(d.children[0].item) == 6.0 &&
      S::is_zero(outside.get(d.children[1].item)), "zero-child adjoint was lost");
}

void check(std::string seq, std::string g, WmbClosurePolicy policy) {
  weights::RuntimeSnapshot owner(seq, g);
  replay::WFinalLogWeightContext ctx(owner, policy);
  const auto n = static_cast<cand_pos_t>(seq.size());
  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  const auto allows = [](const auto &) { return true; };
  const auto ignore = [](const auto &...) {};
  const auto reference = replay::run_w_final_exact_dp_schedule<S>(schedule, ctx, allows, ignore);
  const Key root{NT::W, 1, n, -1, -1};
  const auto beta_ref = replay::run_w_final_exact_regenerated_outside(schedule, ctx, reference, root);
  const auto max_ref = replay::run_w_final_exact_viterbi_schedule_with_context<semirings::LogMaxSemiring>(
      schedule, ctx, allows);
  for (const auto &key : schedule) {
    if (key.nonterminal == NT::VM && !deductions::exact_v_family_active(ctx, key.i, key.j) &&
        !S::is_zero(reference.get(key))) ++unused_live_vm;
    for (const auto &d : deductions::generate_exact_basic_deductions_for_item(key, ctx))
      for (std::size_t c = 0; c < d.child_count; ++c)
        if (d.children[c].item.nonterminal == NT::VM)
          require(d.rule == grammar::RuleId::V_VM &&
              d.parent.i == d.children[c].item.i && d.parent.j == d.children[c].item.j,
              "VM gained a consumer outside V_VM: re-audit omission policy");
  }
  for (const Options options : {Options{true, false, false}, Options{false, true, false},
           Options{false, false, true}, Options::numeric()}) {
    const auto chart = replay::run_w_final_exact_dp_schedule<S>(schedule, ctx, allows, ignore, options);
    const auto maximum = replay::run_w_final_exact_viterbi_schedule_with_context<semirings::LogMaxSemiring>(
        schedule, ctx, allows, options);
    for (const auto &key : schedule) {
      require(chart.contains(key) == reference.contains(key), "optimized chart presence changed");
      if (options.omit_unused_vm && key.nonterminal == NT::VM &&
          !deductions::exact_v_family_active(ctx, key.i, key.j)) {
        require(S::is_zero(chart.get(key)) && S::is_zero(maximum.chart.get(key)), "unused VM not zero");
        continue;
      }
      require(same(reference.get(key), chart.get(key)), "exact sum value bits changed");
      require(same(max_ref.chart.get(key), maximum.chart.get(key)), "exact max value bits changed");
      const auto *a = max_ref.choices.find(key), *b = maximum.choices.find(key);
      require(bool(a) == bool(b), "exact choice presence changed");
      if (a) require(same(a->contribution, b->contribution) && equal(a->deduction, b->deduction),
          "Viterbi choice/tie changed");
      ++compared;
    }
    const auto beta = replay::run_w_final_exact_regenerated_outside(schedule, ctx, chart, root);
    for (const auto &key : schedule)
      require(same(beta_ref.get(key), beta.get(key)), "outside adjoint bits changed");
  }
}
}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  storage_contract(); product_contract();
  std::vector<std::pair<std::string, std::string>> cases{
      {"GGGGGAACCAACCC", "((.....))....."},
      {"UUAUAUGGGUAUUAUUCUCA", "((....))((...))(...)"},
      {"GCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAU", "(......)................................"},
      {"GCAUGCAUGCAUGCAU", "................"},
      {"GCGC", "(())"}, {"AAAA", "xxxx"}};
  std::mt19937 rng(20261020);
  for (int c = 0; c < 12; ++c) {
    std::string seq(28, 'A'), g(28, '.');
    for (char &base : seq) base = "ACGU"[rng() % 4];
    g[2] = '('; g[12] = ')'; seq[2] = 'G'; seq[12] = 'C';
    if (c % 2) { g[17] = '('; g[26] = ')'; seq[17] = 'A'; seq[26] = 'U'; }
    cases.emplace_back(seq, g);
  }
  for (const auto &[seq, g] : cases) {
    require(seq.size() == g.size(), "bad execution witness length");
    for (const auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible})
      check(seq, g, policy);
  }
  require(unused_live_vm > 0, "VM omission witness did not exercise positive unused VM");
  std::cout << "scfg2_exact_numeric_execution=ok cases=" << cases.size() * 2
      << " compared_states=" << compared << " unused_live_vm=" << unused_live_vm << '\n';
}
