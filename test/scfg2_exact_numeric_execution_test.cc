#include "scfg2/replay/w_final_exact_outside.hh"

#include <cstring>
#include <cstdint>
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
  auto optional_provider = [&](const Key &, const auto &...extra) {
    require(sizeof...(extra) == 0, "legacy variadic provider received an extra chart argument");
    return std::vector<D>{d};
  };
  replay::update_w_final_exact_dp_item_with_provider_and_contribution<S>(
      parent, chart, weight, optional_provider, allows, observe, before_item, after_item, reference);
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

void selected_numeric_streams(WmbClosurePolicy policy) {
  std::string seq(40, 'G'), g(40, '.');
  for (std::size_t i = 0; i < seq.size(); ++i) seq[i] = "GCAU"[i % 4];
  g[0] = '('; g[7] = ')';
  weights::RuntimeSnapshot owner(seq, g);
  replay::WFinalLogWeightContext ctx(owner, policy);
  engine::Chart<double> chart(S::zero(), nullptr, 40);
  std::mt19937 rng(20261021);
  for (const auto nt : {NT::V, NT::WMB, NT::WM, NT::WMv, NT::WMp, NT::WI, NT::WIP, NT::VP, NT::BE, NT::WMBP, NT::WMBW})
    for (int j = 1; j <= 40; ++j)
      for (int i = 1; i <= j; ++i)
        if (rng() % 3 == 0) chart.set({nt, i, j, -1, -1}, -double(rng() % 100));
  // A finite but arbitrarily small mass is not a zero. The following slots
  // also distinguish START from SPLIT and VM BASE from a binary left prefix.
  chart.set({NT::V, 10, 35, -1, -1}, -1e300);
  chart.set({NT::WMp, 11, 34, -1, -1}, -1e300);
  chart.set({NT::VP, 2, 35, -1, -1}, 0.0);
  chart.set({NT::VP, 3, 35, -1, -1}, S::zero());
  replay::ExactNumericSplitScratch scratch;
  std::size_t kept = 0, omitted = 0;
  std::array<std::size_t, 3> wmbp_kept{}, wmbp_omitted{};
  for (const auto nt : {NT::WM, NT::VM, NT::WI, NT::WIP, NT::VPR, NT::WMBP}) {
    const auto kept_before = kept, omitted_before = omitted;
    const int begin = nt == NT::WMBP ? 1 : nt == NT::VPR ? 2 : 9;
    for (int i = begin; i < begin + 6; ++i) {
      const Key parent{nt, i, 35, -1, -1};
      const auto reference = deductions::generate_exact_basic_deductions_for_item(parent, ctx);
      scratch.prepare<S>(parent, chart, ctx);
      const auto views = scratch.views();
      std::vector<D> expected;
      for (const auto &d : reference) {
        if (nt == NT::WMBP && views.for_rule(d.rule, nullptr)) {
          std::size_t vp_children = 0;
          for (std::size_t c = 0; c < d.child_count; ++c) {
            const auto &child = d.children[c].item;
            if (child.nonterminal != NT::VP) continue;
            ++vp_children;
            require(equal(child, deductions::exact_wmbp_vp_child(parent,
                d.rule == grammar::RuleId::WMBP_SPLIT_BE_WI_VP ? d.split.p : d.split.k)) &&
                child.i > parent.i, "WMBP VP-child selection contract changed");
          }
          require(vp_children == 1, "WMBP family lost its unique VP child");
        }
        bool live = true;
        if (views.for_rule(d.rule, nullptr))
          for (std::size_t c = 0; c < d.child_count; ++c)
            if (nt != NT::WMBP || d.children[c].item.nonterminal == NT::VP)
              live &= !S::is_zero(chart.get(d.children[c].item));
        if (nt == NT::WMBP && views.for_rule(d.rule, nullptr)) {
          const std::size_t family = d.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP ? 0 :
              d.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP ? 1 : 2;
          ++(live ? wmbp_kept : wmbp_omitted)[family];
        }
        if (live) expected.push_back(d); else ++omitted;
      }
      const auto selected = deductions::generate_exact_basic_deductions_for_item(parent, ctx, nullptr, &views);
      require(selected.size() == expected.size(), "numeric split selection count changed");
      for (std::size_t c = 0; c < expected.size(); ++c) {
        require(equal(expected[c], selected[c]), "numeric selected child fields/order changed");
        ++kept;
      }
    }
    require(kept > kept_before && omitted > omitted_before,
        "numeric split witness missed live/dead alternatives in a selected family");
  }
  for (std::size_t family = 0; family < wmbp_kept.size(); ++family)
    require(wmbp_kept[family] && wmbp_omitted[family],
        "WMBP selection witness missed a live/dead shared or WI rule block");
  require(kept > 0 && omitted > 0, "numeric split selection witness has no live/dead branches");
  const Key wm{NT::WM, 10, 35, -1, -1}, vm{NT::VM, 10, 35, -1, -1};
  scratch.prepare<S>(wm, chart, ctx);
  require(std::find(scratch.selected[0].begin(), scratch.selected[0].end(), 10) != scratch.selected[0].end(),
      "finite tiny V mass omitted from START");
  scratch.prepare<S>(vm, chart, ctx);
  require(std::find(scratch.selected[6].begin(), scratch.selected[6].end(), 11) != scratch.selected[6].end(),
      "finite tiny WMp mass omitted from BASE");
  require(std::find(scratch.selected[4].begin(), scratch.selected[4].end(), 11) == scratch.selected[4].end(),
      "empty VM binary prefix reopened");
}

void ordered_generator_reference() {
  const std::string seq = "GCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAUGCAU";
  const std::string g = "(......)................................";
  for (const auto policy : {WmbClosurePolicy::LegacyCompatible, WmbClosurePolicy::RequireWeakClosure}) {
    weights::RuntimeSnapshot owner(seq, g);
    replay::WFinalLogWeightContext ctx(owner, policy);
    const auto schedule = engine::build_pkfree_basic_exact_schedule(seq.size(), ctx);
    std::uint64_t hash = 14695981039346656037ull;
    const auto add = [&](auto value) { hash = (hash ^ std::uint64_t(value)) * 1099511628211ull; };
    const auto key = [&](const Key &k) { add(k.nonterminal); add(k.i); add(k.j); add(k.ip); add(k.jp); };
    std::size_t count = 0, singleton = 0, shared = 0, wi = 0;
    for (const auto &parent : schedule) {
      key(parent);
      const auto buffered = deductions::generate_exact_basic_deductions_for_item(parent, ctx);
      for (const D &d : buffered) {
        ++count; key(d.parent); add(d.rule); add(d.split.k); add(d.split.l); add(d.split.p); add(d.split.q);
        add(d.child_count);
        for (const auto &child : d.children) key(child.item);
        singleton += parent.nonterminal == NT::W && parent.i == parent.j;
        shared += d.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP ||
                  d.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP;
        wi += d.rule == grammar::RuleId::WMBP_SPLIT_BE_WI_VP;
      }
    }
    // Frozen from a272d73, before item generation used the shared append
    // body. This detects changes in the complete rule/block/field order.
    const bool legacy = policy == WmbClosurePolicy::LegacyCompatible;
    require(count == (legacy ? 120285u : 119913u) &&
        hash == (legacy ? 0xdce54ca6063761aeull : 0x29d3540a4c09a458ull),
        "ordered deduction fingerprint differs from frozen reference");
    require(singleton && shared && wi, "stream witness missed singleton or shared/wi rule blocks");
  }
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
           Options{false, false, true}, Options{false, false, false, true},
           Options{false, false, false, false, true}, Options::numeric()}) {
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
  storage_contract(); product_contract(); ordered_generator_reference();
  for (const auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible})
    selected_numeric_streams(policy);
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
