#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/be_core_left_retained_events.hh"
#include "scfg2/beam/be_core_left_prefix_agenda.hh"
#include "scfg2/beam/be_core_left_rule_index.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/engine/outside.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

bool close(double lhs, double rhs) {
  const double scale = std::max({1.0, std::fabs(lhs), std::fabs(rhs)});
  return std::fabs(lhs - rhs) <= 1e-10 * scale;
}

struct Counts {
  std::array<std::size_t, 2> sources{};
  std::array<std::size_t, 2> emitted{};
  std::size_t child_events = 0;
  std::size_t ancestor_probes = 0;
  std::size_t max_ancestors_per_child = 0;
  std::size_t prefix_internal_probes = 0;
  std::size_t prefix_left_probes = 0;
  std::size_t max_left_prefixes_per_core = 0;

  bool operator==(const Counts &other) const {
    return sources == other.sources && emitted == other.emitted &&
        child_events == other.child_events &&
        ancestor_probes == other.ancestor_probes &&
        max_ancestors_per_child == other.max_ancestors_per_child &&
        prefix_internal_probes == other.prefix_internal_probes &&
        prefix_left_probes == other.prefix_left_probes &&
        max_left_prefixes_per_core == other.max_left_prefixes_per_core;
  }
};

template <class Semiring>
Counts check_events(cparty::scfg2::weights::RuntimeSnapshot &owner) {
  using namespace cparty::scfg2;
  replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const auto schedule =
      engine::build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
  const auto chart = replay::run_w_final_exact_dp_schedule<Semiring>(
      schedule,
      ctx,
      [](const Deduction &) { return true; },
      [](const ItemKey &, const Deduction &,
         typename Semiring::value_type, const auto &, double) {});
  engine::AllLocalWeightModel<replay::WFinalExactInsideContext> weight(&ctx);

  using Totals = std::unordered_map<ItemKey,
                                    typename Semiring::value_type,
                                    deductions::ItemKeyHash,
                                    deductions::ItemKeyEq>;
  using SourceId = std::tuple<grammar::RuleId,
                              cand_pos_t,
                              cand_pos_t,
                              cand_pos_t,
                              cand_pos_t,
                              cand_pos_t,
                              cand_pos_t>;
  Totals driven;
  Totals prefix_driven;
  std::set<SourceId> driven_ids;
  std::set<SourceId> prefix_driven_ids;
  auto events = beam::BeCoreLeftRetainedEvents::build(n, ctx);
  auto prefix_agenda = beam::BeCoreLeftPrefixAgenda::build(n, ctx);
  auto rule_index = beam::BeCoreLeftRuleIndex::build(n, ctx);
  Totals regenerated;
  std::set<SourceId> regenerated_ids;
  std::size_t retained_core_count = 0, retained_wip_count = 0;
  const auto consume = [&](const auto &source,
                           Totals &totals,
                           std::set<SourceId> &ids) {
    const auto &parent = source.parent;
    const auto term = replay::runtime_deduction_contribution<Semiring>(
        chart, source, weight.weight_for(source));
    if (term == Semiring::zero()) return;
    assert(ids.emplace(source.rule,
                       parent.i,
                       parent.j,
                       parent.ip,
                       parent.jp,
                       source.split.k,
                       source.split.l)
               .second);
    auto [found, inserted] = totals.emplace(parent, Semiring::zero());
    (void)inserted;
    found->second = Semiring::plus(found->second, term);
  };
  for (cand_pos_t j = 1; j <= n; ++j) {
    events.prepare_column(j);
    prefix_agenda.prepare_column(j, ctx);
    std::vector<ItemKey> ordered_expected;
    rule_index.for_each_column_source(j, ctx, [&](const auto &source) {
      ordered_expected.push_back(source.parent);
    });
    std::sort(ordered_expected.begin(), ordered_expected.end(), [](const auto &a, const auto &b) {
      return std::tie(a.i, a.ip, a.jp) > std::tie(b.i, b.ip, b.jp);
    });
    ordered_expected.erase(std::unique(ordered_expected.begin(), ordered_expected.end(),
                                      deductions::ItemKeyEq{}), ordered_expected.end());
    const auto ordered_actual = prefix_agenda.parents();
    assert(ordered_actual.size() == ordered_expected.size());
    assert(std::equal(ordered_actual.begin(), ordered_actual.end(), ordered_expected.begin(),
                      deductions::ItemKeyEq{}));
    assert(prefix_agenda.stats().unordered_parent_columns == 0);
    for (const auto &parent : events.parents()) {
      for (const auto &source : events.deductions_for(parent)) {
        consume(source, driven, driven_ids);
      }
    }
    for (const auto &parent : prefix_agenda.parents()) {
      for (const auto &source : prefix_agenda.deductions_for(parent)) {
        consume(source, prefix_driven, prefix_driven_ids);
      }
    }
    std::vector<ItemKey> retained;
    for (const auto &item : schedule) {
      if (item.nonterminal == grammar::NonTerminal::BE && item.j == j &&
          chart.get(item) != Semiring::zero()) {
        retained.push_back(item);
      }
    }
    events.observe_committed_column(retained, ctx);
    prefix_agenda.observe_committed_be_column(retained, ctx);
    rule_index.observe_committed_be_column(retained);
    retained_core_count += retained.size();
    std::vector<cand_pos_t> retained_wip_starts;
    for (const auto &item : schedule) {
      if (item.nonterminal == grammar::NonTerminal::WIP && item.j == j &&
          chart.get(item) != Semiring::zero()) {
        retained_wip_starts.push_back(item.i);
      }
    }
    prefix_agenda.observe_retained_wip_column(
        retained_wip_starts, j, ctx);
    rule_index.observe_retained_wip_column(retained_wip_starts, j, ctx);
    retained_wip_count += retained_wip_starts.size();
  }
  // Regenerate in reverse column order from final retained identities. The
  // index stores no pending full deductions or numerical child snapshots.
  for (cand_pos_t j = n; j >= 1; --j) {
    rule_index.for_each_column_source(j, ctx, [&](const auto &source) {
      consume(source, regenerated, regenerated_ids);
    });
  }
  assert(rule_index.core_count() == retained_core_count);
  assert(rule_index.event_count() <=
         static_cast<std::size_t>(n) + retained_wip_count);
  assert(regenerated_ids == driven_ids);
  for (const auto &[parent, value] : driven)
    assert(close(value, regenerated.at(parent)));
  const auto rebuilt = beam::BeCoreLeftRuleIndex::from_retained_chart<Semiring>(
      n, chart, ctx);
  std::set<SourceId> rebuilt_ids;
  Totals rebuilt_totals;
  for (cand_pos_t j = n; j >= 1; --j)
    rebuilt.for_each_column_source(j, ctx, [&](const auto &source) {
      consume(source, rebuilt_totals, rebuilt_ids);
    });
  assert(rebuilt_ids == driven_ids);
  assert(rebuilt.core_count() == rule_index.core_count());
  assert(rebuilt.event_count() == rule_index.event_count());

  // A sparse final chart, including removed anchored cores and left WIP.
  // Only surviving parents/children enter the regenerated outside surface.
  engine::Chart<double> retained(0.0);
  for (const auto &[item, value] : chart.values()) {
    if (value > 0 && deductions::ItemKeyHash{}(item) % 3 != 0)
      retained.set(item, value);
  }
  const auto sparse =
      beam::BeCoreLeftRuleIndex::from_retained_chart<semirings::SumProductSemiring>(
          n, retained, ctx);
  std::set<SourceId> sparse_ids, sparse_expected;
  Totals sparse_totals, sparse_reference;
  const auto admitted = [&](const auto &source) {
    if (!retained.contains(source.parent)) return false;
    for (std::size_t c = 0; c < source.child_count; ++c)
      if (!retained.contains(source.children[c].item)) return false;
    return true;
  };
  for (cand_pos_t j = n; j >= 1; --j)
    sparse.for_each_column_source(j, ctx, [&](const auto &source) {
      if (admitted(source)) consume(source, sparse_totals, sparse_ids);
    });
  for (const auto &[parent, value] : retained.values()) {
    (void)value;
    if (parent.nonterminal != grammar::NonTerminal::BE) continue;
    for (const auto rule : {grammar::RuleId::BE_INTERNAL_LOOP,
                            grammar::RuleId::BE_WIP_BASEPAIR})
      for (const auto &source : deductions::generate_exact_basic_deductions(
               parent, rule, ctx))
        if (engine::legacy_span_schedule_supports(source) && admitted(source))
          consume(source, sparse_reference, sparse_expected);
  }
  assert(sparse_ids == sparse_expected);

  // The same retained index drives an outside consumer. Child values are
  // independent tape leaves in this local objective, so its derivatives
  // provide an arithmetic oracle independent of the rule-level reverse code.
  using L = semirings::LogSumExpSemiring;
  using T = semirings::TapeLogSumExpSemiring;
  engine::Chart<double> alpha(L::zero()), beta(L::zero());
  semirings::LogSumExpTape tape;
  T::Scope scope(tape);
  engine::Chart<T::value_type> leaves(T::zero());
  for (const auto &[item, value] : chart.values()) {
    if (value <= 0 || (item.nonterminal != grammar::NonTerminal::BE &&
                      item.nonterminal != grammar::NonTerminal::WIP)) continue;
    alpha.set(item, std::log(value));
    leaves.set(item, T::from_weight(value));
  }
  auto objective = T::zero();
  for (cand_pos_t j = n; j >= 1; --j) {
    rule_index.for_each_column_source(j, ctx, [&](const auto &source) {
      objective = T::plus(objective, replay::runtime_deduction_contribution<T>(
          leaves, source, weight.weight_for(source)));
      engine::propagate_outside_deduction<L>(
          source, alpha, L::one(), &beta,
          [&](const auto &d) { return weight.weight_for(d); });
    });
  }
  if (!T::is_zero(objective)) {
    const auto posterior = tape.reverse(objective.node);
    for (const auto &[item, value] : leaves.values()) {
      const auto expected = posterior[value.node];
      const auto actual = std::exp(beta.get(item) + alpha.get(item) -
                                   objective.log_value);
      assert(close(actual, expected));
    }
  }
  assert(events.pending_empty());
  assert(prefix_agenda.pending_empty());
  assert(events.stats().duplicate_sources == 0);
  assert(prefix_agenda.stats().duplicate_sources == 0);
  assert(prefix_driven_ids == driven_ids);
  assert(prefix_driven.size() == driven.size());
  for (const auto &[parent, value] : driven) {
    const auto found = prefix_driven.find(parent);
    assert(found != prefix_driven.end());
    assert(close(value, found->second));
  }

  const std::array<grammar::RuleId, 2> rules{
      grammar::RuleId::BE_INTERNAL_LOOP,
      grammar::RuleId::BE_WIP_BASEPAIR};
  std::set<SourceId> exact_ids;
  std::array<std::size_t, 2> counts{};
  for (const auto &parent : schedule) {
    if (parent.nonterminal != grammar::NonTerminal::BE) continue;
    auto aggregate = Semiring::zero();
    for (std::size_t idx = 0; idx < rules.size(); ++idx) {
      for (const auto &source : deductions::generate_exact_basic_deductions(
               parent, rules[idx], ctx)) {
        if (!engine::legacy_span_schedule_supports(source)) continue;
        const double local = weight.weight_for(source);
        if (local == 0.0) continue;
        const auto term = replay::runtime_deduction_contribution<Semiring>(
            chart, source, local);
        if (term == Semiring::zero()) continue;
        aggregate = Semiring::plus(aggregate, term);
        exact_ids.emplace(source.rule,
                          parent.i,
                          parent.j,
                          parent.ip,
                          parent.jp,
                          source.split.k,
                          source.split.l);
        ++counts[idx];
      }
    }
    const auto found = driven.find(parent);
    const auto actual = found == driven.end()
        ? Semiring::zero()
        : found->second;
    if (!close(aggregate, actual)) {
      std::cerr << "BE core/left mismatch parent=" << parent.i << ','
                << parent.j << ',' << parent.ip << ',' << parent.jp
                << " exact=" << aggregate << " driven=" << actual << '\n';
      assert(false);
    }
  }
  assert(driven_ids == exact_ids);
  return Counts{counts,
                events.stats().emitted_sources,
                events.stats().child_events,
                events.stats().ancestor_probes,
                events.stats().max_ancestors_per_child,
                prefix_agenda.stats().internal_parent_probes,
                prefix_agenda.stats().left_parent_probes,
                prefix_agenda.stats().max_left_prefixes_per_core};
}

}  // namespace

int main(int argc, char **argv) {
  assert(argc == 1 || argc == 3 || argc == 4);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 4 ? argv[3] : "params/rna_Turner04.par",
         1);
  std::string sequence = argc >= 3
      ? argv[1]
      : "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  std::string scaffold = argc >= 3
      ? argv[2]
      : "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  const auto inside =
      check_events<cparty::scfg2::semirings::SumProductSemiring>(owner);
  const auto viterbi =
      check_events<cparty::scfg2::semirings::MaxProductSemiring>(owner);
  assert(inside == viterbi);
  std::string siblings = "(";
  for (int i = 0; i < 12; ++i) siblings += "(....)..";
  siblings += ")";
  for (auto g : {siblings, std::string(".(..(....)..(..(....)..(....)..)..(....)..)..")}) {
    std::string seq(g.size(), 'A');
    for (std::size_t i = 0; i < g.size(); ++i)
      seq[i] = g[i] == '(' ? 'G' : g[i] == ')' ? 'C' : (i % 2 ? 'U' : 'A');
    cparty::scfg2::weights::RuntimeSnapshot extra(
      seq, g);
    (void)check_events<cparty::scfg2::semirings::SumProductSemiring>(extra);
  }
  std::cout << "be_internal\tbe_left\temitted_internal\temitted_left"
               "\tchild_events\tancestor_probes"
               "\tmax_ancestors_per_child\tprefix_internal_probes"
               "\tprefix_left_probes\tmax_left_prefixes_per_core\n"
            << inside.sources[0] << '\t' << inside.sources[1] << '\t'
            << inside.emitted[0] << '\t' << inside.emitted[1] << '\t'
            << inside.child_events << '\t' << inside.ancestor_probes << '\t'
            << inside.max_ancestors_per_child << '\t'
            << inside.prefix_internal_probes << '\t'
            << inside.prefix_left_probes << '\t'
            << inside.max_left_prefixes_per_core << '\n';
  return EXIT_SUCCESS;
}
