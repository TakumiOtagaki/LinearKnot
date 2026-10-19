#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_tail_cache.hh"
#include "scfg2/beam/vp_case3_border_event_plan.hh"
#include "scfg2/beam/vp_case3_left_cache.hh"
#include "scfg2/beam/vp_case3_right_endpoint_producer_registry.hh"
#include "scfg2/beam/vp_case3_tail_attachment_segments.hh"
#include "scfg2/beam/vp_case3_tail_border_profile.hh"
#include "scfg2/beam/vp_case3_tail_column_plan.hh"
#include "scfg2/beam/vp_case3_tail_schedule_profile.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <utility>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

struct TestLogMaxSemiring {
  using value_type = double;
  static value_type zero() {
    return -std::numeric_limits<value_type>::infinity();
  }
  static value_type one() { return 0.0; }
  static value_type from_weight(double weight) { return std::log(weight); }
  static bool is_zero(value_type value) { return value == zero(); }
  static double priority(value_type value) { return value; }
  static value_type plus(value_type lhs, value_type rhs) {
    return std::max(lhs, rhs);
  }
  static value_type times(value_type lhs, value_type rhs) {
    return lhs + rhs;
  }
  static bool better(value_type lhs, value_type rhs) { return lhs > rhs; }
};

struct KeepAllResult {
  cparty::scfg2::beam::ExactVpCase3TailCache cache;
  cparty::scfg2::beam::VpCase3TailScheduleProfile schedule_profile;
  std::size_t planned_source_count = 0;
  std::size_t attachment_segment_count = 0;
  std::size_t border_event_endpoint_pair_count = 0;
  std::size_t child_event_endpoint_pair_count = 0;
  std::size_t right_event_endpoint_pair_count = 0;
  std::size_t registry_event_endpoint_pair_count = 0;
};

using SourceIdentity =
    std::tuple<cparty::scfg2::grammar::NonTerminal,
               cand_pos_t,
               cand_pos_t,
               cand_pos_t,
               cand_pos_t,
               cand_pos_t,
               cand_pos_t>;

template <class Plan>
std::set<SourceIdentity> source_identities(const Plan &plan) {
  std::set<SourceIdentity> out;
  for (const auto &source : plan.sources()) {
    out.emplace(source.deduction.parent.nonterminal,
                source.deduction.parent.i,
                source.tail.j,
                source.deduction.split.p,
                source.tail.q,
                source.tail.k,
                source.tail.l);
  }
  return out;
}

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

void check_source_observation_boundary() {
  namespace beam = cparty::scfg2::beam;
  namespace engine = cparty::scfg2::engine;
  namespace grammar = cparty::scfg2::grammar;
  using Semiring = cparty::scfg2::semirings::SumProductSemiring;
  engine::Chart<double> chart(0.0);
  chart.set(ItemKey{grammar::NonTerminal::WI, 2, 3, -1, -1}, 2.0);
  chart.set(ItemKey{grammar::NonTerminal::WI, 6, 6, -1, -1}, 3.0);
  chart.set(ItemKey{grammar::NonTerminal::WI, 9, 9, -1, -1}, 4.0);
  const beam::VpCase3TailKey key{10, 5, 7, 8};
  const beam::VpCase3TailKey missing{11, 5, 7, 8};
  beam::ExactVpCase3TailCache cache;
  std::size_t visits = 0;
  cache.for_each_tail_key([&](const auto &) { ++visits; });
  assert(visits == 0);
  assert(cache.source_attachments(key).empty());
  assert(cache.source_attachments(missing).empty());
  assert(cache.materialization_count() == 0);
  assert(cache.attachment_count() == 0);

  assert(cache.materialize_tail<Semiring>(key, chart));
  cache.for_each_tail_key([&](const auto &observed) {
    assert(observed == key);
    ++visits;
  });
  assert(visits == 1);
  assert(cache.source_attachments(key).empty());
  assert(cache.source_attachments(missing).empty());
  assert(cache.materialization_count() == 1);
  assert(cache.attachment_count() == 0);
  assert(!cache.materialize_tail<Semiring>(key, chart));
  assert(cache.tail_count() == 1);

  Deduction source;
  source.parent = ItemKey{grammar::NonTerminal::VP, 1, 10, -1, -1};
  source.rule = grammar::RuleId::VP_WI_CASE3;
  source.split = grammar::RuleSplit{7, 8, 4, 5};
  auto absent = source;
  absent.parent.j = missing.j;
  assert(!cache.reconstruct_materialized_left_source<Semiring>(
      absent, chart, 5.0));
  assert(cache.source_attachments(missing).empty());
  assert(cache.attachment_count() == 0);
  const auto value = cache.reconstruct_materialized_left_source<Semiring>(
      source, chart, 5.0);
  assert(value.has_value() && close(*value, 120.0));
  assert(cache.source_attachments(key).size() == 1);
  assert(beam::same_attachment_source(cache.source_attachments(key)[0],
      beam::VpCase3SourceAttachment{source.parent, source.split}));
  assert(cache.attachment_count() == 1);
  assert(cache.full_source_attachments_are_unique());

  auto copied = cache;
  auto second = source;
  second.parent.i = 2;
  assert(copied.reconstruct_materialized_left_source<Semiring>(
      second, chart, 5.0).has_value());
  assert(cache.source_attachments(key).size() == 1);
  assert(copied.source_attachments(key).size() == 2);
  assert(beam::same_attachment_source(copied.source_attachments(key)[1],
      beam::VpCase3SourceAttachment{second.parent, second.split}));
  assert(copied.full_source_attachments_are_unique());
  assert(copied.reconstruct_materialized_left_source<Semiring>(
      source, chart, 5.0).has_value());
  assert(copied.source_attachments(key).size() == 3);
  assert(!copied.full_source_attachments_are_unique());
  assert(cache.full_source_attachments_are_unique());
  auto moved = std::move(copied);
  assert(moved.source_attachments(key).size() == 3);
  assert(moved.attachment_count() == 3);
  assert(moved.materialization_count() == 1);
  assert(moved.max_attachments_per_tail() == 3);
  assert(close(*moved.tail_factor(key), 12.0));
  moved = {};
  assert(moved.tail_count() == 0 && moved.attachment_count() == 0);
  assert(moved.source_attachments(key).empty());
}

void check_coupled_value_factorization() {
  namespace beam = cparty::scfg2::beam;
  namespace engine = cparty::scfg2::engine;
  namespace grammar = cparty::scfg2::grammar;
  using Semiring = cparty::scfg2::semirings::CoupledLogSumMaxSemiring;
  using Value = Semiring::value_type;

  engine::Chart<Value> chart(Semiring::zero());
  chart.set(ItemKey{grammar::NonTerminal::WI, 6, 6, -1, -1},
            Value{std::log(5.0), std::log(3.0)});
  chart.set(ItemKey{grammar::NonTerminal::WI, 9, 9, -1, -1},
            Value{std::log(7.0), std::log(4.0)});

  beam::BasicExactVpCase3TailCache<Value> tails;
  const beam::VpCase3TailKey tail{10, 5, 7, 8};
  assert(tails.materialize_tail<Semiring>(tail, chart));
  const auto tail_factor = tails.tail_factor(tail);
  assert(tail_factor.has_value());
  assert(close(tail_factor->log_sum, std::log(35.0)));
  assert(close(tail_factor->log_max, std::log(12.0)));

  beam::BasicExactVpCase3LeftCache<Value> left;
  const beam::VpCase3LeftKey left_key{tail, 1, 4};
  assert(left.materialize_scored(
             {beam::BasicVpCase3LeftRecord<Value>{
                 left_key,
                 Value{std::log(70.0), std::log(24.0)},
                 beam::vp_case3_carrier_bit(grammar::NonTerminal::VP)}},
             [](const beam::VpCase3LeftKey &) { return true; }) == 1);
  cparty::scfg2::deductions::Deduction source;
  source.parent = ItemKey{grammar::NonTerminal::VP, 1, 10, -1, -1};
  source.rule = grammar::RuleId::VP_WI_CASE3;
  source.split = grammar::RuleSplit{7, 8, 4, 5};
  const auto reconstructed = left.reconstruct_source<Semiring>(source, 2.0);
  assert(reconstructed.has_value());
  assert(close(reconstructed->log_sum, std::log(140.0)));
  assert(close(reconstructed->log_max, std::log(48.0)));
}

void check_log_domain_factorization() {
  namespace beam = cparty::scfg2::beam;
  namespace deductions = cparty::scfg2::deductions;
  namespace engine = cparty::scfg2::engine;
  namespace grammar = cparty::scfg2::grammar;

  engine::Chart<double> chart(TestLogMaxSemiring::zero());
  chart.set(ItemKey{grammar::NonTerminal::WI, 2, 3, -1, -1},
            std::log(2.0));
  chart.set(ItemKey{grammar::NonTerminal::WI, 6, 6, -1, -1},
            std::log(3.0));
  chart.set(ItemKey{grammar::NonTerminal::WI, 9, 9, -1, -1},
            std::log(4.0));

  beam::ExactVpCase3TailCache cache;
  const beam::VpCase3TailKey tail{10, 5, 7, 8};
  assert(cache.materialize_tail<TestLogMaxSemiring>(tail, chart));
  assert(close(*cache.tail_factor(tail), std::log(12.0)));

  deductions::Deduction source;
  source.parent =
      ItemKey{grammar::NonTerminal::VP, 1, 10, -1, -1};
  source.rule = grammar::RuleId::VP_WI_CASE3;
  source.split = grammar::RuleSplit{7, 8, 4, 5};
  const auto reconstructed =
      cache.reconstruct_materialized_left_source<TestLogMaxSemiring>(
          source, chart, 5.0);
  assert(reconstructed.has_value());
  assert(close(*reconstructed, std::log(120.0)));

  beam::ExactVpCase3LeftCache left_cache;
  const beam::VpCase3LeftKey left_key{tail, 1, 4};
  assert(left_cache.materialize_scored(
             {beam::VpCase3LeftRecord{
                 left_key,
                 std::log(24.0),
                 beam::vp_case3_carrier_bit(grammar::NonTerminal::VP)}},
             [](const beam::VpCase3LeftKey &) { return true; }) == 1);
  const auto promoted =
      left_cache.reconstruct_source<TestLogMaxSemiring>(source, 5.0);
  assert(promoted.has_value());
  assert(close(*promoted, std::log(120.0)));
}

template <class Semiring>
KeepAllResult check_keep_all_factorization(
    const std::string &sequence,
    const std::string &scaffold) {
  std::string baseline_sequence = sequence;
  std::string baseline_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot baseline_owner(
      baseline_sequence, baseline_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext baseline_ctx(
      baseline_owner);
  const auto baseline_schedule =
      cparty::scfg2::engine::build_pkfree_basic_exact_schedule(
          static_cast<cand_pos_t>(sequence.size()), baseline_ctx);
  const auto baseline =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<
          Semiring>(
          baseline_schedule,
          baseline_ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &,
             const Deduction &,
             double,
             const cparty::scfg2::engine::Chart<double> &,
             double) {});

  std::string factorized_sequence = sequence;
  std::string factorized_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot factorized_owner(
      factorized_sequence, factorized_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext factorized_ctx(
      factorized_owner);
  const auto factorized_schedule =
      cparty::scfg2::engine::build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), factorized_ctx);
  assert(baseline_schedule.size() == factorized_schedule.size());

  cparty::scfg2::beam::ExactVpCase3TailCache tail_cache;
  const auto reference_tail_plan =
      cparty::scfg2::beam::ExactVpCase3TailColumnPlan::discover(
          static_cast<cand_pos_t>(sequence.size()), factorized_ctx);
  auto border_event_plan =
      cparty::scfg2::beam::ExactVpCase3BorderEventPlan::discover(
          static_cast<cand_pos_t>(sequence.size()), factorized_ctx);
  assert(border_event_plan.border_mismatch_count() == 0);
  assert(border_event_plan.duplicate_source_count() == 0);
  assert(border_event_plan.source_count() == reference_tail_plan.source_count());
  assert(source_identities(border_event_plan) ==
         source_identities(reference_tail_plan));
  const auto border_event_endpoint_pair_count =
      border_event_plan.endpoint_pair_count();
  const auto child_event_index =
      cparty::scfg2::beam::ExactVpCase3BorderEventChildIndex::build(
          static_cast<cand_pos_t>(sequence.size()), factorized_ctx);
  std::vector<cparty::scfg2::beam::VpCase3PlannedSource> child_sources;
  std::size_t child_endpoint_pair_count = 0;
  for (const auto &item : factorized_schedule) {
    const auto query =
        child_event_index.query_completed_wi(item, factorized_ctx);
    assert(query.border_mismatch_count == 0);
    assert(query.duplicate_source_count == 0);
    child_endpoint_pair_count += query.endpoint_pair_count;
    child_sources.insert(child_sources.end(),
                         query.sources.begin(),
                         query.sources.end());
  }
  const auto child_event_plan =
      cparty::scfg2::beam::ExactVpCase3TailColumnPlan::from_sources(
          std::move(child_sources));
  assert(child_endpoint_pair_count <= border_event_endpoint_pair_count);
  assert(source_identities(child_event_plan) ==
         source_identities(reference_tail_plan));
  std::vector<ItemKey> left_wi_items;
  for (const auto &item : factorized_schedule) {
    if (item.nonterminal == cparty::scfg2::grammar::NonTerminal::WI) {
      left_wi_items.push_back(item);
    }
  }
  std::vector<cparty::scfg2::beam::VpCase3PlannedSource> right_sources;
  std::size_t right_endpoint_pair_count = 0;
  for (cand_pos_t j = 1;
       j <= static_cast<cand_pos_t>(sequence.size());
       ++j) {
    const auto query = child_event_index.query_target_endpoint(
        j, left_wi_items, factorized_ctx);
    assert(query.border_mismatch_count == 0);
    assert(query.duplicate_source_count == 0);
    assert(query.readiness_violation_count == 0);
    right_endpoint_pair_count += query.endpoint_pair_count;
    right_sources.insert(right_sources.end(),
                         query.sources.begin(),
                         query.sources.end());
  }
  const auto right_event_plan =
      cparty::scfg2::beam::ExactVpCase3TailColumnPlan::from_sources(
          std::move(right_sources));
  assert(source_identities(right_event_plan) ==
         source_identities(reference_tail_plan));
  std::map<cand_pos_t, std::vector<ItemKey>> left_wi_by_right;
  for (const auto &item : left_wi_items) {
    left_wi_by_right[item.j].push_back(item);
  }
  cparty::scfg2::beam::ExactVpCase3RightEndpointProducerRegistry registry;
  std::vector<cparty::scfg2::beam::VpCase3PlannedSource> registry_sources;
  std::size_t registry_endpoint_pair_count = 0;
  for (cand_pos_t j = 1;
       j <= static_cast<cand_pos_t>(sequence.size());
       ++j) {
    const auto query = registry.query_target_endpoint(
        j, child_event_index, factorized_ctx);
    assert(query.border_mismatch_count == 0);
    assert(query.duplicate_source_count == 0);
    assert(query.readiness_violation_count == 0);
    registry_endpoint_pair_count += query.endpoint_pair_count;
    registry_sources.insert(registry_sources.end(),
                            query.sources.begin(),
                            query.sources.end());
    const auto completed = left_wi_by_right.find(j);
    if (completed == left_wi_by_right.end()) continue;
    for (const auto &item : completed->second) {
      registry.observe_completed_wi(
          item, 1.0, child_event_index, factorized_ctx);
    }
  }
  const auto registry_event_plan =
      cparty::scfg2::beam::ExactVpCase3TailColumnPlan::from_sources(
          std::move(registry_sources));
  assert(source_identities(registry_event_plan) ==
         source_identities(reference_tail_plan));
  assert(registry_endpoint_pair_count == right_endpoint_pair_count);
  const auto tail_plan = std::move(border_event_plan).into_tail_column_plan();
  const auto attachment_segments =
      cparty::scfg2::beam::build_vp_case3_tail_attachment_segments(tail_plan);
  const auto border_profile =
      cparty::scfg2::beam::profile_vp_case3_tail_borders(tail_plan);
  assert(cparty::scfg2::beam::vp_case3_tail_segments_reconstruct_sources(
      tail_plan, attachment_segments));
  assert(border_profile.exact_segment_count == attachment_segments.size());
  assert(border_profile.constant_p_segment_count ==
         border_profile.contiguous_i_run_count);
  assert(border_profile.exact_segment_count ==
         border_profile.constant_p_segment_count);
  cand_pos_t materialized_j = -1;
  const auto factorized =
      cparty::scfg2::replay::
          run_w_final_exact_dp_schedule_with_item_and_contribution<Semiring>(
              factorized_schedule,
              factorized_ctx,
              [](const Deduction &) { return true; },
              [](const ItemKey &,
                 const Deduction &,
                 double,
                 const cparty::scfg2::engine::Chart<double> &,
                 double) {},
              [&](const ItemKey &item,
                  const cparty::scfg2::engine::Chart<double> &chart) {
                if (item.j == materialized_j) return;
                materialized_j = item.j;
                (void)tail_plan.template materialize_column<Semiring>(
                    item.j, chart, tail_cache);
              },
              [&](const Deduction &deduction,
                  const cparty::scfg2::engine::Chart<double> &chart,
                  double local) {
                const auto reconstructed =
                    tail_cache.template
                        reconstruct_materialized_left_source<Semiring>(
                        deduction, chart, local);
                if (reconstructed.has_value()) return *reconstructed;
                assert(cparty::scfg2::beam::vp_case3_partition(deduction) !=
                       cparty::scfg2::beam::VpCase3Partition::LeftNonempty);
                const double direct =
                    cparty::scfg2::replay::
                        runtime_deduction_contribution<Semiring>(
                            chart, deduction, local);
                return direct;
              });

  for (const auto &item : baseline_schedule) {
    assert(factorized.contains(item));
    assert(close(baseline.get(item), factorized.get(item)));
  }
  for (const auto &item : factorized_schedule) {
    assert(baseline.contains(item));
  }

  assert(tail_cache.full_source_attachments_are_unique());
  assert(tail_cache.materialization_count() == tail_cache.tail_count());
  assert(tail_cache.reuse_count() == tail_cache.attachment_count());
  if (tail_cache.tail_count() > 0) {
    assert(tail_cache.attachment_count() >= tail_cache.tail_count());
  assert(tail_cache.max_attachments_per_tail() >= 1);
  }
  assert(tail_plan.tail_count() == tail_cache.tail_count());
  assert(tail_plan.source_count() == tail_cache.attachment_count());
  const auto schedule_profile =
      cparty::scfg2::beam::profile_vp_case3_tail_schedule(
          tail_cache, factorized_schedule);
  assert(schedule_profile.missing_schedule_item_count == 0);
  assert(schedule_profile.topological_tail_count == tail_cache.tail_count());
  assert(schedule_profile.committable_column_count ==
         schedule_profile.column_count);
  assert(schedule_profile.noncommittable_column_count == 0);
  return KeepAllResult{
      std::move(tail_cache),
      schedule_profile,
      tail_plan.source_count(),
      attachment_segments.size(),
      border_event_endpoint_pair_count,
      child_endpoint_pair_count,
      right_endpoint_pair_count,
      registry_endpoint_pair_count};
}

}  // namespace

int main(int argc, char **argv) {
  check_source_observation_boundary();
  assert(argc == 1 || argc == 3 || argc == 4);
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc == 4 ? argv[3] : "params/rna_Turner04.par",
         1);
  const std::string sequence = argc >= 3
      ? argv[1]
      : "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  const std::string scaffold = argc >= 3
      ? argv[2]
      : "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  assert(sequence.size() == scaffold.size());
  check_log_domain_factorization();
  check_coupled_value_factorization();
  const auto inside_cache =
      check_keep_all_factorization<
          cparty::scfg2::semirings::SumProductSemiring>(sequence, scaffold);
  const auto viterbi_cache =
      check_keep_all_factorization<
          cparty::scfg2::semirings::MaxProductSemiring>(sequence, scaffold);
  assert(inside_cache.cache.tail_count() == viterbi_cache.cache.tail_count());
  assert(inside_cache.cache.attachment_count() ==
         viterbi_cache.cache.attachment_count());
  assert(inside_cache.planned_source_count ==
         viterbi_cache.planned_source_count);
  assert(inside_cache.attachment_segment_count ==
         viterbi_cache.attachment_segment_count);
  assert(inside_cache.border_event_endpoint_pair_count ==
         viterbi_cache.border_event_endpoint_pair_count);
  assert(inside_cache.child_event_endpoint_pair_count ==
         viterbi_cache.child_event_endpoint_pair_count);
  assert(inside_cache.right_event_endpoint_pair_count ==
         viterbi_cache.right_event_endpoint_pair_count);
  assert(inside_cache.registry_event_endpoint_pair_count ==
         viterbi_cache.registry_event_endpoint_pair_count);
  assert(inside_cache.schedule_profile.committable_column_count ==
         viterbi_cache.schedule_profile.committable_column_count);
  std::cout << "n\ttails\tsources\tsegments\tevent_pairs\tchild_event_pairs"
               "\tright_event_pairs"
               "\tregistry_event_pairs"
               "\tattachments\tmaterialized\treused"
               "\tcolumns\tcommittable\tnoncommittable\tparity\n"
            << sequence.size() << '\t' << inside_cache.cache.tail_count()
            << '\t' << inside_cache.planned_source_count << '\t'
            << inside_cache.attachment_segment_count << '\t'
            << inside_cache.border_event_endpoint_pair_count << '\t'
            << inside_cache.child_event_endpoint_pair_count << '\t'
            << inside_cache.right_event_endpoint_pair_count << '\t'
            << inside_cache.registry_event_endpoint_pair_count << '\t'
            << inside_cache.cache.attachment_count() << '\t'
            << inside_cache.cache.materialization_count() << '\t'
            << inside_cache.cache.reuse_count() << '\t'
            << inside_cache.schedule_profile.column_count << '\t'
            << inside_cache.schedule_profile.committable_column_count << '\t'
            << inside_cache.schedule_profile.noncommittable_column_count
            << "\t1\n";
  return EXIT_SUCCESS;
}
