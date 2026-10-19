#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/be_frame_reuse_profile.hh"
#include "scfg2/beam/be_frame_cache.hh"
#include "scfg2/beam/be_left_right_frame_profile.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>

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

struct ProfileRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  cparty::scfg2::beam::BeFrameReuseProfile profile{};
  cparty::scfg2::beam::BeRightShellFrameProfile frame_profile{};
  std::size_t readiness_checks = 0;
};

struct CacheRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t cores = 0;
  std::size_t core_reuses = 0;
  std::size_t right_tails = 0;
  std::size_t right_tail_reuses = 0;
  std::size_t reconstructed_sources = 0;
};

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

void check_log_domain_cache() {
  namespace beam = cparty::scfg2::beam;
  namespace engine = cparty::scfg2::engine;
  namespace grammar = cparty::scfg2::grammar;

  engine::Chart<double> chart(TestLogMaxSemiring::zero());
  chart.set(ItemKey{grammar::NonTerminal::BE, 3, 6, 4, 5},
            std::log(2.0));
  chart.set(ItemKey{grammar::NonTerminal::WIP, 7, 9, -1, -1},
            std::log(3.0));
  chart.set(ItemKey{grammar::NonTerminal::WIP, 2, 2, -1, -1},
            std::log(4.0));

  beam::ExactBeFrameCache cache;
  const beam::BeRightTailKey tail{10, beam::BeCoreKey{3, 6, 4, 5}};
  assert(close(cache.materialize_right_tail<TestLogMaxSemiring>(tail, chart),
               std::log(6.0)));

  Deduction source;
  source.parent = ItemKey{grammar::NonTerminal::BE, 1, 10, 4, 5};
  source.rule = grammar::RuleId::BE_WIP_WIP;
  source.split = grammar::RuleSplit{3, 6, -1, -1};
  source.children[0].item =
      ItemKey{grammar::NonTerminal::WIP, 2, 2, -1, -1};
  source.child_count = 3;
  const auto reconstructed = cache.reconstruct_source<TestLogMaxSemiring>(
      source, chart, 5.0);
  assert(reconstructed.has_value());
  assert(close(*reconstructed, std::log(120.0)));
}

void check_coupled_value_cache() {
  namespace beam = cparty::scfg2::beam;
  namespace engine = cparty::scfg2::engine;
  namespace grammar = cparty::scfg2::grammar;
  using Semiring = cparty::scfg2::semirings::CoupledLogSumMaxSemiring;
  using Value = Semiring::value_type;

  engine::Chart<Value> chart(Semiring::zero());
  chart.set(ItemKey{grammar::NonTerminal::BE, 3, 6, 4, 5},
            Value{std::log(5.0), std::log(2.0)});
  chart.set(ItemKey{grammar::NonTerminal::WIP, 7, 9, -1, -1},
            Value{std::log(7.0), std::log(3.0)});
  chart.set(ItemKey{grammar::NonTerminal::WIP, 2, 2, -1, -1},
            Value{std::log(11.0), std::log(4.0)});

  beam::BasicExactBeFrameCache<Value> cache;
  Deduction source;
  source.parent = ItemKey{grammar::NonTerminal::BE, 1, 10, 4, 5};
  source.rule = grammar::RuleId::BE_WIP_WIP;
  source.split = grammar::RuleSplit{3, 6, -1, -1};
  source.children[0].item =
      ItemKey{grammar::NonTerminal::WIP, 2, 2, -1, -1};
  source.child_count = 3;
  const auto reconstructed = cache.reconstruct_source<Semiring>(
      source, chart, 13.0);
  assert(reconstructed.has_value());
  assert(close(reconstructed->log_sum, std::log(5.0 * 7.0 * 11.0 * 13.0)));
  assert(close(reconstructed->log_max, std::log(2.0 * 3.0 * 4.0 * 13.0)));
}

template <class Semiring>
ProfileRun run_profile(const std::string &sequence,
                       const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  std::unordered_map<ItemKey,
                     std::size_t,
                     cparty::scfg2::deductions::ItemKeyHash,
                     cparty::scfg2::deductions::ItemKeyEq>
      schedule_position;
  for (std::size_t idx = 0; idx < schedule.size(); ++idx) {
    schedule_position.emplace(schedule[idx], idx);
  }

  cparty::scfg2::beam::BeFrameReuseProfiler profiler;
  cparty::scfg2::beam::BeRightShellFrameProfiler frame_profiler;
  std::size_t readiness_checks = 0;
  const auto observed = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule<Semiring>(
          schedule,
          ctx,
          [](const Deduction &) { return true; },
          [&](const ItemKey &parent,
              const Deduction &deduction,
              double contribution,
              const cparty::scfg2::engine::Chart<double> &chart,
              double local) {
            const bool recursive =
                deduction.rule ==
                    cparty::scfg2::grammar::RuleId::BE_INTERNAL_LOOP ||
                deduction.rule ==
                    cparty::scfg2::grammar::RuleId::BE_WIP_WIP ||
                deduction.rule ==
                    cparty::scfg2::grammar::RuleId::BE_WIP_BASEPAIR ||
                deduction.rule ==
                    cparty::scfg2::grammar::RuleId::BE_BASEPAIR_WIP;
            if (!recursive) {
              return;
            }
            const auto parent_position = schedule_position.find(parent);
            assert(parent_position != schedule_position.end());
            double factorized = local;
            for (std::size_t child_idx = 0;
                 child_idx < deduction.child_count;
                 ++child_idx) {
              const auto &child = deduction.children[child_idx].item;
              const auto child_position = schedule_position.find(child);
              assert(child_position != schedule_position.end());
              assert(child_position->second < parent_position->second);
              factorized *= chart.get(child);
              ++readiness_checks;
            }
            assert(close(contribution, factorized));
            profiler.observe(deduction);
            frame_profiler.observe(deduction, chart, local, contribution);
          });

  std::string baseline_sequence = sequence;
  std::string baseline_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot baseline_owner(
      baseline_sequence, baseline_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext baseline_ctx(baseline_owner);
  const auto baseline =
      cparty::scfg2::replay::run_w_final_exact_dp_schedule<Semiring>(
          schedule,
          baseline_ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          });
  for (const auto &item : schedule) {
    assert(close(observed.get(item), baseline.get(item)));
  }
  assert(frame_profiler.all_factor_identities_hold());
  return ProfileRun{
      observed, profiler.profile(), frame_profiler.profile(), readiness_checks};
}

template <class Semiring>
CacheRun run_cache(const std::string &sequence,
                   const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  cparty::scfg2::beam::ExactBeFrameCache cache;
  auto chart = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_contribution<Semiring>(
          schedule,
          ctx,
          [](const Deduction &) { return true; },
          [](const ItemKey &, const Deduction &, double, const auto &, double) {
          },
          [&](const Deduction &deduction,
              const auto &current_chart,
              double local) {
            const auto reconstructed =
                cache.template reconstruct_source<Semiring>(
                    deduction, current_chart, local);
            if (reconstructed.has_value()) return *reconstructed;
            return cparty::scfg2::replay::
                runtime_deduction_contribution<Semiring>(
                    current_chart, deduction, local);
          });
  return CacheRun{std::move(chart),
                  cache.core_factor_count(),
                  cache.core_reuse_count(),
                  cache.right_tail_factor_count(),
                  cache.right_tail_reuse_count(),
                  cache.reconstructed_source_count()};
}

}  // namespace

int main(int argc, char **argv) {
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
  check_log_domain_cache();
  check_coupled_value_cache();

  const auto inside =
      run_profile<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi =
      run_profile<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  const auto inside_cache =
      run_cache<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi_cache =
      run_cache<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  const auto &profile = inside.profile;
  assert(profile.source_count == viterbi.profile.source_count);
  assert(profile.core_count == viterbi.profile.core_count);
  assert(profile.right_tail_count == viterbi.profile.right_tail_count);
  assert(profile.frame_count == viterbi.profile.frame_count);
  assert(profile.recursive_source_count ==
         viterbi.profile.recursive_source_count);
  assert(profile.shared_core_count == viterbi.profile.shared_core_count);
  assert(profile.shared_right_tail_count ==
         viterbi.profile.shared_right_tail_count);
  assert(inside.frame_profile.source_count ==
         viterbi.frame_profile.source_count);
  assert(inside.frame_profile.frame_count ==
         viterbi.frame_profile.frame_count);
  assert(inside.frame_profile.tail_count ==
         viterbi.frame_profile.tail_count);
  assert(inside.frame_profile.source_count ==
         profile.right_shell_source_count);
  assert(inside.frame_profile.tail_count == profile.shared_right_tail_count);
  assert(inside.readiness_checks >= 3 * profile.source_count);
  assert(viterbi.readiness_checks == inside.readiness_checks);
  for (const auto &[item, value] : inside.chart.values()) {
    assert(close(value, inside_cache.chart.get(item)));
  }
  for (const auto &[item, value] : viterbi.chart.values()) {
    assert(close(value, viterbi_cache.chart.get(item)));
  }
  assert(inside_cache.cores == profile.shared_core_count);
  assert(inside_cache.core_reuses + inside_cache.cores ==
         profile.recursive_source_count);
  assert(inside_cache.right_tails == profile.shared_right_tail_count);
  assert(inside_cache.right_tail_reuses + inside_cache.right_tails ==
         profile.right_shell_source_count);
  assert(inside_cache.reconstructed_sources ==
         profile.recursive_source_count);
  assert(viterbi_cache.cores == inside_cache.cores);
  assert(viterbi_cache.core_reuses == inside_cache.core_reuses);
  assert(viterbi_cache.right_tails == inside_cache.right_tails);
  assert(viterbi_cache.right_tail_reuses ==
         inside_cache.right_tail_reuses);
  if (argc == 1) assert(profile.source_count > 0);

  std::cout
      << "n\tsources\tcores\treused_cores\tmax_sources_per_core"
         "\tmax_right_tails_per_core\tright_tails\treused_right_tails"
         "\tmax_sources_per_right_tail\tframes\tcore_columns"
         "\tmulti_core_columns\tmax_cores_per_column\tright_tail_columns"
         "\tmulti_right_tail_columns\tmax_right_tails_per_column"
         "\tframe_columns\tmulti_frame_columns\tmax_frames_per_column"
         "\trecursive_sources\tshared_cores\treused_shared_cores"
         "\tmax_recursive_sources_per_core\tright_shell_sources"
         "\tshared_right_tails\treused_shared_right_tails"
         "\tmax_right_shell_sources_per_tail\tcache_cores"
         "\tcache_core_reuses\tcache_right_tails"
         "\tcache_right_tail_reuses\tright_only_frames"
         "\tleft_right_frames\tfinal_frames\tmulti_frame_tails"
         "\tmax_frames_per_tail\tmax_final_frames_per_column"
         "\treadiness_checks\tparity\n"
      << sequence.size() << '\t' << profile.source_count << '\t'
      << profile.core_count << '\t' << profile.reused_core_count << '\t'
      << profile.maximum_sources_per_core << '\t'
      << profile.maximum_right_tails_per_core << '\t'
      << profile.right_tail_count << '\t'
      << profile.reused_right_tail_count << '\t'
      << profile.maximum_sources_per_right_tail << '\t'
      << profile.frame_count << '\t' << profile.core_column_count << '\t'
      << profile.multi_core_column_count << '\t'
      << profile.maximum_cores_per_column << '\t'
      << profile.right_tail_column_count << '\t'
      << profile.multi_right_tail_column_count << '\t'
      << profile.maximum_right_tails_per_column << '\t'
      << profile.frame_column_count << '\t'
      << profile.multi_frame_column_count << '\t'
      << profile.maximum_frames_per_column << '\t'
      << profile.recursive_source_count << '\t'
      << profile.shared_core_count << '\t'
      << profile.reused_shared_core_count << '\t'
      << profile.maximum_recursive_sources_per_core << '\t'
      << profile.right_shell_source_count << '\t'
      << profile.shared_right_tail_count << '\t'
      << profile.reused_shared_right_tail_count << '\t'
      << profile.maximum_right_shell_sources_per_tail << '\t'
      << inside_cache.cores << '\t' << inside_cache.core_reuses << '\t'
      << inside_cache.right_tails << '\t'
      << inside_cache.right_tail_reuses << '\t'
      << inside.frame_profile.right_only_count << '\t'
      << inside.frame_profile.left_right_count << '\t'
      << inside.frame_profile.frame_count << '\t'
      << inside.frame_profile.multi_frame_tail_count << '\t'
      << inside.frame_profile.maximum_frames_per_tail << '\t'
      << inside.frame_profile.maximum_frames_per_column << '\t'
      << inside.readiness_checks << "\t1\n";
  return EXIT_SUCCESS;
}
