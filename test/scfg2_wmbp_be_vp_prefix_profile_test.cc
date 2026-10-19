#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/wmbp_be_vp_prefix_cache.hh"
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

bool close(double lhs, double rhs) {
  return std::fabs(lhs - rhs) <=
         std::max(1e-14, std::max(std::fabs(lhs), std::fabs(rhs)) * 1e-12);
}

void check_log_domain_cache() {
  namespace beam = cparty::scfg2::beam;
  namespace engine = cparty::scfg2::engine;
  namespace grammar = cparty::scfg2::grammar;
  const ItemKey be{grammar::NonTerminal::BE, 1, 4, 2, 3};
  const ItemKey left{grammar::NonTerminal::WMBP, 5, 7, -1, -1};
  const ItemKey right{grammar::NonTerminal::VP, 8, 10, -1, -1};
  engine::Chart<double> chart(TestLogMaxSemiring::zero());
  chart.set(be, std::log(2.0));
  chart.set(left, std::log(3.0));
  chart.set(right, std::log(4.0));
  Deduction source;
  source.parent = ItemKey{grammar::NonTerminal::WMBP, 1, 10, -1, -1};
  source.rule = grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP;
  source.children[0].item = be;
  source.children[1].item = left;
  source.children[2].item = right;
  source.child_count = 3;
  beam::ExactWmbpBeVpPrefixCache cache;
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
  const ItemKey be{grammar::NonTerminal::BE, 1, 4, 2, 3};
  const ItemKey left{grammar::NonTerminal::WMBP, 5, 7, -1, -1};
  const ItemKey right{grammar::NonTerminal::VP, 8, 10, -1, -1};
  engine::Chart<Value> chart(Semiring::zero());
  chart.set(be, Value{std::log(5.0), std::log(2.0)});
  chart.set(left, Value{std::log(7.0), std::log(3.0)});
  chart.set(right, Value{std::log(11.0), std::log(4.0)});
  Deduction source;
  source.parent = ItemKey{grammar::NonTerminal::WMBP, 1, 10, -1, -1};
  source.rule = grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP;
  source.children[0].item = be;
  source.children[1].item = left;
  source.children[2].item = right;
  source.child_count = 3;
  beam::BasicExactWmbpBeVpPrefixCache<Value> cache;
  const auto reconstructed = cache.reconstruct_source<Semiring>(
      source, chart, 13.0);
  assert(reconstructed.has_value());
  assert(close(reconstructed->log_sum, std::log(5.0 * 7.0 * 11.0 * 13.0)));
  assert(close(reconstructed->log_max, std::log(2.0 * 3.0 * 4.0 * 13.0)));
}

struct Run {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  cparty::scfg2::beam::WmbpBeVpPrefixProfile profile{};
  std::size_t readiness_checks = 0;
};

struct CacheRun {
  cparty::scfg2::engine::Chart<double> chart{0.0};
  std::size_t prefixes = 0;
  std::size_t reuses = 0;
  std::size_t inconsistent_reuses = 0;
  std::size_t reconstructed_sources = 0;
};

template <class Semiring>
Run run(const std::string &sequence, const std::string &scaffold) {
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

  cparty::scfg2::beam::WmbpBeVpPrefixProfiler profiler;
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
            if (!cparty::scfg2::beam::is_wmbp_be_vp_rule(deduction.rule)) {
              return;
            }
            assert(deduction.child_count == 2 || deduction.child_count == 3);
            assert(deduction.children[0].item.nonterminal ==
                   cparty::scfg2::grammar::NonTerminal::BE);
            assert(deduction.children[deduction.child_count - 1]
                       .item.nonterminal ==
                   cparty::scfg2::grammar::NonTerminal::VP);
            const auto parent_it = schedule_position.find(parent);
            assert(parent_it != schedule_position.end());
            for (std::size_t child_idx = 0;
                 child_idx < deduction.child_count;
                 ++child_idx) {
              const auto child_it =
                  schedule_position.find(deduction.children[child_idx].item);
              assert(child_it != schedule_position.end());
              assert(child_it->second < parent_it->second);
              ++readiness_checks;
            }
            profiler.observe(deduction, chart, local, contribution);
          });
  assert(profiler.all_factor_identities_hold());

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
  return Run{observed, profiler.profile(), readiness_checks};
}

template <class Semiring>
CacheRun run_cache(const std::string &sequence, const std::string &scaffold) {
  std::string run_sequence = sequence;
  std::string run_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      run_sequence, run_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  cparty::scfg2::beam::ExactWmbpBeVpPrefixCache cache;
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
                  cache.prefix_count(),
                  cache.reuse_count(),
                  cache.inconsistent_reuse_count(),
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
      run<cparty::scfg2::semirings::SumProductSemiring>(sequence, scaffold);
  const auto viterbi =
      run<cparty::scfg2::semirings::MaxProductSemiring>(sequence, scaffold);
  const auto inside_cache =
      run_cache<cparty::scfg2::semirings::SumProductSemiring>(
          sequence, scaffold);
  const auto viterbi_cache =
      run_cache<cparty::scfg2::semirings::MaxProductSemiring>(
          sequence, scaffold);
  const auto &profile = inside.profile;
  assert(profile.source_count == viterbi.profile.source_count);
  assert(profile.prefix_count == viterbi.profile.prefix_count);
  assert(profile.right_vp_count == viterbi.profile.right_vp_count);
  assert(profile.local_variant_prefix_count ==
         viterbi.profile.local_variant_prefix_count);
  assert(profile.wmbp_count + profile.wmbw_count + profile.wi_empty_count +
             profile.wi_count ==
         profile.source_count);
  assert(inside.readiness_checks == viterbi.readiness_checks);
  assert(inside.readiness_checks >= 2 * profile.source_count);
  for (const auto &[item, value] : inside.chart.values()) {
    assert(close(value, inside_cache.chart.get(item)));
  }
  for (const auto &[item, value] : viterbi.chart.values()) {
    assert(close(value, viterbi_cache.chart.get(item)));
  }
  assert(inside_cache.prefixes == profile.prefix_count);
  assert(inside_cache.reuses == profile.reused_prefix_count);
  assert(inside_cache.inconsistent_reuses == 0);
  assert(inside_cache.reconstructed_sources == profile.source_count);
  assert(viterbi_cache.prefixes == inside_cache.prefixes);
  assert(viterbi_cache.reuses == inside_cache.reuses);
  assert(viterbi_cache.inconsistent_reuses == 0);
  if (argc == 1) assert(profile.source_count > 0);

  std::cout << "n\tsources\twmbp\twmbw\twi_empty\twi\tprefixes"
               "\treused_prefixes\tright_vps\treused_right_vps"
               "\tmulti_source_prefixes\tmax_sources_per_prefix"
               "\tlocal_variant_prefixes\tmax_locals_per_prefix"
               "\tmax_prefixes_per_right_vp\tcolumns"
               "\tmax_sources_per_column\tcache_prefixes\tcache_reuses"
               "\treadiness_checks\tparity\n"
            << sequence.size() << '\t' << profile.source_count << '\t'
            << profile.wmbp_count << '\t' << profile.wmbw_count << '\t'
            << profile.wi_empty_count << '\t' << profile.wi_count << '\t'
            << profile.prefix_count << '\t' << profile.reused_prefix_count
            << '\t' << profile.right_vp_count << '\t'
            << profile.reused_right_vp_count << '\t'
            << profile.multi_source_prefix_count << '\t'
            << profile.maximum_sources_per_prefix << '\t'
            << profile.local_variant_prefix_count << '\t'
            << profile.maximum_locals_per_prefix << '\t'
            << profile.maximum_prefixes_per_right_vp << '\t'
            << profile.column_count << '\t'
            << profile.maximum_sources_per_column << '\t'
            << inside_cache.prefixes << '\t' << inside_cache.reuses << '\t'
            << inside.readiness_checks << "\t1\n";
}
