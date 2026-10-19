#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_left_selection.hh"
#include "scfg2/beam/vp_case3_tail_selection.hh"
#include "scfg2/beam/vp_case3_wi_frontier.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {

struct OnlineProfile {
  double root = 0.0;
  std::size_t producer_candidates = 0;
  std::size_t producer_selected = 0;
  std::size_t endpoint_pairs = 0;
  std::size_t sources = 0;
  std::size_t tails = 0;
  std::size_t tails_selected = 0;
  std::size_t left_candidates = 0;
  std::size_t left_selected = 0;
  std::size_t rejected_sources = 0;
  std::size_t active_expiry_probes = 0;
  std::size_t tail_materialization_probes = 0;
  std::size_t prepared_score_probes = 0;
  std::size_t max_active_segments = 0;
};

template <class Semiring, class Context>
OnlineProfile run_online(
    cand_pos_t n,
    std::size_t producer_width,
    std::size_t tail_width,
    std::size_t left_width,
    cparty::scfg2::beam::VpCase3WiFrontierScope frontier_scope,
    const Context &ctx,
    const std::vector<cparty::scfg2::deductions::ItemKey> &schedule,
    const cparty::scfg2::beam::ExactVpCase3BorderEventChildIndex &index) {
  cparty::scfg2::beam::PendingVpCase3LeftCandidates pending(left_width);
  cparty::scfg2::beam::ExactVpCase3TailCache tails;
  cparty::scfg2::beam::ExactVpCase3LeftCache left;
  cparty::scfg2::beam::TopWidthVpCase3LeftSelection left_selection(
      left_width);
  cparty::scfg2::beam::TopWidthVpCase3TailSelection tail_selection(
      tail_width);
  cparty::scfg2::beam::VpCase3WiProducerFrontier producer(
      producer_width, frontier_scope);
  cand_pos_t materialized_j = -1;
  std::size_t rejected_sources = 0;
  const auto chart = cparty::scfg2::replay::
      run_w_final_exact_dp_schedule_with_item_observers_and_contribution<
          Semiring>(
          schedule,
          ctx,
          [](const auto &) { return true; },
          [](const auto &, const auto &, double, const auto &, double) {},
          [&](const auto &item, const auto &current_chart) {
            if (item.j == materialized_j) return;
            materialized_j = item.j;
            producer.flush_before(
                item.j, current_chart, pending, index, ctx);
            pending.prepare_column(item.j, current_chart, tails);
            tail_selection.commit_column(item.j, tails);
            const auto scored = pending.score_prepared_column(
                item.j,
                tails,
                [&](const auto &key) { return tail_selection.allows(key); });
            left_selection.commit_column(item.j, scored);
            (void)left.materialize_scored(
                scored, [&](const auto &key) {
                  return left_selection.allows(key);
                });
          },
          [&](const auto &item, const auto &current_chart) {
            producer.observe_completed_child(item, current_chart, index, ctx);
          },
          [&](const auto &deduction,
              const auto &current_chart,
              double local) {
            if (!left_selection.allows(deduction)) {
              ++rejected_sources;
              return Semiring::zero();
            }
            const auto reconstructed =
                left.reconstruct_source(deduction, local);
            if (reconstructed.has_value()) return *reconstructed;
            return cparty::scfg2::replay::
                runtime_deduction_contribution<Semiring>(
                    current_chart, deduction, local);
          });
  return OnlineProfile{
      chart.get(cparty::scfg2::deductions::ItemKey{
          cparty::scfg2::grammar::NonTerminal::W, 1, n, -1, -1}),
      producer.candidate_count(),
      producer.selected_count(),
      pending.queried_endpoint_pair_count(),
      pending.queried_source_count(),
      tails.tail_count(),
      tail_selection.selected_count(),
      left_selection.candidate_count(),
      left_selection.selected_count(),
      rejected_sources,
      pending.active_expiry_probe_count(),
      pending.tail_materialization_probe_count(),
      pending.prepared_score_probe_count(),
      pending.max_active_segment_count()};
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 4 || argc > 10) {
    std::cerr << "usage: scfg2_vp_case3_online_profile_tool <seq> <G> "
                 "<producer_width|all> [repetitions] [inside|max] "
                 "[param_file] [eligible|all-wi] [tail_width|all] "
                 "[left_width|all]\n";
    return 2;
  }
  try {
    std::string sequence = argv[1];
    std::string scaffold = argv[2];
    if (sequence.size() != scaffold.size()) {
      throw std::runtime_error("sequence and scaffold lengths differ");
    }
    const auto parse_width = [](const std::string &arg) {
      return arg == "all" ? std::numeric_limits<std::size_t>::max()
                          : static_cast<std::size_t>(std::stoull(arg));
    };
    const std::string producer_width_arg = argv[3];
    const std::size_t producer_width = parse_width(producer_width_arg);
    const std::size_t repetitions =
        argc >= 5 ? static_cast<std::size_t>(std::stoull(argv[4])) : 1;
    if (repetitions == 0) throw std::runtime_error("repetitions must be positive");
    const std::string mode = argc >= 6 ? argv[5] : "inside";
    if (mode != "inside" && mode != "max") {
      throw std::runtime_error("mode must be inside or max");
    }
    const std::string scope_arg = argc >= 8 ? argv[7] : "eligible";
    if (scope_arg != "eligible" && scope_arg != "all-wi") {
      throw std::runtime_error("scope must be eligible or all-wi");
    }
    const auto frontier_scope = scope_arg == "eligible"
        ? cparty::scfg2::beam::VpCase3WiFrontierScope::EligibleCase3Producers
        : cparty::scfg2::beam::VpCase3WiFrontierScope::AllWiItems;
    const std::string tail_width_arg = argc >= 9 ? argv[8] : "all";
    const std::string left_width_arg = argc >= 10 ? argv[9] : "all";
    const std::size_t tail_width = parse_width(tail_width_arg);
    const std::size_t left_width = parse_width(left_width_arg);
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
           argc >= 7 ? argv[6] : "params/rna_Turner04.par",
           1);

    const auto setup_begin = std::chrono::steady_clock::now();
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
    const cand_pos_t n = static_cast<cand_pos_t>(sequence.size());
    const auto schedule = cparty::scfg2::engine::
        build_pkfree_basic_exact_right_boundary_schedule(n, ctx);
    const auto index =
        cparty::scfg2::beam::ExactVpCase3BorderEventChildIndex::build(n, ctx);
    const auto setup_end = std::chrono::steady_clock::now();

    OnlineProfile profile;
    const auto dp_begin = std::chrono::steady_clock::now();
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition) {
      profile = mode == "inside"
          ? run_online<cparty::scfg2::semirings::SumProductSemiring>(
                n,
                producer_width,
                tail_width,
                left_width,
                frontier_scope,
                ctx,
                schedule,
                index)
          : run_online<cparty::scfg2::semirings::MaxProductSemiring>(
                n,
                producer_width,
                tail_width,
                left_width,
                frontier_scope,
                ctx,
                schedule,
                index);
    }
    const auto dp_end = std::chrono::steady_clock::now();
    const double setup_ms =
        std::chrono::duration<double, std::milli>(setup_end - setup_begin)
            .count();
    const double dp_ms =
        std::chrono::duration<double, std::milli>(dp_end - dp_begin).count() /
        static_cast<double>(repetitions);
    std::cout << std::setprecision(17)
              << "n\tmode\tscope\tproducer_width\ttail_width\tleft_width"
                 "\trepetitions\tsetup_ms\tdp_ms_per_run\troot"
                 "\tproducer_candidates\tproducer_selected\tendpoint_pairs"
                 "\tsources\ttails\ttails_selected\tleft_candidates"
                 "\tleft_selected\trejected_sources\tactive_expiry_probes"
                 "\ttail_materialization_probes\tprepared_score_probes"
                 "\tmax_active_segments\n"
              << n << '\t' << mode << '\t' << scope_arg << '\t'
              << producer_width_arg << '\t' << tail_width_arg << '\t'
              << left_width_arg << '\t'
              << repetitions << '\t' << setup_ms << '\t' << dp_ms << '\t'
              << profile.root << '\t' << profile.producer_candidates << '\t'
              << profile.producer_selected << '\t' << profile.endpoint_pairs
              << '\t' << profile.sources << '\t' << profile.tails << '\t'
              << profile.tails_selected << '\t' << profile.left_candidates
              << '\t' << profile.left_selected << '\t'
              << profile.rejected_sources << '\t'
              << profile.active_expiry_probes << '\t'
              << profile.tail_materialization_probes << '\t'
              << profile.prepared_score_probes << '\t'
              << profile.max_active_segments
              << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "scfg2_vp_case3_online_profile_tool failed: "
              << error.what() << '\n';
    return 1;
  }
}
