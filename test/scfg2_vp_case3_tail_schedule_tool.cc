#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_tail_cache.hh"
#include "scfg2/beam/vp_case3_tail_schedule_profile.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "usage: scfg2_vp_case3_tail_schedule_tool <seq> <G> "
                 "[param_file] [span|right]\n";
    return 2;
  }
  try {
    std::string sequence = argv[1];
    std::string scaffold = argv[2];
    if (sequence.size() != scaffold.size()) {
      throw std::runtime_error("sequence and scaffold lengths differ");
    }
    if (argc >= 4) {
      setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[3], 1);
    } else if (!std::getenv("CPARTY_FIXED_ENERGY_PARAM_FILE")) {
      setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    }

    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
    const bool right_boundary = argc == 5 && std::string(argv[4]) == "right";
    if (argc == 5 && !right_boundary && std::string(argv[4]) != "span") {
      throw std::runtime_error("schedule must be span or right");
    }
    const auto schedule = right_boundary
        ? cparty::scfg2::engine::
              build_pkfree_basic_exact_right_boundary_schedule(
                  static_cast<cand_pos_t>(sequence.size()), ctx)
        : cparty::scfg2::engine::build_pkfree_basic_exact_schedule(
              static_cast<cand_pos_t>(sequence.size()), ctx);
    cparty::scfg2::beam::ExactVpCase3TailCache cache;
    (void)cparty::scfg2::replay::
        run_w_final_exact_dp_schedule_with_contribution<
            cparty::scfg2::semirings::SumProductSemiring>(
            schedule,
            ctx,
            [](const cparty::scfg2::deductions::Deduction &) { return true; },
            [](const auto &, const auto &, double, const auto &, double) {},
            [&](const cparty::scfg2::deductions::Deduction &deduction,
                const cparty::scfg2::engine::Chart<double> &chart,
                double local) {
              const auto reconstructed =
                  cache.observe_left_source(deduction, chart, local);
              if (reconstructed.has_value()) return *reconstructed;
              return cparty::scfg2::replay::runtime_deduction_contribution<
                  cparty::scfg2::semirings::SumProductSemiring>(
                  chart, deduction, local);
            });
    const auto profile =
        cparty::scfg2::beam::profile_vp_case3_tail_schedule(cache, schedule);
    std::cout << "n\ttails\tattachments\tmax_attachments\tmaterialized"
                 "\treused\tcolumns\ttopological_tails"
                 "\tcommittable_columns\tnoncommittable_columns"
                 "\tmissing_schedule_items\tmax_overlap\n"
              << sequence.size() << '\t' << cache.tail_count() << '\t'
              << cache.attachment_count() << '\t'
              << cache.max_attachments_per_tail() << '\t'
              << cache.materialization_count() << '\t' << cache.reuse_count()
              << '\t' << profile.column_count << '\t'
              << profile.topological_tail_count << '\t'
              << profile.committable_column_count << '\t'
              << profile.noncommittable_column_count << '\t'
              << profile.missing_schedule_item_count << '\t'
              << profile.maximum_column_overlap << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "scfg2_vp_case3_tail_schedule_tool failed: " << error.what()
              << '\n';
    return 1;
  }
}
