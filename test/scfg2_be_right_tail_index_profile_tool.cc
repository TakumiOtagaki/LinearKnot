#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/be_right_tail_child_index.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "usage: scfg2_be_right_tail_index_profile_tool <seq> <G> "
                 "[repetitions] [param_file]\n";
    return 2;
  }
  try {
    std::string sequence = argv[1];
    std::string scaffold = argv[2];
    if (sequence.size() != scaffold.size()) {
      throw std::runtime_error("sequence and scaffold lengths differ");
    }
    const std::size_t repetitions =
        argc >= 4 ? static_cast<std::size_t>(std::stoull(argv[3])) : 10;
    if (repetitions == 0) {
      throw std::runtime_error("repetitions must be positive");
    }
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
           argc >= 5 ? argv[4] : "params/rna_Turner04.par",
           1);
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
    const cand_pos_t n = static_cast<cand_pos_t>(sequence.size());
    const auto schedule = cparty::scfg2::engine::
        build_pkfree_basic_exact_right_boundary_schedule(n, ctx);

    std::size_t exact_count = 0;
    const auto exact_begin = std::chrono::steady_clock::now();
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition) {
      const auto plan =
          cparty::scfg2::beam::ExactBeRightTailColumnPlan::build(schedule, ctx);
      const cparty::scfg2::beam::ExactBeRightTailChildIndex index(plan);
      exact_count += index.record_count();
    }
    const auto exact_end = std::chrono::steady_clock::now();

    std::size_t coordinate_count = 0;
    cparty::scfg2::beam::ExactBeRightTailChildIndex coordinate_last =
        cparty::scfg2::beam::ExactBeRightTailChildIndex::
            build_from_scaffold_coordinate_oracle(n, ctx);
    const auto coordinate_begin = std::chrono::steady_clock::now();
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition) {
      coordinate_last = cparty::scfg2::beam::ExactBeRightTailChildIndex::
          build_from_scaffold_coordinate_oracle(n, ctx);
      coordinate_count += coordinate_last.record_count();
    }
    const auto coordinate_end = std::chrono::steady_clock::now();

    std::size_t tree_count = 0;
    cparty::scfg2::beam::ExactBeRightTailChildIndex tree_last =
        cparty::scfg2::beam::ExactBeRightTailChildIndex::
            build_from_scaffold_pair_tree(n, ctx);
    const auto tree_begin = std::chrono::steady_clock::now();
    for (std::size_t repetition = 0; repetition < repetitions; ++repetition) {
      tree_last =
          cparty::scfg2::beam::ExactBeRightTailChildIndex::
              build_from_scaffold_pair_tree(n, ctx);
      tree_count += tree_last.record_count();
    }
    const auto tree_end = std::chrono::steady_clock::now();

    if (!tree_last.has_same_keys(coordinate_last) ||
        exact_count != coordinate_count || exact_count != tree_count) {
      throw std::runtime_error("index builders disagree");
    }
    const auto milliseconds_per = [&](auto begin, auto end) {
      return std::chrono::duration<double, std::milli>(end - begin).count() /
             static_cast<double>(repetitions);
    };
    std::cout << std::setprecision(10)
              << "n\trepetitions\ttails\texact_plan_ms\tcoordinate_ms"
                 "\tpair_tree_ms\tcoordinate_core_visits"
                 "\tcoordinate_anchor_visits\tpair_tree_ancestor_visits"
                 "\tpair_tree_anchor_visits\tparity\n"
              << n << '\t' << repetitions << '\t'
              << tree_last.record_count() << '\t'
              << milliseconds_per(exact_begin, exact_end) << '\t'
              << milliseconds_per(coordinate_begin, coordinate_end) << '\t'
              << milliseconds_per(tree_begin, tree_end) << '\t'
              << coordinate_last.coordinate_core_visit_count() << '\t'
              << coordinate_last.coordinate_anchor_visit_count() << '\t'
              << tree_last.core_outer_visit_count() << '\t'
              << tree_last.anchor_visit_count() << "\t1\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "scfg2_be_right_tail_index_profile_tool failed: "
              << error.what() << '\n';
    return 1;
  }
}
