#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/vp_case3_tail_border_profile.hh"
#include "scfg2/beam/vp_case3_border_event_plan.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 4) {
    std::cerr << "usage: scfg2_vp_case3_tail_border_profile_tool <seq> <G> "
                 "[param_file]\n";
    return 2;
  }
  try {
    std::string sequence = argv[1];
    std::string scaffold = argv[2];
    if (sequence.size() != scaffold.size()) {
      throw std::runtime_error("sequence and scaffold lengths differ");
    }
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
           argc == 4 ? argv[3] : "params/rna_Turner04.par",
           1);
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
    const auto plan =
        cparty::scfg2::beam::build_vp_case3_tail_column_plan_from_border_events(
            static_cast<cand_pos_t>(sequence.size()), ctx);
    const auto profile =
        cparty::scfg2::beam::profile_vp_case3_tail_borders(plan);
    std::cout << "n\tsources\ttails\tgeometries\ti_runs\tp_segments"
                 "\texact_segments"
                 "\tdiscontinuous_tails\tmulti_p_tails\tmax_runs"
                 "\tmax_geometries\n"
              << sequence.size() << '\t' << profile.source_count << '\t'
              << profile.tail_count << '\t' << profile.geometry_count << '\t'
              << profile.contiguous_i_run_count << '\t'
              << profile.constant_p_segment_count << '\t'
              << profile.exact_segment_count << '\t'
              << profile.discontinuous_tail_count << '\t'
              << profile.multi_p_tail_count << '\t'
              << profile.maximum_runs_per_tail << '\t'
              << profile.maximum_geometries_per_tail << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "scfg2_vp_case3_tail_border_profile_tool failed: "
              << error.what() << '\n';
    return 1;
  }
}
