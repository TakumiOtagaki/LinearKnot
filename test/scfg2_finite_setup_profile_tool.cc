#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/be_right_tail_retained_core_index.hh"
#include "scfg2/beam/vp_case12_last_child_index.hh"
#include "scfg2/beam/vp_case3_border_event_plan.hh"
#include "scfg2/beam/vp_case3_nonleft_last_child_index.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

int main(int argc, char **argv) {
  if (argc < 4 || argc > 6) {
    std::cerr << "usage: " << argv[0]
              << " SEQUENCE SCAFFOLD WIDTH [PARAM_FILE] [REPETITIONS]\n";
    return EXIT_FAILURE;
  }
  const std::string sequence = argv[1];
  const std::string scaffold = argv[2];
  const std::size_t width = std::stoull(argv[3]);
  const std::size_t repetitions =
      argc >= 6 ? std::stoull(argv[5]) : 1;
  if (sequence.size() != scaffold.size() || repetitions == 0) {
    std::cerr << "invalid sequence/scaffold length or repetitions\n";
    return EXIT_FAILURE;
  }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc >= 5 ? argv[4] : "params/rna_Turner04.par",
         1);

  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  double minimum_ms = std::numeric_limits<double>::infinity();
  std::size_t case3_left_paths = 0;
  std::size_t case3_nonleft_paths = 0;
  std::size_t case12_empty_summary_entries = 0;
  std::size_t be_parent_links = 0;

  for (std::size_t repetition = 0; repetition < repetitions; ++repetition) {
    const auto started = std::chrono::steady_clock::now();
    const auto case3_left =
        cparty::scfg2::beam::ExactVpCase3BorderEventChildIndex::build(
            n, ctx, 30);
    const auto case3_nonleft =
        cparty::scfg2::beam::VpCase3NonleftChildEventIndex::build(
            n, ctx, 30);
    const auto case12_nonempty =
        cparty::scfg2::beam::VpCase12DirectNonemptyLeftQuery::build(n, ctx);
    const auto case12_empty =
        cparty::scfg2::beam::VpCase12DirectEmptyLeftQuery::build(
            n, ctx, width);
    const auto be =
        cparty::scfg2::beam::BeRightTailRetainedCoreIndex::build(n, ctx);
    const auto stopped = std::chrono::steady_clock::now();
    minimum_ms = std::min(
        minimum_ms,
        std::chrono::duration<double, std::milli>(stopped - started).count());
    case3_left_paths = case3_left.path_entry_count();
    case3_nonleft_paths = case3_nonleft.path_entry_count();
    case12_empty_summary_entries =
        case12_empty.case1_summary_entry_count();
    be_parent_links = be.scaffold_parent_link_count();
    (void)case12_nonempty;
  }

  std::cout << std::setprecision(17)
            << "n\twidth\trepetitions\tsetup_min_ms\tcase3_left_paths"
               "\tcase3_nonleft_paths\tcase12_empty_summary_entries"
               "\tbe_parent_links\n"
            << n << '\t' << width << '\t' << repetitions << '\t'
            << minimum_ms << '\t' << case3_left_paths << '\t'
            << case3_nonleft_paths << '\t'
            << case12_empty_summary_entries << '\t' << be_parent_links
            << '\n';
  return EXIT_SUCCESS;
}
