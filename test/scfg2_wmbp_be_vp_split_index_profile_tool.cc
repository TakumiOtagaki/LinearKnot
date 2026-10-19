#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/wmbp_be_vp_split_child_index.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    std::cerr << "usage: scfg2_wmbp_be_vp_split_index_profile_tool "
                 "SEQUENCE SCAFFOLD [PARAM_FILE] [REPEATS]\n";
    return EXIT_FAILURE;
  }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE",
         argc >= 4 ? argv[3] : "params/rna_Turner04.par",
         1);
  const std::string sequence = argv[1];
  const std::string scaffold = argv[2];
  const int repeats = argc == 5 ? std::max(1, std::atoi(argv[4])) : 20;
  assert(sequence.size() == scaffold.size());

  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto schedule = cparty::scfg2::engine::
      build_pkfree_basic_exact_right_boundary_schedule(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  std::vector<cparty::scfg2::deductions::ItemKey> parents;
  for (const auto &item : schedule) {
    if (cparty::scfg2::deductions::exact_wmbp_item_active(item, ctx)) {
      parents.push_back(item);
    }
  }

  const auto build_begin = std::chrono::steady_clock::now();
  const auto index = cparty::scfg2::beam::
      ExactWmbpBeVpSplitChildIndex::build(
          static_cast<cand_pos_t>(sequence.size()), ctx);
  const auto build_end = std::chrono::steady_clock::now();

  std::size_t coordinate_visits = 0;
  std::size_t coordinate_geometry_visits = 0;
  std::size_t indexed_visits = 0;
  std::size_t bucket_visits = 0;
  for (const auto &parent : parents) {
    const auto coordinate = cparty::scfg2::deductions::
        scan_exact_wmbp_be_vp_split_geometries(parent, ctx, true, true);
    const auto indexed = cparty::scfg2::beam::
        scan_exact_wmbp_be_vp_split_geometries_from_index(
            parent, ctx, index, true, true);
    assert(cparty::scfg2::beam::same_wmbp_be_vp_split_geometry_scans(
        coordinate, indexed.geometry));
    coordinate_visits += coordinate.l_visits;
    coordinate_geometry_visits += coordinate.geometry_visits;
    indexed_visits += indexed.geometry.l_visits;
    bucket_visits += indexed.parent_bucket_visits;
  }

  std::size_t coordinate_checksum = 0;
  const auto coordinate_begin = std::chrono::steady_clock::now();
  for (int repeat = 0; repeat < repeats; ++repeat) {
    for (const auto &parent : parents) {
      const auto scan = cparty::scfg2::deductions::
          scan_exact_wmbp_be_vp_split_geometries(parent, ctx, true, true);
      coordinate_checksum += scan.shared.size() + scan.wi.size();
    }
  }
  const auto coordinate_end = std::chrono::steady_clock::now();

  std::size_t indexed_checksum = 0;
  const auto indexed_begin = std::chrono::steady_clock::now();
  for (int repeat = 0; repeat < repeats; ++repeat) {
    for (const auto &parent : parents) {
      const auto scan = cparty::scfg2::beam::
          scan_exact_wmbp_be_vp_split_geometries_from_index(
              parent, ctx, index, true, true);
      indexed_checksum +=
          scan.geometry.shared.size() + scan.geometry.wi.size();
    }
  }
  const auto indexed_end = std::chrono::steady_clock::now();
  assert(coordinate_checksum == indexed_checksum);

  const auto milliseconds = [](auto begin, auto end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
  };
  const double coordinate_ms =
      milliseconds(coordinate_begin, coordinate_end) / repeats;
  const double indexed_ms = milliseconds(indexed_begin, indexed_end) / repeats;
  std::cout << std::fixed << std::setprecision(6)
            << "n\tparents\tcoordinate_visits\tcoordinate_geometry_visits"
               "\tindexed_visits"
               "\tbucket_visits\tindex_build_ms\tcoordinate_ms"
               "\tindexed_ms\tindexed_over_coordinate\tparity\n"
            << sequence.size() << '\t' << parents.size() << '\t'
            << coordinate_visits << '\t' << coordinate_geometry_visits << '\t'
            << indexed_visits << '\t'
            << bucket_visits << '\t' << milliseconds(build_begin, build_end)
            << '\t' << coordinate_ms << '\t' << indexed_ms << '\t'
            << (coordinate_ms == 0.0 ? 0.0 : indexed_ms / coordinate_ms)
            << "\t1\n";
}
