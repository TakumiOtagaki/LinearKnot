#include "scfg2/api/hotspot_viterbi.hh"

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>

namespace {

std::size_t parse_size(const char *text, const char *name) {
  char *end = nullptr;
  const auto value = std::strtoull(text, &end, 10);
  if (!text[0] || !end || *end != '\0' || value == 0) {
    std::cerr << name << " must be a positive integer\n";
    std::exit(EXIT_FAILURE);
  }
  return static_cast<std::size_t>(value);
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 2 || argc > 9) {
    std::cerr
        << "usage: scfg2_hotspot_viterbi_tool SEQUENCE"
        << " [PARAM_FILE] [BEAM_WIDTH] [ANCHOR_CAP] [HOTSPOT_WIDTH]"
        << " [PASSES] [MAX_RESULTS] [DIVERSITY:none|deterministic]\n";
    return EXIT_FAILURE;
  }

  const std::string sequence = argv[1];
  const std::string parameter_file =
      argc >= 3 ? argv[2] : "params/rna_Turner04.par";
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", parameter_file.c_str(), 1);

  cparty::scfg2::api::BeamConfig beam;
  cparty::scfg2::proposal::WavefrontHotspotConfig proposal;
  if (argc >= 4) beam.beam_size = parse_size(argv[3], "beam width");
  if (argc >= 5) proposal.anchor_cap = parse_size(argv[4], "anchor cap");
  if (argc >= 6) {
    proposal.beam_width = parse_size(argv[5], "hotspot width");
  }
  if (argc >= 7) proposal.passes = parse_size(argv[6], "passes");
  if (argc >= 8) {
    proposal.max_results = parse_size(argv[7], "max results");
  }
  if (argc >= 9) {
    const std::string mode = argv[8];
    if (mode == "none") {
      proposal.diversity_mode =
          cparty::scfg2::proposal::HotspotDiversityMode::None;
    } else if (mode == "deterministic") {
      proposal.diversity_mode =
          cparty::scfg2::proposal::HotspotDiversityMode::
              DeterministicPassDiversification;
    } else {
      std::cerr << "diversity must be none or deterministic\n";
      return EXIT_FAILURE;
    }
  }

  const auto prediction = cparty::scfg2::api::predict_hotspot_viterbi(
      sequence, proposal, beam);
  std::size_t case3_truncated_records = 0;
  for (const auto &record : prediction.records) {
    if (record.viterbi.work.case3_truncated_endpoints > 0) {
      ++case3_truncated_records;
    }
  }
  std::cout << std::setprecision(12)
            << "SUMMARY"
            << "\tok=" << prediction.ok
            << "\treachable=" << prediction.reachable
            << "\tlength=" << sequence.size()
            << "\tbeam_width=" << beam.beam_size
            << "\tanchor_cap=" << proposal.anchor_cap
            << "\thotspot_width=" << proposal.beam_width
            << "\tpasses=" << proposal.passes
            << "\tmax_results=" << proposal.max_results
            << "\trecords=" << prediction.records.size()
            << "\tgroups=" << prediction.structures.size()
            << "\tproposal_seconds=" << prediction.proposal_elapsed_seconds
            << "\tviterbi_seconds=" << prediction.viterbi_elapsed_seconds
            << "\ttotal_seconds=" << prediction.total_elapsed_seconds
            << "\tcase3_truncated_records=" << case3_truncated_records
            << "\terror=" << prediction.error << '\n';

  std::cout << "RECORD_HEADER"
            << "\tindex\treachable\tscaffold_pairs\theuristic_score\tpass\torigin"
            << "\tenergy\tseconds\tscaffold_G\tgenerated_Gprime\tunion_R\n";
  for (std::size_t index = 0; index < prediction.records.size(); ++index) {
    const auto &record = prediction.records[index];
    std::cout << "RECORD"
              << '\t' << index
              << '\t' << record.viterbi.reachable
              << '\t' << record.proposal.scaffold_pairs.size()
              << '\t' << record.proposal.heuristic_score
              << '\t' << record.proposal.pass
              << '\t' << record.proposal.origin
              << '\t'
              << (record.viterbi.reachable
                      ? record.viterbi.energy
                      : std::numeric_limits<double>::infinity())
              << '\t' << record.elapsed_seconds
              << '\t' << record.scaffold_structure
              << '\t' << record.viterbi.generated_structure
              << '\t' << record.viterbi.union_structure << '\n';
    const auto &work = record.viterbi.work;
    std::cout << "WORK"
              << '\t' << index
              << "\titem_visits=" << work.total_item_visits
              << "\tpositive_items=" << work.total_positive_items
              << "\tfrontier_candidates=" << work.total_frontier_candidates
              << "\tfrontier_selected=" << work.total_frontier_selected
              << "\tvp_candidates=" << work.vp_candidates
              << "\tvp_selected=" << work.vp_selected
              << "\twmbp_candidates=" << work.wmbp_candidates
              << "\twmbp_selected=" << work.wmbp_selected
              << "\tcase3_candidates=" << work.case3_candidates
              << "\tcase3_selected=" << work.case3_selected
              << "\tbinary_join_probes=" << work.binary_join_probes
              << "\textension_probes=" << work.extension_probes
              << "\tv_parent_probes=" << work.v_parent_probes
              << "\tvp_loop_parent_probes=" << work.vp_loop_parent_probes
              << "\tcase3_truncated_endpoints="
              << work.case3_truncated_endpoints
              << '\n';
    for (const auto &family : work.families) {
      if (family.item_visits == 0 && family.positive_items == 0 &&
          family.frontier_candidates == 0 &&
          family.frontier_selected == 0) {
        continue;
      }
      std::cout << "FAMILY"
                << '\t' << index
                << '\t' << family.nonterminal
                << '\t' << family.item_visits
                << '\t' << family.positive_items
                << '\t' << family.frontier_candidates
                << '\t' << family.frontier_selected << '\n';
    }
  }

  std::cout << "GROUP_HEADER"
            << "\tindex\trecord_count\tminimum_energy\tmaximum_energy"
            << "\tenergy_spread\tselected_record\tselected\n";
  for (std::size_t index = 0; index < prediction.structures.size(); ++index) {
    const auto &group = prediction.structures[index];
    std::cout << "GROUP"
              << '\t' << index
              << '\t' << group.record_indices.size()
              << '\t' << group.minimum_energy
              << '\t' << group.maximum_energy
              << '\t' << group.energy_spread
              << '\t' << group.selected_record
              << '\t'
              << (prediction.reachable &&
                  index == prediction.selected_structure)
              << '\n';
  }
  return prediction.ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
