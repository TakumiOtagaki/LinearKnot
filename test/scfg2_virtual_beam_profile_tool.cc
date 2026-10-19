#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/source_selection.hh"
#include "scfg2/beam/virtual_profile.hh"
#include "scfg2/beam/vp_case3_factor.hh"
#include "scfg2/deductions/generate_exact_basic.hh"
#include "scfg2/engine/basic_schedule.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using cparty::scfg2::beam::VirtualBeamBucketKey;
using cparty::scfg2::beam::VirtualBeamProfile;
using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;

void print_usage() {
  std::cerr << "usage: scfg2_virtual_beam_profile_tool <seq> <G> "
               "[widths=10,20,50,100,200] [param_file]\n";
}

std::vector<std::size_t> parse_widths(const std::string &text) {
  std::vector<std::size_t> widths;
  std::istringstream input(text);
  std::string token;
  while (std::getline(input, token, ',')) {
    if (token.empty()) throw std::runtime_error("empty beam width");
    std::size_t consumed = 0;
    const unsigned long parsed = std::stoul(token, &consumed);
    if (consumed != token.size() || parsed == 0) {
      throw std::runtime_error("beam widths must be positive integers");
    }
    widths.push_back(static_cast<std::size_t>(parsed));
  }
  if (widths.empty()) throw std::runtime_error("no beam widths supplied");
  std::sort(widths.begin(), widths.end());
  widths.erase(std::unique(widths.begin(), widths.end()), widths.end());
  return widths;
}

VirtualBeamProfile build_profile(const std::string &sequence,
                                 const std::string &scaffold) {
  cparty::scfg2::weights::RuntimeSnapshot viterbi_owner(
      sequence, scaffold);
  const auto viterbi =
      cparty::scfg2::replay::run_w_final_exact_viterbi_traceback(
          viterbi_owner);
  std::vector<Deduction> selected_sources;
  cparty::scfg2::beam::collect_selected_source_deductions(
      viterbi.trace, &selected_sources);

  cparty::scfg2::weights::RuntimeSnapshot inside_owner(
      sequence, scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(inside_owner);
  const auto schedule = cparty::scfg2::engine::build_pkfree_basic_exact_schedule(
      static_cast<cand_pos_t>(sequence.size()), ctx);

  VirtualBeamProfile profile;
  (void)cparty::scfg2::replay::run_w_final_exact_dp_schedule<
      cparty::scfg2::semirings::SumProductSemiring>(
      schedule,
      ctx,
      [](const Deduction &) { return true; },
      [&](const ItemKey &,
          const Deduction &deduction,
          double contribution,
          const cparty::scfg2::engine::Chart<double> &chart,
          double) {
        const auto factors =
            cparty::scfg2::beam::vp_case3_left_factor_values(
                deduction, chart);
        cparty::scfg2::beam::VirtualBeamPriorityOverrides priorities;
        if (factors.available) {
          priorities.vp_case3_left_tail = factors.private_tail;
          priorities.vp_case3_left = factors.parent_facing_left;
        }
        profile.observe(
            deduction,
            contribution,
            cparty::scfg2::beam::source_deduction_is_selected(
                deduction, selected_sources),
            priorities);
      });
  return profile;
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3 || argc > 5) {
    print_usage();
    return 2;
  }

  try {
    const std::string sequence = argv[1];
    const std::string scaffold = argv[2];
    if (sequence.size() != scaffold.size()) {
      throw std::runtime_error("sequence and scaffold lengths differ");
    }
    const auto widths = parse_widths(argc >= 4 ? argv[3] : "10,20,50,100,200");
    if (argc == 5) {
      setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[4], 1);
    } else if (!std::getenv("CPARTY_FIXED_ENERGY_PARAM_FILE")) {
      setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    }

    const auto profile = build_profile(sequence, scaffold);
    std::vector<std::pair<VirtualBeamBucketKey,
                          const cparty::scfg2::beam::VirtualBeamBucketSummary *>>
        rows;
    rows.reserve(profile.bucket_count());
    for (const auto &[key, summary] : profile.buckets()) {
      if (!summary.full_source_identities_are_unique() ||
          !summary.proposed_identities_have_consistent_ranking_values()) {
        throw std::runtime_error(
            "virtual beam identity consistency check failed");
      }
      rows.emplace_back(key, &summary);
    }
    std::sort(rows.begin(), rows.end(), [](const auto &lhs, const auto &rhs) {
      if (lhs.first.right_boundary != rhs.first.right_boundary) {
        return lhs.first.right_boundary < rhs.first.right_boundary;
      }
      return static_cast<unsigned char>(lhs.first.family) <
             static_cast<unsigned char>(rhs.first.family);
    });

    std::cout << "n\tfamily\tright_boundary\twidth\tpriority_kind"
                 "\tcandidate_count"
                 "\tproposed_unique_count\tmax_sources_per_proposed_identity"
                 "\tparent_facing_count\tselected_source_count\tkept_count"
                 "\tselected_in_top\tkept_unique_count"
                 "\tselected_in_top_unique\tpriority_max\tsource_sum"
                 "\ttop_source_sum\tsource_coverage\ttop_unique_source_sum"
                 "\tunique_source_coverage\n";
    std::cout << std::setprecision(17);
    for (const auto &[key, summary] : rows) {
      for (const std::size_t width : widths) {
        const std::size_t kept = std::min(width, summary->candidate_count);
        const std::size_t kept_unique =
            std::min(width, summary->proposed_unique_identity_count());
        const double top_sum = summary->top_source_contribution_sum(width);
        const double top_unique_sum =
            summary->top_unique_source_contribution_sum(width);
        const double coverage = summary->source_contribution_sum > 0.0
                                    ? top_sum / summary->source_contribution_sum
                                    : 0.0;
        const double unique_coverage = summary->source_contribution_sum > 0.0
                                           ? top_unique_sum /
                                                 summary->source_contribution_sum
                                           : 0.0;
        std::cout << sequence.size() << '\t'
                  << cparty::scfg2::beam::virtual_beam_family_name(key.family)
                  << '\t' << key.right_boundary << '\t' << width << '\t'
                  << (summary->uses_factor_priority ? "factor" : "source")
                  << '\t'
                  << summary->candidate_count << '\t'
                  << summary->proposed_unique_identity_count() << '\t'
                  << summary->max_sources_per_proposed_identity() << '\t'
                  << summary->parent_facing_count << '\t'
                  << summary->selected_source_count << '\t' << kept << '\t'
                  << summary->selected_source_count_in_top(width) << '\t'
                  << kept_unique << '\t'
                  << summary->selected_source_count_in_top_unique(width) << '\t'
                  << summary->max_ranking_value << '\t'
                  << summary->source_contribution_sum << '\t' << top_sum << '\t'
                  << coverage << '\t' << top_unique_sum << '\t'
                  << unique_coverage << '\n';
      }
    }
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "scfg2_virtual_beam_profile_tool failed: " << error.what()
              << '\n';
    return 1;
  }
}
