#include "scfg2/weights/runtime_snapshot.hh"
// Stage clocks are diagnostic only; repeated throughput runs have them off.
#include "scfg2_pilot_guided_inside.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

void equal_outer(const Outer &a, const Outer &b) {
  if (a.values.values().size() != b.values.values().size() || a.links.size() != b.links.size())
    throw std::runtime_error("optimized/reference support mismatch");
  for (const auto &[key, value] : a.values.values())
    if (value != b.values.get(key)) throw std::runtime_error("optimized/reference score mismatch");
  for (const auto &[key, link] : a.links) {
    const auto &other = b.links.at(key);
    if (!sc::deductions::ItemKeyEq{}(link.parent, other.parent) || link.edge != other.edge)
      throw std::runtime_error("optimized/reference predecessor mismatch");
  }
}

// Outside the timed region: compare complete guides across frozen binaries,
// including each chosen predecessor, rather than checking only the root.
void write_outer(const std::string &path, const Outer &outer) {
  std::vector<std::pair<Item, double>> ordered(
      outer.values.values().begin(), outer.values.values().end());
  std::sort(ordered.begin(), ordered.end(), [](const auto &a, const auto &b) {
    if (a.first.nonterminal != b.first.nonterminal)
      return a.first.nonterminal < b.first.nonterminal;
    return less_key(a.first, b.first);
  });
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write outer snapshot");
  out << "family\ti\tj\tip\tjp\tvalue_hex\tparent_family\tparent_i\tparent_j"
         "\tparent_ip\tparent_jp\tedge_hex\n" << std::hexfloat;
  for (const auto &[key, value] : ordered) {
    out << static_cast<int>(key.nonterminal) << '\t' << key.i << '\t' << key.j
        << '\t' << key.ip << '\t' << key.jp << '\t' << value;
    const auto found = outer.links.find(key);
    if (found == outer.links.end()) out << "\t-1\t-1\t-1\t-1\t-1\t0x0p+0\n";
    else {
      const auto &[parent, edge] = found->second;
      out << '\t' << static_cast<int>(parent.nonterminal) << '\t' << parent.i
          << '\t' << parent.j << '\t' << parent.ip << '\t' << parent.jp
          << '\t' << edge << '\n';
    }
  }
}

int main(int argc, char **argv) {
  try {
    if (argc != 7) throw std::invalid_argument(
        "usage: cost-probe INPUT given|empty OUTPUT_PREFIX PARAM_FILE PILOT_WIDTHS REPEATS");
    std::ifstream input(argv[1]);
    std::string sequence, scaffold, line;
    while (std::getline(input, line)) {
      if (line.empty() || line.front() == '>') continue;
      if (sequence.empty()) sequence = line;
      else { scaffold = line; break; }
    }
    if (std::string(argv[2]) == "empty") scaffold.assign(sequence.size(), '.');
    else if (std::string(argv[2]) != "given") throw std::invalid_argument("invalid mode");
    if (sequence.empty() || sequence.size() != scaffold.size())
      throw std::invalid_argument("invalid input");
    const int repeats = std::stoi(argv[6]);
    if (repeats < 1) throw std::invalid_argument("positive repeats required");
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[4], 1);
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    sc::replay::WFinalLogWeightContext ctx(owner);
    const auto n = static_cast<cand_pos_t>(sequence.size());
    const Item root{NT::W, 1, n, -1, -1};
    const auto forced = sc::replay::scaffold_only_derivation_items(ctx, n,
        sc::replay::WFinalVpWmbpOnlineWidths::uniform(100));
    std::ofstream summary(std::string(argv[3]) + ".tsv");
    std::ofstream families(std::string(argv[3]) + ".families.tsv");
    if (!summary || !families) throw std::runtime_error("cannot write cost report");
    summary << std::setprecision(17);
    families << std::setprecision(17);
    summary << "n\tp\trepeats\tmedian_seconds\tmin_seconds\tmax_seconds\tprofile_seconds"
               "\tpreparation_seconds\tproposal_seconds\tgeneration_seconds\tpropagation_seconds"
               "\tfinalization_seconds\tcleanup_seconds\tparents\tdeductions\toffers\tweight_evaluations\tchild_alpha_lookups"
               "\tleaf_deductions\tzero_sibling_deductions\treference_median_seconds\tspeedup"
               "\treference_weight_evaluations\treference_child_alpha_lookups\n";
    families << "p\tfamily\tgeneration_seconds\tpropagation_seconds\n";
    std::stringstream widths(argv[5]);
    std::string token;
    while (std::getline(widths, token, ',')) {
      const auto p = std::stoull(token);
      if (!p) throw std::invalid_argument("positive pilot width required");
      auto widths = sc::replay::WFinalVpWmbpOnlineWidths::uniform(p);
      widths.forced_items = forced;
      const auto pilot = sc::replay::run_w_final_vp_wmbp_online<Sum>(ctx, n, widths);
      const auto profile_start = Clock::now();
      const auto profiled = bounded_outer(pilot.chart, root, ctx, p, 4 * p, false, true);
      const double profile_seconds = seconds(profile_start);
      const auto reference = bounded_outer<decltype(ctx), false>(pilot.chart, root, ctx, p, 4 * p);
      equal_outer(profiled.outer, reference.outer);
      write_outer(std::string(argv[3]) + ".p" + std::to_string(p) + ".outer.tsv", profiled.outer);
      std::vector<double> timings, reference_timings;
      for (int repeat = 0; repeat < repeats; ++repeat) {
        const auto optimized = [&]() {
          const auto start = Clock::now();
          const auto measured = bounded_outer(pilot.chart, root, ctx, p, 4 * p);
          timings.push_back(seconds(start));
        };
        const auto original = [&]() {
          const auto start = Clock::now();
          const auto measured = bounded_outer<decltype(ctx), false>(pilot.chart, root, ctx, p, 4 * p);
          reference_timings.push_back(seconds(start));
        };
        if (repeat % 2) { original(); optimized(); }
        else { optimized(); original(); }
      }
      std::sort(timings.begin(), timings.end());
      std::sort(reference_timings.begin(), reference_timings.end());
      const auto &w = profiled.work;
      summary << n << '\t' << p << '\t' << repeats << '\t' << timings[timings.size() / 2]
              << '\t' << timings.front() << '\t' << timings.back() << '\t' << profile_seconds
              << '\t' << w.preparation_seconds << '\t' << w.proposal_seconds
              << '\t' << w.generation_seconds << '\t' << w.propagation_seconds
              << '\t' << w.finalization_seconds << '\t' << profile_seconds - profiled.outer.elapsed
              << '\t' << profiled.outer.visited
              << '\t' << profiled.outer.generated << '\t' << profiled.outer.selected
              << '\t' << w.weight_evaluations << '\t' << w.child_alpha_lookups
              << '\t' << w.leaf_deductions << '\t' << w.zero_sibling_deductions
              << '\t' << reference_timings[reference_timings.size() / 2]
              << '\t' << reference_timings[reference_timings.size() / 2] / timings[timings.size() / 2]
              << '\t' << reference.work.weight_evaluations << '\t' << reference.work.child_alpha_lookups << '\n';
      for (std::size_t f = 0; f < outer_probe_families; ++f)
        families << p << '\t' << sc::grammar::nonterminal_name(static_cast<NT>(f))
                 << '\t' << w.family_generation[f] << '\t' << w.family_propagation[f] << '\n';
      std::cerr << "p=" << p << " outer_median=" << timings[timings.size() / 2] << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
