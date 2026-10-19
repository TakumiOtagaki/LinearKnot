#include "scfg2/weights/runtime_snapshot.hh"
// Actual finite-beam inside comparison. Exact inside is an optional scorer
// only; no exact alpha or probability outside is supplied to the guide.
#include "scfg2_pilot_guided_inside.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

std::vector<std::size_t> widths_from(const std::string &csv) {
  std::stringstream stream(csv);
  std::string token;
  std::vector<std::size_t> widths;
  while (std::getline(stream, token, ',')) {
    if (token.empty() || token.front() == '-') throw std::invalid_argument("invalid width");
    std::size_t parsed = 0;
    const auto width = std::stoull(token, &parsed);
    if (parsed != token.size()) throw std::invalid_argument("invalid width");
    sqrt_pilot_width(width);
    widths.push_back(width);
  }
  if (widths.empty()) throw std::invalid_argument("empty widths");
  return widths;
}
std::ofstream output_file(const std::string &path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write " + path);
  out << std::setprecision(17);
  return out;
}

std::vector<PilotExponent> exponents_from(const std::string &csv) {
  std::stringstream stream(csv);
  std::string token;
  std::vector<PilotExponent> result;
  while (std::getline(stream, token, ',')) {
    PilotExponent q;
    if (token == "0.5" || token == "1/2") q = PilotExponent::Half;
    else if (token == "2/3") q = PilotExponent::TwoThirds;
    else throw std::invalid_argument("supported pilot exponents: 0.5,2/3");
    if (std::find(result.begin(), result.end(), q) != result.end())
      throw std::invalid_argument("duplicate pilot exponent");
    result.push_back(q);
  }
  if (result.empty()) throw std::invalid_argument("empty pilot exponents");
  return result;
}

int main(int argc, char **argv) {
  try {
    if (argc < 7 || argc > 11) throw std::invalid_argument(
        "usage: probe INPUT given|empty OUTPUT_PREFIX WIDTHS exact|no-exact PARAM_FILE [OUTER_FACTOR [weak|legacy [Q_LIST [REPEATS]]]]");
    std::ifstream input(argv[1]);
    if (!input) throw std::runtime_error("cannot open input");
    std::string sequence, scaffold, line;
    while (std::getline(input, line)) {
      if (!line.empty() && line.back() == '\r') line.pop_back();
      if (line.empty() || line.front() == '>') continue;
      if (sequence.empty()) sequence = line;
      else { scaffold = line; break; }
    }
    if (std::string(argv[2]) == "empty") scaffold.assign(sequence.size(), '.');
    else if (std::string(argv[2]) != "given") throw std::invalid_argument("invalid scaffold mode");
    if (sequence.empty() || sequence.size() != scaffold.size())
      throw std::invalid_argument("expected sequence and scaffold of equal nonzero length");
    const auto widths = widths_from(argv[4]);
    const bool score_exact = std::string(argv[5]) == "exact";
    if (!score_exact && std::string(argv[5]) != "no-exact") throw std::invalid_argument("invalid exact mode");
    const auto factor = argc >= 8 ? std::stoull(argv[7]) : 4;
    const std::string policy = argc >= 9 ? argv[8] : "weak";
    if (policy != "weak" && policy != "legacy") throw std::invalid_argument("invalid closure policy");
    const auto exponents = exponents_from(argc >= 10 ? argv[9] : "2/3");
    const int repeats = argc >= 11 ? std::stoi(argv[10]) : 1;
    if (repeats < 1) throw std::invalid_argument("positive repeats required");
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[6], 1);
    const auto prep_start = Clock::now();
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    sc::replay::WFinalLogWeightContext ctx(owner,
        policy == "weak" ? sc::WmbClosurePolicy::RequireWeakClosure : sc::WmbClosurePolicy::LegacyCompatible);
    const double context_seconds = seconds(prep_start);
    const auto n = static_cast<cand_pos_t>(sequence.size());
    const Item root{NT::W, 1, n, -1, -1};
    const auto force_start = Clock::now();
    const auto forced = sc::replay::scaffold_only_derivation_items(ctx, n,
        sc::replay::WFinalVpWmbpOnlineWidths::uniform(*std::max_element(widths.begin(), widths.end())));
    const double force_seconds = seconds(force_start);
    auto summary = output_file(std::string(argv[3]) + ".summary.tsv");
    auto families = output_file(std::string(argv[3]) + ".families.tsv");
    summary << "n\twidth\tpilot_width\touter_factor\tclosure\tbaseline_log_z\tpilot_log_z\tguided_log_z"
               "\texact_log_z\tbaseline_log_gap\tguided_log_gap\tbaseline_mass_ratio\tguided_mass_ratio"
               "\tcontext_seconds\tforce_seconds\tbaseline_seconds\tpilot_seconds\touter_seconds"
               "\tmain_seconds\tguided_seconds\tbaseline_total_seconds\tguided_total_seconds"
               "\tbaseline_items\tpilot_items\tmain_items\touter_items\tcovered_candidates"
               "\tfallback_candidates\tnovel_candidates\tproposals\tmax_proposals\tgenerated_deductions"
               "\tpeak_guide_states\tvalidation_checks\tpilot_exponent\trepeat\tphase_seconds_sum"
               "\tweight_evaluations\tchild_alpha_lookups\tzero_sibling_deductions\n";
    families << "width\tpilot_width\tfamily\tcovered_candidates\tfallback_candidates\tnovel_candidates\tpilot_exponent\trepeat\n";
    struct Row {
      std::size_t b;
      double baseline, guided, pilot, baseline_seconds, pilot_seconds, outer_seconds, main_seconds;
      std::size_t baseline_items, pilot_items, main_items, outer_items;
      PilotRankingCounts ranks;
      bounded_outer_detail::Work work;
      std::size_t generated, peak_states;
      PilotExponent exponent;
      int repeat;
      double pipeline_seconds;
    };
    std::vector<Row> rows;
    for (const auto b : widths) {
      auto main_widths = sc::replay::WFinalVpWmbpOnlineWidths::uniform(b);
      main_widths.forced_items = forced;
      for (int repeat = 0; repeat < repeats; ++repeat) {
        const auto baseline_start = Clock::now();
        auto baseline = sc::replay::run_w_final_vp_wmbp_online<Sum>(ctx, n, main_widths);
        const double baseline_seconds = seconds(baseline_start);
        // Rotate order to avoid always giving the largest pilot the warmest run.
        for (std::size_t offset = 0; offset < exponents.size(); ++offset) {
          const auto q = exponents[(offset + repeat) % exponents.size()];
          const auto pipeline_start = Clock::now();
          auto guided = run_pilot_guided_inside(ctx, n, main_widths, b, factor, q);
          const double pipeline_seconds = seconds(pipeline_start);
          const double z0 = baseline.chart.get(root), z1 = guided.inside.chart.get(root);
          if (!std::isfinite(z0) || !std::isfinite(z1)) throw std::runtime_error("unreachable finite root");
          rows.push_back({b, z0, z1, guided.pilot.get(root), baseline_seconds,
              guided.pilot_seconds, guided.guide.outer.elapsed, guided.main_seconds,
              baseline.chart.values().size(), guided.pilot.values().size(),
              guided.inside.chart.values().size(), guided.guide.outer.values.values().size(),
              guided.ranks, guided.guide.work, guided.guide.outer.generated, guided.guide.outer.peak_guide_states,
              q, repeat, pipeline_seconds});
          for (std::size_t f = 0; f < outer_probe_families; ++f)
            families << b << '\t' << guided.pilot_width << '\t'
                     << sc::grammar::nonterminal_name(static_cast<NT>(f)) << '\t'
                     << guided.ranks.covered[f] << '\t' << guided.ranks.fallback[f] << '\t'
                     << guided.ranks.novel[f] << '\t' << pilot_exponent_name(q) << '\t' << repeat << '\n';
          std::cerr << "b=" << b << " q=" << pilot_exponent_name(q) << " repeat=" << repeat
                    << " p=" << guided.pilot_width << " baseline=" << z0
                    << " guided=" << z1 << " pilot_s=" << guided.pilot_seconds
                    << " outer_s=" << guided.guide.outer.elapsed << " main_s=" << guided.main_seconds << '\n';
        }
      }
    }
    // Diagnostic reference follows all timed beam runs and never enters them.
    double exact_z = std::numeric_limits<double>::quiet_NaN();
    if (score_exact) {
      const auto exact = sc::replay::run_w_final_exact_dp_schedule<Sum>(
          sc::engine::build_pkfree_basic_exact_schedule(n, ctx), ctx,
          [](const auto &) { return true; },
          [](const auto &, const auto &, auto, const auto &, auto) {});
      exact_z = exact.get(root);
      if (!std::isfinite(exact_z)) throw std::runtime_error("unreachable exact root");
    }
    for (const auto &row : rows) {
      if (score_exact && (row.baseline > exact_z + 1e-8 || row.guided > exact_z + 1e-8 ||
                          row.pilot > exact_z + 1e-8))
        throw std::runtime_error("finite root exceeds exact root");
      const auto phase_seconds = row.pilot_seconds + row.outer_seconds + row.main_seconds;
      const auto guided_seconds = row.pipeline_seconds;
      summary << n << '\t' << row.b << '\t' << pilot_width_for(row.b, row.exponent) << '\t' << factor << '\t' << policy
              << '\t' << row.baseline << '\t' << row.pilot << '\t' << row.guided << '\t' << exact_z
              << '\t' << exact_z - row.baseline << '\t' << exact_z - row.guided
              << '\t' << std::exp(row.baseline - exact_z) << '\t' << std::exp(row.guided - exact_z)
              << '\t' << context_seconds << '\t' << force_seconds << '\t' << row.baseline_seconds
              << '\t' << row.pilot_seconds << '\t' << row.outer_seconds << '\t' << row.main_seconds
              << '\t' << guided_seconds << '\t' << context_seconds + force_seconds + row.baseline_seconds
              << '\t' << context_seconds + force_seconds + guided_seconds
              << '\t' << row.baseline_items << '\t' << row.pilot_items << '\t' << row.main_items
              << '\t' << row.outer_items << '\t' << row.ranks.total(row.ranks.covered)
              << '\t' << row.ranks.total(row.ranks.fallback) << '\t' << row.ranks.total(row.ranks.novel)
              << '\t' << row.work.proposals << '\t' << row.work.max_proposals << '\t' << row.generated
              << '\t' << row.peak_states << '\t' << row.work.checks
              << '\t' << pilot_exponent_name(row.exponent) << '\t' << row.repeat << '\t' << phase_seconds
              << '\t' << row.work.weight_evaluations << '\t' << row.work.child_alpha_lookups
              << '\t' << row.work.zero_sibling_deductions << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
