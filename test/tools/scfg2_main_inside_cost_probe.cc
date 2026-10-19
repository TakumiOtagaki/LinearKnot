#include "scfg2/weights/runtime_snapshot.hh"
// Main-only timings and canonical chart/counter snapshots. Guide preparation
// is separate. Define SCFG2_COST_GPROF with -pg to sample only the main calls.
#include "scfg2_pilot_guided_inside.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#ifdef SCFG2_COST_GPROF
extern "C" void moncontrol(int);
void profiling(bool enabled) { moncontrol(enabled); }
#else
void profiling(bool) {}
#endif

using Run = sc::replay::WFinalVpWmbpOnlineRun<double>;

std::string counters(const Run &run) {
  std::ostringstream out;
  const auto family = [&](const char *name, const auto &counts) {
    for (std::size_t f = 0; f < counts.size(); ++f)
      out << name << '.' << sc::grammar::nonterminal_name(static_cast<NT>(f))
          << '\t' << counts[f] << '\n';
  };
  family("visits", run.item_visits);
  family("positive", run.positive_items);
  family("frontier", run.frontier_selected);
#define COUNT(name) out << #name << '\t' << run.name << '\n'
  COUNT(vp_selected); COUNT(wmbp_selected); COUNT(case3_producer_selected);
  COUNT(case3_tail_selected); COUNT(case3_left_selected); COUNT(be_right_tail_selected);
  COUNT(be_right_shell_selected); COUNT(wmbp_prefix_selected);
  COUNT(vm_base_prefix_selected); COUNT(wm_start_prefix_selected);
  COUNT(wip_start_prefix_selected); COUNT(vpl_prefix_selected); COUNT(vpr_base_state_selected);
  COUNT(wmb_btail_selected); COUNT(wmb_retained_split_selected);
  COUNT(wmbp_shared_emitted_sources); COUNT(wmbp_wi_prefix_candidates);
  COUNT(wmbp_wi_emitted_sources); COUNT(vp_stack_event_sources); COUNT(vp_internal_event_sources);
  COUNT(vp_case3_nonleft_emitted_sources); COUNT(vp_case12_emitted_sources);
  COUNT(be_stack_sources); COUNT(be_core_left_sources); COUNT(be_right_shell_sources);
  COUNT(v_hairpin_sources); COUNT(v_internal_sources); COUNT(v_vm_sources);
#undef COUNT
  return out.str();
}

void write_chart(const std::string &path, const Chart &chart) {
  std::vector<std::pair<Item, double>> ordered(chart.values().begin(), chart.values().end());
  std::sort(ordered.begin(), ordered.end(), [](const auto &a, const auto &b) {
    if (a.first.nonterminal != b.first.nonterminal)
      return a.first.nonterminal < b.first.nonterminal;
    return less_key(a.first, b.first);
  });
  std::ofstream out(path);
  if (!out) throw std::runtime_error("cannot write chart snapshot");
  out << "family\ti\tj\tip\tjp\talpha_hex\n" << std::hexfloat;
  for (const auto &[key, value] : ordered)
    out << static_cast<int>(key.nonterminal) << '\t' << key.i << '\t' << key.j
        << '\t' << key.ip << '\t' << key.jp << '\t' << value << '\n';
}

int main(int argc, char **argv) {
  profiling(false);
  try {
    if (argc != 9) throw std::invalid_argument(
        "usage: main-cost INPUT given|empty PREFIX WIDTH baseline|0.5|2/3 PARAM REPEATS weak|legacy");
    std::ifstream input(argv[1]);
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
      throw std::invalid_argument("invalid input");
    const auto b = std::stoull(argv[4]);
    sqrt_pilot_width(b);
    const std::string ranking = argv[5], policy = argv[8];
    if (ranking != "baseline" && ranking != "0.5" && ranking != "2/3")
      throw std::invalid_argument("invalid ranking");
    if (policy != "weak" && policy != "legacy") throw std::invalid_argument("invalid closure policy");
    const int repeats = std::stoi(argv[7]);
    if (repeats < 1) throw std::invalid_argument("positive repeats required");
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[6], 1);
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    sc::replay::WFinalLogWeightContext ctx(owner,
        policy == "weak" ? sc::WmbClosurePolicy::RequireWeakClosure : sc::WmbClosurePolicy::LegacyCompatible);
    const auto n = static_cast<cand_pos_t>(sequence.size());
    const Item root{NT::W, 1, n, -1, -1};
    auto widths = sc::replay::WFinalVpWmbpOnlineWidths::uniform(b);
    widths.forced_items = sc::replay::scaffold_only_derivation_items(ctx, n, widths);
    Chart pilot(zero), guide(zero);
    std::size_t p = 0;
    if (ranking != "baseline") {
      p = pilot_width_for(b, ranking == "0.5" ? PilotExponent::Half : PilotExponent::TwoThirds);
      auto pilot_widths = sc::replay::WFinalVpWmbpOnlineWidths::uniform(p);
      pilot_widths.forced_items = widths.forced_items;
      pilot = sc::replay::run_w_final_vp_wmbp_online<Sum>(ctx, n, pilot_widths).chart;
      guide = bounded_outer(pilot, root, ctx, p, 4 * p).outer.values;
    }
    PilotRankingCounts counts;
    sc::replay::WFinalPublicItemRanking rank;
    if (ranking != "baseline") rank = PilotOuterRanking{guide, pilot, counts};
    const std::string prefix = argv[3];
    std::ofstream times(prefix + ".times.tsv");
    if (!times) throw std::runtime_error("cannot write timings");
    times << std::setprecision(17)
          << "n\twidth\tranking\tpilot_width\trepeat\tmain_seconds\tlog_z\titems\n";
    Chart first(zero);
    std::string first_counts;
    for (int repeat = 0; repeat < repeats; ++repeat) {
      profiling(true);
      const auto start = Clock::now();
      auto run = sc::replay::run_w_final_vp_wmbp_online_ranked_inside(ctx, n, widths, rank);
      const double elapsed = seconds(start);
      profiling(false);
      if (!std::isfinite(run.chart.get(root))) throw std::runtime_error("unreachable root");
      const auto observed_counts = counters(run);
      if (repeat == 0) {
        write_chart(prefix + ".chart.tsv", run.chart);
        std::ofstream(prefix + ".counts.tsv") << observed_counts;
        first = run.chart;
        first_counts = observed_counts;
      } else {
        if (run.chart.values().size() != first.values().size() || observed_counts != first_counts)
          throw std::runtime_error("repeated chart/counter support changed");
        for (const auto &[key, value] : first.values())
          if (value != run.chart.get(key)) throw std::runtime_error("repeated chart changed bitwise");
      }
      times << n << '\t' << b << '\t' << ranking << '\t' << p << '\t' << repeat
            << '\t' << elapsed << '\t' << run.chart.get(root) << '\t' << run.chart.values().size() << '\n';
    }
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    profiling(false);
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
