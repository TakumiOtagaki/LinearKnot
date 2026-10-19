// Observation-only numeric forward probe. TSV goes to stdout, result to stderr.
#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/beam_detail.hh"
#include "scfg2/api/beam.hh"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#if defined(__GLIBC__)
#include <malloc.h>
#endif

namespace sc = cparty::scfg2;
namespace {
struct ProcessMemory {
  std::size_t rss = 0, hwm = 0, arena = 0, allocated = 0, free = 0, mapped = 0;
};
ProcessMemory process_memory() {
  ProcessMemory out;
#if defined(__linux__)
  std::ifstream status("/proc/self/status");
  std::string line;
  while (std::getline(status, line)) {
    std::istringstream row(line);
    std::string key;
    std::size_t value;
    if (!(row >> key >> value)) continue;
    if (key == "VmRSS:") out.rss = value;
    if (key == "VmHWM:") out.hwm = value;
  }
#endif
#if defined(__GLIBC__)
#if __GLIBC_PREREQ(2, 33)
  const auto info = mallinfo2();
  out.arena = info.arena;
  out.allocated = info.uordblks;
  out.free = info.fordblks;
  out.mapped = info.hblkhd;
#endif
#endif
  return out;
}
void observe(const sc::beam::MemoryObservation &observation) {
  const auto memory = process_memory();
  const auto prefix = [&](const char *kind, const std::string &name) {
    std::cout << sc::beam::memory_phase_name(observation.phase) << '\t'
        << observation.event << '\t' << observation.width << '\t'
        << observation.column << '\t' << kind << '\t' << name;
  };
  const auto suffix = [&] {
    std::cout << '\t' << memory.rss << '\t' << memory.hwm << '\t'
        << memory.arena << '\t' << memory.allocated << '\t' << memory.free
        << '\t' << memory.mapped << '\n';
  };
  prefix("process", "process");
  std::cout << '\t' << observation.chart_items << "\t0\t0\t0\t0";
  suffix();
  for (const auto &container : observation.containers) {
    prefix("container", container.name);
    std::cout << '\t' << container.entries << '\t' << container.buckets
        << '\t' << container.value_bytes << '\t' << container.vector_size_bytes
        << '\t' << container.vector_capacity_bytes;
    suffix();
  }
  for (const auto &family : observation.families) {
    prefix("family", std::string(sc::grammar::nonterminal_name(family.nonterminal)));
    std::cout << '\t' << family.live_items << "\t0\t0\t0\t0";
    suffix();
  }
  std::cout.flush();
}
} // namespace

int main(int argc, char **argv) {
  if (argc < 3 || argc > 8) {
    std::cerr << "usage: " << argv[0]
        << " INPUT WIDTH [baseline|outer] [HAIRPIN_CAP] [CASE3_DEPTH] [COLUMN_INTERVAL] [inside|mfe]\n";
    return 1;
  }
  try {
    std::ifstream input(argv[1]);
    std::string sequence, scaffold, line;
    while (std::getline(input, line)) {
      if (line.empty() || line[0] == '>') continue;
      if (sequence.empty()) sequence = line;
      else { scaffold = line; break; }
    }
    if (sequence.empty() || scaffold.size() != sequence.size())
      throw std::invalid_argument("expected sequence and equally long scaffold");
    sc::api::BeamConfig config;
    config.beam_size = std::stoull(argv[2]);
    if (argc > 3) {
      const std::string mode = argv[3];
      if (mode == "baseline") config.ranking = sc::api::BeamRanking::Baseline;
      else if (mode != "outer") throw std::invalid_argument("invalid ranking");
    }
    if (argc > 4) config.max_hairpin_unpaired = std::stoull(argv[4]);
    if (argc > 5) config.max_case3_ancestor_depth = std::stoull(argv[5]);
    config.memory_column_interval = argc > 6 ? std::stoull(argv[6]) : 256;
    std::string error;
    if (!sc::api::detail::validate_beam_config(config, &error))
      throw std::invalid_argument(error);
    config.memory_observer = observe;
    std::cout << "phase\tevent\twidth\tcolumn\tkind\tname\tentries\tbuckets"
        "\tvalue_bytes\tvector_size_bytes\tvector_capacity_bytes\trss_kib"
        "\thwm_kib\tarena_bytes\tallocated_bytes\tfree_bytes\tmapped_bytes\n";
    sc::beam::MemoryObservation start;
    start.event = "process_begin";
    observe(start);
    if (argc > 7) {
      const std::string operation = argv[7];
      if (operation == "mfe") {
        const auto result = sc::api::beam_viterbi(sequence, scaffold, config);
        if (!result.ok) throw std::runtime_error(result.error);
        std::cerr << std::setprecision(17) << "mfe_log_weight=" << result.log_weight
            << " mfe_energy=" << result.energy << '\n';
        return result.reachable ? 0 : 2;
      }
      if (operation != "inside") throw std::invalid_argument("invalid operation");
    }
    cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
    sc::replay::WFinalLogWeightContext context(owner, config.wmb_closure_policy);
    const auto n = static_cast<cand_pos_t>(sequence.size());
    const auto prepared = sc::api::detail::prepare_beam_ranking(context, n, config);
    const auto run = sc::replay::run_w_final_vp_wmbp_online<
        sc::semirings::LogSumExpSemiring>(context, n, prepared.widths, prepared.ranking);
    observe(sc::beam::chart_memory(sc::beam::MemoryPhase::Main, "returned",
        config.beam_size, n, run.chart));
    const double log_z = run.chart.get({sc::grammar::NonTerminal::W, 1, n, -1, -1});
    std::cerr << std::setprecision(17) << "log_partition_weight=" << log_z
        << " live_chart_items=" << run.chart.values().size() << '\n';
    return std::isfinite(log_z) ? 0 : 2;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
