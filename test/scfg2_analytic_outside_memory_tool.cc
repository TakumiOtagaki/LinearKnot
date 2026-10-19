// Same-process-work measurement in separate invocations, so ru_maxrss is
// comparable across analytic/taped outside and profiles. Input is FASTA-like:
// optional >header, then one sequence line and one scaffold line.
#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/beam.hh"
#include "scfg2/engine/outside_marginals.hh"
#include "scfg2/profiles/beam_profiles.hh"

#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sys/resource.h>

namespace {
template <class Semiring, class Graph>
bool outside(std::string &sequence, std::string &scaffold,
             const cparty::scfg2::replay::WFinalVpWmbpOnlineWidths &widths,
             double &log_z, std::size_t &states, std::size_t &terms,
             std::size_t &factors, std::size_t &bytes, std::size_t &tape_nodes,
             const char *dump) {
  using namespace cparty::scfg2;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  Graph graph;
  const auto run = [&] {
    typename Semiring::Scope scope(graph);
    return replay::run_w_final_vp_wmbp_online<Semiring>(ctx, n, widths);
  }();
  const auto root = run.chart.get(deductions::ItemKey{grammar::NonTerminal::W, 1, n, -1, -1});
  log_z = root.log_value;
  if constexpr (std::is_same_v<Graph, semirings::LogSumExpTape>) {
    tape_nodes = graph.nodes().size();
  } else {
    states = graph.states().size(); terms = graph.terms().size();
    factors = graph.factors().size(); bytes = graph.storage_bytes();
  }
  if (!std::isfinite(log_z)) return false;
  const auto posterior = [&] {
    if constexpr (std::is_same_v<Graph, semirings::LogSumExpTape>) return graph.reverse(root.node);
    else return graph.posteriors(root.node);
  }();
  const auto marginal = engine::collect_outside_marginals(run.chart, posterior, n, ctx);
  if (!marginal.ok) std::cerr << marginal.error << '\n';
  if (dump) {
    std::map<std::pair<cand_pos_t, cand_pos_t>, std::array<double, 2>> pairs;
    for (const auto &[item, value] : run.chart.values()) {
      if (value.node >= posterior.size() || posterior[value.node] < 1e-15) continue;
      const bool round = item.nonterminal == grammar::NonTerminal::V;
      const bool square = item.nonterminal == grammar::NonTerminal::VP;
      if ((!round && !square) || ctx.pair_partner(item.i) == item.j) continue;
      pairs[{item.i, item.j}][round ? 0 : 1] += posterior[value.node];
    }
    std::ofstream output(dump);
    output << std::setprecision(17);
    for (const auto &[key, p] : pairs) output << key.first << '\t' << key.second
        << '\t' << p[0] << '\t' << p[1] << '\n';
  }
  return marginal.ok;
}
}  // namespace

int main(int argc, char **argv) {
  using namespace cparty::scfg2;
  if (argc != 4 && argc != 5) {
    std::cerr << "usage: " << argv[0] << " INPUT WIDTH inside|outside|graph-outside|taped-outside|profiles|taped-profiles [RAW_OUTPUT_DUMP]\n";
    return 1;
  }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  std::ifstream input(argv[1]);
  std::string sequence, scaffold, line;
  while (std::getline(input, line)) {
    if (line.empty() || line[0] == '>') continue;
    if (sequence.empty()) sequence = line;
    else { scaffold = line; break; }
  }
  if (sequence.empty() || scaffold.size() != sequence.size()) return 1;
  const std::size_t width = std::stoull(argv[2]);
  const std::string mode = argv[3];
  const auto start = std::chrono::steady_clock::now();
  double log_z = 0;
  std::size_t states = 0, terms = 0, factors = 0, bytes = 0, tape_nodes = 0;
  bool ok = false;
  const auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(width);
  if (mode == "profiles" || mode == "taped-profiles") {
    const auto result = mode == "profiles"
        ? profiles::beam_profiles(sequence, scaffold, widths, {1, 2, 4, 7})
        : profiles::beam_profiles_taped(sequence, scaffold, widths, {1, 2, 4, 7});
    log_z = result.log_partition_weight;
    states = result.reverse_records; terms = result.reverse_terms;
    factors = result.reverse_factors; bytes = result.reverse_storage_bytes;
    ok = result.ok;
    if (!ok) std::cerr << result.error << '\n';
    if (argc == 5 && ok) {
      std::ofstream output(argv[4]);
      output << std::setprecision(17);
      for (std::size_t i = 1; i < result.context.size(); ++i) {
        output << "context\t" << i;
        for (const double probability : result.context[i]) output << '\t' << probability;
        output << '\n';
      }
      for (std::size_t w = 0; w < result.accessibility.size(); ++w)
        for (std::size_t i = 1; i < result.accessibility[w].size(); ++i)
          output << "accessibility\t" << result.window_lengths[w] << '\t' << i
                 << '\t' << result.accessibility[w][i] << '\n';
    }
  } else if (mode == "inside") {
    const auto result = api::beam_inside(sequence, scaffold, api::BeamConfig{width, 30, 30});
    log_z = result.log_partition_weight; ok = result.ok && result.reachable;
  } else if (mode == "outside") {
    const auto result = api::beam_outside(sequence, scaffold, api::BeamConfig{width, 30, 30});
    log_z = result.log_partition_weight;
    states = result.reverse_records; terms = result.reverse_terms;
    factors = result.reverse_factors; bytes = result.reverse_storage_bytes;
    ok = result.ok && result.reachable;
    if (!ok) std::cerr << result.error << '\n';
    if (argc == 5 && ok) {
      std::ofstream output(argv[4]);
      output << std::setprecision(17);
      std::vector<cand_pos_t> fixed(sequence.size() + 1, -1), opening;
      for (std::size_t i = 0; i < scaffold.size(); ++i) {
        if (scaffold[i] == '(') opening.push_back(i + 1);
        if (scaffold[i] == ')') {
          const auto left = opening.back(); opening.pop_back();
          fixed[left] = i + 1;
        }
      }
      for (const auto &p : result.pair_probabilities) {
        if (std::max(p.round_probability, p.square_probability) < 1e-15) continue;
        if (fixed[p.i] == p.j) continue;
        output << p.i << '\t' << p.j << '\t' << p.round_probability
               << '\t' << p.square_probability << '\n';
      }
    }
  } else if (mode == "graph-outside") {
    ok = outside<semirings::FactorizedLogSumExpSemiring, semirings::FactorizedInsideGraph>(
        sequence, scaffold, widths, log_z, states, terms, factors, bytes, tape_nodes,
        argc == 5 ? argv[4] : nullptr);
  } else if (mode == "taped-outside") {
    ok = outside<semirings::TapeLogSumExpSemiring, semirings::LogSumExpTape>(
        sequence, scaffold, widths, log_z, states, terms, factors, bytes, tape_nodes,
        argc == 5 ? argv[4] : nullptr);
  } else return 1;
  const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
  std::cout << std::setprecision(17) << "n=" << sequence.size() << " b=" << width
            << " mode=" << mode << " ok=" << ok << " log_z=" << log_z
            << " seconds=" << seconds << " peak_kb=" << usage.ru_maxrss
            << " states=" << states << " terms=" << terms << " factors=" << factors
            << " graph_bytes=" << bytes << " tape_nodes=" << tape_nodes << '\n';
  return ok ? 0 : 1;
}
