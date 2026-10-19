// Observation-only audit of a single V pair on an exact or frozen beam forest.
// It computes P(0/1/2+ occurrences), independently of the outside derivative,
// and projects a positive derivation with repeated pair-state use.
#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"

#include <array>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <unordered_map>

int main(int argc, char **argv) {
  using namespace cparty::scfg2;
  using S = semirings::FactorizedLogSumExpSemiring;
  using Id = std::uint32_t;
  if (argc != 5) {
    std::cerr << "usage: " << argv[0]
              << " INPUT WIDTH|exact|exact-legacy-order I J\n";
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
  if (sequence.empty() || sequence.size() != scaffold.size()) return 1;
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const std::string mode = argv[2];
  const bool legacy_order = mode == "exact-legacy-order";
  const bool exact = mode == "exact" || legacy_order;
  const auto width = exact ? 0 : std::stoull(argv[2]);
  const auto i = static_cast<cand_pos_t>(std::stoi(argv[3]));
  const auto j = static_cast<cand_pos_t>(std::stoi(argv[4]));
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  semirings::FactorizedInsideGraph graph;
  const auto chart = [&] {
    S::Scope scope(graph);
    if (exact)
      return replay::run_w_final_exact_dp_schedule<S>(
          engine::build_pkfree_basic_exact_schedule(n, ctx), ctx,
          [&](const deductions::Deduction &ded) {
            // Observation-only reconstruction of legacy's i-descending
            // availability boundary. Its leftward BE read sees initialized
            // zero, because the BE opening column has not been computed yet.
            const bool shared =
                ded.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP ||
                ded.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP;
            return !legacy_order || !shared || ded.children[0].item.i >= ded.parent.i;
          },
          [](const deductions::ItemKey &, const deductions::Deduction &,
             S::value_type, const engine::Chart<S::value_type> &, double) {});
    return replay::run_w_final_vp_wmbp_online<S>(
        ctx, n, replay::WFinalVpWmbpOnlineWidths::uniform(width)).chart;
  }();
  const auto root = chart.get({grammar::NonTerminal::W, 1, n, -1, -1});
  const auto target = chart.get({grammar::NonTerminal::V, i, j, -1, -1});
  if (S::is_zero(root) || S::is_zero(target)) return 1;
  std::cout << std::setprecision(17);
  const auto post = graph.posteriors(root.node);
  std::cout << "n=" << n << " mode=" << (exact ? mode : "beam")
            << " b=" << width << " log_partition=" << S::priority(root)
            << " energy=" << owner.partition_free_energy_from_log_weight(S::priority(root))
            << " V(" << i << ',' << j
            << ") node=" << target.node << " expected_count=" << post[target.node] << '\n';

  const auto &states = graph.states();
  const auto &terms = graph.terms();
  const auto &factors = graph.factors();
  std::vector<std::array<double, 3>> distribution(states.size());
  std::vector<unsigned char> maximum(states.size(), 0);
  std::vector<Id> chosen(states.size(), semirings::invalid_outside_state);
  for (std::size_t node = 0; node < states.size(); ++node) {
    if (node == target.node) {
      distribution[node] = {0.0, 1.0, 0.0};
      maximum[node] = 1;
      continue;
    }
    const auto &state = states[node];
    double best_log = -std::numeric_limits<double>::infinity();
    for (std::size_t t = state.first_term;
         t < std::size_t(state.first_term) + state.term_count; ++t) {
      const auto &term = terms[t];
      std::array<double, 3> product{1.0, 0.0, 0.0};
      unsigned count = 0;
      for (std::size_t f = term.first_factor;
           f < std::size_t(term.first_factor) + term.factor_count; ++f) {
        const auto child = factors[f];
        std::array<double, 3> next{};
        for (unsigned a = 0; a < 3; ++a)
          for (unsigned b = 0; b < 3; ++b)
            next[std::min(2u, a + b)] += product[a] * distribution[child][b];
        product = next;
        count = std::min(3u, count + maximum[child]);
      }
      const double probability = std::exp(term.log_value - state.log_value);
      for (unsigned k = 0; k < 3; ++k) distribution[node][k] += probability * product[k];
      if (count > maximum[node] || (count == maximum[node] && term.log_value > best_log)) {
        maximum[node] = count;
        chosen[node] = static_cast<Id>(t);
        best_log = term.log_value;
      }
    }
  }
  const auto &p = distribution[root.node];
  if (std::fabs(p[0] + p[1] + p[2] - 1.0) > 1e-7 ||
      (maximum[root.node] <= 2 &&
       std::fabs(p[1] + 2.0 * p[2] - post[target.node]) > 1e-7)) {
    std::cerr << "occurrence distribution/derivative consistency failed\n";
    return 2;
  }
  std::cout << "p_zero=" << p[0] << " p_once=" << p[1]
            << " p_two_or_more=" << p[2] << " p_present=" << p[1] + p[2]
            << " max_count_capped3=" << unsigned(maximum[root.node]) << '\n';

  std::unordered_map<Id, deductions::ItemKey> labels;
  for (const auto &[item, value] : chart.values()) labels.emplace(value.node, item);
  const auto describe = [&](Id node) {
    const auto found = labels.find(node);
    if (found == labels.end()) return std::string("private#") + std::to_string(node);
    const auto &k = found->second;
    return std::string(grammar::nonterminal_name(k.nonterminal)) + '(' +
           std::to_string(k.i) + ',' + std::to_string(k.j) + ',' +
           std::to_string(k.ip) + ',' + std::to_string(k.jp) + ')';
  };
  if (maximum[root.node] < 2) return 0;
  Id fork = root.node;
  std::vector<Id> ancestors;
  while (fork != target.node) {
    ancestors.push_back(fork);
    const auto &term = terms[chosen[fork]];
    Id next = semirings::invalid_outside_state;
    for (std::size_t f = term.first_factor;
         f < std::size_t(term.first_factor) + term.factor_count; ++f)
      if (maximum[factors[f]] >= 2) { next = factors[f]; break; }
    if (next == semirings::invalid_outside_state) break;
    fork = next;
  }
  std::cout << "last_ancestors";
  const auto begin = ancestors.size() > 12 ? ancestors.size() - 12 : 0;
  for (std::size_t k = begin; k < ancestors.size(); ++k)
    std::cout << " -> " << describe(ancestors[k]);
  std::cout << '\n';
  const auto &fork_term = terms[chosen[fork]];
  std::cout << "fork " << describe(fork) << " term=" << chosen[fork] << '\n';
  for (std::size_t f = fork_term.first_factor;
       f < std::size_t(fork_term.first_factor) + fork_term.factor_count; ++f) {
    auto node = factors[f];
    std::cout << "factor " << describe(node) << " count=" << unsigned(maximum[node]);
    std::vector<Id> path;
    while (maximum[node] > 0) {
      if (labels.count(node) || node == target.node) path.push_back(node);
      if (node == target.node) break;
      const auto &term = terms[chosen[node]];
      Id next = semirings::invalid_outside_state;
      for (std::size_t q = term.first_factor;
           q < std::size_t(term.first_factor) + term.factor_count; ++q)
        if (maximum[factors[q]] > 0) { next = factors[q]; break; }
      if (next == semirings::invalid_outside_state) return 2;
      node = next;
    }
    for (const auto step : path) std::cout << " -> " << describe(step);
    std::cout << '\n';
  }
}
