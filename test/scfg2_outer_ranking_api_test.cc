#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/beam.hh"
#include "scfg2/api/profiles.hh"
#include "scfg2/api/sampling.hh"
#include "scfg2/api/subopt.hh"
#include "scfg2/api/mea.hh"
#include "scfg2/api/centroid.hh"
#include "scfg2/profiles/beam_profiles.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace cparty::scfg2;
std::size_t cases = 0, taped_window_differences = 0;
double maximum_taped_window_difference = 0.0;

void require(bool ok, const std::string &message) {
  if (!ok) throw std::runtime_error(message);
}
bool close(double a, double b, double tolerance = 2e-9) {
  return a == b || std::fabs(a - b) <= tolerance;
}

// Independently condition the actual forward forest on every requested
// window being unpaired. This does not use loop emissions, chain-window
// assembly, outside, or production sampling/decoding rewards.
void check_conditioned_windows(std::string sequence, std::string scaffold,
                              const api::BeamConfig &config,
                              const api::StructuralProfileResult &profile) {
  weights::RuntimeSnapshot owner(sequence, scaffold);
  replay::WFinalLogWeightContext ctx(owner, config.wmb_closure_policy);
  const cand_pos_t n = sequence.size();
  const auto prepared = beam::prepare_ranking(ctx, n,
      replay::WFinalVpWmbpOnlineWidths::uniform(config.beam_size), config.beam_size,
      config.ranking, config.outer_pilot_exponent);
  using F = semirings::FactorizedLogSumExpSemiring;
  using S = semirings::LogSumExpSemiring;
  semirings::FactorizedInsideGraph graph;
  const auto run = [&] {
    F::Scope scope(graph);
    return replay::run_w_final_vp_wmbp_online<F>(ctx, n, prepared.widths, prepared.ranking);
  }();
  const auto root = run.chart.get({grammar::NonTerminal::W, 1, n, -1, -1});
  require(close(root.log_value, profile.log_partition_weight), "conditioned forest uses another beam");
  std::vector<std::pair<cand_pos_t, cand_pos_t>> pairs(graph.states().size());
  for (const auto &[key, value] : run.chart.values()) {
    using NT = grammar::NonTerminal;
    if (key.nonterminal != NT::V && key.nonterminal != NT::VP) continue;
    if (ctx.pair_partner(key.i) == key.j) continue;
    auto &pair = pairs.at(value.node);
    require(pair.first == 0 || pair == std::make_pair(key.i, key.j), "ambiguous test pair identity");
    pair = {key.i, key.j};
  }
  std::vector<double> filtered(graph.states().size());
  for (std::size_t w = 0; w < profile.access_lengths.size(); ++w) {
    const cand_pos_t length = profile.access_lengths[w];
    for (cand_pos_t left = 1; left + length - 1 <= n; ++left) {
      const cand_pos_t right = left + length - 1;
      bool fixed = false;
      for (cand_pos_t k = left; k <= right; ++k) fixed |= ctx.pair_partner(k) > 0;
      double expected = 0.0;
      if (!fixed) {
        for (std::size_t id = 0; id < graph.states().size(); ++id) {
          const auto pair = pairs[id];
          const bool forbidden = (pair.first >= left && pair.first <= right) ||
                                 (pair.second >= left && pair.second <= right);
          double total = S::zero();
          if (!forbidden) {
            const auto &state = graph.states()[id];
            for (std::size_t t = state.first_term; t < state.first_term + state.term_count; ++t) {
              const auto &term = graph.terms()[t];
              double value = term.local_log_weight;
              for (std::size_t f = term.first_factor; f < term.first_factor + term.factor_count; ++f) {
                const auto child = graph.factors()[f];
                require(child < id, "conditioned reference is not a DAG");
                value = S::times(value, filtered[child]);
              }
              total = S::plus(total, value);
            }
          }
          filtered[id] = total;
        }
        expected = std::exp(filtered[root.node] - root.log_value);
      }
      if (!close(expected, profile.accessibility[w][left]))
        std::cerr << std::setprecision(17) << "conditioned window mismatch case=" << cases
                  << " length=" << length << " start=" << left
                  << " conditioned=" << expected << " profile=" << profile.accessibility[w][left] << '\n';
      require(close(expected, profile.accessibility[w][left]),
              "accessibility differs from the independently conditioned forest");
    }
  }
}

std::vector<int> predicted_partners(const std::string &scaffold,
                                    const std::vector<api::BasePair> &added) {
  std::vector<int> partners(scaffold.size() + 1), stack;
  for (std::size_t i = 0; i < scaffold.size(); ++i) {
    if (scaffold[i] == '(') stack.push_back(i + 1);
    if (scaffold[i] == ')') {
      const int left = stack.back(); stack.pop_back();
      partners[left] = i + 1; partners[i + 1] = left;
    }
  }
  for (const auto &pair : added) {
    require(!partners[pair.i] && !partners[pair.j], "prediction reuses a scaffold endpoint");
    partners[pair.i] = pair.j; partners[pair.j] = pair.i;
  }
  return partners;
}

double check(const std::string &sequence, const std::string &scaffold,
             std::size_t b, api::BeamRanking mode, api::OuterPilotExponent exponent,
             WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure) {
  api::BeamConfig config{b};
  config.ranking = mode;
  config.outer_pilot_exponent = exponent;
  config.wmb_closure_policy = policy;
  api::BeamSession session(sequence, scaffold, config);
  // Viterbi first prepares the session guide. The PF owner is destroyed
  // before the same immutable guide is reused by the other consumers.
  const auto best = session.viterbi();
  const auto inside = session.inside();
  const auto cached = session.inside();
  const auto outside = session.outside();
  require(best.ok && best.reachable, "guided Viterbi: " + best.error);
  require(inside.ok && inside.reachable, "guided inside: " + inside.error);
  require(outside.ok && outside.reachable, "guided outside: " + outside.error);
  require(inside.log_partition_weight == cached.log_partition_weight &&
          inside.elapsed_seconds == cached.elapsed_seconds &&
          inside.log_partition_weight == outside.log_partition_weight &&
          inside.work.total_item_visits == outside.work.total_item_visits,
          "session consumers use different inside computations");
  auto copy = session;
  const auto copied = copy.outside();
  require(copied.ok && copied.log_partition_weight == inside.log_partition_weight,
          "copied session lost its numeric guide");

  const auto profile = api::beam_structural_profiles(sequence, scaffold, {1, 3}, config);
  const auto reference = profiles::beam_profiles_taped(sequence, scaffold,
      replay::WFinalVpWmbpOnlineWidths::uniform(b), {1, 3}, policy, mode, exponent);
  require(profile.ok, "guided profiles: " + profile.error);
  require(reference.ok, "guided profile tape: " + reference.error);
  require(close(profile.log_partition_weight, inside.log_partition_weight) &&
          close(reference.log_partition_weight, inside.log_partition_weight),
          "profile consumers use another beam");
  check_conditioned_windows(sequence, scaffold, config, profile);
  for (std::size_t i = 1; i < profile.contexts.size(); ++i) {
    for (std::size_t k = 0; k < profile.contexts[i].size(); ++k)
      require(close(profile.contexts[i][k], reference.context[i][k]),
              "guided profile differs from arithmetic tape");
    require(close(1.0 - profile.contexts[i][0], outside.unpaired_probabilities[i]),
            "guided pair and context posteriors differ");
  }
  // The historical tape profile assembler can alias private/public chain
  // identities. Observe its windows, but establish their physical meaning
  // with the independently conditioned forest above, not that assembler.
  for (std::size_t w = 0; w < profile.accessibility.size(); ++w)
    for (std::size_t i = 1; i + profile.access_lengths[w] - 1 <= sequence.size(); ++i) {
      const double difference = std::fabs(profile.accessibility[w][i] - reference.accessibility[w][i]);
      maximum_taped_window_difference = std::max(maximum_taped_window_difference, difference);
      taped_window_differences += difference > 2e-9;
    }

  const auto samples = api::beam_sample(sequence, scaffold, 8, 73, config);
  require(samples.ok, "guided sampling: " + samples.error);
  require(close(samples.log_partition_weight, inside.log_partition_weight),
          "guided sampling uses another beam");
  for (const auto &pairs : samples.generated_pairs) {
    const auto p = session.probability(pairs);
    require(p.ok && p.target_included && p.probability > 0.0 && p.probability <= 1.0 &&
            close(p.log_partition_weight, inside.log_partition_weight),
            "sample is absent from the matching guided target forest: " + p.error);
  }
  api::SuboptConfig subopt; subopt.top_k = 3;
  const auto ranked = api::beam_subopt(sequence, scaffold, subopt, config);
  require(ranked.ok && ranked.reachable, "guided subopt: " + ranked.error);
  require(close(ranked.log_partition_weight, inside.log_partition_weight),
          "guided subopt uses another beam");
  for (const auto &structure : ranked.structures) {
    const auto p = session.probability(structure.generated_pairs);
    require(p.ok && p.target_included && close(p.log_target_weight, structure.log_weight),
            "guided k-best and target weights differ: " + p.error);
  }

  constexpr double gamma = 1.5;
  const auto mea = api::beam_mea(sequence, scaffold, gamma, config);
  const auto centroid = api::beam_centroid(sequence, scaffold, gamma, config);
  require(mea.ok && mea.reachable && centroid.ok && centroid.reachable,
          "guided decoding: " + mea.error + " / " + centroid.error);
  require(close(mea.log_partition_weight, inside.log_partition_weight) &&
          close(centroid.log_partition_weight, inside.log_partition_weight),
          "guided decoding uses another beam");
  const auto mea_partners = predicted_partners(scaffold, mea.generated_pairs);
  double score = 0.0;
  for (const auto &p : outside.pair_probabilities)
    if (mea_partners[p.i] == p.j) score += 2.0 * gamma * p.union_probability;
  for (std::size_t i = 1; i < mea_partners.size(); ++i)
    if (!mea_partners[i]) score += outside.unpaired_probabilities[i];
  require(close(score, mea.expected_accuracy), "guided MEA score uses another posterior");
  const auto cp = predicted_partners(scaffold, centroid.generated_pairs);
  double tp = 0.0, tn = 0.0;
  std::vector<std::vector<double>> probabilities(sequence.size() + 1,
                                               std::vector<double>(sequence.size() + 1));
  for (const auto &p : outside.pair_probabilities) probabilities[p.i][p.j] = p.union_probability;
  for (std::size_t i = 1; i <= sequence.size(); ++i)
    for (std::size_t j = i + 1; j <= sequence.size(); ++j)
      if (cp[i] == static_cast<int>(j)) tp += probabilities[i][j];
      else tn += 1.0 - probabilities[i][j];
  require(close(gamma * tp + tn, centroid.expected_accuracy, 1e-7) &&
          close(tp, centroid.expected_true_positive_pairs) &&
          close(tn, centroid.expected_true_negative_pairs, 1e-7),
          "guided centroid score uses another posterior");
  ++cases;
  return inside.log_partition_weight;
}

void invalid_configs() {
  const std::string seq = "AAAAAA", g = "......";
  for (int which = 0; which < 3; ++which) {
    api::BeamConfig config;
    if (which == 0) config.beam_size = 0;
    if (which == 1) config.ranking = static_cast<api::BeamRanking>(99);
    if (which == 2) config.outer_pilot_exponent = static_cast<api::OuterPilotExponent>(99);
    api::BeamSession session(seq, g, config);
    require(!session.inside().ok && !session.outside().ok && !session.viterbi().ok &&
            !session.probability(std::vector<api::BasePair>{}).ok &&
            !api::beam_structural_profiles(seq, g, {1}, config).ok &&
            !api::beam_sample(seq, g, 1, 1, config).ok &&
            !api::beam_subopt(seq, g, {}, config).ok &&
            !api::beam_mea(seq, g, 1., config).ok &&
            !api::beam_centroid(seq, g, 1., config).ok,
            "a consumer accepted an invalid ranking configuration");
  }
}
} // namespace

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
    require(api::BeamConfig{}.ranking == api::BeamRanking::Outer &&
            api::BeamConfig{}.outer_pilot_exponent == api::OuterPilotExponent::TwoThirds,
            "public API default is not outer with q=2/3");
    invalid_configs();
    std::ifstream fixture("test/fixtures/scfg2_outer_ranking_20261014/bpRNA_RFAM_27342.txt");
    std::string name, sequence, scaffold;
    require(bool(std::getline(fixture, name) && std::getline(fixture, sequence) &&
                 std::getline(fixture, scaffold)), "missing outer integration fixture");
    for (auto mode : {api::BeamRanking::Baseline, api::BeamRanking::Outer}) {
      for (auto q : {api::OuterPilotExponent::Half, api::OuterPilotExponent::TwoThirds}) {
        check("GGAAACGAAACC", "((...)(...))", 3, mode, q);
        check("GCGCGCGCGCGCGC", "(.(......)...)", 3, mode, q);
        check("GGAAACGAAACC", "((...)(...))", 3, mode, q, WmbClosurePolicy::LegacyCompatible);
        const double z = check(sequence, scaffold, 20, mode, q);
        if (q == api::OuterPilotExponent::TwoThirds) {
          // Frozen direct pilot/outer/main composition, refreshed after
          // multiloop closing patch0009 (caps30/30, b20, q=2/3).
          const double expected = mode == api::BeamRanking::Baseline
              ? 56.550132047788956 : 60.715913333152933;
          require(close(z, expected, 1e-11), "public integration changed the frozen research root");
          if (mode == api::BeamRanking::Outer) {
            api::BeamSession implicit(sequence, scaffold, 20);
            require(implicit.config().ranking == mode && implicit.inside().log_partition_weight == z,
                    "implicit finite-beam default differs from explicit outer");
          }
        }
      }
    }
    std::cout << "scfg2_outer_ranking_api=ok cases=" << cases
              << " modes=2 exponents=2 default=outer"
              << " observed_tape_window_differences=" << taped_window_differences
              << " max_tape_window_difference=" << std::setprecision(17)
              << maximum_taped_window_difference << '\n';
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
