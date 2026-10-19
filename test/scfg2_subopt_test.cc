#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/subopt.hh"
#include "scfg2/api/detail.hh"
#include "scfg2/engine/runtime_child_item.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stdexcept>

namespace {
using namespace cparty::scfg2;
using Pair = std::pair<int, int>;
struct Witness {
  double log_weight;
  std::vector<Pair> pairs;
};
void require(bool ok, const std::string &message) { if (!ok) throw std::runtime_error(message); }

std::vector<Pair> pairs_of(const api::SuboptStructure &structure) {
  std::vector<Pair> pairs;
  for (const auto &p : structure.generated_pairs) pairs.emplace_back(p.i, p.j);
  return pairs;
}

// Independent exhaustive deduction DP: no factorized forest or k-best code.
// Each child Cartesian product is fully expanded, even when weights tie.
std::vector<Witness> exhaustive(std::string sequence, std::string scaffold) {
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalLogWeightContext ctx(owner);
  engine::AllLocalWeightModel<decltype(ctx)> weights(&ctx);
  std::unordered_map<deductions::ItemKey, std::vector<Witness>,
      deductions::ItemKeyHash, deductions::ItemKeyEq> chart;
  const auto n = static_cast<cand_pos_t>(sequence.size());
  for (const auto &item : engine::build_pkfree_basic_exact_schedule(n, ctx)) {
    std::vector<Witness> outputs;
    for (const auto &d : deductions::generate_exact_basic_deductions_for_item(item, ctx)) {
      if (!engine::legacy_span_schedule_supports(d)) continue;
      const auto local = weights.log_weight_for(d);
      if (!std::isfinite(local)) continue;
      std::vector<engine::TracebackPair> emitted;
      engine::append_deduction_generated_pairs(d, ctx, &emitted);
      Witness initial{local, {}};
      for (const auto &p : emitted) initial.pairs.emplace_back(p.i, p.j);
      std::vector<Witness> partial{initial};
      for (std::size_t c = 0; c < d.child_count; ++c) {
        std::vector<Witness> joined;
        const auto it = chart.find(engine::RuntimeDeductionChildItem{}(d, c));
        if (it != chart.end()) for (const auto &a : partial) for (const auto &b : it->second) {
          Witness output{a.log_weight + b.log_weight, a.pairs};
          output.pairs.insert(output.pairs.end(), b.pairs.begin(), b.pairs.end());
          std::sort(output.pairs.begin(), output.pairs.end());
          output.pairs.erase(std::unique(output.pairs.begin(), output.pairs.end()), output.pairs.end());
          joined.push_back(std::move(output));
        }
        partial = std::move(joined);
      }
      outputs.insert(outputs.end(), partial.begin(), partial.end());
      require(outputs.size() < 100000, "small-n exhaustive fixture is too large");
    }
    chart[item] = std::move(outputs);
  }
  auto outputs = chart[{grammar::NonTerminal::W, 1, n, -1, -1}];
  std::sort(outputs.begin(), outputs.end(), [](const auto &a, const auto &b) {
    return a.log_weight > b.log_weight;
  });
  return outputs;
}

void check_order_and_mass(const api::SuboptResult &result) {
  require(result.ok && result.reachable, "subopt failed: " + result.error);
  require(!result.structures.empty(), "reachable forest returned no structures");
  std::set<std::vector<Pair>> unique;
  double log_z = -std::numeric_limits<double>::infinity();
  double previous = std::numeric_limits<double>::infinity();
  for (const auto &structure : result.structures) {
    require(structure.log_weight <= previous + 1e-12, "subopt energy order is wrong");
    require(unique.insert(pairs_of(structure)).second, "duplicate physical structure");
    previous = structure.log_weight;
    log_z = semirings::LogSumExpSemiring::plus(log_z, structure.log_weight);
  }
  require(std::abs(result.mfe_energy - result.structures.front().energy) < 1e-10,
          "MFE was computed on a different support");
  if (result.stop_reason == api::SuboptStopReason::Exhausted)
    require(std::abs(log_z - result.log_partition_weight) < 1e-9,
            "exhaustive subopt weights do not sum to inside");
}

void check_fixture(const std::string &sequence, const std::string &scaffold) {
  const auto want = exhaustive(sequence, scaffold);
  require(!want.empty(), "empty exhaustive fixture");
  std::map<std::vector<Pair>, double> expected;
  for (const auto &w : want)
    require(expected.emplace(w.pairs, w.log_weight).second,
            "independent deduction enumeration found grammar ambiguity");
  api::SuboptConfig all;
  all.top_k = want.size() + 1;
  const auto exact = api::exact_subopt(sequence, scaffold, all);
  check_order_and_mass(exact);
  require(exact.stop_reason == api::SuboptStopReason::Exhausted &&
          exact.structures.size() == want.size(), "exact missed an exhaustive structure");
  for (std::size_t k = 0; k < want.size(); ++k) {
    const auto &got = exact.structures[k];
    require(std::abs(got.log_weight - want[k].log_weight) < 1e-9,
            "exact k-best rank differs from independent exhaustive sort");
    require(expected.count(pairs_of(got)) &&
            std::abs(expected.at(pairs_of(got)) - got.log_weight) < 1e-9,
            "exact reconstructed a wrong pair map/weight");
    std::vector<api::BasePair> rendered;
    std::string error;
    require(api::detail::parse_dot_bracket_pairs(got.generated_structure, &rendered, &error),
            "invalid generated brackets: " + error);
    std::vector<Pair> rendered_pairs;
    for (const auto &p : rendered) rendered_pairs.emplace_back(p.i, p.j);
    std::sort(rendered_pairs.begin(), rendered_pairs.end());
    require(rendered_pairs == pairs_of(got), "rendered structure and pairs differ");
  }
  api::ExactSession session(sequence, scaffold);
  const auto mfe = session.viterbi();
  require(mfe.ok && mfe.reachable && std::abs(mfe.energy - exact.mfe_energy) < 1e-9,
          "exact top-1 differs from shared Viterbi");

  api::BeamConfig keep_all;
  keep_all.beam_size = keep_all.max_hairpin_unpaired = keep_all.max_case3_ancestor_depth =
      std::numeric_limits<std::size_t>::max();
  for (auto beam : {keep_all, api::BeamConfig{1}, api::BeamConfig{3}}) {
    const auto result = api::beam_subopt(sequence, scaffold, all, beam);
    check_order_and_mass(result);
    require(result.stop_reason == api::SuboptStopReason::Exhausted,
            "small beam fixture did not exhaust");
    if (beam.beam_size == keep_all.beam_size)
      require(result.structures.size() == want.size(), "keep-all beam differs from exact support");
    for (const auto &structure : result.structures)
      require(expected.count(pairs_of(structure)) &&
              std::abs(expected.at(pairs_of(structure)) - structure.log_weight) < 1e-9,
              "beam returned a structure/energy absent from exact grammar");
    const auto inside = api::beam_inside(sequence, scaffold, beam);
    require(inside.ok && std::abs(inside.log_partition_weight - result.log_partition_weight) < 1e-9,
            "beam subopt uses different retained support from inside");
    api::SuboptConfig limited;
    limited.top_k = std::min<std::size_t>(3, result.structures.size());
    const auto top = api::beam_subopt(sequence, scaffold, limited, beam);
    require(top.ok && top.structures.size() == limited.top_k &&
            top.stop_reason == api::SuboptStopReason::TopK, "top-k did not stop at K");
    for (std::size_t k = 0; k < top.structures.size(); ++k)
      require(top.structures[k].union_structure == result.structures[k].union_structure,
              "top-k prefix changed with K");
    limited.top_k = all.top_k;
    limited.delta_energy = 0.0;
    const auto band = api::beam_subopt(sequence, scaffold, limited, beam);
    require(band.ok && !band.structures.empty(), "zero-width energy band failed");
    std::size_t ties = 0;
    for (const auto &structure : result.structures)
      if (structure.energy - result.mfe_energy <= 1e-9) ++ties;
    require(band.structures.size() == ties, "zero band lost a tied optimum or included a higher energy");
  }
  all.delta_energy = 0.0;
  const auto band = api::exact_subopt(sequence, scaffold, all);
  require(band.ok && !band.structures.empty(), "exact energy band failed");
}

void check_larger_and_log_domain() {
  const std::string sequence = "AAAGGGAAAGGGCCCAAACCCAAA";
  const std::string g = "...(((......))).........";
  api::SuboptConfig config;
  config.top_k = 20;
  api::ExactSession session(sequence, g);
  const auto exact = api::exact_subopt(sequence, g, config);
  check_order_and_mass(exact);
  for (const auto &result : {exact, api::beam_subopt(sequence, g, config, api::BeamConfig{3}),
                                   api::beam_subopt(sequence, g, config, api::BeamConfig{20})}) {
    check_order_and_mass(result);
    for (const auto &structure : result.structures) {
      const auto evaluated = session.evaluate_energy(structure.generated_pairs);
      require(evaluated.ok && evaluated.reachable &&
              std::abs(evaluated.energy - structure.energy) < 1e-9,
              "larger suboptimal disagrees with constrained energy evaluation");
    }
  }
  std::string long_sequence(150, 'G'), long_scaffold(150, '(');
  long_sequence += "AAA" + std::string(150, 'C');
  long_scaffold += "..." + std::string(150, ')');
  const auto long_result = api::beam_subopt(long_sequence, long_scaffold, config, api::BeamConfig{3});
  check_order_and_mass(long_result);
  require(long_result.structures.size() == 1 && std::isfinite(long_result.mfe_energy) &&
          long_result.structures.front().union_structure == long_scaffold,
          "long log-domain scaffold lost its sole structure");
}
}

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  try {
    check_fixture("AAAAAA", "......");
    check_fixture("GGGAAACCC", ".........");
    check_fixture("GCGCGCGCGCGCGC", "(.(......)...)");  // corrected WMB closure
    check_fixture("GGGAAACCCAGGGAAACCC", "(((...)))..........");
    std::mt19937 rng(20261011);
    for (int c = 0; c < 12; ++c) {
      const int n = 12 + c % 5;
      std::string sequence(n, 'A'), scaffold(n, '.');
      for (auto &base : sequence) base = "ACGU"[rng() % 4];
      if (c % 2) {
        scaffold[2] = '('; scaffold[n - 4] = ')';
        sequence[2] = 'G'; sequence[n - 4] = 'C';
      }
      check_fixture(sequence, scaffold);
    }
    check_larger_and_log_domain();
    require(!api::exact_subopt("AAAAAA", "......", api::SuboptConfig{0, {}}).ok,
            "zero top_k accepted");
    require(!api::beam_subopt("AAAAAA", "......", {}, api::BeamConfig{0}).ok,
            "zero beam accepted");
    for (const double delta : {-1.0, std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity()})
      require(!api::exact_subopt("AAAAAA", "......", api::SuboptConfig{10, delta}).ok,
              "invalid delta accepted");
    require(!api::exact_subopt("AAAAAA", "[....]").ok, "invalid scaffold accepted");
    api::SuboptConfig legacy_config;
    legacy_config.top_k = 10000;
    const auto legacy = api::exact_subopt("GCGCGCGCGCGCGC", "(.(......)...)", legacy_config,
                                         api::WmbClosurePolicy::LegacyCompatible);
    require(!legacy.ok && legacy.structures.empty() && legacy.error.find("ambiguity") != std::string::npos,
            "legacy ambiguity was silently deduplicated");
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n'; return EXIT_FAILURE;
  }
  std::cout << "scfg2_subopt=ok (16 exhaustive fixtures, exact/keep-all/finite beam)\n";
}
