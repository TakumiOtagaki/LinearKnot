#include "scfg2/weights/runtime_snapshot.hh"
// Freeze the former exact arithmetic-tape BPP/profile path as a small oracle.
// The public APIs must preserve both lanes, chain windows and both closures.
#include "scfg2/api/exact.hh"
#include "scfg2/api/profiles.hh"
#include "scfg2/engine/outside_marginals.hh"
#include "scfg2/profiles/exact_profiles.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

namespace {
using namespace cparty::scfg2;

struct Reference {
  engine::OutsideMarginalSummary marginal;
  profiles::ProfileResult profile;
};

Reference taped_reference(std::string seq, std::string g,
                          const std::vector<int> &windows,
                          WmbClosurePolicy policy) {
  using T = semirings::TapeLogSumExpSemiring;
  using V = T::value_type;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      seq, g);
  replay::WFinalLogWeightContext ctx(owner, policy);
  const auto n = static_cast<cand_pos_t>(seq.size());
  struct Recorded {
    deductions::Deduction deduction;
    profiles::UnpairedEmission emission;
    V contribution;
    weights::LogWeight local;
  };
  std::vector<Recorded> recorded;
  semirings::LogSumExpTape tape;
  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  const auto chart = [&] {
    T::Scope scope(tape);
    return replay::run_w_final_exact_dp_schedule<T>(schedule, ctx,
        [](const auto &) { return true; },
        [&](const auto &, const auto &d, V contribution, const auto &, auto local) {
          const auto emission = profiles::unpaired_emission(d,
              [&](cand_pos_t p) { return ctx.is_unpaired(p); });
          if (emission.count) recorded.push_back({d, emission, contribution, local});
        });
  }();
  const deductions::ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
  const auto value = chart.get(root);
  const auto adjoints = tape.reverse(value.node);
  const auto adjoint = [&](std::size_t node) {
    return node < adjoints.size() ? adjoints[node] : 0.;
  };
  // Deliberately retain the old batch assembler as the differential reference.
  profiles::ProfileAccumulator acc(n);
  for (const auto &entry : recorded) {
    const auto probability = adjoint(entry.contribution.node);
    if (!entry.emission.chain_step) {
      acc.add_whole_runs(entry.emission, probability);
      continue;
    }
    const auto &parent = entry.deduction.parent;
    auto &step = acc.chain_step(static_cast<int>(parent.nonterminal), parent.i,
                                parent.j, entry.emission.runs[0].context);
    step.log_entry = entry.contribution.log_value;
    step.log_step = weights::log_weight(entry.local);
    step.has_prev = entry.deduction.child_count > 0;
    step.emit_probability = probability;
    step.log_state = chart.get(parent).log_value;
    step.state_adjoint = adjoint(chart.get(parent).node);
  }
  return {engine::collect_outside_marginals(chart, adjoints, n, ctx),
          acc.finish(value.log_value, std::vector<cand_pos_t>(windows.begin(), windows.end()))};
}

bool close(double a, double b) {
  return a == b || std::fabs(a - b) <= 1e-10 * std::max({1., std::fabs(a), std::fabs(b)});
}

bool check(const std::string &seq, const std::string &g, WmbClosurePolicy policy) {
  const std::vector<int> windows{1, 2, 3, 5, 13, static_cast<int>(seq.size()),
                                static_cast<int>(seq.size()) + 1, 3};
  const auto ref = taped_reference(seq, g, windows, policy);
  const auto outside = api::outside(seq, g, policy);
  const auto profile = api::exact_structural_profiles(seq, g, windows, policy);
  if (!ref.marginal.ok || !ref.profile.ok || !outside.ok || !profile.ok ||
      outside.log_partition_weight != ref.profile.log_partition_weight ||
      profile.log_partition_weight != ref.profile.log_partition_weight ||
      outside.pair_probabilities.size() != ref.marginal.pairs.size()) return false;
  for (std::size_t p = 0; p < ref.marginal.pairs.size(); ++p) {
    const auto &a = ref.marginal.pairs[p];
    const auto &b = outside.pair_probabilities[p];
    if (a.i != b.i || a.j != b.j || !close(a.round_probability, b.round_probability) ||
        !close(a.square_probability, b.square_probability) ||
        !close(a.union_probability, b.union_probability)) return false;
  }
  for (std::size_t i = 1; i <= seq.size(); ++i) {
    if (!close(ref.marginal.unpaired_probability[i], outside.unpaired_probabilities[i])) return false;
    for (std::size_t c = 0; c < profiles::loop_context_count; ++c)
      if (!close(ref.profile.context[i][c], profile.contexts[i][c])) return false;
    for (std::size_t w = 0; w < windows.size(); ++w)
      if (!close(ref.profile.accessibility[w][i], profile.accessibility[w][i])) return false;
  }
  return true;
}
}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  std::vector<std::pair<std::string, std::string>> cases{
      {"AAAAAA", "......"}, {"GAGAAACAC", "........."},
      {"GCGCGCGCGCGC", "..(....)...."},
      {"GCGCGCGCGCGCGC", "(.(......)...)"},
      {"GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
       "....(...)((....)..)((.(.((...))))(.(....)((...))))"},
      {"G" + std::string(40, 'A') + "C", "(" + std::string(40, '.') + ")"},
      {"GCG" + std::string(31, 'A') + "CGC", "(((" + std::string(31, '.') + ")))"}};
  std::mt19937 rng(20261015);
  for (int k = 0; k < 16; ++k) {
    const int n = 16 + rng() % 20;
    std::string seq(n, 'A'), g(n, '.');
    for (auto &ch : seq) ch = "ACGU"[rng() % 4];
    if (k % 2) {
      seq[2] = 'G'; seq[n - 3] = 'C'; g[2] = '('; g[n - 3] = ')';
    }
    cases.push_back({seq, g});
  }
  std::size_t count = 0;
  for (const auto policy : {WmbClosurePolicy::RequireWeakClosure, WmbClosurePolicy::LegacyCompatible})
    for (const auto &[seq, g] : cases) {
      if (!check(seq, g, policy)) {
        std::cerr << "exact regenerated mismatch sequence=" << seq << " scaffold=" << g
                  << " policy=" << static_cast<int>(policy) << '\n';
        return EXIT_FAILURE;
      }
      ++count;
    }
  std::cout << "scfg2_exact_regenerated_outside=ok cases=" << count << '\n';
}
