#include "runtime_observations.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/semirings/logsumexp_tape.hh"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

using cparty::scfg2::deductions::Deduction;
using cparty::scfg2::deductions::ItemKey;
using cparty::scfg2::engine::TracebackPair;
using cparty::scfg2::engine::TracebackPairLane;

struct DerivationWitness {
  std::string signature;
  double weight = 0.0;
};

struct OutputGroup {
  std::vector<TracebackPair> pairs;
  double summed_weight = 0.0;
  double best_weight = 0.0;
  std::uint64_t derivation_count = 0;
  std::vector<DerivationWitness> witnesses;
};

using GroupMap = std::map<std::string, OutputGroup>;
using GroupChart = std::unordered_map<
    ItemKey,
    GroupMap,
    cparty::scfg2::deductions::ItemKeyHash,
    cparty::scfg2::deductions::ItemKeyEq>;

using PairProbabilityMap = std::map<std::pair<cand_pos_t, cand_pos_t>, double>;

void poison_storage_only_carriers(cparty::scfg2::testing::ObservedRuntimeSnapshot &owner, cand_pos_t n) {
  using cparty::scfg2::storage::ExactItemKey;
  using cparty::scfg2::storage::NonTerminal;
  constexpr double kPoison = 1.0e100;
  const std::vector<NonTerminal> exact_carriers = {
      NonTerminal::WIP_VP_L,
      NonTerminal::WIP_VP_R,
      NonTerminal::WIP_BE_L,
      NonTerminal::WIP_BE_R,
      NonTerminal::VP_CASE3_BOUNDARY,
      NonTerminal::VP_CASE3_R,
      NonTerminal::VP_CASE3_MR,
      NonTerminal::VP_CASE3_LEFT,
      NonTerminal::VP_C12_R2,
      NonTerminal::VP_C12_2,
      NonTerminal::VP_INTERNAL_LOCAL,
      NonTerminal::WMB_BTAIL,
      NonTerminal::WMB_BSTEP,
      NonTerminal::WMB_NONEMPTY,
      NonTerminal::WMB_EMPTY_LOCAL,
      NonTerminal::WM,
      NonTerminal::WMp,
  };
  for (cand_pos_t i = 1; i <= n; ++i) {
    for (cand_pos_t j = i + 1; j <= n; ++j) {
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::WIP_VP_L, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::WIP_VP_R, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::WIP_BE_L, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::WIP_BE_R, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::VP_CASE3_BOUNDARY, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::VP_CASE3_R, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::VP_CASE3_MR, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::VP_CASE3_LEFT, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::VP_C12_R2, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::VP_C12_2, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::VP_INTERNAL_LOCAL, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::WMB_BTAIL, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::WMB_BSTEP, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::WMB_NONEMPTY, i, j, kPoison);
      owner.diagnostics().set(cparty::scfg2::testing::DiagnosticTable::WMB_EMPTY_LOCAL, i, j, kPoison);
      for (const auto nt : exact_carriers) {
        ExactItemKey item;
        item.nt = nt;
        item.i = i;
        item.j = j;
        cparty::scfg2::storage::exact_store::set_total(item, kPoison);
      }
    }
  }
}

std::string pair_set_key(std::vector<TracebackPair> pairs) {
  std::sort(pairs.begin(), pairs.end(), [](const auto &lhs, const auto &rhs) {
    if (lhs.i != rhs.i) return lhs.i < rhs.i;
    if (lhs.j != rhs.j) return lhs.j < rhs.j;
    return static_cast<unsigned char>(lhs.lane) <
           static_cast<unsigned char>(rhs.lane);
  });
  std::string key;
  for (const auto &pair : pairs) {
    key += std::to_string(pair.i) + ":" + std::to_string(pair.j) + ":" +
           (pair.lane == TracebackPairLane::Round ? "R;" : "S;");
  }
  return key;
}

std::uint64_t saturated_product(std::uint64_t lhs, std::uint64_t rhs) {
  if (lhs == 0 || rhs == 0) return 0;
  if (lhs > std::numeric_limits<std::uint64_t>::max() / rhs) {
    return std::numeric_limits<std::uint64_t>::max();
  }
  return lhs * rhs;
}

std::uint64_t saturated_sum(std::uint64_t lhs, std::uint64_t rhs) {
  if (lhs > std::numeric_limits<std::uint64_t>::max() - rhs) {
    return std::numeric_limits<std::uint64_t>::max();
  }
  return lhs + rhs;
}

void merge_group(GroupMap *groups, OutputGroup candidate) {
  const std::string key = pair_set_key(candidate.pairs);
  auto [it, inserted] = groups->try_emplace(key, std::move(candidate));
  if (inserted) return;
  for (const auto &witness : candidate.witnesses) {
    if (it->second.witnesses.size() >= 2) break;
    const bool already_present = std::any_of(
        it->second.witnesses.begin(),
        it->second.witnesses.end(),
        [&](const auto &existing) {
          return existing.signature == witness.signature;
        });
    if (!already_present) {
      it->second.witnesses.push_back(witness);
    }
  }
  it->second.derivation_count = saturated_sum(
      it->second.derivation_count, candidate.derivation_count);
  it->second.summed_weight += candidate.summed_weight;
  if (candidate.best_weight > it->second.best_weight) {
    it->second.best_weight = candidate.best_weight;
    it->second.pairs = std::move(candidate.pairs);
  }
}

std::string deduction_signature(const Deduction &deduction) {
  return std::string(cparty::scfg2::grammar::nonterminal_name(
                         deduction.parent.nonterminal)) +
         "(" + std::to_string(deduction.parent.i) + "," +
         std::to_string(deduction.parent.j) + ")/" +
         std::string(cparty::scfg2::grammar::rule_name(deduction.rule)) +
         "[" + std::to_string(deduction.split.k) + "," +
         std::to_string(deduction.split.l) + "," +
         std::to_string(deduction.split.p) + "," +
         std::to_string(deduction.split.q) + "]";
}

ItemKey runtime_child_item(const Deduction &deduction, std::size_t child_idx) {
  ItemKey child = deduction.children[child_idx].item;
  if (deduction.rule == cparty::scfg2::grammar::RuleId::WMBP_DIRECT_VP &&
      child.nonterminal == cparty::scfg2::grammar::NonTerminal::VP) {
    child.nonterminal = cparty::scfg2::grammar::NonTerminal::VP;
  }
  return child;
}

template <class Context>
GroupMap enumerate_root_output_groups(cparty::scfg2::testing::ObservedRuntimeSnapshot &owner,
                                      Context &ctx,
                                      bool poison_carriers = false,
                                      bool enforce_completed_wmb_weak_closure = false) {
  using namespace cparty::scfg2;
  owner.reset_diagnostics();
  storage::exact_store::reset();
  // Deliberately do not run the production DP first. Storage-only carriers
  // are observer-derived provenance, not grammar inputs; output grouping must
  // be recoverable from deductions, local weights, and child groups alone.

  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  if (poison_carriers) poison_storage_only_carriers(owner, n);
  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  engine::AllLocalWeightModel<Context> weight(&ctx);
  GroupChart chart;

  for (const auto &item : schedule) {
    GroupMap item_groups;
    const auto candidates = deductions::generate_exact_basic_deductions_for_item(
        item, ctx);
    for (const auto &deduction : candidates) {
      // SI defines WMB/ZP as a completed pseudoknotted closed region.  Keep
      // this observation-only switch separate from the legacy-compatible
      // production grammar so the missing implementation-domain guard can be
      // regression-tested without silently changing the semantic oracle.
      if (enforce_completed_wmb_weak_closure &&
          deduction.parent.nonterminal == grammar::NonTerminal::WMB &&
          !ctx.weakly_closed(deduction.parent.i, deduction.parent.j)) {
        continue;
      }
      if (!engine::legacy_span_schedule_supports(deduction)) continue;
      const double local = weight.weight_for(deduction);
      if (!(local > 0.0)) continue;

      OutputGroup initial;
      initial.summed_weight = local;
      initial.best_weight = local;
      initial.derivation_count = 1;
      initial.witnesses.push_back(
          DerivationWitness{deduction_signature(deduction), local});
      engine::append_deduction_generated_pairs(deduction, ctx, &initial.pairs);
      GroupMap partial;
      merge_group(&partial, std::move(initial));

      for (std::size_t child_idx = 0;
           child_idx < deduction.child_count && !partial.empty();
           ++child_idx) {
        const auto child_it = chart.find(runtime_child_item(deduction, child_idx));
        if (child_it == chart.end() || child_it->second.empty()) {
          partial.clear();
          break;
        }
        GroupMap joined;
        for (const auto &[_, prefix] : partial) {
          for (const auto &[__, child] : child_it->second) {
            OutputGroup combined = prefix;
            for (const auto &pair : child.pairs) {
              engine::append_unique_traceback_pair(&combined.pairs, pair);
            }
            std::vector<DerivationWitness> joined_witnesses;
            for (const auto &prefix_witness : combined.witnesses) {
              for (const auto &child_witness : child.witnesses) {
                if (joined_witnesses.size() >= 2) break;
                joined_witnesses.push_back(DerivationWitness{
                    prefix_witness.signature + "{" +
                        child_witness.signature + "}",
                    prefix_witness.weight * child_witness.weight,
                });
              }
              if (joined_witnesses.size() >= 2) break;
            }
            combined.witnesses = std::move(joined_witnesses);
            combined.best_weight *= child.best_weight;
            combined.summed_weight *= child.summed_weight;
            combined.derivation_count = saturated_product(
                combined.derivation_count, child.derivation_count);
            merge_group(&joined, std::move(combined));
          }
        }
        partial = std::move(joined);
      }
      for (auto &[_, group] : partial) merge_group(&item_groups, std::move(group));
    }
    chart[item] = std::move(item_groups);
  }

  const ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
  const auto root_it = chart.find(root);
  return root_it == chart.end() ? GroupMap{} : root_it->second;
}

template <class Context>
PairProbabilityMap tape_pair_probabilities(cparty::scfg2::testing::ObservedRuntimeSnapshot &owner,
                                           Context &ctx,
                                           double *log_z) {
  using namespace cparty::scfg2;
  const cand_pos_t n = static_cast<cand_pos_t>(owner.length());
  const auto schedule = engine::build_pkfree_basic_exact_schedule(n, ctx);
  semirings::LogSumExpTape tape;
  const auto run = [&] {
    semirings::TapeLogSumExpSemiring::Scope scope(tape);
    return replay::run_w_final_exact_dp_schedule<
        semirings::TapeLogSumExpSemiring>(
        schedule,
        ctx,
        [](const Deduction &) { return true; },
        [](const ItemKey &, const Deduction &, auto, const auto &, double) {});
  }();
  const ItemKey root{grammar::NonTerminal::W, 1, n, -1, -1};
  const auto root_value = run.get(root);
  if (log_z) *log_z = root_value.log_value;
  const auto adjoints = tape.reverse(root_value.node);
  PairProbabilityMap out;
  for (const auto &[item, value] : run.values()) {
    const auto nt = item.nonterminal;
    const bool owns_pair =
        nt == grammar::NonTerminal::V ||
        nt == grammar::NonTerminal::VP;
    if (!owns_pair ||
        value.node == semirings::invalid_tape_node ||
        ctx.pair_partner(item.i) == item.j ||
        ctx.pair_partner(item.j) == item.i) {
      continue;
    }
    out[{item.i, item.j}] += adjoints[value.node];
  }
  return out;
}

PairProbabilityMap grouped_pair_probabilities(const GroupMap &groups,
                                              double *partition_weight) {
  double total = 0.0;
  PairProbabilityMap pair_weights;
  for (const auto &[_, group] : groups) {
    total += group.summed_weight;
    for (const auto &pair : group.pairs) {
      pair_weights[{pair.i, pair.j}] += group.summed_weight;
    }
  }
  if (partition_weight) *partition_weight = total;
  for (auto &[_, weight] : pair_weights) weight /= total;
  return pair_weights;
}

bool check_tape_pair_probabilities(cparty::scfg2::testing::ObservedRuntimeSnapshot &owner,
                                   cparty::scfg2::replay::
                                       WFinalExactInsideContext &ctx,
                                   const GroupMap &groups,
                                   const std::string &sequence,
                                   const std::string &scaffold) {
  double grouped_z = 0.0;
  const auto expected = grouped_pair_probabilities(groups, &grouped_z);
  double tape_log_z = -std::numeric_limits<double>::infinity();
  const auto observed = tape_pair_probabilities(owner, ctx, &tape_log_z);
  const double grouped_log_z = std::log(grouped_z);
  if (std::fabs(grouped_log_z - tape_log_z) > 1e-9) {
    std::cerr << "grouped/tape log-Z mismatch sequence=" << sequence
              << " scaffold=" << scaffold << " grouped=" << grouped_log_z
              << " tape=" << tape_log_z << '\n';
    return false;
  }
  PairProbabilityMap all = expected;
  all.insert(observed.begin(), observed.end());
  for (const auto &[pair, _] : all) {
    const auto expected_it = expected.find(pair);
    const auto observed_it = observed.find(pair);
    const double expected_probability =
        expected_it == expected.end() ? 0.0 : expected_it->second;
    const double observed_probability =
        observed_it == observed.end() ? 0.0 : observed_it->second;
    if (std::fabs(expected_probability - observed_probability) > 1e-9) {
      std::cerr << "grouped/tape pair probability mismatch sequence="
                << sequence << " scaffold=" << scaffold << " pair=("
                << pair.first << ',' << pair.second << ") grouped="
                << expected_probability << " tape=" << observed_probability
                << '\n';
      return false;
    }
  }
  for (cand_pos_t position = 1;
       position <= static_cast<cand_pos_t>(sequence.size());
       ++position) {
    double enumerated_unpaired_weight = 0.0;
    for (const auto &[_, group] : groups) {
      const bool paired = std::any_of(
          group.pairs.begin(), group.pairs.end(), [&](const auto &pair) {
            return pair.i == position || pair.j == position;
          });
      if (!paired) enumerated_unpaired_weight += group.summed_weight;
    }
    double taped_paired_probability = 0.0;
    for (const auto &[pair, probability] : observed) {
      if (pair.first == position || pair.second == position) {
        taped_paired_probability += probability;
      }
    }
    // Fixed G belongs to the union even though output groups contain only G'.
    if (ctx.pair_partner(position) > 0) {
      taped_paired_probability = 1.0;
      enumerated_unpaired_weight = 0.0;
    }
    const double enumerated_unpaired =
        enumerated_unpaired_weight / grouped_z;
    const double taped_unpaired = 1.0 - taped_paired_probability;
    if (std::fabs(enumerated_unpaired - taped_unpaired) > 1e-9) {
      std::cerr << "grouped/tape unpaired probability mismatch sequence="
                << sequence << " scaffold=" << scaffold
                << " position=" << position
                << " grouped=" << enumerated_unpaired
                << " tape=" << taped_unpaired << '\n';
      return false;
    }
  }
  return true;
}

bool same_group_map(const GroupMap &lhs, const GroupMap &rhs) {
  if (lhs.size() != rhs.size()) return false;
  for (const auto &[key, left] : lhs) {
    const auto right_it = rhs.find(key);
    if (right_it == rhs.end()) return false;
    const auto &right = right_it->second;
    if (left.derivation_count != right.derivation_count ||
        left.summed_weight != right.summed_weight ||
        left.best_weight != right.best_weight ||
        !cparty::scfg2::engine::same_traceback_pair_set(left.pairs,
                                                        right.pairs)) {
      return false;
    }
  }
  return true;
}

struct KnownAmbiguityExpectation {
  std::string target;
  std::uint64_t derivation_count = 0;
  std::vector<double> energies;
};

bool check_case(
    const std::string &sequence,
    const std::string &input_scaffold,
    bool require_scaffold_crossing,
    const KnownAmbiguityExpectation *known_ambiguity = nullptr) {
  using namespace cparty::scfg2;
  std::string owner_sequence = sequence;
  std::string scaffold = input_scaffold;
  cparty::scfg2::testing::ObservedRuntimeSnapshot owner(
      owner_sequence, scaffold);
  const auto policy = known_ambiguity
      ? cparty::scfg2::WmbClosurePolicy::LegacyCompatible
      : cparty::scfg2::WmbClosurePolicy::RequireWeakClosure;
  replay::WFinalExactInsideContext ctx(owner, policy);
  const GroupMap groups = enumerate_root_output_groups(owner, ctx);
  if (groups.empty()) {
    std::cerr << "no root output groups for " << sequence
              << " scaffold=" << input_scaffold << '\n';
    return false;
  }
  if (!check_tape_pair_probabilities(
          owner, ctx, groups, sequence, input_scaffold)) {
    return false;
  }
  if (require_scaffold_crossing) {
    const GroupMap poisoned_groups =
        enumerate_root_output_groups(owner, ctx, true);
    if (!same_group_map(groups, poisoned_groups)) {
      std::cerr << "storage-only carrier poisoning changed output groups for "
                << sequence << " scaffold=" << input_scaffold << '\n';
      return false;
    }
  }

  std::size_t ambiguous_groups = 0;
  std::size_t nonempty_groups = 0;
  std::size_t scaffold_crossing_groups = 0;
  for (const auto &[_, group] : groups) {
    if (group.derivation_count > 1) ++ambiguous_groups;
    if (!group.pairs.empty()) ++nonempty_groups;
    bool crosses_scaffold = false;
    for (const auto &generated : group.pairs) {
      for (cand_pos_t left = 1;
           left <= static_cast<cand_pos_t>(sequence.size());
           ++left) {
        const cand_pos_t right = ctx.pair_round(left);
        if (right <= left) continue;
        const bool crosses_from_inside =
            left < generated.i && generated.i < right && right < generated.j;
        const bool crosses_from_outside =
            generated.i < left && left < generated.j && generated.j < right;
        if (crosses_from_inside || crosses_from_outside) {
          crosses_scaffold = true;
          break;
        }
      }
      if (crosses_scaffold) break;
    }
    if (crosses_scaffold) ++scaffold_crossing_groups;
    const auto constrained =
        replay::run_w_final_exact_output_constrained_viterbi_traceback(
            owner, group.pairs, policy);
    const ItemKey root{grammar::NonTerminal::W,
                       1,
                       static_cast<cand_pos_t>(sequence.size()),
                       -1,
                       -1};
    const double constrained_weight = constrained.run.chart.get(root);
    if (!(constrained_weight > 0.0) ||
        !engine::same_traceback_pair_set(constrained.generated_pairs,
                                         group.pairs)) {
      std::cerr << "reachable output rejected for " << sequence
                << " scaffold=" << input_scaffold
                << " target=" << pair_set_key(group.pairs) << '\n';
      return false;
    }
    const double rt_kcal = owner.parameters().kT / 1000.0;
    const double expected_energy = -rt_kcal * std::log(group.best_weight);
    const double constrained_energy = -rt_kcal * std::log(constrained_weight);
    if (std::fabs(expected_energy - constrained_energy) > 1e-10) {
      std::cerr << "grouped/constrained energy mismatch for " << sequence
                << " scaffold=" << input_scaffold
                << " target=" << pair_set_key(group.pairs)
                << " expected=" << expected_energy
                << " constrained=" << constrained_energy << '\n';
      return false;
    }
  }

  if (require_scaffold_crossing && scaffold_crossing_groups == 0) {
    std::cerr << "no scaffold-crossing output group for " << sequence
              << " scaffold=" << input_scaffold << '\n';
    return false;
  }
  if (ambiguous_groups != 0 && known_ambiguity == nullptr) {
    std::cerr << "same-output derivation multiplicity detected for "
              << sequence << " scaffold=" << input_scaffold
              << " ambiguous_groups=" << ambiguous_groups << '\n';
    for (const auto &[key, group] : groups) {
      if (group.derivation_count <= 1) continue;
      std::cerr << "ambiguous_target=" << key
                << " derivations=" << group.derivation_count << '\n';
      for (const auto &witness : group.witnesses) {
        const double energy =
            -(owner.parameters().kT / 1000.0) * std::log(witness.weight);
        std::cerr << "witness_weight=" << witness.weight
                  << " witness_energy=" << energy
                  << " witness=" << witness.signature << '\n';
      }
    }
    return false;
  }
  if (known_ambiguity != nullptr) {
    const auto it = groups.find(known_ambiguity->target);
    if (ambiguous_groups != 1 || it == groups.end() ||
        it->second.derivation_count != known_ambiguity->derivation_count ||
        it->second.witnesses.size() != known_ambiguity->energies.size()) {
      std::cerr << "known ambiguity shape changed for " << sequence
                << " scaffold=" << input_scaffold
                << " ambiguous_groups=" << ambiguous_groups << '\n';
      return false;
    }
    std::vector<double> observed_energies;
    double witnessed_sum = 0.0;
    for (const auto &witness : it->second.witnesses) {
      observed_energies.push_back(
          -(owner.parameters().kT / 1000.0) * std::log(witness.weight));
      witnessed_sum += witness.weight;
    }
    std::sort(observed_energies.begin(), observed_energies.end());
    auto expected_energies = known_ambiguity->energies;
    std::sort(expected_energies.begin(), expected_energies.end());
    for (std::size_t idx = 0; idx < observed_energies.size(); ++idx) {
      if (std::fabs(observed_energies[idx] - expected_energies[idx]) > 1e-10) {
        std::cerr << "known ambiguity energy changed for " << sequence
                  << " scaffold=" << input_scaffold
                  << " observed=" << observed_energies[idx]
                  << " expected=" << expected_energies[idx] << '\n';
        return false;
      }
    }
    const double sum_tolerance =
        std::numeric_limits<double>::epsilon() * witnessed_sum * 8.0;
    if (std::fabs(it->second.summed_weight - witnessed_sum) > sum_tolerance) {
      std::cerr << "known ambiguity grouped sum changed for " << sequence
                << " scaffold=" << input_scaffold << '\n';
      return false;
    }

    if (ctx.weakly_closed(2, 7) || ctx.weakly_closed(8, 13) ||
        !ctx.weakly_closed(2, 13)) {
      std::cerr << "known ambiguity weak-closure geometry changed for "
                << sequence << " scaffold=" << input_scaffold << '\n';
      return false;
    }
    replay::WFinalExactInsideContext closed_ctx(
        owner, cparty::scfg2::WmbClosurePolicy::RequireWeakClosure);
    const GroupMap paper_domain_groups =
        enumerate_root_output_groups(owner, closed_ctx);
    const auto paper_it = paper_domain_groups.find(known_ambiguity->target);
    const std::size_t paper_ambiguous_groups = std::count_if(
        paper_domain_groups.begin(),
        paper_domain_groups.end(),
        [](const auto &entry) { return entry.second.derivation_count > 1; });
    if (paper_it == paper_domain_groups.end() ||
        paper_it->second.derivation_count != 1 ||
        paper_it->second.best_weight != it->second.best_weight ||
        paper_ambiguous_groups != 0) {
      std::cerr << "paper-domain WMB closure did not canonicalize known "
                   "ambiguity for "
                << sequence << " scaffold=" << input_scaffold
                << " ambiguous_groups=" << paper_ambiguous_groups << '\n';
      return false;
    }
  }

  std::cout << "sequence=" << sequence << " scaffold=" << input_scaffold
            << " policy=" << cparty::scfg2::wmb_closure_policy_name(policy)
            << " groups=" << groups.size()
            << " nonempty_groups=" << nonempty_groups
            << " scaffold_crossing_groups=" << scaffold_crossing_groups
            << " carrier_poison_checked=" << require_scaffold_crossing
            << " ambiguous_groups=" << ambiguous_groups << '\n';
  return true;
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  struct SmallCase {
    std::string sequence;
    std::string scaffold;
    bool require_scaffold_crossing = false;
    KnownAmbiguityExpectation known_ambiguity;
  };
  const std::vector<SmallCase> cases = {
      {"GCGCGCGCGC", "..........", false, {}},
      {"GCAUCGCAUG", "..........", false, {}},
      {"GGGAAACCCU", "..........", false, {}},
      {"GCGCGCGCGCGC", "............", false, {}},
      {"GGGAAACCCUUU", "............", false, {}},
      {"GCGCGCGCGCGC", "((....))....", true, {}},
      {"GCGCGCGCGCGC", "..(....)....", true, {}},
      {"GCGCGCGCGCGC", "(..)..(....)", true, {}},
      {"GCGCGCGCGCGC", "(((....)))..", true, {}},
      {"GCGCGCGCGCGC", "((.x..))....", true, {}},
      {"GCGCGCGCGCGCGC", "((....))(....)", true, {}},
      {"GCGCGCGCGCGCGC", "(.(......)...)", true, {}},
      {"GCGCGCGCGCGCGC",
       "(.(......)...)",
       true,
       // External legacy term accounting after closing-neighbor patch0009.
       {"2:7:S;8:13:S;", 2, {32.78, 25.14}}},
  };
  for (const auto &test_case : cases) {
    if (!check_case(test_case.sequence,
                    test_case.scaffold,
                    test_case.require_scaffold_crossing,
                    test_case.known_ambiguity.target.empty()
                        ? nullptr
                        : &test_case.known_ambiguity)) {
      return EXIT_FAILURE;
    }
  }
  std::cout << "scfg2_constrained_small_n_grouping=ok\n";
  return EXIT_SUCCESS;
}
