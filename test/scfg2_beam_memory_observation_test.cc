#include "scfg2/api/beam.hh"

#include <cstdlib>
#include <iostream>
#include <numeric>
#include <stdexcept>

namespace sc = cparty::scfg2;
namespace {
void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
void check(sc::api::BeamRanking ranking, sc::api::OuterPilotExponent exponent,
           sc::WmbClosurePolicy closure) {
  const std::string sequence = "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  const std::string scaffold = "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  sc::api::BeamConfig config{20, 30, 30};
  config.ranking = ranking;
  config.outer_pilot_exponent = exponent;
  config.wmb_closure_policy = closure;
  sc::api::BeamSession plain(sequence, scaffold, config);
  const auto expected = plain.inside();
  const auto expected_outside = plain.outside();
  const auto expected_best = plain.viterbi();
  std::vector<sc::beam::MemoryObservation> observations;
  config.memory_column_interval = 7;
  config.memory_observer = [&](const auto &observation) { observations.push_back(observation); };
  sc::api::BeamSession observed(sequence, scaffold, config);
  const auto actual = observed.inside();
  require(expected.ok && actual.ok && expected.reachable && actual.reachable,
          "memory observation inside failed");
  require(actual.log_partition_weight == expected.log_partition_weight &&
              actual.work.total_item_visits == expected.work.total_item_visits &&
              actual.work.total_frontier_selected == expected.work.total_frontier_selected,
          "memory observation changed numeric forward or work");
  const auto calls = observations.size();
  require(observed.inside().log_partition_weight == actual.log_partition_weight &&
              observations.size() == calls, "cached inside repeated diagnostic callbacks");
  bool scaffold_seen = false, pilot_seen = false, outer_seen = false, main_seen = false;
  std::size_t retained = 0;
  for (const auto &observation : observations) {
    scaffold_seen |= observation.phase == sc::beam::MemoryPhase::Scaffold;
    pilot_seen |= observation.phase == sc::beam::MemoryPhase::Pilot;
    outer_seen |= observation.phase == sc::beam::MemoryPhase::Outer;
    main_seen |= observation.phase == sc::beam::MemoryPhase::Main;
    if (observation.phase != sc::beam::MemoryPhase::Main ||
        std::string(observation.event) != "forward_end") continue;
    for (const auto &family : observation.families) retained += family.live_items;
    require(retained == observation.chart_items && retained > 0,
            "family live counts do not match the retained chart");
    for (const auto &container : observation.containers) {
      require(container.name.find("VP_DIRECT") == std::string::npos &&
                  container.name.find("VP_CLOSED") == std::string::npos &&
                  container.name.find("vp_direct") == std::string::npos &&
                  container.name.find("vp_closed") == std::string::npos,
              "retired VP carrier storage is still observed");
      if (container.name == "vp_loop_sources" ||
          container.name == "v_loop_sources")
        require(container.entries == 0, "completed loop candidates survived forward");
      if (container.name == "wmb_btail_history")
        require(container.vector_capacity_bytes == container.vector_size_bytes,
                "BTAIL retained quadratic candidate capacity");
    }
  }
  require(scaffold_seen && main_seen && retained <= actual.work.total_positive_items,
          "forcing/main observation missing or live counts mistaken for visits");
  require(pilot_seen == (ranking == sc::api::BeamRanking::Outer) &&
              outer_seen == pilot_seen, "ranking phase attribution is incorrect");
  const auto outside = observed.outside();
  const auto best = observed.viterbi();
  require(outside.ok && expected_outside.ok && best.ok && expected_best.ok &&
              outside.log_partition_weight == expected_outside.log_partition_weight &&
              outside.unpaired_probabilities == expected_outside.unpaired_probabilities &&
              outside.pair_probabilities.size() == expected_outside.pair_probabilities.size() &&
              best.log_weight == expected_best.log_weight &&
              best.union_structure == expected_best.union_structure,
          "memory observation changed a downstream consumer");
  for (std::size_t i = 0; i < outside.pair_probabilities.size(); ++i) {
    const auto &a = outside.pair_probabilities[i];
    const auto &b = expected_outside.pair_probabilities[i];
    require(a.i == b.i && a.j == b.j && a.union_probability == b.union_probability &&
                a.round_probability == b.round_probability && a.square_probability == b.square_probability,
            "memory observation changed BPP");
  }
}
} // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  for (const auto ranking : {sc::api::BeamRanking::Baseline, sc::api::BeamRanking::Outer})
    for (const auto exponent : {sc::api::OuterPilotExponent::Half, sc::api::OuterPilotExponent::TwoThirds})
      for (const auto closure : {sc::WmbClosurePolicy::LegacyCompatible,
                                sc::WmbClosurePolicy::RequireWeakClosure})
        check(ranking, exponent, closure);
  std::cout << "scfg2_beam_memory_observation=ok\n";
}
