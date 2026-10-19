#include "scfg2/api/hotspot_viterbi.hh"

#include "scfg2/api/detail.hh"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <map>
#include <tuple>
#include <utility>
#include <vector>

namespace cparty::scfg2::api {
namespace {

using PairKey = std::vector<std::pair<int, int>>;
using Clock = std::chrono::steady_clock;

double elapsed_seconds(Clock::time_point begin, Clock::time_point end) {
  return std::chrono::duration<double>(end - begin).count();
}

PairKey pair_key(std::vector<BasePair> pairs) {
  PairKey out;
  out.reserve(pairs.size());
  for (const auto &raw : pairs) {
    out.emplace_back(std::min(raw.i, raw.j), std::max(raw.i, raw.j));
  }
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<BasePair> from_pair_key(const PairKey &key) {
  std::vector<BasePair> out;
  out.reserve(key.size());
  for (const auto &pair : key) {
    out.push_back(BasePair{pair.first, pair.second});
  }
  return out;
}

bool parse_pair_key(const std::string &structure,
                    PairKey *key,
                    std::string *error) {
  std::vector<BasePair> pairs;
  if (!detail::parse_dot_bracket_pairs(structure, &pairs, error)) return false;
  *key = pair_key(std::move(pairs));
  return true;
}

bool record_better(const ScaffoldViterbiRecord &lhs,
                   std::size_t lhs_index,
                   const ScaffoldViterbiRecord &rhs,
                   std::size_t rhs_index) {
  if (lhs.viterbi.energy != rhs.viterbi.energy) {
    return lhs.viterbi.energy < rhs.viterbi.energy;
  }
  PairKey lhs_scaffold;
  PairKey rhs_scaffold;
  std::string ignored;
  parse_pair_key(lhs.scaffold_structure, &lhs_scaffold, &ignored);
  parse_pair_key(rhs.scaffold_structure, &rhs_scaffold, &ignored);
  return std::tie(
             lhs_scaffold,
             lhs.proposal.pass,
             lhs.proposal.origin,
             lhs_index) <
         std::tie(
             rhs_scaffold,
             rhs.proposal.pass,
             rhs.proposal.origin,
             rhs_index);
}

bool render_scaffold(int n,
                     const std::vector<BasePair> &pairs,
                     std::string *structure,
                     std::string *error) {
  if (!detail::validate_pair_map(pairs, n, error)) return false;
  *structure = std::string(static_cast<std::size_t>(n), '.');
  auto sorted = pairs;
  std::sort(sorted.begin(), sorted.end(), [](const auto &lhs, const auto &rhs) {
    return lhs.i != rhs.i ? lhs.i < rhs.i : lhs.j > rhs.j;
  });
  for (std::size_t lhs = 0; lhs < sorted.size(); ++lhs) {
    const auto &a = sorted[lhs];
    for (std::size_t rhs = lhs + 1; rhs < sorted.size(); ++rhs) {
      const auto &b = sorted[rhs];
      if ((a.i < b.i && b.i < a.j && a.j < b.j) ||
          (b.i < a.i && a.i < b.j && b.j < a.j)) {
        if (error) *error = "hotspot scaffold proposal contains crossing pairs";
        return false;
      }
    }
    (*structure)[static_cast<std::size_t>(a.i - 1)] = '(';
    (*structure)[static_cast<std::size_t>(a.j - 1)] = ')';
  }
  return true;
}

}  // namespace

HotspotViterbiPrediction select_minimum_energy_records(
    std::vector<ScaffoldViterbiRecord> records) {
  HotspotViterbiPrediction out;
  out.records = std::move(records);
  std::map<PairKey, std::vector<std::size_t>> grouped;

  for (std::size_t index = 0; index < out.records.size(); ++index) {
    const auto &record = out.records[index];
    if (!record.viterbi.ok) {
      out.error = "scaffold Viterbi record " + std::to_string(index) +
                  " failed: " + record.viterbi.error;
      return out;
    }
    if (!record.viterbi.reachable) continue;
    if (!std::isfinite(record.viterbi.energy)) {
      out.error = "reachable scaffold Viterbi record has non-finite energy";
      return out;
    }
    PairKey union_key;
    PairKey scaffold_key;
    if (!parse_pair_key(
            record.viterbi.union_structure, &union_key, &out.error) ||
        !parse_pair_key(
            record.scaffold_structure, &scaffold_key, &out.error)) {
      return out;
    }
    grouped[std::move(union_key)].push_back(index);
  }

  out.structures.reserve(grouped.size());
  for (const auto &entry : grouped) {
    CompletedStructureGroup group;
    group.completed_pairs = from_pair_key(entry.first);
    group.record_indices = entry.second;
    group.selected_record = entry.second.front();
    for (std::size_t index : entry.second) {
      const auto &record = out.records[index];
      group.minimum_energy =
          std::min(group.minimum_energy, record.viterbi.energy);
      group.maximum_energy =
          std::max(group.maximum_energy, record.viterbi.energy);
      if (record_better(
              record,
              index,
              out.records[group.selected_record],
              group.selected_record)) {
        group.selected_record = index;
      }
    }
    group.energy_spread = group.maximum_energy - group.minimum_energy;
    group.representative_union_structure =
        out.records[group.selected_record].viterbi.union_structure;
    out.structures.push_back(std::move(group));
  }

  out.ok = true;
  out.reachable = !out.structures.empty();
  if (!out.reachable) return out;

  out.selected_structure = 0;
  for (std::size_t index = 1; index < out.structures.size(); ++index) {
    const auto &candidate = out.structures[index];
    const auto &selected = out.structures[out.selected_structure];
    if (candidate.minimum_energy < selected.minimum_energy ||
        (candidate.minimum_energy == selected.minimum_energy &&
         pair_key(candidate.completed_pairs) <
             pair_key(selected.completed_pairs))) {
      out.selected_structure = index;
    }
  }
  return out;
}

HotspotViterbiPrediction predict_hotspot_viterbi(
    const std::string &sequence,
    proposal::WavefrontHotspotConfig proposal_config,
    BeamConfig beam_config) {
  const auto total_begin = Clock::now();
  const auto proposal_begin = Clock::now();
  const auto proposed =
      proposal::propose_hotspots_wavefront(sequence, proposal_config);
  const double proposal_elapsed =
      elapsed_seconds(proposal_begin, Clock::now());
  if (!proposed.ok) {
    HotspotViterbiPrediction out;
    out.error = proposed.error;
    out.proposal_elapsed_seconds = proposal_elapsed;
    out.total_elapsed_seconds = elapsed_seconds(total_begin, Clock::now());
    return out;
  }

  std::vector<ScaffoldViterbiRecord> records;
  records.reserve(proposed.proposals.size());
  for (const auto &candidate : proposed.proposals) {
    std::string scaffold;
    std::string error;
    if (!render_scaffold(
            static_cast<int>(sequence.size()),
            candidate.scaffold_pairs,
            &scaffold,
            &error)) {
      HotspotViterbiPrediction out;
      out.proposal_stats = proposed.stats;
      out.error = error;
      out.proposal_elapsed_seconds = proposal_elapsed;
      out.total_elapsed_seconds = elapsed_seconds(total_begin, Clock::now());
      return out;
    }
    const auto viterbi_begin = Clock::now();
    BeamSession session(sequence, scaffold, beam_config);
    auto viterbi = session.viterbi();
    const double viterbi_elapsed =
        elapsed_seconds(viterbi_begin, Clock::now());
    records.push_back(ScaffoldViterbiRecord{
        candidate, scaffold, std::move(viterbi), viterbi_elapsed});
  }

  auto out = select_minimum_energy_records(std::move(records));
  out.proposal_stats = proposed.stats;
  out.proposal_elapsed_seconds = proposal_elapsed;
  for (const auto &record : out.records) {
    out.viterbi_elapsed_seconds += record.elapsed_seconds;
  }
  out.total_elapsed_seconds = elapsed_seconds(total_begin, Clock::now());
  return out;
}

}  // namespace cparty::scfg2::api
