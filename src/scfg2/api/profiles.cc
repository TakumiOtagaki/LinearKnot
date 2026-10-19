#include "scfg2/api/profiles.hh"

#include "scfg2/api/detail.hh"
#include "scfg2/api/beam_detail.hh"
#include "scfg2/profiles/beam_profiles.hh"
#include "scfg2/profiles/exact_profiles.hh"

#include <chrono>
#include <exception>
#include <utility>

namespace cparty::scfg2::api {
namespace {

bool validate(const std::string &sequence,
              const std::string &scaffold,
              const std::vector<int> &access_lengths,
              std::string *error) {
  if (!detail::validate_session_inputs(sequence, scaffold, error)) return false;
  for (const int w : access_lengths) {
    if (w < 1) {
      *error = "access lengths must be positive";
      return false;
    }
  }
  return true;
}

StructuralProfileResult convert(profiles::ProfileResult in,
                                const std::vector<int> &access_lengths) {
  StructuralProfileResult out;
  out.ok = in.ok;
  out.error = in.error;
  out.log_partition_weight = in.log_partition_weight;
  out.ensemble_free_energy = in.ensemble_free_energy;
  out.rt_kcal_mol = in.rt_kcal_mol;
  out.contexts = std::move(in.context);
  out.access_lengths = access_lengths;
  out.accessibility = std::move(in.accessibility);
  out.accessibility_by_loop = std::move(in.accessibility_by_loop);
  return out;
}

std::vector<cand_pos_t> window_lengths(const std::vector<int> &lengths) {
  return std::vector<cand_pos_t>(lengths.begin(), lengths.end());
}

template <class Fn>
StructuralProfileResult timed(const std::string &sequence,
                              const std::string &scaffold,
                              const std::vector<int> &access_lengths,
                              WmbClosurePolicy policy, Fn &&compute) {
  const auto started = std::chrono::steady_clock::now();
  StructuralProfileResult out;
  if (validate(sequence, scaffold, access_lengths, &out.error)) {
    try {
      out = convert(compute(), access_lengths);
    } catch (const std::exception &e) {
      out = StructuralProfileResult{};
      out.error = e.what();
    }
  }
  out.wmb_closure_policy = policy;
  out.elapsed_seconds = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - started).count();
  return out;
}

}  // namespace

StructuralProfileResult exact_structural_profiles(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::vector<int> &access_lengths, WmbClosurePolicy policy,
    bool by_loop) {
  return timed(sequence, scaffold_structure, access_lengths, policy, [&] {
    return profiles::exact_profiles(sequence, scaffold_structure,
                                    window_lengths(access_lengths), policy, by_loop);
  });
}

StructuralProfileResult beam_structural_profiles(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::vector<int> &access_lengths,
    BeamConfig config, bool by_loop) {
  std::string config_error;
  if (!detail::validate_beam_config(config, &config_error)) {
    StructuralProfileResult out;
    out.wmb_closure_policy = config.wmb_closure_policy;
    out.error = config_error;
    return out;
  }
  auto out = timed(sequence, scaffold_structure, access_lengths, config.wmb_closure_policy, [&] {
    return profiles::beam_profiles(sequence, scaffold_structure, detail::beam_widths(config),
        window_lengths(access_lengths), config.wmb_closure_policy,
        config.ranking, config.outer_pilot_exponent, by_loop);
  });
  out.beam_size = config.beam_size;
  return out;
}

}  // namespace cparty::scfg2::api
