#ifndef CPARTY_SCFG2_API_BEAM_DETAIL_HH
#define CPARTY_SCFG2_API_BEAM_DETAIL_HH

#include "scfg2/api/beam.hh"
#include "scfg2/beam/prepare_ranking.hh"

namespace cparty::scfg2::api::detail {

inline bool validate_beam_config(const BeamConfig &config, std::string *error) {
  if (config.beam_size == 0) *error = "beam size must be positive";
  else if (config.ranking != BeamRanking::Baseline && config.ranking != BeamRanking::Outer)
    *error = "invalid beam ranking mode";
  else if (config.outer_pilot_exponent != OuterPilotExponent::Half &&
           config.outer_pilot_exponent != OuterPilotExponent::TwoThirds)
    *error = "invalid outer pilot exponent";
  else return true;
  return false;
}

inline replay::WFinalVpWmbpOnlineWidths beam_widths(const BeamConfig &config) {
  auto widths = replay::WFinalVpWmbpOnlineWidths::uniform(config.beam_size);
  widths.max_hairpin_unpaired = config.max_hairpin_unpaired;
  widths.max_case3_ancestor_depth = config.max_case3_ancestor_depth;
  widths.memory_observer = config.memory_observer;
  widths.memory_column_interval = config.memory_column_interval;
  return widths;
}

template <class Context>
beam::PreparedRanking prepare_beam_ranking(Context &ctx, cand_pos_t n,
                                            const BeamConfig &config) {
  return beam::prepare_ranking(ctx, n, beam_widths(config), config.beam_size,
                               config.ranking, config.outer_pilot_exponent);
}

template <class Semiring, class Context>
auto run_beam_inside(Context &ctx, cand_pos_t n, const BeamConfig &config) {
  const auto prepared = prepare_beam_ranking(ctx, n, config);
  return replay::run_w_final_vp_wmbp_online<Semiring>(
      ctx, n, prepared.widths, prepared.ranking);
}

} // namespace cparty::scfg2::api::detail
#endif
