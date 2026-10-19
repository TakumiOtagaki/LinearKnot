#ifndef CPARTY_SCFG2_WEIGHTS_INTERNAL_LOOP_HH
#define CPARTY_SCFG2_WEIGHTS_INTERNAL_LOOP_HH

#include <algorithm>
#include <cmath>

extern "C" {
#include "ViennaRNA/loops/internal.h"
}

namespace cparty::scfg2::weights {

// Temperature-scaled energies in Vienna's 0.01 kcal/mol units, before PF
// rounding. The runtime boundary supplies these from the loaded parameter
// set; the weight layer does not read mutable global energy parameters.
struct InternalLoopExtrapolation {
  double bulge_energy_30 = 0.0;
  double internal_energy_30 = 0.0;
  double ninio_per_base = 0.0;
  double max_ninio = 0.0;
};

inline double internal_loop_boltzmann_weight(
    int u1, int u2, int type, int type2,
    short si1, short sj1, short sp1, short sq1,
    vrna_exp_param_t *params, const InternalLoopExtrapolation &extended) {
  const int size = u1 + u2;
  // Keep all existing small-loop special cases and their exact PF values.
  if (size <= MAXLOOP)
    return exp_E_IntLoop(u1, u2, type, type2, si1, sj1, sp1, sq1, params);

  if (params->model_details.noGUclosure &&
      (type == 3 || type == 4 || type2 == 3 || type2 == 4)) return 0.0;

  const auto round_energy = [&](double energy) {
    return params->model_details.pf_smooth ? energy : std::trunc(energy);
  };
  const auto boltzmann = [&](double energy) {
    return std::exp(-energy * 10.0 / params->kT);
  };
  const bool bulge = u1 == 0 || u2 == 0;
  // Round the complete initiation energy, as Vienna's parameter-table
  // expansion does; rounding only the logarithmic increment is different
  // away from 37 C when the base energy has a fractional part.
  const double initiation = round_energy(
      (bulge ? extended.bulge_energy_30 : extended.internal_energy_30) +
      params->lxc * std::log(static_cast<double>(size) / 30.0));
  double weight = boltzmann(initiation);
  if (bulge) {
    if (type > 2) weight *= params->expTermAU;
    if (type2 > 2) weight *= params->expTermAU;
    return weight;
  }

  // A long 1xn loop has its own mismatch table. The other small special
  // cases (1x1, 1x2, 2x2, 2x3) cannot occur beyond MAXLOOP.
  if (std::min(u1, u2) == 1) {
    weight *= params->expmismatch1nI[type][si1][sj1] *
              params->expmismatch1nI[type2][sq1][sp1];
  } else {
    weight *= params->expmismatchI[type][si1][sj1] *
              params->expmismatchI[type2][sq1][sp1];
  }
  const double asymmetry = std::min(
      extended.max_ninio,
      std::abs(u1 - u2) * round_energy(extended.ninio_per_base));
  return weight * boltzmann(asymmetry);
}

}  // namespace cparty::scfg2::weights

#endif
