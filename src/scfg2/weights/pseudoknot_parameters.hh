#ifndef CPARTY_SCFG2_WEIGHTS_PSEUDOKNOT_PARAMETERS_HH
#define CPARTY_SCFG2_WEIGHTS_PSEUDOKNOT_PARAMETERS_HH

namespace cparty::scfg2::weights {

// Model data, in 10 cal/mol units. Preserve the CParty/DP09 parameterization;
// reorganizing its storage does not introduce a different energy model.
struct PseudoknotEnergies {
  int exterior = -138;
  int in_multiloop = 1007;
  int in_pseudoloop = 1500;
  int band = 246;
  int unpaired = 6;
  int nested_region = 96;
  int multiloop_closing = 339;
  int multiloop_branch = 3;
  int multiloop_unpaired = 2;
  int band_multiloop_closing = 341;
  int band_multiloop_branch = 56;
  int band_multiloop_unpaired = 12;
  double stack_multiplier = 0.89;
  double internal_multiplier = 0.74;
};

// Serial configuration boundary, sampled by each RuntimeSnapshot. Changes do
// not alter existing snapshots. Public session model settings must stay fixed.
inline PseudoknotEnergies &pseudoknot_energy_defaults() {
  static PseudoknotEnergies model;
  return model;
}

struct PseudoknotFactors {
  double ps = 0, psm = 0, psp = 0, pb = 0, pup = 0, pps = 0;
  double a = 0, b = 0, c = 0, ap = 0, bp = 0, cp = 0;
};

}  // namespace cparty::scfg2::weights
#endif
