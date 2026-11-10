#ifndef CPARTY_SCFG2_WEIGHTS_PSEUDOKNOT_PARAMETERS_HH
#define CPARTY_SCFG2_WEIGHTS_PSEUDOKNOT_PARAMETERS_HH

namespace cparty::scfg2::weights {

// Model data, in 10 cal/mol units. Preserve the CParty/DP09 parameterization;
// reorganizing its storage does not introduce a different energy model.
// Values are stored as double so that continuous parameter studies can set
// non-integer energies; the defaults are the integer DP09 values.
struct PseudoknotEnergies {
  double exterior = -138;
  double in_multiloop = 1007;
  double in_pseudoloop = 1500;
  double band = 246;
  double unpaired = 6;
  double nested_region = 96;
  double multiloop_closing = 339;
  double multiloop_branch = 3;
  double multiloop_unpaired = 2;
  double band_multiloop_closing = 341;
  double band_multiloop_branch = 56;
  double band_multiloop_unpaired = 12;
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
