#ifndef CPARTY_SCFG2_WEIGHTS_RUNTIME_SNAPSHOT_HH
#define CPARTY_SCFG2_WEIGHTS_RUNTIME_SNAPSHOT_HH

#include "scfg2/types.hh"
#include "scfg2/weights/internal_loop.hh"
#include "scfg2/weights/pseudoknot_parameters.hh"

#include <cstdlib>
#include <memory>
#include <string>

extern "C" {
#include "ViennaRNA/params/basic.h"
}

namespace cparty::scfg2::weights {

struct RuntimeOptions {
  int dangles = 2;
};

// Owns normalized sequence, fixed scaffold and a parameter snapshot. Folding,
// outside, sampling and decoding belong to their existing shared consumers.
class RuntimeSnapshot {
 public:
  using LogPkParameters = PseudoknotFactors;
  RuntimeSnapshot(const std::string &sequence, const std::string &scaffold,
                  RuntimeOptions options = {});
  RuntimeSnapshot(const RuntimeSnapshot &) = delete;
  RuntimeSnapshot &operator=(const RuntimeSnapshot &) = delete;

  cand_pos_t length() const { return static_cast<cand_pos_t>(sequence_.size()); }
  const std::string &sequence() const { return sequence_; }
  const std::string &scaffold() const { return scaffold_; }
  short encoded_base(cand_pos_t pos) const { return encoded_.get()[pos]; }
  short adjacent_base(cand_pos_t pos) const { return adjacent_.get()[pos]; }
  bool bases_can_pair(cand_pos_t i, cand_pos_t j) const;
  // Vienna's exp* fields contain log weights, never raw Boltzmann factors.
  const vrna_exp_param_t &parameters() const { return *parameters_; }
  const vrna_exp_param_t &log_pf_parameters() const { return *parameters_; }
  const PseudoknotEnergies &pseudoknot_model() const { return model_; }
  const PseudoknotFactors &log_pk_parameters() const { return log_pk_; }
  const InternalLoopExtrapolation &loop_extrapolation() const { return extrapolation_; }
  double partition_free_energy_from_log_weight(double log_weight) const;

 private:
  void initialize_factors();
  std::string sequence_, scaffold_;
  std::unique_ptr<vrna_exp_param_t, decltype(&std::free)> parameters_{nullptr, &std::free};
  std::unique_ptr<short[], decltype(&std::free)> encoded_{nullptr, &std::free}, adjacent_{nullptr, &std::free};
  PseudoknotEnergies model_;
  PseudoknotFactors log_pk_;
  InternalLoopExtrapolation extrapolation_;
};

}  // namespace cparty::scfg2::weights
#endif
