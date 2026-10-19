#ifndef CPARTY_SCFG2_WEIGHTS_RUNTIME_PARAMETERS_HH
#define CPARTY_SCFG2_WEIGHTS_RUNTIME_PARAMETERS_HH

#include <cstdlib>
#include <stdexcept>
#include <string>

extern "C" {
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/params/io.h"
}

namespace cparty::scfg2::weights {

// Serial compatibility boundary. Load before snapshotting either PF tables
// or long-loop energies; do not cache a pathname independently of Vienna's
// mutable global parameter state.
inline vrna_exp_param_t *runtime_pf_parameters(int dangles) {
  const char *path = std::getenv("CPARTY_FIXED_ENERGY_PARAM_FILE");
  if (path && (!*path ||
               !vrna_params_load(path, VRNA_PARAMETER_FORMAT_DEFAULT))) {
    throw std::invalid_argument(
        "cannot load CPARTY_FIXED_ENERGY_PARAM_FILE: " + std::string(path));
  }
  return cparty_log_pf_parameters_for_globals(dangles);
}

}  // namespace cparty::scfg2::weights

#endif
