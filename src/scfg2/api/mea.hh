#ifndef CPARTY_SCFG2_API_MEA_HH
#define CPARTY_SCFG2_API_MEA_HH

#include "scfg2/api/posterior_decoding.hh"

namespace cparty::scfg2::api {

using MeaResult = PosteriorDecodingResult;

// Maximize 2*gamma*sum_selected_pairs p(i,j) + sum_unpaired_bases q(i)
// for R = G union G', including the constant reward of fixed G pairs.
// The posterior and admissible outputs come from one inside forest. Gamma
// is finite and nonnegative. Beam decoding never reranks or reprunes it.
MeaResult beam_mea(const std::string &sequence,
                   const std::string &scaffold_structure,
                   double gamma = 1.0, BeamConfig config = {});
MeaResult exact_mea(const std::string &sequence,
                    const std::string &scaffold_structure,
                    double gamma = 1.0,
                    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);

}  // namespace cparty::scfg2::api

#endif
