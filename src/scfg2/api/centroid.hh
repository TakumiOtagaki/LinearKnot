#ifndef CPARTY_SCFG2_API_CENTROID_HH
#define CPARTY_SCFG2_API_CENTROID_HH

#include "scfg2/api/posterior_decoding.hh"

namespace cparty::scfg2::api {

using CentroidResult = PosteriorDecodingResult;

// Hamada gamma-centroid: maximize E[gamma*TP_pair + TN_pair].
// Generated pairs receive (gamma+1)*p(i,j)-1. Fixed G is preserved;
// each fixed pair has probability one. TN ranges over all positions i<j.
// expected_accuracy reports the full expected gain, including its constant.
// Gamma is finite and nonnegative. Use the same posterior and frozen support
// as MEA; band geometry is inherited from the shared grammar.
CentroidResult beam_centroid(const std::string &sequence,
                             const std::string &scaffold_structure,
                             double gamma = 1.0, BeamConfig config = {});
CentroidResult exact_centroid(const std::string &sequence,
                              const std::string &scaffold_structure,
                              double gamma = 1.0,
                    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure);

}  // namespace cparty::scfg2::api

#endif
