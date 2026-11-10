#include "scfg2/weights/runtime_snapshot.hh"

#include "rna_pair_table.hh"
#include "scfg2/geometry/scaffold_input.hh"
#include "scfg2/weights/runtime_parameters.hh"
#include "scfg2/weights/sequence.hh"

#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/params/default.h"
}

namespace cparty::scfg2::weights {
namespace {
double temperature_energy(double g, double h, double ratio) {
  return h - (h - g) * ratio;
}
void validate_input(const std::string &sequence, const std::string &scaffold) {
  if (sequence.size() != scaffold.size() ||
      sequence.size() > static_cast<std::size_t>(std::numeric_limits<cand_pos_t>::max()))
    throw std::invalid_argument("runtime sequence/scaffold length mismatch");
  std::string error;
  if (!geometry::validate_scaffold_input(scaffold, &error)) throw std::invalid_argument(error);
}
}  // namespace

RuntimeSnapshot::RuntimeSnapshot(const std::string &sequence, const std::string &scaffold,
                                 RuntimeOptions options)
    : sequence_(normalize_rna_sequence(sequence)), scaffold_(scaffold),
      model_(pseudoknot_energy_defaults()) {
  validate_input(sequence, scaffold);
  parameters_.reset(runtime_pf_parameters(options.dangles));
  if (!parameters_) throw std::runtime_error("cannot create ViennaRNA parameter snapshot");
  parameters_->model_details.dangles = options.dangles;
  make_pair_matrix();
  encoded_.reset(encode_sequence(sequence_.c_str(), 0));
  adjacent_.reset(encode_sequence(sequence_.c_str(), 1));
  if (!encoded_ || !adjacent_) throw std::runtime_error("cannot encode RNA sequence");
  const double ratio = (parameters_->temperature + K0) / Tmeasure;
  extrapolation_ = {temperature_energy(bulge37[30], bulgedH[30], ratio),
      temperature_energy(internal_loop37[30], internal_loopdH[30], ratio),
      temperature_energy(ninio37, niniodH, ratio), static_cast<double>(MAX_NINIO)};
  initialize_factors();
}

void RuntimeSnapshot::initialize_factors() {
  const auto &p = *parameters_;
  const double ratio = (p.temperature + K0) / Tmeasure;
  const auto factor = [&](double g, double h) {
    const double energy = temperature_energy(g, h, ratio);
    return -(p.model_details.pf_smooth ? energy : std::trunc(energy)) * 10. / p.kT;
  };
  const auto pk = [&](double g) { return factor(g, 3. * g); };
  log_pk_ = {pk(model_.exterior), pk(model_.in_multiloop), pk(model_.in_pseudoloop),
      pk(model_.band), pk(model_.unpaired), pk(model_.nested_region),
      factor(model_.multiloop_closing, ML_closingdH), factor(model_.multiloop_branch, ML_interndH),
      factor(model_.multiloop_unpaired, ML_BASEdH), pk(model_.band_multiloop_closing),
      pk(model_.band_multiloop_branch), pk(model_.band_multiloop_unpaired)};
}

bool RuntimeSnapshot::bases_can_pair(cand_pos_t i, cand_pos_t j) const {
  return cparty::rna::pair_type_for(encoded_base(i), encoded_base(j)) != 0;
}

double RuntimeSnapshot::partition_free_energy_from_log_weight(double log_weight) const {
  return -log_weight * parameters_->kT / 1000.;
}

}  // namespace cparty::scfg2::weights
