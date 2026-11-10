#ifndef LINEARKNOT_EXPERIMENTAL_CONTRAFOLD_SNAPSHOT_HH
#define LINEARKNOT_EXPERIMENTAL_CONTRAFOLD_SNAPSHOT_HH

// Experimental derivative/parameter view of the existing local-weight source.
// The generator, executor, reverse pass and output-ownership index are reused.
#include "scfg2/api/detail.hh"
#include "scfg2/engine/generated_target_index.hh"
#include "scfg2/engine/outside.hh"
#include "scfg2/engine/run_schedule.hh"
#include "contrafold_context.hh"

#include <array>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>

namespace linearknot_experiment {
namespace sc = cparty::scfg2;
using Vector = std::array<double,3>;
using Context = sc::replay::ContrafoldContext<true>;
using Model = sc::weights::PseudoknotEnergies;
using Deduction = sc::deductions::Deduction;
using Semiring = sc::semirings::LogSumExpSemiring;

struct Restore {
  Model saved=sc::weights::pseudoknot_energy_defaults();
  ~Restore(){sc::weights::pseudoknot_energy_defaults()=saved;}
};

struct Snapshot {
  std::unique_ptr<sc::weights::RuntimeSnapshot> owner;
  std::unique_ptr<Context> context;
  std::unique_ptr<sc::engine::AllLocalWeightModel<Context>> local;
  Snapshot(const std::string &sequence,const std::string &g,const Model &model) {
    Restore restore;sc::weights::pseudoknot_energy_defaults()=model;
    owner=std::make_unique<sc::weights::RuntimeSnapshot>(sequence,g);
    if(std::abs(owner->parameters().temperature-37.)>1e-9)
      throw std::invalid_argument("learning view requires 37 C");
    context=std::make_unique<Context>(*owner,sc::WmbClosurePolicy::RequireWeakClosure);
    local=std::make_unique<sc::engine::AllLocalWeightModel<Context>>(context.get());
  }
};

} // namespace linearknot_experiment
#endif
