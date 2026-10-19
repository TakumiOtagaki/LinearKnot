#ifndef CPARTY_TEST_RUNTIME_OBSERVATIONS_HH
#define CPARTY_TEST_RUNTIME_OBSERVATIONS_HH

#include "diagnostic_tables.hh"
#include "scfg2/weights/runtime_snapshot.hh"

namespace cparty::scfg2::testing {

// An observation witness owns its projections separately from the weight
// snapshot. The production context sees only the RuntimeSnapshot base.
class ObservedRuntimeSnapshot : public weights::RuntimeSnapshot {
 public:
  ObservedRuntimeSnapshot(const std::string &sequence, const std::string &scaffold,
                          weights::RuntimeOptions options = {})
      : RuntimeSnapshot(sequence, scaffold, options), diagnostics_(length(), true) {
    std::fill(diagnostics_.prefix().begin(), diagnostics_.prefix().end(), 1.);
    auto &wi = diagnostics_.table(DiagnosticTable::WI);
    std::fill(wi.begin(), wi.end(), 1.);
  }
  DiagnosticTables &diagnostics() { return diagnostics_; }
  void reset_diagnostics() {
    diagnostics_.clear();
    // These observations are raw values, with a terminal W seed of one.
    diagnostics_.prefix()[0] = 1.;
  }

 private:
  DiagnosticTables diagnostics_;
};

}  // namespace cparty::scfg2::testing
#endif
