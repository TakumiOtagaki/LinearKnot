#include "scfg2/geometry/static_scaffold_geometry.hh"
#include "scfg2/weights/runtime_snapshot.hh"
#include "runtime_observations.hh"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <type_traits>
#include <utility>

namespace w = cparty::scfg2::weights;
namespace st = cparty::scfg2::testing;

template <class T, class = void> struct exposes_diagnostics : std::false_type {};
template <class T> struct exposes_diagnostics<T,
    std::void_t<decltype(std::declval<T &>().diagnostics())>> : std::true_type {};
static_assert(!exposes_diagnostics<w::RuntimeSnapshot>::value,
              "weight snapshots must not expose test-only diagnostic storage");

static void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
static void close(double a, double b) {
  require(std::isfinite(a) && std::isfinite(b) && std::abs(a - b) < 1e-10,
          "runtime snapshot numeric mismatch");
}

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", SCFG2_TEST_PARAM_FILE, 1);
    unsetenv("CPARTY_PF_OVERRIDE_SCALE");
    unsetenv("CPARTY_PF_EXACT_W_PREFIX_SCALE_SEED");
    w::RuntimeSnapshot log("gggaaatcc", "(((...)))");
    require(log.sequence() == "GGGAAAUCC" && log.scaffold() == "(((...)))",
            "runtime sequence normalization changed the scaffold");
    st::DiagnosticTables disabled(log.length(), false);
    for (std::size_t id = 0; id < static_cast<std::size_t>(st::DiagnosticTable::Count); ++id)
      require(disabled.table(static_cast<st::DiagnosticTable>(id)).empty(),
              "disabled observations allocated quadratic diagnostic tables");
    require(log.encoded_base(1) == log.encoded_base(3) && log.bases_can_pair(1, 9),
            "encoded bases cannot be used for pairability");
    close(log.partition_free_energy_from_log_weight(0.), 0.);
    // Obsolete numerical-scaling controls must not affect any snapshot or
    // the physical length penalties. In particular no MFE seed is needed.
    setenv("CPARTY_PF_OVERRIDE_SCALE", "7.25", 1);
    setenv("CPARTY_PF_EXACT_W_PREFIX_SCALE_SEED", "1", 1);
    st::ObservedRuntimeSnapshot observed("gggaaatcc", "(((...)))");
    for (double w::PseudoknotFactors::*member : {&w::PseudoknotFactors::ps,
         &w::PseudoknotFactors::psm, &w::PseudoknotFactors::psp,
         &w::PseudoknotFactors::pb, &w::PseudoknotFactors::pup,
         &w::PseudoknotFactors::pps, &w::PseudoknotFactors::a,
         &w::PseudoknotFactors::b, &w::PseudoknotFactors::c,
         &w::PseudoknotFactors::ap, &w::PseudoknotFactors::bp, &w::PseudoknotFactors::cp})
      close(observed.log_pk_parameters().*member, log.log_pk_parameters().*member);
    close(observed.log_pf_parameters().exphairpin[3], log.log_pf_parameters().exphairpin[3]);
    close(observed.partition_free_energy_from_log_weight(1234.),
          log.partition_free_energy_from_log_weight(1234.));
    close(observed.partition_free_energy_from_log_weight(1234.),
          -1234. * observed.parameters().kT / 1000.);
    unsetenv("CPARTY_PF_OVERRIDE_SCALE");
    unsetenv("CPARTY_PF_EXACT_W_PREFIX_SCALE_SEED");

    auto &tables = observed.diagnostics();
    for (std::size_t id = 0; id < static_cast<std::size_t>(st::DiagnosticTable::Count); ++id)
      require(!tables.table(static_cast<st::DiagnosticTable>(id)).empty(),
              "explicit observations did not materialize diagnostic tables");
    close(tables.prefix()[0], 1.);
    close(tables.get(st::DiagnosticTable::WI, 1, observed.length()), 1.);
    std::size_t value = 1;
    for (int i = 1; i <= observed.length(); ++i) for (int j = i + 1; j <= observed.length(); ++j)
      tables.set(st::DiagnosticTable::VP, i, j, static_cast<double>(value++));
    value = 1;
    for (int i = 1; i <= observed.length(); ++i) for (int j = i + 1; j <= observed.length(); ++j)
      close(tables.get(st::DiagnosticTable::VP, i, j), static_cast<double>(value++));
    tables.set(st::DiagnosticTable::VP_CASE3_R, 4, 4, 3.);
    close(tables.get(st::DiagnosticTable::VP_CASE3_R, 4, 4), 3.);
    close(tables.get(st::DiagnosticTable::WI, 5, 4), 1.);
    observed.reset_diagnostics();
    close(tables.get(st::DiagnosticTable::VP_CASE3_R, 4, 4), 0.);
    close(tables.get(st::DiagnosticTable::VP, 1, 9), 0.);
    close(tables.prefix()[0], 1.);

    const auto saved = w::pseudoknot_energy_defaults();
    w::pseudoknot_energy_defaults().exterior += 100;
    w::pseudoknot_energy_defaults().stack_multiplier = 0.5;
    w::RuntimeSnapshot later("GGGAAACCC", "(((...)))");
    w::pseudoknot_energy_defaults() = saved;
    require(log.pseudoknot_model().stack_multiplier == 0.89 &&
            later.pseudoknot_model().stack_multiplier == 0.5 &&
            std::abs(later.log_pk_parameters().ps - log.log_pk_parameters().ps) > 0.1,
            "model settings were not isolated between snapshots");

    bool rejected = false;
    try { w::RuntimeSnapshot invalid("AAA", ".."); }
    catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "mismatched input lengths accepted");
    rejected = false;
    try { w::RuntimeSnapshot invalid("AAAA", "[..]"); }
    catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "crossing input notation accepted as scaffold");

    const std::string long_sequence(70000, 'A'), long_scaffold(70000, '.');
    w::RuntimeSnapshot long_log(long_sequence, long_scaffold);
    require(long_log.length() == 70000 && long_log.encoded_base(70000) == long_log.encoded_base(1),
            "runtime sequence has an old 16-bit coordinate limit");
    rejected = false;
    try { st::ObservedRuntimeSnapshot invalid(long_sequence, long_scaffold); }
    catch (const std::length_error &) { rejected = true; }
    require(rejected, "overflowing diagnostic offsets accepted");
    std::cout << "scfg2_runtime_snapshot=ok\n";
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
