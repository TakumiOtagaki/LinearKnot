#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void print_usage() {
  std::cerr
      << "usage: scfg2_restricted_partition_function_tool [--details] [--legacy-compatible] <seq> <structure> [param_file]\n";
}

}  // namespace

int main(int argc, char **argv) {
  if (argc < 3) {
    print_usage();
    return 2;
  }

  int argi = 1;
  bool details = false;
  auto policy = cparty::scfg2::WmbClosurePolicy::RequireWeakClosure;
  while (argi < argc) {
    const std::string arg = argv[argi];
    if (arg == "--details") {
      details = true;
      ++argi;
      continue;
    }
    if (arg == "--legacy-compatible") {
      policy = cparty::scfg2::WmbClosurePolicy::LegacyCompatible;
      ++argi;
      continue;
    }
    break;
  }

  const int remaining = argc - argi;
  if (remaining != 2 && remaining != 3) {
    print_usage();
    return 2;
  }

  try {
    const std::string seq = argv[argi++];
    const std::string restricted = argv[argi++];
    if (remaining == 3) {
      setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", argv[argi], 1);
    }

        cparty::scfg2::weights::RuntimeSnapshot pf(
      seq, restricted);
    namespace scfg2 = cparty::scfg2;
    scfg2::replay::WFinalLogWeightContext ctx(pf, policy);
    const auto n = static_cast<cand_pos_t>(pf.length());
    const auto chart = scfg2::replay::run_w_final_exact_dp_schedule<
        scfg2::semirings::LogSumExpSemiring>(
            scfg2::engine::build_pkfree_basic_exact_schedule(n, ctx), ctx,
            [](const auto &) { return true; },
            [](const auto &, const auto &, const auto &, const auto &, auto) {});
    const double log_z = chart.get({scfg2::grammar::NonTerminal::W, 1, n, -1, -1});
    if (!std::isfinite(log_z)) {
      throw std::runtime_error("scfg2 constrained PF returned non-finite log total");
    }
    const double pf_energy = pf.partition_free_energy_from_log_weight(log_z);
    const double rt_kcal = pf.parameters().kT / 1000.0;

    std::cout << std::setprecision(17);
    if (details) {
      std::cout << "pf_energy=" << pf_energy
                << "\tlog_partition_scaled=" << log_z
                << "\tpf_scale=1"
                << "\trt_kcal=" << rt_kcal
                << "\tlength=" << pf.length()
                << "\twmb_closure_policy=" << scfg2::wmb_closure_policy_name(policy)
                << "\n";
    } else {
      std::cout << pf_energy << "\n";
    }
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "scfg2_restricted_partition_function_tool failed: " << e.what() << "\n";
    return 1;
  }
}
