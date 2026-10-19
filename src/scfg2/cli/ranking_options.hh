#ifndef CPARTY_SCFG2_CLI_RANKING_OPTIONS_HH
#define CPARTY_SCFG2_CLI_RANKING_OPTIONS_HH

#include "scfg2/api/beam.hh"
#include <cerrno>
#include <cstdlib>
#include <limits>
#include <ostream>
#include <string>

namespace cparty::scfg2::cli {

struct RankingOptions {
  api::BeamRanking mode = api::BeamRanking::Outer;
  api::OuterPilotExponent exponent = api::OuterPilotExponent::TwoThirds;
  bool explicit_mode = false, explicit_exponent = false;

  static bool recognizes(const std::string &arg) {
    return arg == "--ranking" || arg == "--outer-pilot-exponent";
  }
  bool parse(const std::string &arg, const char *raw, std::string *error) {
    if (!raw) { *error = "missing value for " + arg; return false; }
    const std::string value(raw);
    if (arg == "--ranking") {
      explicit_mode = true;
      if (value == "outer") mode = api::BeamRanking::Outer;
      else if (value == "baseline") mode = api::BeamRanking::Baseline;
      else { *error = "--ranking must be outer or baseline"; return false; }
    } else {
      explicit_exponent = true;
      if (value == "0.5") exponent = api::OuterPilotExponent::Half;
      else if (value == "2/3") exponent = api::OuterPilotExponent::TwoThirds;
      else { *error = "--outer-pilot-exponent must be 0.5 or 2/3"; return false; }
    }
    return true;
  }
  bool validate(bool exact, std::string *error) const {
    if (exact && (explicit_mode || explicit_exponent)) {
      *error = "ranking options cannot be used with --exact";
      return false;
    }
    if (mode == api::BeamRanking::Baseline && explicit_exponent) {
      *error = "--outer-pilot-exponent requires --ranking outer";
      return false;
    }
    return true;
  }
  void apply(api::BeamConfig &config) const {
    config.ranking = mode;
    config.outer_pilot_exponent = exponent;
  }
  const char *name() const { return mode == api::BeamRanking::Outer ? "outer" : "baseline"; }
  const char *exponent_name() const { return exponent == api::OuterPilotExponent::Half ? "0.5" : "2/3"; }
  static void usage(std::ostream &os) {
    os << "  --ranking outer|baseline beam ranking (default outer)\n"
          "  --outer-pilot-exponent Q pilot exponent: 0.5 or 2/3 (default 2/3)\n";
  }
};

inline bool positive_beam_argument(const char *raw, std::size_t *out) {
  if (!raw || *raw < '0' || *raw > '9') return false;
  char *end = nullptr;
  errno = 0;
  const auto value = std::strtoull(raw, &end, 10);
  if (errno == ERANGE || *end || !value || value > std::numeric_limits<std::size_t>::max())
    return false;
  *out = static_cast<std::size_t>(value);
  return true;
}

} // namespace cparty::scfg2::cli
#endif
