#include "scfg2/cli/commands.hh"
// linearknot profile: per-base structural-context probabilities for the
// pseudoknot ensemble of G' given a scaffold G (LinearKnot).
//
// Input: records of an optional ">name" line, a sequence line and a scaffold
// dot-bracket line.  Output follows CapR / LinearCapR: ">name", then one line
// per context with n values. PseudoLoop labels gaps between crossing bands;
// multiloop wings within a band are Multiloop.

#include "scfg2/api/profiles.hh"
#include "scfg2/cli/cli_common.hh"
#include "scfg2/cli/ranking_options.hh"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace {

using cparty::scfg2::cli::Record;

struct Options {
  std::string input;
  std::string output;
  std::string sequence;
  std::string scaffold;
  std::size_t beam = 100;
  cparty::scfg2::cli::RankingOptions ranking;
  bool exact = false;
  std::size_t max_hairpin_unpaired = 30;
  std::size_t max_case3_ancestor_depth = 30;
};

void usage(std::ostream &os) {
  os << "Usage: linearknot profile [options] (-i FILE | --sequence SEQ --scaffold G)\n"
        "  -i, --input FILE         records: [>name] / sequence / scaffold G\n"
        "  -o, --output FILE        output file (default: stdout)\n"
        "  -b, --beam N             beam size (default 100)\n"
        "  --exact                  exact cubic-time computation (no beam)\n"
        "  --max-hairpin N          beam hairpin length cap (default 30)\n"
        "  --max-case3-depth N      beam CASE3 ancestor cap (default 30)\n";
  cparty::scfg2::cli::RankingOptions::usage(os);
}

bool parse_args(int argc, char **argv, Options *opt) {
  for (int a = 1; a < argc; ++a) {
    const std::string arg = argv[a];
    const auto value = [&]() -> const char * {
      return a + 1 < argc ? argv[++a] : nullptr;
    };
    const char *v = nullptr;
    if (cparty::scfg2::cli::RankingOptions::recognizes(arg)) {
      std::string error;
      if (!opt->ranking.parse(arg, value(), &error)) {
        std::cerr << "error: " << error << '\n';
        return false;
      }
    } else if (arg == "-h" || arg == "--help") {
      usage(std::cout);
      std::exit(EXIT_SUCCESS);
    } else if (arg == "-i" || arg == "--input") {
      if (!(v = value())) return false;
      opt->input = v;
    } else if (arg == "-o" || arg == "--output") {
      if (!(v = value())) return false;
      opt->output = v;
    } else if (arg == "--sequence") {
      if (!(v = value())) return false;
      opt->sequence = v;
    } else if (arg == "--scaffold") {
      if (!(v = value())) return false;
      opt->scaffold = v;
    } else if (arg == "-b" || arg == "--beam") {
      if (!(v = value())) return false;
      if (!cparty::scfg2::cli::positive_beam_argument(v, &opt->beam)) {
        std::cerr << "error: beam size must be a positive integer\n";
        return false;
      }
    } else if (arg == "--exact") {
      opt->exact = true;
    } else if (arg == "--max-hairpin") {
      if (!(v = value())) return false;
      opt->max_hairpin_unpaired = std::strtoul(v, nullptr, 10);
    } else if (arg == "--max-case3-depth") {
      if (!(v = value())) return false;
      opt->max_case3_ancestor_depth = std::strtoul(v, nullptr, 10);
    } else {
      std::cerr << "unknown option: " << arg << "\n";
      return false;
    }
  }
  std::string ranking_error;
  if (!opt->ranking.validate(opt->exact, &ranking_error)) {
    std::cerr << "error: " << ranking_error << '\n';
    return false;
  }
  if (opt->input.empty() == opt->sequence.empty()) return false;
  if (!opt->sequence.empty() && opt->scaffold.empty()) return false;
  return true;
}

// CapR line order, with PseudoLoop added after Multiloop.
struct Line {
  const char *label;
  std::size_t context;
};
constexpr Line kLines[] = {
    {"Bulge", 2},     {"Exterior", 5},   {"Hairpin", 1}, {"Internal", 3},
    {"Multiloop", 4}, {"PseudoLoop", 6}, {"Stem", 0},
};

}  // namespace

int cparty::scfg2::cli::run_profile(int argc, char **argv) {
  Options opt;
  if (!parse_args(argc, argv, &opt)) {
    usage(std::cerr);
    return EXIT_FAILURE;
  }
  std::vector<Record> records;
  if (!opt.input.empty()) {
    std::string error;
    if (!cparty::scfg2::cli::read_records(opt.input, &records, &error)) {
      std::cerr << "error: " << error << "\n";
      return EXIT_FAILURE;
    }
  } else {
    records.push_back(Record{"seq1", opt.sequence, opt.scaffold});
  }
  std::ofstream file;
  if (!opt.output.empty()) {
    file.open(opt.output);
    if (!file) {
      std::cerr << "error: cannot write " << opt.output << "\n";
      return EXIT_FAILURE;
    }
  }
  std::ostream &os = opt.output.empty() ? std::cout : file;
  os << std::setprecision(std::numeric_limits<double>::max_digits10);
  int failures = 0;
  for (const auto &record : records) {
    cparty::scfg2::api::StructuralProfileResult result;
    if (opt.exact) {
      result = cparty::scfg2::api::exact_structural_profiles(
          record.sequence, record.scaffold, {});
    } else {
      cparty::scfg2::api::BeamConfig config;
      config.beam_size = opt.beam;
      config.max_hairpin_unpaired = opt.max_hairpin_unpaired;
      config.max_case3_ancestor_depth = opt.max_case3_ancestor_depth;
      opt.ranking.apply(config);
      result = cparty::scfg2::api::beam_structural_profiles(
          record.sequence, record.scaffold, {}, config);
    }
    if (!result.ok) {
      std::cerr << "error: " << record.name << ": " << result.error << "\n";
      ++failures;
      continue;
    }
    os << ">" << record.name << "\n";
    for (const auto &line : kLines) {
      os << line.label;
      for (std::size_t pos = 1; pos < result.contexts.size(); ++pos) {
        double p = result.contexts[pos][line.context];
        if (p < 0.0) p = 0.0;
        if (p > 1.0) p = 1.0;
        os << " " << p;
      }
      os << "\n";
    }
    os << "\n";
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
