#include "scfg2/cli/commands.hh"
// Posterior decoding from one fixed-scaffold retained ensemble.
#include "scfg2/api/mea.hh"
#include "scfg2/api/centroid.hh"
#include "scfg2/cli/cli_common.hh"
#include "scfg2/cli/ranking_options.hh"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace {


struct Options {
  std::string input, output, sequence, scaffold;
  double gamma = 1.0;
  cparty::scfg2::api::BeamConfig beam;
  cparty::scfg2::cli::RankingOptions ranking;
  bool exact = false;
  bool stats = false;
};

void usage(std::ostream &os, bool centroid) {
  const char *program = centroid ? "linearknot centroid" : "linearknot mea";
  os << "Usage: " << program << " [options] (-i FILE | --sequence SEQ --scaffold G)\n"
        "  -i, --input FILE         records: [>name] / sequence / scaffold G\n"
        "  -o, --output FILE        output file (default: stdout)\n"
        "  --gamma G               nonnegative pair reward multiplier (default 1)\n"
        "  -b, --beam N             beam size (default 100)\n"
        "  --exact                  uncapped exact ensemble and decoding\n"
        "  --max-hairpin N          beam hairpin length cap (default 30)\n"
        "  --max-case3-depth N      beam CASE3 ancestor cap (default 30)\n"
        "  --stats                  report timing and decoding work to stderr\n";
  cparty::scfg2::cli::RankingOptions::usage(os);
}

std::size_t size_argument(const std::string &value) {
  std::size_t consumed = 0;
  if (value.empty() || value[0] == '-') throw std::invalid_argument("invalid size");
  const auto size = std::stoull(value, &consumed);
  if (consumed != value.size() || size > std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument("invalid size");
  return static_cast<std::size_t>(size);
}

bool parse(int argc, char **argv, Options &opt, bool centroid) {
  try {
    for (int a = 1; a < argc; ++a) {
      const std::string arg = argv[a];
      const auto value = [&]() -> std::string {
        if (++a >= argc) throw std::invalid_argument("missing value for " + arg);
        return argv[a];
      };
      if (cparty::scfg2::cli::RankingOptions::recognizes(arg)) {
        const auto text = value();
        std::string error;
        if (!opt.ranking.parse(arg, text.c_str(), &error)) throw std::invalid_argument(error);
      } else if (arg == "-h" || arg == "--help") { usage(std::cout, centroid); std::exit(EXIT_SUCCESS); }
      else if (arg == "-i" || arg == "--input") opt.input = value();
      else if (arg == "-o" || arg == "--output") opt.output = value();
      else if (arg == "--sequence") opt.sequence = value();
      else if (arg == "--scaffold") opt.scaffold = value();
      else if (arg == "--gamma") {
        const auto text = value();
        std::size_t consumed = 0;
        opt.gamma = std::stod(text, &consumed);
        if (consumed != text.size() || !std::isfinite(opt.gamma) || opt.gamma < 0.0)
          throw std::invalid_argument("gamma must be finite and nonnegative");
      } else if (arg == "-b" || arg == "--beam") {
        opt.beam.beam_size = size_argument(value());
        if (!opt.beam.beam_size) throw std::invalid_argument("beam size must be positive");
      } else if (arg == "--max-hairpin") opt.beam.max_hairpin_unpaired = size_argument(value());
      else if (arg == "--max-case3-depth") opt.beam.max_case3_ancestor_depth = size_argument(value());
      else if (arg == "--exact") opt.exact = true;
      else if (arg == "--stats") opt.stats = true;
      else throw std::invalid_argument("unknown option: " + arg);
    }
    std::string ranking_error;
    if (!opt.ranking.validate(opt.exact, &ranking_error)) throw std::invalid_argument(ranking_error);
    opt.ranking.apply(opt.beam);
    return opt.input.empty() != opt.sequence.empty() &&
           (opt.sequence.empty() || !opt.scaffold.empty());
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << '\n';
    return false;
  }
}
int run_posterior_decoding(int argc, char **argv, bool centroid) {
  Options opt;
  if (!parse(argc, argv, opt, centroid)) { usage(std::cerr, centroid); return EXIT_FAILURE; }
  std::vector<cparty::scfg2::cli::Record> records;
  if (opt.input.empty()) records.push_back({"seq1", opt.sequence, opt.scaffold});
  else {
    std::string error;
    if (!cparty::scfg2::cli::read_records(opt.input, &records, &error)) {
      std::cerr << "error: " << error << '\n';
      return EXIT_FAILURE;
    }
  }
  std::ofstream file;
  if (!opt.output.empty()) {
    file.open(opt.output);
    if (!file) { std::cerr << "error: cannot write " << opt.output << '\n'; return EXIT_FAILURE; }
  }
  std::ostream &os = opt.output.empty() ? std::cout : file;
  int failures = 0;
  for (const auto &record : records) {
    const auto result = centroid ? (opt.exact
        ? cparty::scfg2::api::exact_centroid(record.sequence, record.scaffold, opt.gamma)
        : cparty::scfg2::api::beam_centroid(record.sequence, record.scaffold, opt.gamma, opt.beam))
        : (opt.exact
        ? cparty::scfg2::api::exact_mea(record.sequence, record.scaffold, opt.gamma)
        : cparty::scfg2::api::beam_mea(record.sequence, record.scaffold, opt.gamma, opt.beam));
    if (!result.ok || !result.reachable) {
      std::cerr << "error: " << record.name << ": "
                << (result.error.empty() ? "no retained structure is reachable" : result.error) << '\n';
      ++failures;
      continue;
    }
    os << '>' << record.name << '\n' << record.sequence << '\n' << record.scaffold << '\n'
       << result.union_structure << (centroid ? " (EG=" : " (EA=") << std::setprecision(17)
       << result.expected_accuracy << "; gamma=" << result.gamma;
    if (centroid)
      os << "; ETP=" << result.expected_true_positive_pairs
         << "; ETN=" << result.expected_true_negative_pairs;
    os << ")\n\n";
    if (opt.stats)
      std::cerr << std::setprecision(17) << "record=" << record.name
                << " n=" << record.sequence.size() << " beam=" << result.beam_size
                << " ranking=" << (opt.exact ? "exact" : opt.ranking.name())
                << " outer_pilot_exponent=" << (opt.exact || opt.ranking.mode != cparty::scfg2::api::BeamRanking::Outer
                    ? "-" : opt.ranking.exponent_name())
                << " gamma=" << result.gamma << " log_z=" << result.log_partition_weight
                << " states=" << result.decoding_states << " terms=" << result.decoding_terms
                << " factors=" << result.decoding_factors << " seconds=" << result.elapsed_seconds
                << " decoding_seconds=" << result.decoding_seconds << '\n';
  }
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}

}  // namespace

int cparty::scfg2::cli::run_mea(int argc, char **argv) {
  return run_posterior_decoding(argc, argv, false);
}

int cparty::scfg2::cli::run_centroid(int argc, char **argv) {
  return run_posterior_decoding(argc, argv, true);
}
