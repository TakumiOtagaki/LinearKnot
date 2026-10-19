#include "scfg2/cli/commands.hh"
#include "scfg2/api/subopt.hh"
#include "scfg2/cli/cli_common.hh"
#include "scfg2/cli/ranking_options.hh"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
struct Options {
  std::string input, output, sequence, scaffold;
  cparty::scfg2::api::SuboptConfig subopt;
  cparty::scfg2::api::BeamConfig beam;
  cparty::scfg2::cli::RankingOptions ranking;
  bool exact = false, stats = false;
};

void usage(std::ostream &os) {
  os << "Usage: linearknot subopt [options] (-i FILE | --sequence SEQ --scaffold G)\n"
        "  -i, --input FILE         records: [>name] / sequence / scaffold G\n"
        "  -o, --output FILE        output file (default stdout)\n"
        "  -k, --top-k K            maximum structures, in energy order (default 10)\n"
        "  -e, --delta-energy E     also stop above this forest's MFE + E kcal/mol\n"
        "  -b, --beam N             beam size (default 100)\n"
        "  --exact                  full grammar, without beam or beam length caps\n"
        "  --max-hairpin N          beam hairpin length cap (default 30)\n"
        "  --max-case3-depth N      beam CASE3 ancestor cap (default 30)\n"
        "  --stats                  report scope, MFE, stop reason and work to stderr\n";
  cparty::scfg2::cli::RankingOptions::usage(os);
}

std::size_t size_argument(const std::string &text) {
  if (text.empty() || text.find_first_not_of("0123456789") != std::string::npos)
    throw std::invalid_argument("expected a nonnegative integer");
  std::size_t consumed = 0;
  const auto value = std::stoull(text, &consumed);
  if (consumed != text.size() || value > std::numeric_limits<std::size_t>::max())
    throw std::invalid_argument("integer is too large");
  return static_cast<std::size_t>(value);
}

bool parse(int argc, char **argv, Options &opt) {
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
      } else if (arg == "--help" || arg == "-h") { usage(std::cout); std::exit(EXIT_SUCCESS); }
      else if (arg == "-i" || arg == "--input") opt.input = value();
      else if (arg == "-o" || arg == "--output") opt.output = value();
      else if (arg == "--sequence") opt.sequence = value();
      else if (arg == "--scaffold") opt.scaffold = value();
      else if (arg == "-k" || arg == "--top-k") {
        opt.subopt.top_k = size_argument(value());
        if (!opt.subopt.top_k) throw std::invalid_argument("top_k must be positive");
      } else if (arg == "-e" || arg == "--delta-energy") {
        const auto text = value();
        std::size_t consumed = 0;
        const double delta = std::stod(text, &consumed);
        if (consumed != text.size() || !std::isfinite(delta) || delta < 0.0)
          throw std::invalid_argument("delta energy must be finite and nonnegative");
        opt.subopt.delta_energy = delta;
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
    if (!opt.input.empty()) return opt.sequence.empty() && opt.scaffold.empty();
    return !opt.sequence.empty() && !opt.scaffold.empty();
  } catch (const std::exception &e) {
    std::cerr << "error: " << e.what() << '\n';
    return false;
  }
}
}

int cparty::scfg2::cli::run_subopt(int argc, char **argv) {
  Options opt;
  if (!parse(argc, argv, opt)) { usage(std::cerr); return EXIT_FAILURE; }
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
    bool header_written = false;
    const auto consume = [&](cparty::scfg2::api::SuboptStructure &&structure) {
      if (!header_written) {
        os << '>' << record.name << '\n' << record.sequence << '\n' << record.scaffold << '\n';
        header_written = true;
      }
      os << structure.union_structure << " (" << std::setprecision(17) << structure.energy << ")\n";
      if (!os) throw std::runtime_error("cannot write suboptimal structure");
      return true;
    };
    const auto result = opt.exact
        ? cparty::scfg2::api::exact_subopt_each(record.sequence, record.scaffold, consume, opt.subopt)
        : cparty::scfg2::api::beam_subopt_each(record.sequence, record.scaffold, consume, opt.subopt, opt.beam);
    if (!result.ok || !result.reachable) {
      std::cerr << "error: " << record.name << ": "
                << (result.error.empty() ? "no retained structure is reachable" : result.error) << '\n';
      ++failures;
      continue;
    }
    os << '\n';
    if (opt.stats)
      std::cerr << std::setprecision(17) << "record=" << record.name
                << " scope=" << (opt.exact ? "exact" : "inside_retained_beam")
                << " ranking=" << (opt.exact ? "exact" : opt.ranking.name())
                << " outer_pilot_exponent=" << (opt.exact || opt.ranking.mode != cparty::scfg2::api::BeamRanking::Outer
                    ? "-" : opt.ranking.exponent_name())
                << " beam=" << result.beam_size << " mfe=" << result.mfe_energy
                << " returned=" << result.structures_emitted
                << " stop=" << cparty::scfg2::api::subopt_stop_reason_name(result.stop_reason)
                << " log_z=" << result.log_partition_weight
                << " states=" << result.graph_states << " terms=" << result.graph_terms
                << " factors=" << result.graph_factors
                << " candidates=" << result.candidates_created
                << " cached=" << result.derivations_cached
                << " seconds=" << result.elapsed_seconds << '\n';
  }
  return failures ? EXIT_FAILURE : EXIT_SUCCESS;
}
