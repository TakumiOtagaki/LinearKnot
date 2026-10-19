#include "scfg2/cli/commands.hh"
// linearknot access: accessibility of every window of the requested lengths
// for the pseudoknot ensemble of G' given a scaffold G (LinearKnot).
//
// Input: records of an optional ">name" line, a sequence line and a scaffold
// dot-bracket line.  Output follows LinearRaccess / Raccess: for each
// zero-based start position, "len,ED;" for every requested length, where
// ED = -RT log P(window unpaired).  --probability prints P instead.
// --by-loop appends joint probabilities for six unpaired loop contexts.

#include "scfg2/api/profiles.hh"
#include "scfg2/cli/cli_common.hh"
#include "scfg2/cli/ranking_options.hh"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

using cparty::scfg2::cli::Record;

struct Options {
  std::string input;
  std::string output;
  std::string sequence;
  std::string scaffold;
  std::vector<int> access_lengths;
  std::size_t beam = 100;
  cparty::scfg2::cli::RankingOptions ranking;
  bool exact = false;
  bool probability = false;
  bool by_loop = false;
  std::size_t max_hairpin_unpaired = 30;
  std::size_t max_case3_ancestor_depth = 30;
};

void usage(std::ostream &os) {
  os << "Usage: linearknot access --access-len 1,5,10 [options] (-i FILE | "
        "--sequence SEQ --scaffold G)\n"
        "  -i, --input FILE         records: [>name] / sequence / scaffold G\n"
        "  -o, --output FILE        output file (default: stdout)\n"
        "  --access-len L1,L2,...   window lengths (required)\n"
        "  -b, --beam N             beam size (default 100)\n"
        "  --exact                  exact cubic-time computation (no beam)\n"
        "  --probability            print P instead of ED = -RT log P\n"
        "  --by-loop                append joint loop probabilities E,H,B,I,M,P\n"
        "                           aliases: --byloop, -byloop\n"
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
    } else if (arg == "--access-len") {
      if (!(v = value()) || !cparty::scfg2::cli::parse_lengths(v, &opt->access_lengths)) return false;
    } else if (arg == "-b" || arg == "--beam") {
      if (!(v = value())) return false;
      if (!cparty::scfg2::cli::positive_beam_argument(v, &opt->beam)) {
        std::cerr << "error: beam size must be a positive integer\n";
        return false;
      }
    } else if (arg == "--exact") {
      opt->exact = true;
    } else if (arg == "--probability") {
      opt->probability = true;
    } else if (arg == "--by-loop" || arg == "--byloop" || arg == "-byloop") {
      opt->by_loop = true;
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
  if (opt->access_lengths.empty()) return false;
  std::string ranking_error;
  if (!opt->ranking.validate(opt->exact, &ranking_error)) {
    std::cerr << "error: " << ranking_error << '\n';
    return false;
  }
  if (opt->input.empty() == opt->sequence.empty()) return false;
  if (!opt->sequence.empty() && opt->scaffold.empty()) return false;
  return true;
}

void write_metadata(std::ostream &os, const Options &opt) {
  os << "# linearknot_access_metadata=1\n";
  os << "# engine=linearknot-scfg2\n";
  os << "# ensemble=conditional_on_scaffold_G\n";
  os << "# beam=" << (opt.exact ? std::string("exact") : std::to_string(opt.beam))
     << "\n";
  if (!opt.exact) {
    os << "# beam_ranking=" << opt.ranking.name() << "\n";
    if (opt.ranking.mode == cparty::scfg2::api::BeamRanking::Outer)
      os << "# outer_pilot_exponent=" << opt.ranking.exponent_name() << "\n";
    os << "# max_hairpin=" << opt.max_hairpin_unpaired << "\n";
    os << "# max_case3_depth=" << opt.max_case3_ancestor_depth << "\n";
  }
  os << "# output_mode=" << (opt.probability ? "probability" : "energy") << "\n";
  if (opt.by_loop) {
    os << "# byloop_value=joint_probability\n";
    os << "# byloop_contexts=E:Exterior,H:Hairpin,B:Bulge,I:Internal,M:Multi,P:PseudoLoop\n";
  }
  os << "# access_lens=";
  for (std::size_t i = 0; i < opt.access_lengths.size(); ++i) {
    os << (i ? "," : "") << opt.access_lengths[i];
  }
  os << "\n";
}

void write_record(std::ostream &os, const Record &record, const Options &opt,
                  const cparty::scfg2::api::StructuralProfileResult &result) {
  const int n = static_cast<int>(record.sequence.size());
  os << ">" << record.name << "\n";
  os << "# scaffold=" << record.scaffold << "\n";
  os << "# ensemble_free_energy=" << result.ensemble_free_energy << "\n";
  os << "# folding_rt=" << result.rt_kcal_mol << "\n";
  os << "# elapsed_seconds=" << result.elapsed_seconds << "\n";
  if (opt.probability) os << "# value=probability\n";
  for (int pos = 0; pos < n; ++pos) {
    os << pos;
    bool first = true;
    for (std::size_t w = 0; w < opt.access_lengths.size(); ++w) {
      const int len = opt.access_lengths[w];
      if (pos + len > n) continue;
      double p = result.accessibility[w][static_cast<std::size_t>(pos + 1)];
      p = std::min(1.0, std::max(0.0, p));
      os << (first ? "\t" : "") << len << ",";
      if (opt.probability) {
        os << p;
      } else if (p > 0.0) {
        os << -result.rt_kcal_mol * std::log(p);
      } else {
        os << "inf";
      }
      os << ";";
      first = false;
    }
    if (opt.by_loop && !first) {
      // Match LinearRaccess's extra-column layout. The total follows the
      // selected mode; loop contributions are always joint probabilities.
      constexpr std::array<std::size_t, 6> order{4, 0, 1, 2, 3, 5};
      constexpr std::array<char, 6> codes{'E', 'H', 'B', 'I', 'M', 'P'};
      os << "\tbyloop:";
      for (std::size_t w = 0; w < opt.access_lengths.size(); ++w) {
        const int len = opt.access_lengths[w];
        if (pos + len > n) continue;
        os << len;
        const auto &row = result.accessibility_by_loop[w][static_cast<std::size_t>(pos + 1)];
        for (std::size_t c = 0; c < order.size(); ++c)
          os << "," << codes[c] << "=" << std::min(1.0, std::max(0.0, row[order[c]]));
        os << ";";
      }
    }
    os << "\n";
  }
  os << "\n";
}

}  // namespace

int cparty::scfg2::cli::run_access(int argc, char **argv) {
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
  write_metadata(os, opt);
  int failures = 0;
  for (const auto &record : records) {
    cparty::scfg2::api::StructuralProfileResult result;
    if (opt.exact) {
      result = cparty::scfg2::api::exact_structural_profiles(
          record.sequence, record.scaffold, opt.access_lengths,
          cparty::scfg2::WmbClosurePolicy::RequireWeakClosure, opt.by_loop);
    } else {
      cparty::scfg2::api::BeamConfig config;
      config.beam_size = opt.beam;
      config.max_hairpin_unpaired = opt.max_hairpin_unpaired;
      config.max_case3_ancestor_depth = opt.max_case3_ancestor_depth;
      opt.ranking.apply(config);
      result = cparty::scfg2::api::beam_structural_profiles(
          record.sequence, record.scaffold, opt.access_lengths, config, opt.by_loop);
    }
    if (!result.ok) {
      std::cerr << "error: " << record.name << ": " << result.error << "\n";
      ++failures;
      continue;
    }
    write_record(os, record, opt, result);
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
