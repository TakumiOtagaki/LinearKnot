#include "scfg2/cli/commands.hh"
// linearknot fold: MFE structure, ensemble free energy and base-pair
// probabilities for the pseudoknot ensemble of G' given a scaffold G.
//
// Input: records of an optional ">name" line, a sequence line and a scaffold
// dot-bracket line.  For each record it prints the sequence, G, the MFE
// union G u G' with its energy, and the ensemble free energy.  --bpp FILE
// writes "i j p" lines (one-based, i < j) per record, LinearPartition style;
// pairs of G have probability 1.

#include "scfg2/api/beam.hh"
#include "scfg2/api/exact.hh"
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
  std::string sequence;
  std::string scaffold;
  std::string bpp;
  double bpp_cutoff = 0.0;
  std::size_t beam = 100;
  cparty::scfg2::cli::RankingOptions ranking;
  bool exact = false;
  std::size_t max_hairpin_unpaired = 30;
  std::size_t max_case3_ancestor_depth = 30;
};

void usage(std::ostream &os) {
  os << "Usage: linearknot fold [options] (-i FILE | --sequence SEQ --scaffold G)\n"
        "  -i, --input FILE         records: [>name] / sequence / scaffold G\n"
        "  -b, --beam N             beam size (default 100)\n"
        "  --exact                  exact cubic-time computation (no beam)\n"
        "  --bpp FILE               write base-pair probabilities to FILE\n"
        "  --bpp-cutoff P           only write pairs with probability >= P\n"
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
    } else if (arg == "--sequence") {
      if (!(v = value())) return false;
      opt->sequence = v;
    } else if (arg == "--scaffold") {
      if (!(v = value())) return false;
      opt->scaffold = v;
    } else if (arg == "--bpp") {
      if (!(v = value())) return false;
      opt->bpp = v;
    } else if (arg == "--bpp-cutoff") {
      if (!(v = value())) return false;
      opt->bpp_cutoff = std::strtod(v, nullptr);
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

struct Prediction {
  bool ok = false;
  std::string error;
  std::string union_structure;
  double mfe = std::numeric_limits<double>::quiet_NaN();
  double ensemble_free_energy = std::numeric_limits<double>::quiet_NaN();
  std::vector<cparty::scfg2::api::PairProbability> pairs;
};

Prediction predict(const Record &record, const Options &opt, bool want_bpp) {
  namespace api = cparty::scfg2::api;
  Prediction out;
  if (opt.exact) {
    api::ExactSession session(record.sequence, record.scaffold);
    const auto best = session.viterbi();
    if (!best.ok) {
      out.error = best.error;
      return out;
    }
    out.union_structure = best.union_structure;
    out.mfe = best.energy;
    if (want_bpp) {
      const auto outside = session.outside();
      if (!outside.ok) {
        out.error = outside.error;
        return out;
      }
      out.ensemble_free_energy = outside.ensemble_free_energy;
      out.pairs = outside.pair_probabilities;
    } else {
      const auto inside = session.inside();
      if (!inside.ok) {
        out.error = inside.error;
        return out;
      }
      out.ensemble_free_energy = inside.ensemble_free_energy;
    }
  } else {
    api::BeamConfig config;
    config.beam_size = opt.beam;
    config.max_hairpin_unpaired = opt.max_hairpin_unpaired;
    config.max_case3_ancestor_depth = opt.max_case3_ancestor_depth;
    opt.ranking.apply(config);
    api::BeamSession session(record.sequence, record.scaffold, config);
    const auto best = session.viterbi();
    if (!best.ok) {
      out.error = best.error;
      return out;
    }
    if (!best.reachable) {
      out.error = "the beam disconnected the root; increase the beam size";
      return out;
    }
    out.union_structure = best.union_structure;
    out.mfe = best.energy;
    if (want_bpp) {
      const auto outside = session.outside();
      if (!outside.ok) {
        out.error = outside.error;
        return out;
      }
      out.ensemble_free_energy = outside.ensemble_free_energy;
      out.pairs = outside.pair_probabilities;
    } else {
      const auto inside = session.inside();
      if (!inside.ok) {
        out.error = inside.error;
        return out;
      }
      out.ensemble_free_energy = inside.ensemble_free_energy;
    }
  }
  out.ok = true;
  return out;
}

}  // namespace

int cparty::scfg2::cli::run_fold(int argc, char **argv) {
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
  std::ofstream bpp;
  if (!opt.bpp.empty()) {
    bpp.open(opt.bpp);
    if (!bpp) {
      std::cerr << "error: cannot write " << opt.bpp << "\n";
      return EXIT_FAILURE;
    }
    bpp << std::setprecision(std::numeric_limits<double>::max_digits10);
  }
  int failures = 0;
  for (const auto &record : records) {
    const auto result = predict(record, opt, !opt.bpp.empty());
    if (!result.ok) {
      std::cerr << "error: " << record.name << ": " << result.error << "\n";
      ++failures;
      continue;
    }
    std::cout << ">" << record.name << "\n"
              << record.sequence << "\n"
              << record.scaffold << "\n"
              << std::fixed << std::setprecision(2) << result.union_structure
              << " (" << result.mfe << ")\n"
              << "ensemble free energy: " << std::setprecision(4)
              << result.ensemble_free_energy << " kcal/mol\n"
              << std::defaultfloat;
    if (!opt.bpp.empty()) {
      bpp << ">" << record.name << "\n";
      for (const auto &pair : result.pairs) {
        if (pair.union_probability < opt.bpp_cutoff) continue;
        bpp << pair.i << " " << pair.j << " " << pair.union_probability << "\n";
      }
      bpp << "\n";
    }
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
