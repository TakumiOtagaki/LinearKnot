#include "scfg2/cli/commands.hh"
// linearknot sample: structures G u G' drawn from the pseudoknot ensemble of G'
// given a scaffold G (LinearKnot), on the retained beam DAG or exactly.
//
// Input: records of an optional ">name" line, a sequence line and a scaffold
// dot-bracket line.  For each record it prints ">name", the sequence, then
// one sampled union structure per line ('[]' marks G' pairs that cross G).
// --unique prints "count<TAB>structure" sorted by count instead.

#include "scfg2/api/sampling.hh"
#include "scfg2/cli/cli_common.hh"
#include "scfg2/cli/ranking_options.hh"

#include <algorithm>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {

using cparty::scfg2::cli::Record;

struct Options {
  std::string input;
  std::string output;
  std::string sequence;
  std::string scaffold;
  std::size_t count = 10;
  std::uint64_t seed = 1;
  bool unique = false;
  std::size_t beam = 100;
  cparty::scfg2::cli::RankingOptions ranking;
  bool exact = false;
  std::size_t max_hairpin_unpaired = 30;
  std::size_t max_case3_ancestor_depth = 30;
};

void usage(std::ostream &os) {
  os << "Usage: linearknot sample [options] (-i FILE | --sequence SEQ --scaffold G)\n"
        "  -i, --input FILE         records: [>name] / sequence / scaffold G\n"
        "  -o, --output FILE        output file (default: stdout)\n"
        "  -k, --count N            samples per record (default 10)\n"
        "  --seed S                 random seed (default 1)\n"
        "  --unique                 print count<TAB>structure, most frequent first\n"
        "  -b, --beam N             beam size (default 100)\n"
        "  --exact                  exact cubic-time ensemble (no beam)\n"
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
    } else if (arg == "-k" || arg == "--count") {
      if (!(v = value())) return false;
      opt->count = std::strtoul(v, nullptr, 10);
    } else if (arg == "--seed") {
      if (!(v = value())) return false;
      opt->seed = std::strtoull(v, nullptr, 10);
    } else if (arg == "--unique") {
      opt->unique = true;
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

}  // namespace

int cparty::scfg2::cli::run_sample(int argc, char **argv) {
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
  int failures = 0;
  for (std::size_t r = 0; r < records.size(); ++r) {
    const auto &record = records[r];
    const std::uint64_t seed = opt.seed + r;
    bool header_written = false;
    const auto header = [&] {
      if (!header_written) {
        os << ">" << record.name << "\n" << record.sequence << "\n";
        header_written = true;
      }
    };
    std::map<std::string, std::size_t> counts;
    const auto consume = [&](cparty::scfg2::api::SampledStructure &&structure) {
      if (opt.unique) ++counts[std::move(structure.union_structure)];
      else {
        header();
        os << structure.union_structure << "\n";
        if (!os) throw std::runtime_error("cannot write sampled structure");
      }
      return true;
    };
    cparty::scfg2::api::SampleStreamResult result;
    if (opt.exact) {
      result = cparty::scfg2::api::exact_sample_each(record.sequence, record.scaffold,
                                                   opt.count, seed, consume);
    } else {
      cparty::scfg2::api::BeamConfig config;
      config.beam_size = opt.beam;
      config.max_hairpin_unpaired = opt.max_hairpin_unpaired;
      config.max_case3_ancestor_depth = opt.max_case3_ancestor_depth;
      opt.ranking.apply(config);
      result = cparty::scfg2::api::beam_sample_each(record.sequence, record.scaffold,
                                                  opt.count, seed, consume, config);
    }
    if (!result.ok) {
      std::cerr << "error: " << record.name << ": " << result.error << "\n";
      ++failures;
      continue;
    }
    header();  // count=0 still emits the record header.
    if (opt.unique) {
      std::vector<std::pair<std::string, std::size_t>> sorted(counts.begin(), counts.end());
      std::stable_sort(sorted.begin(), sorted.end(), [](const auto &a, const auto &b) {
        return a.second > b.second;
      });
      for (const auto &[structure, count] : sorted)
        os << count << "\t" << structure << "\n";
    }
    os << "\n";
  }
  return failures == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
