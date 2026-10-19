#include "scfg2/api/beam.hh"
#include "scfg2/api/centroid.hh"
#include "scfg2/api/mea.hh"
#include "scfg2/api/profiles.hh"
#include "scfg2/api/sampling.hh"
#include "scfg2/api/subopt.hh"
#include "scfg2/geometry/static_scaffold_geometry.hh"
#include "scfg2/weights/runtime_snapshot.hh"

extern "C" {
#include "ViennaRNA/params/io.h"
}

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unistd.h>
#include <vector>

namespace {
namespace api = cparty::scfg2::api;

void require(bool value, const char *message) {
  if (!value) throw std::runtime_error(message);
}

void close(double actual, double expected) {
  require(std::isfinite(actual) && std::isfinite(expected) &&
          std::abs(actual - expected) < 1e-10, "boundary result mismatch");
}

struct TemporaryParameters {
  std::string path;
  TemporaryParameters() {
    auto pattern = (std::filesystem::temp_directory_path() /
                    "scfg2-parameters-XXXXXX").string();
    std::vector<char> name(pattern.begin(), pattern.end());
    name.push_back('\0');
    const int fd = mkstemp(name.data());
    require(fd >= 0, "cannot create temporary parameter file");
    ::close(fd);
    path = name.data();
  }
  ~TemporaryParameters() { std::filesystem::remove(path); }
};

void default_parameters() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", SCFG2_TEST_PARAM_FILE, 1);
}

void parameter_selection() {
  default_parameters();
  std::string seq = "GGGAAACCC", g = "(((...)))";
  const auto baseline = api::inside(seq, g);
  require(baseline.ok, "baseline inside failed");
  cparty::scfg2::weights::RuntimeSnapshot old_owner(
      seq, g);
  const double old_loop = old_owner.log_pf_parameters().exphairpin[3];

  // Save a complete alternate model, then restore the globals. Only the
  // environment-selected load can supply the additional 10 kcal/mol.
  TemporaryParameters alternate;
  hairpin37[3] += 1000;
  require(vrna_params_save(alternate.path.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT),
          "cannot save alternate parameters");
  require(vrna_params_load(SCFG2_TEST_PARAM_FILE, VRNA_PARAMETER_FORMAT_DEFAULT),
          "cannot restore baseline globals");
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", alternate.path.c_str(), 1);
  const auto changed = api::inside(seq, g);
  require(changed.ok, "alternate inside failed");
  close(changed.ensemble_free_energy, baseline.ensemble_free_energy + 10.);
  close(old_owner.log_pf_parameters().exphairpin[3], old_loop);
  const auto check_partition = [&](const auto &result) {
    require(result.ok, "consumer failed with alternate parameters");
    close(result.log_partition_weight, changed.log_partition_weight);
  };
  check_partition(api::outside(seq, g));
  check_partition(api::beam_inside(seq, g));
  check_partition(api::beam_outside(seq, g));
  check_partition(api::exact_sample(seq, g, 1, 123));
  check_partition(api::beam_sample(seq, g, 1, 123));
  check_partition(api::exact_structural_profiles(seq, g, {1, 3}));
  check_partition(api::beam_structural_profiles(seq, g, {1, 3}));
  check_partition(api::exact_mea(seq, g, 1.));
  check_partition(api::beam_mea(seq, g, 1.));
  check_partition(api::exact_centroid(seq, g, 1.));
  check_partition(api::beam_centroid(seq, g, 1.));
  check_partition(api::exact_subopt(seq, g));
  check_partition(api::beam_subopt(seq, g));

  const auto reject_parameters = [&] {
    const auto check = [](const auto &result) {
      require(!result.ok &&
              result.error.find("CPARTY_FIXED_ENERGY_PARAM_FILE") != std::string::npos,
              "consumer accepted an invalid parameter file");
    };
    check(api::inside(seq, g));
    check(api::outside(seq, g));
    check(api::viterbi(seq, g));
    check(api::beam_inside(seq, g));
    check(api::beam_outside(seq, g));
    check(api::beam_viterbi(seq, g));
    check(api::exact_sample(seq, g, 1, 123));
    check(api::beam_sample(seq, g, 1, 123));
    check(api::exact_structural_profiles(seq, g, {1, 3}));
    check(api::beam_structural_profiles(seq, g, {1, 3}));
    check(api::exact_mea(seq, g, 1.));
    check(api::beam_mea(seq, g, 1.));
    check(api::exact_centroid(seq, g, 1.));
    check(api::beam_centroid(seq, g, 1.));
    check(api::exact_subopt(seq, g));
    check(api::beam_subopt(seq, g));
    api::ExactSession exact(seq, g);
    api::BeamSession beam(seq, g);
    check(exact.evaluate_energy(g));
    check(exact.probability(g));
    check(beam.probability(g));
  };

  // These inputs used to return success without changing any parameters.
  // Rejection must precede both global table and loaded-name mutation.
  const int selected_hairpin = hairpin37[3];
  const std::string selected_name = last_parameter_file();
  TemporaryParameters invalid;
  for (const auto *content : {"x", "\n", "arbitrary parameter text\n",
       "## RNAfold parameter file v1.0\n# Misc\n410 360 50 370\n",
       "## RNAfold parameter file v2.0\n",
       "## RNAfold parameter file v2.0\n#END\n",
       "## RNAfold parameter file v2.0\n# unknown\n1 2 3\n"}) {
    { std::ofstream out(invalid.path); out << content; }
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", invalid.path.c_str(), 1);
    reject_parameters();
    require(hairpin37[3] == selected_hairpin &&
            selected_name == last_parameter_file(),
            "rejected parameter file mutated the active model");
  }

  TemporaryParameters empty;
  for (const auto &path : {alternate.path + ".missing", empty.path, std::string{}}) {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", path.c_str(), 1);
    reject_parameters();
  }
  default_parameters();
  const auto restored = api::inside(seq, g);
  require(restored.ok, "restored inside failed");
  close(restored.log_partition_weight, baseline.log_partition_weight);

  // The content may legitimately equal the active/default model. Valid
  // partial files remain overlays on the current globals, as in ViennaRNA.
  TemporaryParameters partial;
  { std::ofstream out(partial.path);
    out << "## RNAfold parameter file v2.0\n# Misc\n410 360 50 370\n"; }
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", partial.path.c_str(), 1);
  const auto overlay = api::inside(seq, g);
  require(overlay.ok, "valid partial parameter overlay was rejected");
  close(overlay.log_partition_weight, baseline.log_partition_weight);

  TemporaryParameters andronescu;
  require(vrna_params_load_RNA_Andronescu2007() &&
          vrna_params_save(andronescu.path.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT),
          "cannot save Andronescu parameter set");
  require(vrna_params_load(SCFG2_TEST_PARAM_FILE, VRNA_PARAMETER_FORMAT_DEFAULT),
          "cannot restore Turner globals");
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", andronescu.path.c_str(), 1);
  const auto other_model = api::inside(seq, g);
  require(other_model.ok && std::abs(other_model.ensemble_free_energy -
          baseline.ensemble_free_energy) > 0.01,
          "valid Andronescu parameter file was rejected or ignored");
  default_parameters();
}

void sequence_spelling() {
  const std::string g = "(....)";
  const auto reference = api::inside("CUUCGG", g);
  require(reference.ok, "special hairpin baseline failed");
  close(reference.ensemble_free_energy, 3.7);
  for (const auto &seq : {"cuucgg", "CTTCGG", "cTtCgG"}) {
    const auto exact = api::inside(seq, g);
    const auto beam = api::beam_inside(seq, g);
    const auto best = api::viterbi(seq, g);
    require(exact.ok && beam.ok && best.ok, "normalized spelling failed");
    close(exact.log_partition_weight, reference.log_partition_weight);
    close(beam.log_partition_weight, reference.log_partition_weight);
    close(best.energy, reference.ensemble_free_energy);
  }
}

void scaffold_validation() {
  const std::string seq = "GGGAAACCC";
  for (const auto &g : {"[[[...]]]", "{{{...}}}", "<<<...>>>",
                        "([.....)]", "?...?....", "(((...)).", ").......("}) {
    require(!api::inside(seq, g).ok && !api::outside(seq, g).ok &&
            !api::viterbi(seq, g).ok && !api::beam_inside(seq, g).ok &&
            !api::beam_outside(seq, g).ok && !api::beam_viterbi(seq, g).ok &&
            !api::exact_sample(seq, g, 1, 123).ok &&
            !api::beam_sample(seq, g, 1, 123).ok &&
            !api::exact_structural_profiles(seq, g, {1}).ok &&
            !api::beam_structural_profiles(seq, g, {1}).ok &&
            !api::exact_mea(seq, g, 1.).ok && !api::beam_mea(seq, g, 1.).ok &&
            !api::exact_centroid(seq, g, 1.).ok && !api::beam_centroid(seq, g, 1.).ok,
            "consumer accepted an invalid scaffold");
    bool rejected = false;
    try { cparty::scfg2::geometry::StaticScaffoldGeometry geometry(g, 9); }
    catch (const std::invalid_argument &) { rejected = true; }
    require(rejected, "geometry accepted an invalid scaffold");
  }
  require(api::inside(seq, "(((...)))").ok &&
          api::inside(seq, "xx_____XX").ok,
          "valid scaffold symbols were rejected");
  const auto dots = api::viterbi(seq, ".........");
  const auto underscores = api::viterbi(seq, "_________");
  require(dots.ok && underscores.ok && !underscores.generated_pairs.empty(),
          "free underscore scaffold blocked generated pairs");
  close(underscores.energy, dots.energy);
  require(api::beam_viterbi(seq, "_________").ok &&
          api::exact_sample(seq, "_________", 4, 123).ok &&
          api::beam_mea(seq, "_________", 10.).ok,
          "consumer rejected a free scaffold marker during composition");
}

void reversed_pairs() {
  const std::string seq = "GGGAAACCC", g = ".........";
  api::ExactSession exact(seq, g);
  api::BeamSession beam(seq, g);
  for (const auto &forward : {std::vector<api::BasePair>{{1, 9}},
                              std::vector<api::BasePair>{{1, 9}, {2, 8}}}) {
    auto reverse = forward;
    for (auto &pair : reverse) std::swap(pair.i, pair.j);
    const auto a = exact.evaluate_energy(forward), b = exact.evaluate_energy(reverse);
    const auto pa = exact.probability(forward), pb = exact.probability(reverse);
    const auto ba = beam.probability(forward), bb = beam.probability(reverse);
    require(a.ok && b.ok && a.reachable && b.reachable &&
            pa.ok && pb.ok && ba.ok && bb.ok &&
            ba.target_included && bb.target_included,
            "reversed target was lost");
    close(a.energy, b.energy);
    close(pa.log_probability, pb.log_probability);
    close(ba.log_probability, bb.log_probability);
    require(a.union_structure == b.union_structure &&
            bb.target_generated_pairs.size() == forward.size(),
            "reversed target changed output");
    for (const auto &pair : bb.target_generated_pairs)
      require(pair.i < pair.j, "public target output was not normalized");
  }
  for (const auto &bad : {std::vector<api::BasePair>{{1, 1}},
                          std::vector<api::BasePair>{{0, 9}},
                          std::vector<api::BasePair>{{10, 1}},
                          std::vector<api::BasePair>{{1, 9}, {9, 1}}}) {
    require(!exact.evaluate_energy(bad).ok && !beam.probability(bad).ok,
            "invalid pair list passed normalization");
  }
}
}  // namespace

int main() {
  try {
    parameter_selection();
    sequence_spelling();
    scaffold_validation();
    reversed_pairs();
    std::cout << "scfg2_api_input_boundary=ok\n";
    return EXIT_SUCCESS;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return EXIT_FAILURE;
  }
}
