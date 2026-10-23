#ifndef CPARTY_SCFG2_API_PROFILES_HH
#define CPARTY_SCFG2_API_PROFILES_HH

#include "scfg2/api/beam.hh"

#include <array>
#include <limits>
#include <string>
#include <vector>

namespace cparty::scfg2::api {

// Structural contexts of one base in G u G'. PseudoLoop holds unpaired
// bases directly in gaps between crossing bands. Multiloop wings within a
// band are Multi, independently of their pseudoknot energy parameters.
inline constexpr std::size_t structural_context_count = 7;
inline constexpr std::array<const char *, structural_context_count>
    structural_context_names = {"Stem",  "Hairpin",  "Bulge",     "Internal",
                                "Multi", "Exterior", "PseudoLoop"};
inline constexpr std::size_t accessibility_loop_context_count = structural_context_count - 1;
inline constexpr std::array<const char *, accessibility_loop_context_count>
    accessibility_loop_context_names = {"Hairpin", "Bulge", "Internal",
                                        "Multi", "Exterior", "PseudoLoop"};

// PK-CapR profile and PK-Raccess accessibility for the conditional ensemble
// of G' given the scaffold G.  The beam variant normalises by the beam
// partition function, as LinearCapR and LinearRaccess do.
struct StructuralProfileResult {
  bool ok = false;
  WmbClosurePolicy wmb_closure_policy = WmbClosurePolicy::RequireWeakClosure;
  std::string error;
  std::size_t beam_size = 0;  // zero for the exact computation
  double log_partition_weight = -std::numeric_limits<double>::infinity();
  double ensemble_free_energy = std::numeric_limits<double>::quiet_NaN();
  // RT in kcal/mol of the energy model, for ED = -RT log P.
  double rt_kcal_mol = std::numeric_limits<double>::quiet_NaN();
  // contexts[pos][c] with one-based pos (index zero unused); c indexes
  // structural_context_names.
  std::vector<std::array<double, structural_context_count>> contexts;
  // accessibility[w_idx][s] = P(bases s..s+w-1 all unpaired), one-based s,
  // for w = access_lengths[w_idx]; entries past n-w+1 are zero.
  std::vector<int> access_lengths;
  std::vector<std::vector<double>> accessibility;
  // Populated only when by_loop=true. [window][one-based start][context],
  // with contexts ordered by accessibility_loop_context_names. Each entry is
  // P(the entire window is unpaired in that loop), not a conditional ratio
  // or a product of base marginals. Summing contexts recovers accessibility.
  std::vector<std::vector<std::array<double, accessibility_loop_context_count>>>
      accessibility_by_loop;
  double elapsed_seconds = 0.0;
};

StructuralProfileResult exact_structural_profiles(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::vector<int> &access_lengths,
    WmbClosurePolicy policy = WmbClosurePolicy::RequireWeakClosure,
    bool by_loop = false);

StructuralProfileResult beam_structural_profiles(
    const std::string &sequence,
    const std::string &scaffold_structure,
    const std::vector<int> &access_lengths,
    BeamConfig config = {},
    bool by_loop = false);

}  // namespace cparty::scfg2::api

#endif  // CPARTY_SCFG2_API_PROFILES_HH
