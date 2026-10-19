#include "h_globals.hh"
#include "part_func.hh"
#include "sparse_tree.hh"

extern "C" {
#include "ViennaRNA/params/io.h"
}

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>

// A deliberately small observation surface for a frozen legacy checkout.
//
// The public CParty CLI runs legacy MFE before legacy inside.  Some constraint
// sentinels (notably `x`) can therefore crash before hfold_pf is reached.  This
// probe calls only the legacy inside implementation so that such inputs remain
// usable as semantic-oracle witnesses.  Compile this file against the object
// files from the frozen legacy checkout; do not link it to SCFG2.
int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "usage: legacy_inside_oracle_probe <sequence> <scaffold> "
                 "<parameter-file>\n";
    return 2;
  }

  std::string sequence = argv[1];
  std::string scaffold = argv[2];
  if (sequence.size() != scaffold.size()) {
    std::cerr << "sequence and scaffold lengths differ\n";
    return 2;
  }
  if (vrna_params_load(argv[3], VRNA_PARAMETER_FORMAT_DEFAULT) == 0) {
    std::cerr << "failed to load parameter file: " << argv[3] << '\n';
    return 2;
  }

  sparse_tree tree(scaffold, static_cast<cand_pos_t>(sequence.size()));
  W_final_pf inside(sequence,
                    scaffold,
                    false,
                    false,
                    false,
                    2,
                    0.0,
                    0,
                    false);
  const pf_t ensemble_free_energy = inside.hfold_pf(tree);
  std::cout << std::setprecision(17) << ensemble_free_energy << '\n';
  return EXIT_SUCCESS;
}
