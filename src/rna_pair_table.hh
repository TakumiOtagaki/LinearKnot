#ifndef CPARTY_RNA_PAIR_TABLE_HH
#define CPARTY_RNA_PAIR_TABLE_HH

#include "scfg2/types.hh"

namespace cparty::rna {

void ensure_pair_matrix_initialized();
pair_type pair_type_for(short si, short sj);
int reverse_pair_type(int pair_type);

}  // namespace cparty::rna

#endif  // CPARTY_RNA_PAIR_TABLE_HH
