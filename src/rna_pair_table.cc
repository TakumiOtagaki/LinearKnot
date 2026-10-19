#include "rna_pair_table.hh"

#include <string.h>

extern "C" {
#include <ViennaRNA/pair_mat.h>
}

namespace cparty::rna {
namespace {

void ensure_pair_matrix_ready_impl() {
  static thread_local bool ready = false;
  if (ready) return;
  make_pair_matrix();
  ready = true;
}

}  // namespace

void ensure_pair_matrix_initialized() { ensure_pair_matrix_ready_impl(); }

pair_type pair_type_for(short si, short sj) {
  ensure_pair_matrix_ready_impl();
  return static_cast<pair_type>(pair[si][sj]);
}

int reverse_pair_type(int pair_type) {
  ensure_pair_matrix_ready_impl();
  return rtype[pair_type];
}

}  // namespace cparty::rna
