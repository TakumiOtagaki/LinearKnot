#ifndef CPARTY_SCFG2_BEAM_PREFIX_AGENDA_STORAGE_HH
#define CPARTY_SCFG2_BEAM_PREFIX_AGENDA_STORAGE_HH

#include <cstddef>

namespace cparty::scfg2::beam {

// Observation only. Computing a snapshot may scan all agenda buckets; the
// production executor must not request one in its per-column hot path.
struct PrefixAgendaStorage {
  std::size_t pending_records = 0;
  std::size_t largest_pending_bucket = 0;
  std::size_t dedup_keys = 0;
  std::size_t active_records = 0;
  std::size_t saved_frame_records = 0;
  std::size_t saved_frames = 0;
};

}  // namespace cparty::scfg2::beam
#endif
