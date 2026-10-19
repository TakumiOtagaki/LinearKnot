#ifndef CPARTY_SCFG2_ENGINE_TRACEBACK_RENDER_HH
#define CPARTY_SCFG2_ENGINE_TRACEBACK_RENDER_HH

#include "scfg2/engine/traceback_pairs.hh"

#include <string>

namespace cparty::scfg2::engine {

struct TracebackLaneStructures { std::string round; std::string square; };

inline TracebackLaneStructures render_traceback_pairs(
    cand_pos_t n,
    const std::vector<TracebackPair> &pairs) {
  TracebackLaneStructures out{std::string(static_cast<std::size_t>(n), '.'),
                              std::string(static_cast<std::size_t>(n), '.')};
  for (const auto &pair : pairs) {
    auto &structure = pair.lane == TracebackPairLane::Round ? out.round : out.square;
    const char open = pair.lane == TracebackPairLane::Round ? '(' : '[';
    const char close = pair.lane == TracebackPairLane::Round ? ')' : ']';
    if (pair.i > 0 && pair.j <= n && structure[pair.i - 1] == '.' && structure[pair.j - 1] == '.') {
      structure[pair.i - 1] = open;
      structure[pair.j - 1] = close;
    }
  }
  return out;
}

}  // namespace cparty::scfg2::engine

#endif  // CPARTY_SCFG2_ENGINE_TRACEBACK_RENDER_HH
