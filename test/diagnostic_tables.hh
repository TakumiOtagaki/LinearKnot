#ifndef CPARTY_TEST_DIAGNOSTIC_TABLES_HH
#define CPARTY_TEST_DIAGNOSTIC_TABLES_HH

#include "scfg2/types.hh"

#include <algorithm>
#include <array>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cparty::scfg2::testing {

// Test-only materialized projections. These are not the executor's chart;
// production snapshots and contexts neither own nor read this storage.
enum class DiagnosticTable : std::size_t {
  V, VM, WM, WMv, WMp, WI, WIP, WIP_VP_L, WIP_VP_R, WIP_BE_L, WIP_BE_R,
  VP, VP_CASE3_BOUNDARY, VP_CASE3_R, VP_CASE3_MR, VP_CASE3_LEFT,
  VP_C12_R2, VP_C12_2, VP_INTERNAL_LOCAL, VPL, VPR, WMB, WMBP,
  WMB_BTAIL, WMB_BSTEP, WMB_NONEMPTY, WMB_EMPTY_LOCAL, WMBW, BE, Count
};

class DiagnosticTables {
 public:
  explicit DiagnosticTables(cand_pos_t n = 0, bool materialize = false)
      : offsets_(static_cast<std::size_t>(n) + 1),
        prefix_(static_cast<std::size_t>(n) + 1) {
    if (!materialize) return;
    const auto extent = static_cast<std::size_t>(n);
    if (extent * (extent + 1) / 2 > static_cast<std::size_t>(std::numeric_limits<cand_pos_t>::max()))
      throw std::length_error("diagnostic span offsets exceed the coordinate domain");
    for (cand_pos_t i = 1; i <= n; ++i) {
      const std::size_t row = static_cast<std::size_t>(i - 1);
      offsets_[i] = static_cast<cand_pos_t>(row * (2 * static_cast<std::size_t>(n) - i + 2) / 2);
    }
    const std::size_t slots = (static_cast<std::size_t>(n) + 1) * (static_cast<std::size_t>(n) + 2) / 2;
    for (auto &table : spans_) table.resize(slots);
  }

  std::vector<double> &table(DiagnosticTable id) { return spans_[static_cast<std::size_t>(id)]; }
  const std::vector<double> &table(DiagnosticTable id) const { return spans_[static_cast<std::size_t>(id)]; }
  std::vector<double> &prefix() { return prefix_; }

  double get(DiagnosticTable id, cand_pos_t i, cand_pos_t j) const {
    if (i > j) return id == DiagnosticTable::WI ? 1. : 0.;
    if (i == j && id != DiagnosticTable::WI &&
        id != DiagnosticTable::VP_CASE3_R && id != DiagnosticTable::VP_CASE3_MR) return 0.;
    return table(id).at(static_cast<std::size_t>(offsets_.at(i)) + j - i);
  }
  void set(DiagnosticTable id, cand_pos_t i, cand_pos_t j, double value) {
    if (i >= j && !(i == j && (id == DiagnosticTable::VP_CASE3_R || id == DiagnosticTable::VP_CASE3_MR))) return;
    table(id).at(static_cast<std::size_t>(offsets_.at(i)) + j - i) = value;
  }
  void clear() {
    for (auto &values : spans_) std::fill(values.begin(), values.end(), 0.);
    std::fill(prefix_.begin(), prefix_.end(), 0.);
  }

 private:
  std::array<std::vector<double>, static_cast<std::size_t>(DiagnosticTable::Count)> spans_;
  std::vector<cand_pos_t> offsets_;
  std::vector<double> prefix_;
};

}  // namespace cparty::scfg2::testing
#endif
