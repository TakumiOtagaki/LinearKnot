// Independent diagnostic witnesses around the shared bounded guide.
#ifndef SCFG2_BOUNDED_OUTER_PROBE_HH
#define SCFG2_BOUNDED_OUTER_PROBE_HH

#include "scfg2_outer_probe_common.hh"
#include "scfg2/beam/bounded_outer.hh"

namespace bounded_outer_detail {
using namespace sc::beam::outer_detail;
// Dense, nonmonotone priorities and ties; compare every shrinking prefix
// against an independently sorted full history, including widths 1 and 20.
inline void check_prefix_top() {
  for (std::size_t width : {1, 3, 20}) {
    Work work;
    std::vector<Record> records;
    for (int end = 2; end <= 97; ++end)
      records.push_back({{NT::WM, 1, end, -1, -1}, double((end * 37) % 23)});
    PrefixTop index(records, width, work);
    for (int upper = 100; upper >= 0; --upper) {
      const auto &got = index.query(upper);
      std::vector<std::size_t> want;
      for (std::size_t p = 0; p < records.size(); ++p)
        if (records[p].key.j <= upper) want.push_back(p);
      std::sort(want.begin(), want.end(), [&](auto a, auto b) {
        if (records[a].alpha != records[b].alpha) return records[a].alpha > records[b].alpha;
        return records[a].key.j < records[b].key.j;
      });
      if (want.size() > width) want.resize(width);
      if (got != want) throw std::runtime_error("C++ prefix/full-history mismatch");
    }
    bool rejected = false;
    try { index.query(1); } catch (const std::runtime_error &) { rejected = true; }
    if (!rejected || work.comparisons > records.size() * width ||
        work.build_moves > records.size() * (width + 1) ||
        work.undo_moves > records.size() * (width + 1) || work.cursor_steps > records.size())
      throw std::runtime_error("C++ prefix query/work contract violated");
  }
}

} // namespace bounded_outer_detail
using sc::beam::bounded_outer;
#endif
