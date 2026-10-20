#include "scfg2/beam/bounded_outer.hh"

#include <cstring>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
namespace sc = cparty::scfg2;
namespace outer = sc::beam::outer_detail;
using Key = sc::deductions::ItemKey;
using NT = sc::grammar::NonTerminal;

void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(a)) == 0; }

// Point lookup needs only the static scaffold boundary of the context. Keep
// this witness independent of thermodynamic setup and inside computation.
struct Scaffold {
  std::vector<int> partners;
  explicit Scaffold(int n) : partners(n + 1, -2) {
    for (const auto [i, j] : {std::pair{2, 11}, std::pair{3, 10}}) {
      partners[i] = j;
      partners[j] = i;
    }
  }
  int pair_partner(int p) const { return partners.at(p); }
  int pair_square(int) const { return -2; }
};

std::size_t check(std::size_t width) {
  constexpr int n = 12;
  outer::Chart pilot(outer::zero);
  for (std::size_t family = 0; family < outer::outer_probe_families; ++family) {
    const auto nt = static_cast<NT>(family);
    for (int j = 1; j <= n; ++j) {
      for (int i = 1; i <= j; ++i) {
        if (nt == NT::BE) {
          if (j - i < 3) continue;
          for (int ip = i + 1; ip < j - 1; ++ip)
            pilot.set({nt, i, j, ip, j - 1}, double((i + 3 * ip + j) % 7) - 8.25);
        } else {
          pilot.set({nt, i, j, -1, -1}, double((family + i + 3 * j) % 7) - 4.25);
        }
      }
    }
  }
  const Key tiny{NT::V, 1, n, -1, -1};
  const Key negative_zero{NT::WI, 1, n, -1, -1};
  const Key explicit_zero{NT::WM, 4, 8, -1, -1};
  const std::vector<Key> noncanonical{
      {NT::V, 1, n, 2, 7}, {NT::V, 1, n, 2, 8},
      {NT::V, 1, n, -1, 8}, {NT::WI, 2, n, 0, 0}};
  pilot.set(tiny, -1e300);
  pilot.set(negative_zero, -0.0);
  pilot.set(explicit_zero, outer::zero);
  for (std::size_t k = 0; k < noncanonical.size(); ++k)
    pilot.set(noncanonical[k], double(k) + 0.125);

  Scaffold ctx(n);
  outer::Work work;
  const outer::Proposals proposals(pilot, n, width, ctx, work);
  std::size_t queries = 0, missing = 0;
  const auto compare = [&](const Key &key) {
    const auto expected = proposals.alpha.get(key);
    const auto actual = proposals.get_alpha(key);
    require(same(expected, actual), "sorted alpha point lookup differs bitwise from chart");
    ++queries;
    if (expected == outer::zero) ++missing;
  };

  // Every retained record, every rejected pilot record and neighboring full
  // keys are checked against the independent unordered-map representation.
  for (const auto &[key, value] : pilot.values()) {
    (void)value;
    compare(key);
    for (int delta : {-1, 1}) {
      auto neighbor = key; neighbor.i += delta; compare(neighbor);
      neighbor = key; neighbor.j += delta; compare(neighbor);
      neighbor = key; neighbor.ip += delta; compare(neighbor);
      neighbor = key; neighbor.jp += delta; compare(neighbor);
    }
  }
  for (const auto &[key, value] : proposals.alpha.values()) {
    compare(key);
    require(same(value, proposals.get_alpha(key)), "retained alpha value bits changed");
  }
  for (std::size_t family = 0; family <= outer::outer_probe_families; ++family)
    for (int j = -2; j <= n + 2; ++j)
      for (int i = -2; i <= n + 2; ++i)
        compare({static_cast<NT>(family), i, j, -1, -1});
  // The historical hash adds small signed offsets to coordinates. Stay
  // below its overflow boundary while still probing far outside the guide.
  constexpr int far = std::numeric_limits<int>::max() - 32;
  for (int j : {-1, 0, 1, n, n + 1, far})
    compare({static_cast<NT>(255), 1, j, -1, -1});
  compare({NT::V, std::numeric_limits<int>::min(), n, -1, -1});
  compare({NT::BE, 2, 11, far, 10});
  compare({NT::BE, 2, 11, 3, far});
  require(proposals.get_alpha(explicit_zero) == outer::zero,
          "explicit zero acquired a guide record");
  require(missing > 0 && queries > 0, "missing alpha queries were not exercised");
  if (width >= 22) {
    require(same(proposals.get_alpha(tiny), -1e300), "finite tiny alpha was discarded");
    require(same(proposals.get_alpha(negative_zero), -0.0), "negative zero alpha bits changed");
    for (const auto &key : noncanonical)
      require(proposals.alpha.contains(key), "noncanonical witness was not retained");
  }
  return queries;
}
}  // namespace

int main() {
  try {
    std::size_t queries = 0;
    for (std::size_t width : {1, 2, 22}) queries += check(width);
    std::cout << "scfg2_bounded_outer_alpha_lookup=ok widths=3 queries=" << queries << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
