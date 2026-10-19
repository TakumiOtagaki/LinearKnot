// Verify the static index against direct interval definitions, without a
// second tree implementation or privileged access to production internals.
#include "scfg2/geometry/static_scaffold_geometry.hh"

#include <cstdlib>
#include <iostream>
#include <random>
#include <string>
#include <vector>

namespace {

struct IntervalSpecification {
  std::vector<int> partner, parent, up;
  explicit IntervalSpecification(const std::string &g)
      : partner(g.size() + 1, -2), parent(g.size() + 1), up(g.size() + 1) {
    std::vector<int> openings;
    for (std::size_t pos = 1; pos <= g.size(); ++pos) {
      if (g[pos - 1] == 'x' || g[pos - 1] == 'X') partner[pos] = -1;
      if (g[pos - 1] == '(') openings.push_back(static_cast<int>(pos));
      if (g[pos - 1] == ')') {
        const int a = openings.back();
        openings.pop_back();
        partner[a] = static_cast<int>(pos);
        partner[pos] = a;
      }
    }
    for (std::size_t pos = 1; pos <= g.size(); ++pos) {
      // The immediate enclosing arc excludes its own two endpoints.
      for (int a = 1; a < static_cast<int>(pos); ++a)
        if (partner[a] > static_cast<int>(pos)) parent[pos] = a;
      up[pos] = partner[pos] > 0 || pos == 1 ? 0 : up[pos - 1] + 1;
    }
  }
  int bp(int i, int l) const {
    return partner[l] > 0 || parent[l] == 0 ? -2 : parent[l] < i ? -1 : parent[l];
  }
  int Bp(int l, int j) const {
    return partner[l] > 0 || parent[l] == 0 ? -2 : partner[parent[l]] > j ? -1 : partner[parent[l]];
  }
  int b(int i, int l) const {
    if (bp(i, l) < 0) return bp(i, l);
    // Outermost arc enclosing l whose opening lies in [i,l).
    for (int a = i; a < l; ++a) if (partner[a] > l) return a;
    throw std::logic_error("inconsistent interval specification");
  }
  int B(int l, int j) const {
    if (Bp(l, j) < 0) return Bp(l, j);
    // Outermost arc enclosing l whose closing lies in (l,j].
    for (int z = j; z > l; --z) if (partner[z] > 0 && partner[z] < l) return z;
    throw std::logic_error("inconsistent interval specification");
  }
  bool weakly_closed(int i, int j) const {
    for (int pos = i; pos <= j; ++pos)
      if (partner[pos] > 0 && (partner[pos] < i || partner[pos] > j)) return false;
    return true;
  }
};

void fill(std::mt19937 &rng, std::string &g, int lo, int hi, int depth) {
  std::uniform_int_distribution<int> coin(0, 99);
  int pos = lo;
  while (pos + 5 <= hi) {
    if (coin(rng) < 35) {
      ++pos;
      continue;
    }
    const int span = std::uniform_int_distribution<int>(5, hi - pos + 1)(rng);
    const int a = pos, b = pos + span - 1;
    g[a] = '(';
    g[b] = ')';
    if (depth < 5) fill(rng, g, a + 1, b - 1, depth + 1);
    pos = b + 1;
  }
}

std::size_t checked = 0, scaffolds = 0;

void check_forbidden_case_equivalence(const std::string &lowercase) {
  std::string uppercase = lowercase;
  bool has_forbidden = false;
  for (char &ch : uppercase) {
    if (ch == 'x') {
      ch = 'X';
      has_forbidden = true;
    }
  }
  if (!has_forbidden) return;
  const int n = static_cast<int>(lowercase.size());
  cparty::scfg2::geometry::StaticScaffoldGeometry lower(lowercase, n);
  cparty::scfg2::geometry::StaticScaffoldGeometry upper(uppercase, n);
  if (!lower.usable() || !upper.usable()) std::exit(EXIT_FAILURE);
  // This contract is independent of the old in-tree sparse_tree parser:
  // both spellings explicitly mean a pair-forbidden base, not a free dot.
  for (int x = 1; x <= n; ++x) {
    if ((uppercase[x - 1] == 'X' && upper.partner(x) != -1) ||
        upper.partner(x) != lower.partner(x) ||
        upper.parent(x) != lower.parent(x) ||
        upper.unpaired_prefix(x) != lower.unpaired_prefix(x)) {
      std::cerr << "pair-forbidden case mismatch G=" << uppercase
                << " pos=" << x << '\n';
      std::exit(EXIT_FAILURE);
    }
    for (int y = x; y <= n; ++y) {
      if (upper.B(x, y) != lower.B(x, y) ||
          upper.b(x, y) != lower.b(x, y) ||
          upper.bp(x, y) != lower.bp(x, y) ||
          upper.Bp(x, y) != lower.Bp(x, y) ||
          upper.weakly_closed(x, y) != lower.weakly_closed(x, y)) {
        std::cerr << "pair-forbidden geometry mismatch G=" << uppercase
                  << " args=" << x << ',' << y << '\n';
        std::exit(EXIT_FAILURE);
      }
      ++checked;
    }
  }
}

void check_small(const std::string &g) {
  const int n = static_cast<int>(g.size());
  const IntervalSpecification spec(g);
  cparty::scfg2::geometry::StaticScaffoldGeometry index(g, n);
  if (!index.usable()) std::exit(EXIT_FAILURE);
  ++scaffolds;
  for (int x = 1; x <= n; ++x) {
    if (index.partner(x) != spec.partner[x] ||
        index.parent(x) != spec.parent[x] || index.unpaired_prefix(x) != spec.up[x]) {
      std::cerr << "parent/partner/up mismatch G=" << g << " pos=" << x << '\n';
      std::exit(EXIT_FAILURE);
    }
    for (int y = x; y <= n; ++y) {
      if (index.B(x, y) != spec.B(x, y) || index.b(x, y) != spec.b(x, y) ||
          index.bp(x, y) != spec.bp(x, y) || index.Bp(x, y) != spec.Bp(x, y) ||
          index.weakly_closed(x, y) != spec.weakly_closed(x, y)) {
        std::cerr << "static geometry mismatch G=" << g << " args=(" << x << ',' << y
                  << ") B=" << index.B(x, y) << '/' << spec.B(x, y)
                  << " b=" << index.b(x, y) << '/' << spec.b(x, y) << '\n';
        std::exit(EXIT_FAILURE);
      }
      ++checked;
    }
  }
  if (!index.weakly_closed(1, 0)) std::exit(EXIT_FAILURE);
  check_forbidden_case_equivalence(g);
}

void enumerate(std::string &g, int balance, int remaining) {
  if (!remaining) {
    if (!balance) {
      check_small(g);
      std::string forced = g;
      for (std::size_t x = 0; x < forced.size(); ++x)
        if (forced[x] == '.' && x % 2 == 0) forced[x] = 'x';
      if (forced != g) check_small(forced);
    }
    return;
  }
  for (char ch : {'.', '(', ')'}) {
    const int next = balance + (ch == '(' ? 1 : ch == ')' ? -1 : 0);
    if (next < 0 || next > remaining - 1) continue;
    g.push_back(ch);
    enumerate(g, next, remaining - 1);
    g.pop_back();
  }
}

// Independent ancestor-walk oracle for inputs beyond the old tree's fixed
// table/uint16_t limit. It deliberately uses linear navigation and child scans.
void check_long(const std::string &g, std::mt19937 &rng, const char *shape) {
  const int n = static_cast<int>(g.size());
  cparty::scfg2::geometry::StaticScaffoldGeometry index(g, n);
  std::vector<int> partner(n + 1, -2), parent(n + 1, -1), depth(n + 1), stack{0};
  std::vector<std::vector<int>> children(n + 1);
  for (int pos = 1; pos <= n; ++pos) {
    if (g[pos - 1] == ')') {
      partner[pos] = stack.back();
      partner[stack.back()] = pos;
      stack.pop_back();
    }
    parent[pos] = stack.back();
    depth[pos] = depth[parent[pos]] + 1;
    if (g[pos - 1] == '(') {
      children[parent[pos]].push_back(pos);
      stack.push_back(pos);
    }
  }
  auto ancestor = [&](int a, int b) {
    while (a != b) {
      if (depth[a] >= depth[b]) a = parent[a];
      else b = parent[b];
    }
    return a;
  };
  for (int q = 0; q < 2000; ++q) {
    int i = 1 + rng() % n, j = 1 + rng() % n;
    if (i > j) std::swap(i, j);
    int b = -2, B = -2;
    if (parent[j] && partner[j] < 0) {
      if (parent[j] < i) b = -1;
      else {
        const int a = ancestor(i, j);
        b = i == a ? i : -100;
        if (i != a) for (int child : children[a]) {
          if (child < j && partner[child] > j) { b = child; break; }
        }
      }
    }
    if (parent[i] && partner[i] < 0) {
      if (partner[parent[i]] > j) B = -1;
      else {
        const int a = ancestor(i, j);
        B = j == a ? j : -100;
        if (j != a) for (int child : children[a]) {
          if (child < j && partner[child] > i) { B = partner[child]; break; }
        }
      }
    }
    if (index.b(i, j) != b || index.B(i, j) != B ||
        index.lca(i, j) != ancestor(i, j)) {
      std::cerr << "long geometry mismatch shape=" << shape << " args=" << i << ',' << j << '\n';
      std::exit(EXIT_FAILURE);
    }
    ++checked;
  }
  if (!index.usable() || index.storage_bytes() > 256 * static_cast<std::size_t>(n + 1))
    std::exit(EXIT_FAILURE);
  std::cout << "n=" << n << " shape=" << shape << " storage_bytes=" << index.storage_bytes() << '\n';
}

}  // namespace

int main(int argc, char **argv) {
  const int cases = argc > 1 ? std::atoi(argv[1]) : 100;
  std::mt19937 rng(20261010);
  for (int n = 1; n <= 9; ++n) {
    std::string g;
    enumerate(g, 0, n);
  }
  for (int c = 0; c < cases; ++c) {
    const int n = std::uniform_int_distribution<int>(8, 150)(rng);
    std::string g(n, '.');
    fill(rng, g, 0, n - 1, 0);
    check_small(g);
  }
  for (int n : {4096, 16384, 65536, 262144}) {
    check_long(std::string(n, '.'), rng, "unpaired");
    std::string deep(n / 2, '(');
    deep += std::string(n % 2, '.') + std::string(n / 2, ')');
    // Include unpaired bases in the deep case to exercise borders as well as LCA.
    for (int x = 1; x + 1 < n / 2; x += 5) deep[x] = deep[n - x - 1] = '.';
    check_long(deep, rng, "nested");
    std::string wide(n, '.');
    for (int x = 0; x + 2 < n; x += 3) { wide[x] = '('; wide[x + 2] = ')'; }
    check_long(wide, rng, "wide");
  }
  std::cout << "scaffolds=" << scaffolds << " checked=" << checked << "\nscfg2_border_query=ok\n";
  return EXIT_SUCCESS;
}
