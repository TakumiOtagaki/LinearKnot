#include "scfg2/engine/primary_lookup_index.hh"

#include <algorithm>
#include <array>
#include <cstring>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace {
namespace sc = cparty::scfg2;
using Index = sc::engine::PrimaryLookupIndex;
using Key = sc::deductions::ItemKey;
using NT = sc::grammar::NonTerminal;
using Identity = std::tuple<int, cand_pos_t, cand_pos_t, cand_pos_t, cand_pos_t>;
using Reference = std::map<Identity, double>;
constexpr std::size_t families = static_cast<std::size_t>(NT::BE);
std::size_t comparisons = 0;

void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(a)) == 0; }
Identity identity(const Key &key) {
  return {static_cast<int>(key.nonterminal), key.i, key.j, key.ip, key.jp};
}
bool supported(cand_pos_t n, const Key &key) {
  return key.nonterminal < NT::BE && key.i >= 1 && key.i <= key.j &&
         key.j <= n && key.ip == -1 && key.jp == -1;
}
void compare(const Index &index, const Reference &reference, const Key &key,
             double default_value = -123.25) {
  const auto actual = index.get(key, default_value);
  if (!supported(index.length(), key)) {
    require(!actual, "unsupported primary lookup did not request fallback");
  } else {
    require(actual.has_value(), "supported primary lookup requested fallback");
    const auto found = reference.find(identity(key));
    const auto expected = found == reference.end() ? default_value : found->second;
    require(same(*actual, expected), "primary lookup value bits differ from reference map");
  }
  ++comparisons;
}
void set(Index &index, Reference &reference, const Key &key, double value) {
  index.set(key, value);
  if (supported(index.length(), key)) reference[identity(key)] = value;
}
void erase(Index &index, Reference &reference, const Key &key) {
  index.erase(key);
  reference.erase(identity(key));
}
void check_counts(const Index &index, const Reference &reference) {
  require(index.size() == reference.size(), "primary lookup entry count differs");
  std::size_t capacity = 0, active = 0;
  for (std::size_t f = 0; f < families; ++f)
    for (cand_pos_t j = 1; j <= index.length(); ++j) {
      const auto row_capacity = index.row_capacity({static_cast<NT>(f), 1, j, -1, -1});
      capacity += row_capacity;
      if (row_capacity) ++active;
    }
  require(index.row_count() == families * static_cast<std::size_t>(index.length()),
          "primary lookup row-header count differs");
  require(index.active_rows() == active && index.allocated_capacity() == capacity,
          "primary lookup capacity counters are stale");
  require(capacity <= 8 * index.size(), "primary lookup retained candidate-peak capacity");
}

void key_and_value_contract() {
  Index empty(0);
  Reference none;
  const Key first{NT::W, 1, 1, -1, -1};
  compare(empty, none, first);
  empty.set(first, 3.0);
  empty.erase(first);
  require(empty.size() == 0 && empty.row_count() == 0 && empty.memory_bytes() == 0,
          "zero-length primary lookup allocated storage");
  bool rejected = false;
  try { Index invalid(-1); } catch (const std::invalid_argument &) { rejected = true; }
  require(rejected, "negative primary lookup length was accepted");

  Index index(32);
  Reference reference;
  const auto header_bytes = index.memory_bytes();
  require(index.allocated_capacity() == 0 && index.active_rows() == 0,
          "primary lookup allocated slot buffers eagerly");
  const std::array<double, 7> values{
      -std::numeric_limits<double>::infinity(), -0.0, 0.0, -1e300,
      std::numeric_limits<double>::infinity(), std::numeric_limits<double>::denorm_min(), 1e300};
  for (std::size_t f = 0; f < families; ++f) {
    for (std::size_t v = 0; v < values.size(); ++v) {
      const Key key{static_cast<NT>(f), static_cast<cand_pos_t>(v + 1), 32, -1, -1};
      set(index, reference, key, values[v]);
      compare(index, reference, key);
    }
  }
  for (double default_value : values) {
    compare(index, reference, {NT::V, 8, 32, -1, -1}, default_value);
    compare(index, reference, {NT::V, 1, 31, -1, -1}, default_value);
  }
  const auto before_size = index.size(), before_capacity = index.allocated_capacity();
  const auto before_bytes = index.memory_bytes();
  const std::vector<Key> unsupported{
      {NT::BE, 1, 32, 2, 31}, {NT::BE, 1, 32, -1, -1},
      {NT::V, 1, 32, 2, 31}, {NT::V, 1, 32, -1, 31},
      {NT::V, 0, 32, -1, -1}, {NT::V, -1, 32, -1, -1},
      {NT::V, 2, 1, -1, -1}, {NT::V, 1, 0, -1, -1},
      {NT::V, 1, 33, -1, -1}, {static_cast<NT>(255), 1, 32, -1, -1},
      {NT::V, std::numeric_limits<cand_pos_t>::min(), 32, -1, -1},
      {NT::V, 1, std::numeric_limits<cand_pos_t>::max(), -1, -1}};
  for (const auto &key : unsupported) {
    set(index, reference, key, 1.25);
    erase(index, reference, key);
    compare(index, reference, key);
    require(index.row_capacity(key) == 0, "unsupported key exposed a row capacity");
  }
  require(index.size() == before_size && index.allocated_capacity() == before_capacity &&
          index.memory_bytes() == before_bytes && before_bytes > header_bytes,
          "unsupported updates changed primary lookup storage");
  check_counts(index, reference);
}

void collisions_and_wrap() {
  // These starts share the last bucket of the minimum eight-slot table.
  // Four insertions occupy slots 7,0,1,2, exercising wraparound deletion;
  // expected values come exclusively from the independent map below.
  std::array<cand_pos_t, 4> order{11, 16, 27, 32};
  do {
    Index index(128);
    Reference reference;
    for (const auto i : order) set(index, reference, {NT::WI, i, 128, -1, -1}, double(i) + 0.125);
    require(index.row_capacity({NT::WI, 1, 128, -1, -1}) == 8,
            "collision witness grew beyond the minimum table");
    for (const auto erased : order) {
      const Key key{NT::WI, erased, 128, -1, -1};
      erase(index, reference, key);
      for (cand_pos_t i = 1; i <= 128; ++i)
        compare(index, reference, {NT::WI, i, 128, -1, -1});
      erase(index, reference, key); // A second erase must not alter counts.
    }
    require(index.size() == 0 && index.allocated_capacity() == 0 && index.active_rows() == 0,
            "empty collision row kept a slot buffer");
    for (const auto i : order) set(index, reference, {NT::WI, i, 128, -1, -1}, -double(i));
    for (cand_pos_t i = 1; i <= 128; ++i)
      compare(index, reference, {NT::WI, i, 128, -1, -1});
    check_counts(index, reference);
  } while (std::next_permutation(order.begin(), order.end()));

  // Mix colliding and different-home keys inside one wrapped cluster. A
  // deletion must leave keys at their home or move only across their own
  // probe path, rather than shifting every subsequent slot unconditionally.
  order = {2, 4, 11, 16};
  do {
    Index index(128);
    Reference reference;
    for (const auto i : order) set(index, reference, {NT::VP, i, 128, -1, -1}, double(i));
    for (const auto i : order) {
      erase(index, reference, {NT::VP, i, 128, -1, -1});
      for (cand_pos_t query = 1; query <= 128; ++query)
        compare(index, reference, {NT::VP, query, 128, -1, -1});
    }
    check_counts(index, reference);
  } while (std::next_permutation(order.begin(), order.end()));
}

void random_updates_and_ownership() {
  constexpr cand_pos_t n = 257;
  Index index(n);
  Reference reference;
  std::mt19937 random(20261020);
  for (std::size_t step = 0; step < 30000; ++step) {
    // Most operations revisit one busy row; the others route across all
    // families/ends. Missing deletes, repeated updates and reinsertions occur
    // independently of the index's collision/deletion implementation.
    const bool hot = random() % 4 != 0;
    const cand_pos_t j = hot ? n : 1 + random() % n;
    const Key key{hot ? NT::WM : static_cast<NT>(random() % families),
                  1 + static_cast<cand_pos_t>(random() % j), j, -1, -1};
    if (random() % 3) {
      const double value = step % 17 == 0 ? -0.0 : double(static_cast<int>(random() % 2000) - 1000);
      set(index, reference, key, value);
    } else {
      erase(index, reference, key);
    }
    compare(index, reference, key);
    compare(index, reference, {key.nonterminal, 1 + static_cast<cand_pos_t>(random() % j), j, -1, -1});
    if (step % 257 == 0) check_counts(index, reference);
  }
  check_counts(index, reference);
  const Key key{NT::WM, 1, n, -1, -1};
  set(index, reference, key, -0.0);
  auto copy = index;
  copy.set(key, 17.0);
  compare(index, reference, key);
  require(copy.get(key, 0.0) == 17.0, "primary lookup copy aliases the original");
  copy.erase(key);
  compare(index, reference, key);
  copy = index;
  copy = copy;
  check_counts(copy, reference);
  auto moved = std::move(copy);
  require(copy.length() == 0 && copy.size() == 0 && copy.row_count() == 0 &&
          copy.allocated_capacity() == 0 && copy.memory_bytes() == 0,
          "moved primary lookup retained stale storage metadata");
  compare(copy, {}, key);
  copy.set(key, 9.0);
  copy.erase(key);
  compare(moved, reference, key);
  check_counts(moved, reference);

  Index assigned(3);
  assigned.set({NT::W, 1, 3, -1, -1}, 4.0);
  assigned = std::move(moved);
  assigned = std::move(assigned);
  require(moved.length() == 0 && moved.size() == 0 && moved.allocated_capacity() == 0,
          "move-assigned primary lookup retained stale metadata");
  check_counts(assigned, reference);
  for (const auto &[id, value] : reference) {
    const auto &[family, i, j, ip, jp] = id;
    const Key query{static_cast<NT>(family), i, j, ip, jp};
    compare(assigned, reference, query);
    require(same(*assigned.get(query, 0.0), value), "copy/move changed primary value bits");
  }
  copy = assigned;
  check_counts(copy, reference);
  assigned = Index(0);
  require(assigned.size() == 0 && assigned.memory_bytes() == 0,
          "reset primary lookup kept old buffers");
  compare(copy, reference, key);
}

void erased_columns_release_capacity() {
  constexpr cand_pos_t n = 256;
  Index index(n);
  const auto headers = index.memory_bytes();
  const Key first{NT::W, 1, 1, -1, -1};
  index.set(first, -0.0);
  const auto slot_bytes = (index.memory_bytes() - headers) / index.allocated_capacity();
  require(slot_bytes >= sizeof(cand_pos_t) + sizeof(double), "primary lookup omits slot storage bytes");
  index.erase(first);
  require(index.memory_bytes() == headers, "empty row did not release its capacity");

  std::size_t live = 0;
  for (std::size_t f = 0; f < families; ++f) {
    for (cand_pos_t j = 1; j <= n; ++j) {
      const auto nt = static_cast<NT>(f);
      for (cand_pos_t i = 1; i <= j; ++i)
        index.set({nt, i, j, -1, -1}, double(f * n + i));
      // Every row first grows with a full candidate column, then retains
      // zero or one value. Keeping the old peak would produce quadratic
      // allocated capacity despite a linear live-state count.
      for (cand_pos_t i = 1; i < j; ++i) index.erase({nt, i, j, -1, -1});
      if (j % 3 == 0) index.erase({nt, j, j, -1, -1});
      else ++live;
      require(index.size() == live && index.allocated_capacity() <= 8 * live,
              "erase-heavy columns accumulated candidate-peak capacity");
      require(index.memory_bytes() == headers + slot_bytes * index.allocated_capacity(),
              "constant-time primary memory accounting differs from allocated capacity");
    }
  }
  require(index.active_rows() == live, "erase-heavy active-row counter differs");
  for (std::size_t f = 0; f < families; ++f) {
    for (cand_pos_t j = 1; j <= n; ++j) {
      const Key key{static_cast<NT>(f), j, j, -1, -1};
      const auto value = index.get(key, -1.0);
      require(value.has_value() && same(*value, j % 3 == 0 ? -1.0 : double(f * n + j)),
              "erase-heavy row lost its surviving value");
      index.erase(key);
    }
  }
  require(index.size() == 0 && index.active_rows() == 0 && index.allocated_capacity() == 0 &&
          index.memory_bytes() == headers, "retired columns kept allocated slot buffers");
}

// This explicit boundary is the only promotion point. The independent map
// and existing comparison/count helpers cover both representations.
void post_selection_interval_contract() {
  constexpr cand_pos_t n = 32;
  const auto key = [](cand_pos_t i) { return Key{NT::VP, i, n, -1, -1}; };
  std::uint64_t nan_pattern = UINT64_C(0x7ff8000000000517);
  double nan;
  std::memcpy(&nan, &nan_pattern, sizeof(nan));
  const double defaults[]{-0.0, 0.0, -std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::infinity(), 17.25,
      std::numeric_limits<double>::denorm_min(), nan};
  const auto check_row = [&](const Index &index, const Reference &reference) {
    for (double default_value : defaults)
      for (cand_pos_t i = 1; i <= n; ++i)
        compare(index, reference, key(i), default_value);
    for (const auto &other : std::vector<Key>{{NT::BE, 5, n, 6, 9},
             {NT::VP, 5, n, 0, -1}, {NT::VP, 5, n, -1, 0},
             {NT::VP, 0, n, -1, -1}, {NT::VP, 1, n + 1, -1, -1}})
      compare(index, reference, other, nan);
    check_counts(index, reference);
  };
  Index index(n);
  Reference reference;
  const auto headers = index.memory_bytes();
  set(index, reference, key(5), -0.0);
  set(index, reference, key(7), nan);
  set(index, reference, key(8), -std::numeric_limits<double>::infinity());
  set(index, reference, key(9), std::numeric_limits<double>::infinity());
  const auto hash_bytes = index.memory_bytes();
  index.seal_row(NT::VP, n);
  const auto direct_bytes = index.memory_bytes();
  require(direct_bytes < hash_bytes && index.row_capacity(key(1)) == 5,
          "eligible interval did not replace its hash slots");
  check_row(index, reference); // Hole 6 and both sides use exact default bits.
  index.seal_row(NT::VP, n);
  index.seal_row(NT::BE, n);
  index.seal_row(NT::VP, 0);
  index.seal_row(NT::VP, n + 1);
  require(index.memory_bytes() == direct_bytes, "repeated/unsupported seal changed storage");

  auto copied = index;
  const Reference sealed_reference = reference;
  auto moved = std::move(copied);
  require(copied.length() == 0 && copied.memory_bytes() == 0 && copied.size() == 0,
          "moving an interval retained stale source metadata");
  Index assigned(3);
  assigned.set({NT::W, 1, 3, -1, -1}, 4.0);
  assigned = moved;
  check_row(moved, sealed_reference);
  check_row(assigned, sealed_reference);
  set(index, reference, key(7), -7.25); // Existing-key mutation demotes once.
  require(index.memory_bytes() == hash_bytes, "interval mutation did not restore hash storage");
  index.seal_row(NT::VP, n);
  require(index.memory_bytes() == hash_bytes, "mutated interval was promoted a second time");
  set(index, reference, key(1), -1e300); // New key outside the old interval.
  check_row(index, reference);
  check_row(moved, sealed_reference); // Copies retain independent direct owners.
  check_row(assigned, sealed_reference);

  Reference erased = sealed_reference;
  for (const auto i : {7, 5, 9, 8}) {
    erase(assigned, erased, key(i));
    check_row(assigned, erased);
  }
  require(assigned.memory_bytes() == headers && assigned.allocated_capacity() == 0,
          "erasing an interval retained direct/hash capacity");

  Index sparse(n);
  Reference sparse_reference;
  set(sparse, sparse_reference, key(1), -0.0);
  set(sparse, sparse_reference, key(n), 17.25);
  const auto sparse_bytes = sparse.memory_bytes();
  sparse.seal_row(NT::VP, n); // Extent 32 exceeds the current eight hash slots.
  require(sparse.memory_bytes() == sparse_bytes, "wide sparse row became a dense interval");
  erase(sparse, sparse_reference, key(1));
  sparse.seal_row(NT::VP, n);
  require(sparse.memory_bytes() == sparse_bytes, "an ineligible row retried its seal");
  check_row(sparse, sparse_reference);
}

}  // namespace

int main() {
  try {
    key_and_value_contract();
    collisions_and_wrap();
    random_updates_and_ownership();
    erased_columns_release_capacity();
    post_selection_interval_contract();
    std::cout << "scfg2_primary_lookup_index=ok comparisons=" << comparisons << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
