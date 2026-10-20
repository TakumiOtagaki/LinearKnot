#include "scfg2/engine/chart.hh"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <limits>
#include <new>
#include <optional>
#include <random>
#include <utility>
#include <vector>

// Fail one allocation only while an explicitly armed Chart::set is running.
// Outside that call this executable uses the ordinary malloc/free allocator.
namespace allocation_fault {
thread_local long remaining = -1;
thread_local bool fired = false;

void arm(long successful_allocations) {
  remaining = successful_allocations;
  fired = false;
}
void disable() { remaining = -1; }

void *allocate(std::size_t bytes) {
  if (remaining >= 0) {
    if (remaining-- == 0) {
      remaining = -1;
      fired = true;
      throw std::bad_alloc{};
    }
  }
  if (void *memory = std::malloc(bytes ? bytes : 1)) return memory;
  throw std::bad_alloc{};
}
}  // namespace allocation_fault

void *operator new(std::size_t bytes) { return allocation_fault::allocate(bytes); }
void *operator new[](std::size_t bytes) { return allocation_fault::allocate(bytes); }
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete(void *memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory, std::size_t) noexcept { std::free(memory); }

namespace {
namespace sc = cparty::scfg2;
using Chart = sc::engine::Chart<double>;
using Key = sc::deductions::ItemKey;
using NT = sc::grammar::NonTerminal;
constexpr cand_pos_t length = 17;
constexpr double log_zero = -std::numeric_limits<double>::infinity();
std::size_t comparisons = 0, mutations = 0;
std::size_t reference_discard_calls = 0, indexed_discard_calls = 0;

void require(bool ok, const char *message) {
  if (!ok) {
    allocation_fault::disable();
    std::cerr << message << '\n';
    std::exit(EXIT_FAILURE);
  }
}
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(a)) == 0; }
bool equal(const Key &a, const Key &b) { return sc::deductions::ItemKeyEq{}(a, b); }

bool reference_discard(const double &value) {
  ++reference_discard_calls;
  return value == log_zero;
}
bool indexed_discard(const double &value) {
  ++indexed_discard_calls;
  return value == log_zero;
}
bool discard_marker(const double &value) { return value == 99.0; }

std::vector<Key> all_keys() {
  std::vector<Key> keys;
  for (int family = 0; family < static_cast<int>(NT::BE); ++family)
    for (cand_pos_t j = 1; j <= length; ++j)
      for (cand_pos_t i = 1; i <= j; ++i)
        keys.push_back({static_cast<NT>(family), i, j, -1, -1});
  // Distinct secondary identities and malformed/out-of-route primary spans
  // must retain the old unordered_map behavior, without routing aliases.
  for (const Key key : std::vector<Key>{
           {NT::BE, 1, length, 2, 9}, {NT::BE, 1, length, 3, 9},
           {NT::BE, 1, length, 2, 10}, {NT::BE, 1, length, -1, -1},
           {NT::V, 1, length, 2, 9}, {NT::V, 1, length, -1, 0},
           {NT::W, 0, 7, -1, -1}, {NT::WI, -2, 7, -1, -1},
           {NT::WIP, 2, 1, -1, -1}, {NT::WM, 1, length + 1, -1, -1},
           {NT::V, 1, 0, -1, -1}, {NT::V, 1, -2, -1, -1},
           {static_cast<NT>(255), 1, 5, -1, -1}})
    keys.push_back(key);
  return keys;
}

void compare(const Chart &reference, const Chart &indexed,
             const std::vector<Key> &keys) {
  require(reference.size() == indexed.size(), "lookup changed Chart::size");
  const auto &expected = reference.values();
  const auto &actual = indexed.values();
  require(expected.size() == actual.size(), "lookup duplicated or lost map entries");
  require(expected.bucket_count() == actual.bucket_count(),
          "lookup changed canonical map bucket growth");
  require(expected.max_load_factor() == actual.max_load_factor(),
          "lookup changed canonical map load policy");
  auto a = expected.begin(), b = actual.begin();
  while (a != expected.end()) {
    require(b != actual.end() && equal(a->first, b->first) && same(a->second, b->second),
            "lookup changed canonical map iteration or value bits");
    ++a;
    ++b;
  }
  require(b == actual.end(), "lookup added canonical iteration entries");
  auto cursor = expected.begin();
  indexed.for_each_entry([&](const Key &key, double value) {
    require(cursor != expected.end() && equal(key, cursor->first) && same(value, cursor->second),
            "lookup changed for_each_entry iteration");
    ++cursor;
  });
  require(cursor == expected.end(), "lookup iteration omitted canonical entries");
  const auto bytes = indexed.primary_lookup_bytes();
  for (const auto &key : keys) {
    require(reference.contains(key) == indexed.contains(key), "lookup changed presence");
    require(same(reference.get(key), indexed.get(key)), "lookup returned stale or aliased bits");
    ++comparisons;
  }
  require(indexed.primary_lookup_bytes() == bytes,
          "read-only lookup allocated or changed cache storage");
}

void set_pair(Chart &reference, Chart &indexed, const Key &key, double value) {
  reference.set(key, value);
  indexed.set(key, value);
  ++mutations;
}

void differential(double default_value, Chart::discard_fn discard_reference,
                  Chart::discard_fn discard_indexed, std::uint32_t seed) {
  Chart reference(default_value, discard_reference);
  Chart indexed(default_value, discard_indexed, 0, length);
  require(reference.primary_lookup_length() == 0 && reference.primary_lookup_bytes() == 0,
          "default Chart enabled lookup storage");
  require(indexed.primary_lookup_length() == length && indexed.dense_length() == 0,
          "lookup opt-in enabled the wrong storage mode");
  const auto keys = all_keys();
  compare(reference, indexed, keys);
  // Every primary route and all fallback shapes get at least one explicit
  // insertion, regardless of subsequent random coverage.
  for (std::size_t k = 0; k < keys.size(); ++k) {
    const double value = double(int(k % 97) - 48) / 8.0;
    set_pair(reference, indexed, keys[k], value);
    require(indexed.contains(keys[k]) && same(reference.get(keys[k]), indexed.get(keys[k])),
            "lookup route insertion failed");
  }
  require(indexed.primary_lookup_bytes() > 0, "populated lookup has no observed storage");
  compare(reference, indexed, keys);
  const double values[]{log_zero, 0.0, -0.0, 99.0, 1.25, -7.5, -1e300, 1e300,
      std::numeric_limits<double>::denorm_min(), -std::numeric_limits<double>::denorm_min(),
      std::numeric_limits<double>::infinity()};
  std::mt19937 rng(seed);
  for (int step = 0; step < 1200; ++step) {
    const auto &key = keys[rng() % keys.size()];
    set_pair(reference, indexed, key, values[rng() % (sizeof(values) / sizeof(values[0]))]);
    compare(reference, indexed, keys);
  }
  // Erase and reinsert in a different order. Erasure must update the lookup,
  // while reinsertion must retain unordered_map's original iteration rules.
  for (std::size_t k = 0; k < keys.size(); k += 3)
    set_pair(reference, indexed, keys[k], discard_reference == discard_marker ? 99.0 : log_zero);
  compare(reference, indexed, keys);
  for (std::size_t k = keys.size(); k > 0; --k)
    set_pair(reference, indexed, keys[k - 1], -double(k % 53));
  compare(reference, indexed, keys);
}

void zero_and_default_contract() {
  const Key key{NT::W, 1, length, -1, -1};
  for (double default_value : {17.25, log_zero, -0.0}) {
    Chart reference(default_value, discard_marker);
    Chart indexed(default_value, discard_marker, 0, length);
    const auto &view = indexed.values();
    for (double value : {3.0, 99.0, 0.0, -0.0, log_zero, 99.0}) {
      set_pair(reference, indexed, key, value);
      compare(reference, indexed, {key});
      require(view.size() == reference.size(), "canonical map reference stopped following set");
      if (value == 99.0) {
        require(!indexed.contains(key) && same(indexed.get(key), default_value),
                "discard did not restore the exact default value");
      } else {
        require(indexed.contains(key) && same(view.at(key), value),
                "explicit zero/default value lost presence or sign");
      }
    }
  }
}

void copy_and_move_contract() {
  const auto keys = all_keys();
  std::optional<Chart> copied_reference, copied_indexed;
  Chart assigned_reference(123.0);
  Chart assigned_indexed(123.0, nullptr, 0, 3);
  {
    Chart reference(log_zero, reference_discard);
    Chart indexed(log_zero, indexed_discard, 0, length);
    for (std::size_t k = 0; k < keys.size(); ++k)
      set_pair(reference, indexed, keys[k], double(k % 41));
    for (std::size_t k = 0; k < keys.size(); k += 5)
      set_pair(reference, indexed, keys[k], log_zero);
    copied_reference.emplace(reference);
    copied_indexed.emplace(indexed);
    assigned_reference = reference;
    assigned_indexed = indexed;
    require(copied_indexed->primary_lookup_length() == length &&
                assigned_indexed.primary_lookup_length() == length,
            "copy lost lookup configuration");
    // The copies' canonical maps stay unchanged; a pointer-based index copied
    // without rebinding would incorrectly observe these source mutations.
    for (std::size_t k = 0; k < keys.size(); ++k)
      set_pair(reference, indexed, keys[k], -123.0);
    compare(*copied_reference, *copied_indexed, keys);
    compare(assigned_reference, assigned_indexed, keys);
  }
  // Source maps have now been destroyed.
  compare(*copied_reference, *copied_indexed, keys);
  compare(assigned_reference, assigned_indexed, keys);
  Chart moved_reference(std::move(*copied_reference));
  Chart moved_indexed(std::move(*copied_indexed));
  require(moved_indexed.primary_lookup_length() == length &&
              copied_indexed->primary_lookup_length() == 0,
          "move constructor did not transfer/reset lookup ownership");
  compare(moved_reference, moved_indexed, keys);
  set_pair(*copied_reference, *copied_indexed, keys.front(), 13.0);
  compare(*copied_reference, *copied_indexed, keys);

  // Replace a dense destination too: the incoming map/index must not coexist
  // with the destination's former dense values or its different default.
  Chart destination_reference(5.0, nullptr, 3);
  Chart destination_indexed(5.0, nullptr, 3);
  destination_reference = std::move(moved_reference);
  destination_indexed = std::move(moved_indexed);
  require(destination_indexed.dense_length() == 0 &&
              destination_indexed.primary_lookup_length() == length &&
              moved_indexed.primary_lookup_length() == 0,
          "move assignment retained the wrong storage owner");
  compare(destination_reference, destination_indexed, keys);
  set_pair(moved_reference, moved_indexed, keys.back(), -4.0);
  compare(moved_reference, moved_indexed, keys);

  destination_reference = destination_reference;
  destination_indexed = destination_indexed;
  destination_reference = std::move(destination_reference);
  destination_indexed = std::move(destination_indexed);
  compare(destination_reference, destination_indexed, keys);
  std::swap(destination_reference, moved_reference);
  std::swap(destination_indexed, moved_indexed);
  compare(destination_reference, destination_indexed, keys);
  compare(moved_reference, moved_indexed, keys);
}

void allocation_failure_contract() {
  const Key bootstrap{NT::BE, 1, length, 2, 9};
  const Key primary{NT::VP, 3, length, -1, -1};
  bool observed_disabled_cache = false;
  // Discover the cache-row allocation without depending on a particular
  // standard library's node/bucket allocation count. The bootstrap creates
  // the canonical map first; the primary row is still unallocated.
  for (long failure = 0; failure < 16 && !observed_disabled_cache; ++failure) {
    Chart reference(log_zero), indexed(log_zero, nullptr, 0, length);
    set_pair(reference, indexed, bootstrap, 4.0);
    bool threw = false;
    allocation_fault::arm(failure);
    try { indexed.set(primary, -7.25); }
    catch (const std::bad_alloc &) { threw = true; }
    const bool fired = allocation_fault::fired;
    allocation_fault::disable();
    if (indexed.contains(primary)) reference.set(primary, -7.25);
    else require(threw, "successful set lost a key after allocation failure");
    compare(reference, indexed, {bootstrap, primary});
    if (fired && indexed.contains(primary) && indexed.primary_lookup_length() == 0) {
      observed_disabled_cache = true;
      require(same(indexed.get(primary), -7.25), "cache failure exposed stale numeric value");
      set_pair(reference, indexed, primary, -9.0);
      compare(reference, indexed, {bootstrap, primary});
    }
  }
  require(observed_disabled_cache, "allocation witness did not exercise cache disable/fallback");
}
}  // namespace

int main() {
  differential(log_zero, nullptr, nullptr, 20261021);
  differential(log_zero, reference_discard, indexed_discard, 20261022);
  differential(-0.0, discard_marker, discard_marker, 20261023);
  zero_and_default_contract();
  copy_and_move_contract();
  require(reference_discard_calls == indexed_discard_calls,
          "lookup changed discard callback count");
  allocation_failure_contract();
  std::cout << "scfg2_numeric_primary_lookup=ok mutations=" << mutations
            << " compared_keys=" << comparisons << '\n';
}
