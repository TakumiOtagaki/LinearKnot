#include "scfg2/beam/outer_frontier.hh"

#include <cstring>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iterator>
#include <limits>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {
namespace outer = cparty::scfg2::beam::outer_detail;
namespace sc = cparty::scfg2;
using Item = outer::Item;
using NT = outer::NT;
using Link = outer::Outer::Link;
using Pool = outer::MaxPool;

// Frozen c58b311 MaxPool. Keep the original set implementation and every
// canonical-map mutation so differential checks cover bucket/iteration order
// as well as final maxima. Type aliases/qualification are adjusted locally.
struct ReferencePool {
  struct Entry { double value; Link link; bool mandatory; };
  struct Rank { Item key; double value; };
  struct Better {
    bool operator()(const Rank &a, const Rank &b) const {
      if (a.value != b.value) return a.value > b.value;
      return outer::less_key(a.key, b.key);
    }
  };
  using Entries = std::unordered_map<Item, Entry, sc::deductions::ItemKeyHash,
                                    sc::deductions::ItemKeyEq>;
  Entries entries;
  std::set<Rank, Better> novel;
  std::size_t evicted = 0;

  void offer(const Item &key, double value, const Link &link,
              bool mandatory, std::size_t width, bool prefilter = false) {
    const auto old = entries.find(key);
    if (prefilter && old == entries.end() && !mandatory && novel.size() == width &&
        (width == 0 || !Better{}({key, value}, *std::prev(novel.end())))) {
      ++evicted;
      return;
    }
    if (old != entries.end()) {
      if (old->second.mandatory != mandatory) throw std::runtime_error("pool ownership changed");
      if (value < old->second.value ||
          (value == old->second.value && !outer::better_link(link, old->second.link))) return;
      if (value == old->second.value) {
        old->second.link = link;
        return;
      }
      if (!mandatory) novel.erase({key, old->second.value});
      old->second = Entry{value, link, mandatory};
    } else {
      entries.insert_or_assign(key, Entry{value, link, mandatory});
    }
    if (!mandatory) {
      novel.insert({key, value});
      if (novel.size() > width) {
        const auto worst = std::prev(novel.end());
        entries.erase(worst->key);
        novel.erase(worst);
        ++evicted;
      }
    }
  }
};

std::size_t checked = 0;
void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(a)) == 0; }
bool same(const Item &a, const Item &b) { return sc::deductions::ItemKeyEq{}(a, b); }
template <class Actual, class Expected>
void compare_entries(const Actual &actual, const Expected &expected) {
  require(actual.size() == expected.size(), "max pool map size differs");
  require(actual.bucket_count() == expected.bucket_count(), "max pool map bucket count differs");
  auto a = actual.begin();
  auto b = expected.begin();
  for (; a != actual.end(); ++a, ++b) {
    require(same(a->first, b->first), "max pool canonical map iteration order differs");
    require(same(a->second.value, b->second.value) &&
            a->second.mandatory == b->second.mandatory &&
            same(a->second.link.parent, b->second.link.parent) &&
            same(a->second.link.edge, b->second.link.edge), "max pool entry value/link/ownership differs");
  }
}
void compare(const Pool &actual, const ReferencePool &expected) {
  compare_entries(actual.entries, expected.entries);
  require(actual.evicted == expected.evicted, "max pool eviction count differs");
  require(actual.novel.size() == expected.novel.size() &&
          actual.novel.empty() == expected.novel.empty(), "max pool novel count differs");
  auto a = actual.novel.begin();
  auto b = expected.novel.begin();
  for (; a != actual.novel.end(); ++a, ++b)
    require(same(a->key, b->key) && same(a->value, b->value), "max pool sorted novel rank differs");
  // Observe the cached view again, including std::prev(end()) compatibility.
  require(static_cast<std::size_t>(std::distance(actual.novel.begin(), actual.novel.end())) ==
          expected.novel.size(), "max pool repeated sorted view differs");
  if (!actual.novel.empty()) {
    const auto &worst = *std::prev(actual.novel.end());
    const auto &old_worst = *std::prev(expected.novel.end());
    require(same(worst.key, old_worst.key) && same(worst.value, old_worst.value),
            "max pool observed cutoff differs");
  }
  std::set<std::size_t> indices;
  for (const auto &[key, entry] : actual.entries) {
    (void)key;
    if (entry.mandatory) require(entry.heap_index == Pool::no_heap_index,
                                 "mandatory max pool item acquired a width slot");
    else require(entry.heap_index < actual.novel.size() && indices.insert(entry.heap_index).second,
                 "max pool heap index is stale or aliases another entry");
  }
  require(indices.size() == actual.novel.size(), "max pool heap index coverage differs");
  ++checked;
}

Item key_for(int id, bool be) {
  return be ? Item{NT::BE, 1, 200, 2 + id % 19, 100 + id / 19}
            : Item{NT::V, 1 + id, 200, -1, -1};
}
Link link_for(int id, double edge) {
  return {{id % 2 ? NT::VM : NT::W, 1 + id % 17, 300, -1, -1}, edge};
}
void offer(Pool &actual, ReferencePool &expected, const Item &key, double value,
           const Link &link, bool mandatory, std::size_t width, bool prefilter) {
  std::string actual_error, expected_error;
  try {
    if (prefilter) actual.offer(key, value, link, mandatory, width, true);
    else actual.offer(key, value, link, mandatory, width); // Exercise the default.
  } catch (const std::runtime_error &error) { actual_error = error.what(); }
  try { expected.offer(key, value, link, mandatory, width, prefilter); }
  catch (const std::runtime_error &error) { expected_error = error.what(); }
  require(actual_error == expected_error, "max pool offer rejection differs");
  compare(actual, expected);
}

void ownership_and_order(Pool &actual, ReferencePool &expected,
                         std::size_t width, bool prefilter, bool be) {
  const auto fixed = key_for(0, be);
  offer(actual, expected, fixed, -4.0, {{NT::VM, 17, 300, -1, -1}, 1.0}, true, width, prefilter);
  offer(actual, expected, fixed, -4.0, {{NT::W, 2, 300, -1, -1}, 2.0}, true, width, prefilter);
  // Same predecessor and equal score must keep the old edge, while a
  // strictly higher score wins even with a worse predecessor tie key.
  offer(actual, expected, fixed, -4.0, {{NT::W, 2, 300, -1, -1}, 7.0}, true, width, prefilter);
  offer(actual, expected, fixed, -3.0, {{NT::VM, 19, 300, -1, -1}, -0.0}, true, width, prefilter);
  offer(actual, expected, fixed, -4.0, {{NT::W, 1, 300, -1, -1}, 9.0}, true, width, prefilter);
  offer(actual, expected, fixed, 20.0, link_for(0, 3.0), false, width, prefilter);
  require(actual.entries.at(fixed).mandatory && same(actual.entries.at(fixed).value, -3.0),
          "ownership-rejection witness was vacuous");

  offer(actual, expected, key_for(1, be), -1.0, link_for(9, 1.0), false, width, prefilter);
  offer(actual, expected, key_for(2, be), -1.0, link_for(8, 2.0), false, width, prefilter);
  offer(actual, expected, key_for(1, be), -2.0, link_for(0, 3.0), false, width, prefilter);
  offer(actual, expected, key_for(1, be), -1.0, link_for(0, 4.0), false, width, prefilter);
  offer(actual, expected, key_for(3, be), 7.0, link_for(3, 5.0), false, width, prefilter);
  offer(actual, expected, key_for(2, be), 10.0, link_for(2, 6.0), false, width, prefilter);
  offer(actual, expected, key_for(3, be), 15.0, link_for(1, 7.0), false, width, prefilter);
  offer(actual, expected, key_for(4, be), -0.0, link_for(5, -0.0), false, width, prefilter);
  offer(actual, expected, key_for(4, be), 0.0, link_for(0, 0.0), false, width, prefilter);
  offer(actual, expected, key_for(5, be), -1e300, link_for(0, -1e300), false, width, prefilter);
}

void copy_move_and_reuse(Pool &actual, ReferencePool &expected,
                        std::size_t width, bool prefilter, bool be, int salt) {
  auto copy = actual;
  auto reference_copy = expected;
  compare(copy, reference_copy);
  for (int id = 1; id <= 9; ++id)
    offer(copy, reference_copy, key_for(id, be), 50.0 + id + salt,
          link_for(id, double(salt)), false, width, prefilter);
  compare(actual, expected); // Copy updates must not rewrite source metadata.
  for (int id = 9; id >= 1; --id)
    offer(actual, expected, key_for(id, be), 100.0 + id + salt,
          link_for(id, -double(salt)), false, width, prefilter);
  compare(copy, reference_copy);

  Pool assigned;
  ReferencePool reference_assigned;
  for (int id = 100; id < 145; ++id)
    offer(assigned, reference_assigned, key_for(id, be), double(id), link_for(id, 0.0),
          false, 100, false);
  assigned = copy;
  reference_assigned = reference_copy;
  compare(assigned, reference_assigned);
  auto moved = std::move(assigned);
  auto reference_moved = std::move(reference_assigned);
  compare(moved, reference_moved);
  compare(assigned, reference_assigned);
  require(assigned.novel.capacity() == 0 && assigned.novel.view_capacity() == 0,
          "moved max pool kept rank history capacity");
  assigned = std::move(moved);
  reference_assigned = std::move(reference_moved);
  compare(assigned, reference_assigned);
  for (int id = 9; id >= 1; --id)
    offer(assigned, reference_assigned, key_for(id, be), 200.0 + id + salt,
          link_for(id, double(salt)), false, width, prefilter);
  compare(copy, reference_copy);
  compare(actual, expected);

  // The production lifecycle moves the map before releasing novel. Its old
  // Entry pointers are still owned by parents, and clear must not touch them.
  auto parents = std::move(assigned.entries);
  auto reference_parents = std::move(reference_assigned.entries);
  assigned.novel.clear();
  reference_assigned.novel.clear();
  require(assigned.novel.capacity() == 0 && assigned.novel.view_capacity() == 0,
          "completed max pool retained heap or sorted-view capacity");
  compare(assigned, reference_assigned);
  compare_entries(parents, reference_parents);
  for (int id = 1; id <= 7; ++id)
    offer(assigned, reference_assigned, key_for(id, be), double(id), link_for(id, 1.0),
          false, width, prefilter);
  compare_entries(parents, reference_parents);
}

void trace(std::size_t initial_width, bool prefilter, bool be, bool changing_width,
           bool trusted_index = false) {
  Pool actual;
  if (trusted_index) actual.enable_trusted_index();
  ReferencePool expected;
  ownership_and_order(actual, expected, initial_width, prefilter, be);
  std::mt19937 random(20261020 + initial_width + 31 * be + 61 * prefilter);
  for (int step = 0; step < 1800; ++step) {
    std::size_t width = initial_width;
    if (changing_width) {
      constexpr std::size_t widths[]{20, 2, 1, 0, 20};
      width = widths[(step / 100) % 5];
    }
    const int id = 1 + random() % 96;
    const bool mandatory = id > 9 && id % 13 == 0;
    const auto value = double(static_cast<int>(random() % 29) - 14);
    offer(actual, expected, key_for(id, be), value,
          link_for(static_cast<int>(random() % 31), double(step)), mandatory, width, prefilter);
    if (step % 257 == 0)
      copy_move_and_reuse(actual, expected, width, prefilter, be, step);
  }
  require(actual.evicted > 0, "max pool eviction/re-entry trace was vacuous");
}

Item primary_key(int id, bool fixed_span) {
  const cand_pos_t start = 1 + id;
  return {NT::V, start, fixed_span ? start + 200 : 200, -1, -1};
}
void enable_primary(Pool &pool, bool fixed_span) {
  if (fixed_span) pool.enable_trusted_primary_span_index(NT::V, 200);
  else pool.enable_trusted_primary_end_index(NT::V, 200);
}
void primary_trace(bool fixed_span, bool prefilter) {
  Pool actual;
  ReferencePool expected;
  enable_primary(actual, fixed_span);
  offer(actual, expected, primary_key(0, fixed_span), -0.0,
        link_for(7, 1.0), true, 20, prefilter);
  offer(actual, expected, primary_key(0, fixed_span), 0.0,
        link_for(0, -0.0), true, 20, prefilter);
  offer(actual, expected, primary_key(0, fixed_span), 7.0,
        link_for(0, 2.0), false, 20, prefilter); // ownership must still throw
  std::mt19937 random(20261021 + 71 * fixed_span + 31 * prefilter);
  constexpr std::size_t widths[]{20, 2, 1, 0, 20};
  for (int step = 0; step < 1200; ++step) {
    const int id = 1 + random() % 96;
    const auto value = step % 127 == 0 ? -0.0
        : step % 149 == 0 ? -1e300 : double(static_cast<int>(random() % 41) - 20);
    offer(actual, expected, primary_key(id, fixed_span), value,
          link_for(random() % 31, double(step)), id > 12 && id % 13 == 0,
          widths[(step / 100) % 5], prefilter);
  }

  // Copies/moves must never retain pointers into another pool's map. Rebuild
  // a primary mirror before each ownership change and keep the source live.
  auto copy = actual;
  auto reference_copy = expected;
  enable_primary(copy, fixed_span);
  copy = copy;
  copy = std::move(copy);
  for (int id = 1; id <= 12; ++id)
    offer(copy, reference_copy, primary_key(id, fixed_span), 50.0 + id,
          link_for(id, -0.0), false, 20, prefilter);
  compare(actual, expected);
  Pool assigned;
  enable_primary(assigned, !fixed_span);
  assigned = copy;
  auto reference_assigned = reference_copy;
  enable_primary(assigned, fixed_span);
  auto moved = std::move(assigned);
  auto reference_moved = std::move(reference_assigned);
  compare(assigned, reference_assigned);
  compare(moved, reference_moved);
  enable_primary(moved, fixed_span);
  assigned = std::move(moved);
  reference_assigned = std::move(reference_moved);
  compare(moved, reference_moved);
  offer(assigned, reference_assigned, primary_key(1, fixed_span), 100.0,
        link_for(0, 3.0), false, 20, prefilter);
  enable_primary(assigned, fixed_span);
  assigned.disable_trusted_index();
  auto parents = std::move(assigned.entries);
  auto reference_parents = std::move(reference_assigned.entries);
  assigned.novel.clear();
  reference_assigned.novel.clear();
  compare(assigned, reference_assigned);
  compare_entries(parents, reference_parents);
  enable_primary(assigned, fixed_span);
  offer(assigned, reference_assigned, primary_key(1, fixed_span), -7.0,
        link_for(0, 2.0), false, 20, prefilter);
  compare_entries(parents, reference_parents);
}

void primary_full_key_fallback() {
  for (bool fixed_span : {false, true}) {
    const auto canonical = primary_key(64, fixed_span);
    std::vector<Item> outside;
    auto key = canonical; key.nonterminal = NT::W; outside.push_back(key);
    key = canonical; ++key.j; outside.push_back(key);
    key = canonical; key.ip = 3; outside.push_back(key);
    key = canonical; key.jp = 190; outside.push_back(key);
    key = canonical; key.i = 0; outside.push_back(key);
    key = canonical; key.i = -1; outside.push_back(key);
    key = canonical; key.j = key.i - 1; outside.push_back(key);
    key = canonical; key.nonterminal = NT::BE; key.ip = 70; key.jp = 190; outside.push_back(key);
    for (const auto &unsupported : outside) {
      Pool actual;
      ReferencePool expected;
      enable_primary(actual, fixed_span);
      // Mandatory entries exercise complete-key identity without requesting
      // unsupported mixed-family novel ranking from either implementation.
      offer(actual, expected, canonical, -4.0, link_for(8, 0.0), true, 2, false);
      offer(actual, expected, unsupported, 9.0, link_for(3, -0.0), true, 2, false);
      require(actual.entries.size() == 2, "noncanonical full key aliased primary start");
      enable_primary(actual, fixed_span); // mixed pre-existing map must fall back
      offer(actual, expected, canonical, -5.0, link_for(0, 17.0), true, 2, false);
      offer(actual, expected, unsupported, 10.0, link_for(5, 3.0), true, 2, false);
    }
    Pool be;
    ReferencePool expected;
    if (fixed_span) be.enable_trusted_primary_span_index(NT::BE, 200);
    else be.enable_trusted_primary_end_index(NT::BE, 200);
    for (int id = 1; id < 32; ++id)
      offer(be, expected, key_for(id, true), double(id % 9), link_for(id, 0.0), false, 4, true);
  }
  Pool ordinary;
  ReferencePool expected;
  const Item key{NT::V, 8, 200, 3, 100};
  ordinary.entries.emplace(key, Pool::Entry{7.0, link_for(3, 1.0), true});
  expected.entries.emplace(key, ReferencePool::Entry{7.0, link_for(3, 1.0), true});
  offer(ordinary, expected, key, 6.0, link_for(0, 2.0), true, 1, false);
}
}  // namespace

int main() {
  try {
    for (const bool be : {false, true}) {
      for (const bool prefilter : {false, true}) {
        for (std::size_t width : {0, 1, 2, 20}) trace(width, prefilter, be, false);
        trace(20, prefilter, be, true);
      }
      trace(20, true, be, true, true); // generic full-key trusted index
    }
    for (bool fixed_span : {false, true})
      for (bool prefilter : {false, true}) primary_trace(fixed_span, prefilter);
    primary_full_key_fallback();
    std::cout << "scfg2_outer_max_pool=ok compared_states=" << checked << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
