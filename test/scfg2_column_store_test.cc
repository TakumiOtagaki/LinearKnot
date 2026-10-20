#include "scfg2/beam/column_store.hh"
#include "scfg2/deductions/deduction.hh"

#include <cstdlib>
#include <iostream>
#include <limits>
#include <memory>
#include <new>
#include <random>
#include <stdexcept>
#include <unordered_map>

namespace allocation_fault {
enum class Kind { BadAlloc, LengthError };
thread_local long remaining = -1;
thread_local bool fired = false;
thread_local std::size_t calls = 0;
thread_local Kind kind = Kind::BadAlloc;
void arm(long before, Kind value) {
  remaining = before; fired = false; calls = 0; kind = value;
}
void disable() { remaining = -1; }
void *allocate(std::size_t bytes) {
  if (remaining >= 0) ++calls;
  if (remaining >= 0 && remaining-- == 0) {
    remaining = -1; fired = true;
    if (kind == Kind::LengthError) throw std::length_error("column metadata fault");
    throw std::bad_alloc{};
  }
  if (void *memory = std::malloc(bytes ? bytes : 1)) return memory;
  throw std::bad_alloc{};
}
}
void *operator new(std::size_t bytes) { return allocation_fault::allocate(bytes); }
void *operator new[](std::size_t bytes) { return allocation_fault::allocate(bytes); }
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete(void *memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory, std::size_t) noexcept { std::free(memory); }

namespace {
namespace sc = cparty::scfg2;
using NT = sc::grammar::NonTerminal;
using Key = sc::deductions::ItemKey;
using D = sc::deductions::Deduction;
using Value = std::vector<D>;
template <std::size_t Hint>
using StoreWithHint = sc::beam::ColumnStore<Value, NT::WM, Hint>;
using Store = StoreWithHint<0>;
using Reference = std::unordered_map<Key, Value, sc::deductions::ItemKeyHash,
                                   sc::deductions::ItemKeyEq>;
std::size_t comparisons = 0, mutations = 0;

void require(bool ok, const char *message) {
  if (!ok) { allocation_fault::disable(); throw std::runtime_error(message); }
}
bool same(const Key &a, const Key &b) { return sc::deductions::ItemKeyEq{}(a, b); }
bool same(const D &a, const D &b) {
  if (!same(a.parent, b.parent) || a.rule != b.rule || a.child_count != b.child_count ||
      a.split.k != b.split.k || a.split.l != b.split.l ||
      a.split.p != b.split.p || a.split.q != b.split.q) return false;
  for (std::size_t c = 0; c < a.children.size(); ++c)
    if (!same(a.children[c].item, b.children[c].item)) return false;
  return true;
}
bool same(const Value &a, const Value &b) {
  if (a.size() != b.size()) return false;
  for (std::size_t p = 0; p < a.size(); ++p) if (!same(a[p], b[p])) return false;
  return true;
}
D record(Key key, int serial) {
  D d;
  d.parent = key;
  d.rule = sc::grammar::RuleId::WM_SPLIT_V;
  d.split = {serial, serial + 1, serial + 2, serial + 3};
  d.child_count = serial % 4;
  for (std::size_t p = 0; p < d.children.size(); ++p)
    d.children[p].item = {NT::BE, serial + int(p) + 1, serial + int(p) + 9,
                         serial + int(p) + 2, serial + int(p) + 7};
  return d;
}
template <class StoreType>
void compare(const StoreType &store, const Reference &reference, const std::vector<Key> &probes) {
  require(store.size() == reference.size() && store.empty() == reference.empty(),
      "column store presence count differs from full-key oracle");
  const auto check = [&](const Key &key) {
    const auto *actual = store.find(key);
    const auto expected = reference.find(key);
    require((actual == store.end()) == (expected == reference.end()),
        "column store aliases or loses a complete key");
    if (actual) require(same(actual->first, expected->first) &&
        same(actual->second, expected->second), "column payload fields/order differ");
    ++comparisons;
  };
  for (const auto &[key, value] : reference) { (void)value; check(key); }
  for (const auto &key : probes) check(key);
}
std::vector<Key> keys(int column) {
  std::vector<Key> out;
  for (int i = 1; i <= column; ++i) out.push_back({NT::WM, i, column, -1, -1});
  for (int i : {-2, 0, 1, column, column + 1})
    for (int j : {-2, 0, column - 1, column, column + 1}) {
      out.push_back({NT::WM, i, j, -1, -1});
      out.push_back({NT::VM, i, j, -1, -1});
      out.push_back({NT::WM, i, j, 0, -1});
      out.push_back({NT::WM, i, j, -1, 0});
      out.push_back({NT(255), i, j, -1, -1});
    }
  return out;
}
template <class StoreType>
void copy_move_contract(const StoreType &source, const Reference &reference,
                        const std::vector<Key> &probes) {
  StoreType copy(source);
  compare(copy, reference, probes);
  for (const auto &[key, value] : reference) if (!value.empty()) {
    require(copy.find(key)->second.data() != source.find(key)->second.data(),
        "copied source buffer shares mutable ownership");
    break;
  }
  StoreType assigned;
  assigned.reset_column(3);
  assigned[{NT::WM, 1, 3, -1, -1}].push_back(record({}, 2));
  assigned = copy;
  compare(assigned, reference, probes);
  StoreType moved(std::move(copy));
  require(copy.column() == -1 && copy.empty(), "moved column source is not reusable");
  compare(moved, reference, probes);
  copy.reset_column(4);
  copy[{NT::WM, 1, 4, -1, -1}].push_back(record({}, 3));
  assigned = std::move(moved);
  require(moved.column() == -1 && moved.empty(), "move assignment retained source entries");
  compare(assigned, reference, probes);
  assigned = assigned;
  assigned = std::move(assigned);
  compare(assigned, reference, probes);
  moved.swap(assigned);
  require(assigned.empty(), "swap retained previous private values");
  compare(moved, reference, probes);
}
template <std::size_t Hint>
void differential() {
  StoreWithHint<Hint> store;
  Reference reference;
  std::mt19937 rng(20261020);
  for (int column : {1, 12, 64, 7, 129, 257, 31, 500, 2, 256}) {
    store.reset_column(column);
    reference.clear();
    const auto probes = keys(column);
    compare(store, reference, probes);
    for (int step = 0; step < 2500; ++step) {
      const auto &key = probes[rng() % probes.size()];
      const auto d = record(key, step);
      switch (rng() % 5) {
        case 0: {
          auto [actual, inserted] = store.try_emplace(key);
          auto [expected, oracle_inserted] = reference.try_emplace(key);
          require(inserted == oracle_inserted, "try_emplace changed first-arrival status");
          actual->second.push_back(d); expected->second.push_back(d);
          break;
        }
        case 1: store[key].push_back(d); reference[key].push_back(d); break;
        case 2: {
          auto *actual = store.find(key);
          auto expected = reference.find(key);
          if (expected == reference.end()) require(!actual, "missing take acquired a source");
          else {
            auto got = std::move(actual->second);
            auto want = std::move(expected->second);
            require(same(got, want), "source transfer changed candidates or order");
            require(actual == store.find(key), "empty consumed source lost presence");
          }
          break;
        }
        case 3: store[key] = Value{d, d}; reference[key] = Value{d, d}; break;
        default: compare(store, reference, {key}); break;
      }
      ++mutations;
      if (step % 173 == 0) compare(store, reference, probes);
    }
    compare(store, reference, probes);
    copy_move_contract(store, reference, probes);
    store.clear(); reference.clear();
    compare(store, reference, probes);
  }
}
template <std::size_t Hint>
void pointer_stability() {
  StoreWithHint<Hint> store;
  store.reset_column(128);
  const Key primary{NT::WM, 1, 128, -1, -1}, extra{NT::VM, 1, 128, -1, -1};
  auto *a = store.try_emplace(primary).first;
  auto *b = store.try_emplace(extra).first;
  for (int i = 2; i <= 128; ++i) {
    store.try_emplace({NT::WM, i, 128, -1, -1});
    store.try_emplace({NT::VM, i, 128, -1, -1});
  }
  require(store.find(primary) == a && store.find(extra) == b,
      "in-column insertion invalidated an Entry pointer");
}
struct OwnedValue {
  inline static int live = 0;
  inline static bool fail = false;
  std::shared_ptr<int> token;
  OwnedValue() = delete;
  explicit OwnedValue(std::shared_ptr<int> token) : token(std::move(token)) {
    if (fail) throw std::runtime_error("owned constructor");
    ++live;
  }
  OwnedValue(const OwnedValue &other) : token(other.token) { ++live; }
  OwnedValue(OwnedValue &&other) noexcept : token(std::move(other.token)) { ++live; }
  ~OwnedValue() { --live; }
};
void value_ownership() {
  sc::beam::ColumnStore<OwnedValue, NT::VM> store;
  store.reset_column(500);
  require(OwnedValue::live == 0, "metadata constructed semiring Values");
  auto token = std::make_shared<int>(7);
  std::weak_ptr<int> weak(token);
  const Key key{NT::VM, 1, 500, -1, -1};
  OwnedValue::fail = true;
  bool threw = false;
  try { store.try_emplace(key, token); } catch (const std::runtime_error &) { threw = true; }
  OwnedValue::fail = false;
  require(threw && store.empty() && !store.find(key) && OwnedValue::live == 0,
      "failed Value construction left a present optional or touched record");
  store.try_emplace(key, token);
  auto copy = store;
  token.reset();
  require(OwnedValue::live == 2 && weak.use_count() == 2,
      "private Value copy/checkpoint ownership changed");
  store.clear();
  require(OwnedValue::live == 1 && !weak.expired(), "clear destroyed another store's Value");
  auto moved = std::move(copy);
  require(copy.empty() && copy.column() == -1 && OwnedValue::live == 1,
      "Value move created or retained an identity");
  moved.reset_column(12);
  require(OwnedValue::live == 0 && weak.expired(), "next column retained private forest references");
}
template <std::size_t Hint>
void metadata_failure() {
  namespace fault = allocation_fault;
  for (const auto kind : {fault::Kind::BadAlloc, fault::Kind::LengthError})
    for (long failure : {0, 1})
      for (bool warmed : {false, true}) {
        StoreWithHint<Hint> store;
        Reference reference;
        if (warmed) {
          store.reset_column(4);
          store[{NT::WM, 1, 4, -1, -1}].push_back(record({}, 1));
        }
        fault::arm(failure, kind);
        store.reset_column(129);
        const bool fired = fault::fired;
        fault::disable();
        require(fired && store.empty() && store.column() == 129,
            "metadata fault escaped fallback or retained the old column");
        const auto probes = keys(129);
        for (const auto &key : probes) {
          store[key].push_back(record(key, 5));
          reference[key].push_back(record(key, 5));
        }
        compare(store, reference, probes);
        copy_move_contract(store, reference, probes);
        store.reset_column(257); reference.clear();
        store[{NT::WM, 3, 257, -1, -1}].push_back(record({}, 9));
        reference[{NT::WM, 3, 257, -1, -1}].push_back(record({}, 9));
        compare(store, reference, keys(257));
      }
}

void metadata_growth() {
  namespace fault = allocation_fault;
  Store store;
  fault::arm(std::numeric_limits<long>::max(), fault::Kind::BadAlloc);
  for (int column = 0; column <= 1024; ++column) store.reset_column(column);
  const auto allocations = fault::calls;
  fault::disable();
  require(allocations <= 32 && store.empty() && store.column() == 1024,
      "column metadata allocations do not grow geometrically");
}

void capacity_hint_contract() {
  namespace fault = allocation_fault;
  const Key primary{NT::WM, 1, 17, -1, -1};
  Store ordinary;
  ordinary.reset_column(17);
  fault::arm(0, fault::Kind::BadAlloc);
  auto *plain = ordinary.try_emplace(primary).first;
  const bool default_allocated = fault::fired;
  fault::disable();
  require(!default_allocated && plain->second.capacity() == 0,
      "neutral ColumnStore default allocated a source buffer");

  StoreWithHint<8> hinted;
  hinted.reset_column(17);
  auto *empty = hinted.try_emplace(primary).first;
  require(empty->second.empty() && empty->second.capacity() >= 8,
      "explicit source hint did not reserve an empty buffer");
  Value supplied;
  supplied.reserve(3);
  const auto supplied_capacity = supplied.capacity();
  auto *preserved = hinted.try_emplace({NT::WM, 2, 17, -1, -1}, std::move(supplied)).first;
  require(preserved->second.empty() && preserved->second.capacity() == supplied_capacity,
      "source hint replaced a caller-supplied allocation");
  Value nonempty{record(primary, 1)};
  const auto nonempty_capacity = nonempty.capacity();
  auto *populated = hinted.try_emplace({NT::WM, 3, 17, -1, -1}, std::move(nonempty)).first;
  require(populated->second.size() == 1 && populated->second.capacity() == nonempty_capacity,
      "source hint changed an initialized Value");

  // Failure constructing the Value is still a real insertion failure. Only
  // the later optional reserve is swallowed, after a valid empty Value exists.
  const Key fresh{NT::WM, 4, 17, -1, -1};
  const Value seed{record(fresh, 2)};
  fault::arm(0, fault::Kind::BadAlloc);
  bool threw = false;
  try { hinted.try_emplace(fresh, seed); } catch (const std::bad_alloc &) { threw = true; }
  fault::disable();
  require(threw && !hinted.find(fresh), "source hint swallowed Value construction failure");

  for (const auto kind : {fault::Kind::BadAlloc, fault::Kind::LengthError}) {
    // Direct slots, a noncanonical full-key map, and a canonical key in a
    // metadata-fallback column must all tolerate the optional reserve failure.
    for (int route : {0, 1, 2}) {
      bool observed_hint_failure = false;
      for (long failure = 0; failure < 8 && !observed_hint_failure; ++failure) {
        StoreWithHint<8> store;
        store.reset_column(17);
        if (route == 2) {
          fault::arm(0, kind);
          store.reset_column(129);
          const bool metadata_failed = fault::fired;
          fault::disable();
          require(metadata_failed, "source hint witness missed metadata fallback");
        }
        const int column = route == 2 ? 129 : 17;
        const NT family = route == 1 ? NT::VM : NT::WM;
        const Key target{family, 2, column, -1, -1};
        const Key bootstrap{family, 1, column, -1, -1};
        Reference reference;
        if (route != 0) { store.try_emplace(bootstrap); reference.try_emplace(bootstrap); }
        typename StoreWithHint<8>::Entry *entry = nullptr;
        bool insertion_threw = false;
        fault::arm(failure, kind);
        try { entry = store.try_emplace(target).first; }
        catch (const std::bad_alloc &) { insertion_threw = true; }
        catch (const std::length_error &) { insertion_threw = true; }
        const bool fired = fault::fired;
        fault::disable();
        if (insertion_threw) {
          require(!store.find(target), "failed map insertion left a private source");
          continue;
        }
        reference.try_emplace(target);
        compare(store, reference, {bootstrap, target});
        if (!fired) continue;
        observed_hint_failure = true;
        require(entry && entry->second.empty() && entry->second.capacity() == 0,
            "optional source reserve failure altered the inserted Value");
        fault::arm(0, kind);
        auto [again, inserted] = store.try_emplace(target);
        const bool retried = fault::fired;
        fault::disable();
        require(again == entry && !inserted && !retried,
            "duplicate source insertion retried optional allocation");
        entry->second.push_back(record(target, 8));
        reference[target].push_back(record(target, 8));
        compare(store, reference, {bootstrap, target});
      }
      require(observed_hint_failure, "source reserve fault was not exercised on every route");
    }
  }
}
}

int main() {
  try {
    differential<0>(); differential<8>();
    pointer_stability<0>(); pointer_stability<8>();
    value_ownership(); metadata_failure<0>(); metadata_failure<8>();
    metadata_growth(); capacity_hint_contract();
    std::cout << "scfg2_column_store=ok mutations=" << mutations
              << " comparisons=" << comparisons << '\n';
  } catch (const std::exception &error) {
    allocation_fault::disable();
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
