#include "scfg2/beam/bounded_outer.hh"

#include <cstring>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <limits>
#include <new>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

// Fail one optional-index allocation in a narrowly armed call. Normal test
// setup and the frozen sorter keep the ordinary malloc/free allocator.
namespace proposal_allocation_fault {
enum class Kind { BadAlloc, LengthError };
thread_local long remaining = -1;
thread_local bool fired = false;
thread_local std::size_t calls = 0;
thread_local Kind kind = Kind::BadAlloc;
void arm(long before, Kind failure_kind) {
  remaining = before; fired = false; calls = 0; kind = failure_kind;
}
void disable() { remaining = -1; }
void *allocate(std::size_t bytes) {
  if (remaining >= 0) ++calls;
  if (remaining >= 0 && remaining-- == 0) {
    remaining = -1; fired = true;
    if (kind == Kind::LengthError) throw std::length_error("proposal allocation fault");
    throw std::bad_alloc{};
  }
  if (void *memory = std::malloc(bytes ? bytes : 1)) return memory;
  throw std::bad_alloc{};
}
}
void *operator new(std::size_t bytes) { return proposal_allocation_fault::allocate(bytes); }
void *operator new[](std::size_t bytes) { return proposal_allocation_fault::allocate(bytes); }
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete(void *memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory, std::size_t) noexcept { std::free(memory); }

namespace {
namespace sc = cparty::scfg2;
namespace outer = sc::beam::outer_detail;
using Key = sc::deductions::ItemKey;
using NT = sc::grammar::NonTerminal;

void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
bool same(double a, double b) { return std::memcmp(&a, &b, sizeof(a)) == 0; }

void reference_sort_unique(std::vector<int> &input) {
  std::sort(input.begin(), input.end());
  input.erase(std::unique(input.begin(), input.end()), input.end());
}

std::size_t proposal_position_contract() {
  std::size_t compared = 0;
  std::mt19937 rng(20261021);
  for (const int n : {0, 1, 60, 61, 62, 250, 500, 1000000}) {
    outer::ProposalPositionBits bits(n);
    require(bits.word_count() == 0, "proposal bitset allocated before use");
    const auto compare = [&](std::vector<int> input) {
      auto expected = input;
      reference_sort_unique(expected);
      bits.sort_unique(input);
      require(input == expected, "proposal bitset differs from frozen sort/unique");
      require(bits.active_word_count() == 0, "proposal bitset retained dirty words");
      ++compared;
    };
    compare({});
    compare({-2, -1, 0, 1, n, n + 1, -2, n + 1});
    require(bits.word_count() == (static_cast<std::size_t>(n) + 67) / 64,
        "proposal bitset metadata extent differs");
    compare({n + 1, n + 1, -2, -2});
    compare({-3, n + 2, -3, n + 2, -2, -1});
    compare({std::numeric_limits<int>::min(), std::numeric_limits<int>::max(), 0});
    std::uniform_int_distribution<int> position(-2, n + 1);
    for (int repeat = 0; repeat < 1200; ++repeat) {
      std::vector<int> input;
      const auto size = rng() % 700;
      input.reserve(size + 2);
      for (std::size_t p = 0; p < size; ++p) input.push_back(position(rng));
      if (repeat % 7 == 0) input.push_back(-3);
      if (repeat % 11 == 0) input.push_back(n + 2);
      if (repeat % 13 == 0) {
        input.push_back(std::numeric_limits<int>::min());
        input.push_back(std::numeric_limits<int>::max());
      }
      compare(std::move(input));
      compare({-2, n + 1});  // Sparse reuse after a wide or fallback parent.
    }
  }
  return compared;
}

void proposal_allocation_contract() {
  namespace fault = proposal_allocation_fault;
  // Both metadata and touched-ID allocation are optional. The original
  // integer sort can finish without either allocation or altered input.
  for (const auto kind : {fault::Kind::BadAlloc, fault::Kind::LengthError}) {
    for (long failure : {0, 1}) {
      outer::ProposalPositionBits bits(500);
      std::vector<int> input{501, -2, 61, 500, 0, 250, 61, 60, -1};
      auto expected = input;
      reference_sort_unique(expected);
      fault::arm(failure, kind);
      bits.sort_unique(input);
      const bool fired = fault::fired;
      fault::disable();
      require(fired && input == expected && !bits.enabled() && bits.word_count() == 0 &&
          bits.active_word_count() == 0 && bits.touched_capacity() == 0,
          "proposal allocation failure did not release index and preserve sorted output");
      input = {500, -2, 61, 61};
      expected = input;
      reference_sort_unique(expected);
      fault::arm(0, kind);
      bits.sort_unique(input);
      const bool retried = fault::fired;
      fault::disable();
      require(!retried && input == expected, "disabled proposal index retried allocation");
    }
    outer::ProposalPositionBits bits(1000000);
    std::vector<int> warm{-2, -1};
    bits.sort_unique(warm);
    std::vector<int> input;
    for (int p = 255; p >= 0; --p) input.push_back(p * 64 - 2);
    auto expected = input;
    reference_sort_unique(expected);
    fault::arm(0, kind);  // Grow touched IDs after metadata already exists.
    bits.sort_unique(input);
    const bool fired = fault::fired;
    fault::disable();
    require(fired && input == expected && !bits.enabled() && bits.word_count() == 0 &&
        bits.active_word_count() == 0,
        "proposal scratch growth failure left dirty words or changed positions");
  }
}

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

std::size_t rank_rows() {
  std::size_t queries = 0;
  const double values[]{-1e300, -0.0, 0.0, std::numeric_limits<double>::infinity(),
      std::numeric_limits<double>::quiet_NaN(), -7.125};
  for (int n : {1, 7, 8, 9, 63, 64, 65, 66, 127, 128, 129, 250, 500})
    for (std::size_t family = 0; family < static_cast<std::size_t>(NT::BE); ++family) {
      const auto nt = static_cast<NT>(family);
      std::vector<outer::Record> records;
      for (int i = 1; i <= n; ++i)
        if (i % 5 != 2) records.push_back({{nt, i, n, -1, -1}, values[i % 6]});
      outer::PrimaryAlphaRankRow row;
      row.build(records, nt, n);
      require(row.enabled() && row.word_count() <= records.size(),
          "dense canonical alpha row failed indexing or exceeded its record budget");
      for (int i = -2; i <= n + 2; ++i) {
        const auto found = std::lower_bound(records.begin(), records.end(), i,
            [](const outer::Record &record, int start) { return record.key.i < start; });
        const bool present = found != records.end() && found->key.i == i;
        const auto rank = row.rank(i);
        require(bool(rank) == present, "alpha rank membership differs from immutable records");
        if (rank) require(*rank == static_cast<std::size_t>(found - records.begin()) &&
            same(records[*rank].alpha, found->alpha), "alpha rank changed order/value bits");
        ++queries;
      }
      require(!row.rank(std::numeric_limits<int>::min()) &&
          !row.rank(std::numeric_limits<int>::max()), "far alpha rank coordinate aliased a record");

      for (int malformed = 0; malformed < 7; ++malformed) {
        auto invalid = records;
        switch (malformed) {
          case 0: invalid.front().key.nonterminal = NT::BE; break;
          case 1: invalid.front().key.j = n + 1; break;
          case 2: invalid.front().key.ip = 0; break;
          case 3: invalid.front().key.jp = 0; break;
          case 4: invalid.front().key.i = 0; break;
          case 5: invalid.back().key.i = n + 1; break;
          default: invalid.insert(invalid.begin(), invalid.front()); break;
        }
        row.build(invalid, nt, n);
        require(!row.enabled() && row.word_count() == 0 && row.capacity() == 0,
            "mixed/duplicate/noncanonical alpha row retained a rank index");
      }
      row.build(records, NT::BE, n);
      require(!row.enabled(), "secondary-family alpha row enabled primary rank indexing");
    }

  outer::PrimaryAlphaRankRow row;
  const std::vector<outer::Record> sparse{{{NT::V, 1, 1000000000, -1, -1}, -1.0},
      {{NT::V, 1000000000, 1000000000, -1, -1}, -2.0}};
  row.build(sparse, NT::V, 1000000000);
  require(!row.enabled() && row.capacity() == 0,
      "sparse alpha row allocated an n-sized bit plane");
  const std::vector<outer::Record> shifted{{{NT::V, 999999998, 1000000000, -1, -1}, -1e300},
      {{NT::V, 999999999, 1000000000, -1, -1}, -0.0},
      {{NT::V, 1000000000, 1000000000, -1, -1}, 7.0}};
  row.build(shifted, NT::V, 1000000000);
  require(row.enabled() && row.word_count() <= shifted.size(),
      "shifted narrow alpha row used absolute coordinate storage");
  for (std::size_t p = 0; p < shifted.size(); ++p)
    require(row.rank(shifted[p].key.i) == p, "shifted alpha rank lost its coordinate offset");
  const std::vector<outer::Record> descending{shifted.rbegin(), shifted.rend()};
  row.build(descending, NT::V, 1000000000);
  require(!row.enabled(), "out-of-order alpha row enabled rank indexing");
  row.build({}, NT::V, 1000000000);
  require(!row.enabled() && row.capacity() == 0, "empty alpha row retained metadata");
  return queries;
}

void rank_allocation_contract() {
  namespace fault = proposal_allocation_fault;
  const std::vector<outer::Record> records{{{NT::V, 1, 129, -1, -1}, -1e300},
      {{NT::V, 63, 129, -1, -1}, -0.0}, {{NT::V, 64, 129, -1, -1}, 7.0},
      {{NT::V, 65, 129, -1, -1}, 1.0}, {{NT::V, 128, 129, -1, -1}, 3.0}};
  for (const auto kind : {fault::Kind::BadAlloc, fault::Kind::LengthError}) {
    outer::PrimaryAlphaRankRow row;
    row.build(records, NT::V, 129);
    require(row.enabled(), "row allocation witness was not initially indexed");
    fault::arm(0, kind);
    row.build(records, NT::V, 129);
    const bool fired = fault::fired;
    fault::disable();
    require(fired && !row.enabled() && row.capacity() == 0,
        "failed row allocation retained an old/partial alpha index");
    row.build(records, NT::V, 129);
    require(row.rank(64) == 2 && row.rank(128) == 4,
        "row could not recover after optional allocation failure");
  }

  constexpr int n = 129;
  outer::Chart pilot(outer::zero);
  for (int i : {1, 63, 64, 65, 128, 129}) pilot.set({NT::V, i, n, -1, -1}, -double(i));
  for (int i : {63, 64, 65}) pilot.set({NT::WI, i, n - 1, -1, -1}, -double(i));
  // A sparse row, a mixed row and BE must not allocate rank-word arrays.
  pilot.set({NT::WMB, 1, n, -1, -1}, -1.0);
  pilot.set({NT::WMB, n, n, -1, -1}, -2.0);
  pilot.set({NT::WM, 3, 65, -1, -1}, 1.0);
  pilot.set({NT::WM, 3, 65, 7, 9}, 2.0);
  pilot.set({NT::BE, 2, 11, 3, 10}, 3.0);
  Scaffold ctx(n);
  std::size_t allocations = 0;
  {
    outer::Work work;
    fault::arm(std::numeric_limits<long>::max(), fault::Kind::BadAlloc);
    const outer::Proposals proposals(pilot, n, 128, ctx, work);
    allocations = fault::calls;
    fault::disable();
  }
  require(allocations > 3, "proposal allocation counter did not observe optional tail allocations");
  // Metadata allocation and the two eligible word arrays are exactly the last
  // three allocations; all authoritative pilot/guide allocations precede them.
  for (const auto kind : {fault::Kind::BadAlloc, fault::Kind::LengthError})
    for (std::size_t tail = 1; tail <= 3; ++tail) {
      outer::Work work;
      fault::arm(static_cast<long>(allocations - tail), kind);
      const outer::Proposals proposals(pilot, n, 128, ctx, work);
      const bool fired = fault::fired;
      fault::disable();
      require(fired, "proposal index allocation fault did not fire");
      for (const auto &[key, value] : pilot.values()) {
        (void)value;
        require(same(proposals.get_alpha(key), proposals.alpha.get(key)),
            "optional global/row rank allocation failure changed full-key lookup");
      }
      for (const auto nt : {NT::V, NT::WI, NT::WMB, NT::WM, NT::BE})
        for (int j : {65, 128, 129})
          for (int i = 0; i <= n + 1; ++i)
            for (const Key key : {Key{nt, i, j, -1, -1}, Key{nt, i, j, 7, 9}})
              require(same(proposals.get_alpha(key), proposals.alpha.get(key)),
                  "rank fallback differs for canonical/missing/secondary query");
    }
}
}  // namespace

int main() {
  try {
    const auto positions = proposal_position_contract();
    proposal_allocation_contract();
    const auto ranks = rank_rows();
    rank_allocation_contract();
    std::size_t queries = 0;
    for (std::size_t width : {1, 2, 22}) queries += check(width);
    std::cout << "scfg2_bounded_outer_alpha_rank=ok widths=3 queries=" << queries
              << " proposal_positions=" << positions << " rank_queries=" << ranks << '\n';
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
