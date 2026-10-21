#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/v_retained_agenda.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <new>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

// Fail one allocation only inside an explicitly armed operation. Fixture
// construction and observation continue to use ordinary malloc/free.
namespace allocation_fault {
thread_local long remaining = -1;
thread_local bool fired = false;
void arm(long successful_allocations) { remaining = successful_allocations; fired = false; }
void disable() { remaining = -1; }
void *allocate(std::size_t bytes) {
  if (remaining >= 0 && remaining-- == 0) {
    remaining = -1;
    fired = true;
    throw std::bad_alloc{};
  }
  if (void *memory = std::malloc(bytes ? bytes : 1)) return memory;
  throw std::bad_alloc{};
}
} // namespace allocation_fault
void *operator new(std::size_t bytes) { return allocation_fault::allocate(bytes); }
void *operator new[](std::size_t bytes) { return allocation_fault::allocate(bytes); }
void operator delete(void *memory) noexcept { std::free(memory); }
void operator delete[](void *memory) noexcept { std::free(memory); }
void operator delete(void *memory, std::size_t) noexcept { std::free(memory); }
void operator delete[](void *memory, std::size_t) noexcept { std::free(memory); }

namespace {
using namespace cparty::scfg2;
using Key = deductions::ItemKey;
using Rule = grammar::RuleId;
using Source = deductions::Deduction;
void require(bool ok, const char *message) {
  if (!ok) throw std::runtime_error(message);
}
bool same_source(const Source &a, const Source &b) {
  if (!deductions::ItemKeyEq{}(a.parent, b.parent) ||
      std::tie(a.rule, a.split.k, a.split.l, a.split.p, a.split.q, a.child_count) !=
      std::tie(b.rule, b.split.k, b.split.l, b.split.p, b.split.q, b.child_count)) return false;
  for (std::size_t c = 0; c < a.children.size(); ++c)
    if (!deductions::ItemKeyEq{}(a.children[c].item, b.children[c].item)) return false;
  return true;
}
template <class Agenda>
void absent_aliases(Agenda &agenda, const Key &parent) {
  const auto before = agenda.pending_parent_count();
  for (const auto &other : std::vector<Key>{
           {grammar::NonTerminal::W, parent.i, parent.j, -1, -1},
           {grammar::NonTerminal::BE, parent.i, parent.j, -1, -1},
           {grammar::NonTerminal::V, parent.i, parent.j, 0, -1},
           {grammar::NonTerminal::V, parent.i, parent.j, -1, 0},
           {grammar::NonTerminal::V, 0, parent.j, -1, -1},
           {grammar::NonTerminal::V, parent.j + 1, parent.j, -1, -1},
           {grammar::NonTerminal::V, parent.i, -1, -1, -1}})
    require(agenda.take_for(other).empty(), "packed V parent index accepted a different full key");
  require(agenda.pending_parent_count() == before, "absent V alias consumed a queued parent");
}

// Count the unchanged rectangular event-probe protocol in terms of wing
// lengths, independently of the production parent-coordinate traversal.
template <class Context>
std::size_t internal_probes(cand_pos_t j, const std::vector<cand_pos_t> &starts,
                            cand_pos_t n, const Context &ctx) {
  std::size_t total = 0;
  for (const auto i : starts) {
    if (i <= 1) continue;
    const auto max_left = std::min<cand_pos_t>({MAXLOOP, ctx.unpaired_prefix(i - 1), i - 2});
    for (cand_pos_t left = 0; left <= max_left; ++left)
      total += std::max<cand_pos_t>(0, std::min<cand_pos_t>(n - j, MAXLOOP + 1 - left));
  }
  return total;
}

// Independent full-parent enumeration supplies support. First child-arrival
// order supplies source order; no offset bitmap or forward window is reused.
void repeated_internal_events(bool numeric_pending = false) {
  std::string sequence = std::string(32, 'G') + std::string(32, 'C');
  std::string scaffold(sequence.size(), '.');
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalLogWeightContext ctx(owner);
  beam::VRetainedAgenda<replay::WFinalLogWeightContext> agenda;
  using Pair = std::pair<cand_pos_t, cand_pos_t>;
  std::map<Pair, std::size_t> arrival;
  const cand_pos_t n = sequence.size();
  if (numeric_pending)
    require(agenda.enable_numeric_pending(n), "repeated-event ring did not enable");
  std::size_t probes = 0;
  for (cand_pos_t j = 6; j < n; ++j) {
    std::vector<cand_pos_t> starts;
    for (cand_pos_t i = j - ctx.turn() - 1; i >= 1; --i) {
      if (!ctx.can_pair(i, j)) continue;
      arrival.emplace(Pair{i, j}, arrival.size());
      starts.push_back(i);
    }
    if (starts.empty()) continue;
    starts.push_back(starts.front());
    probes += internal_probes(j, starts, n, ctx);
    agenda.observe_committed_v_column(j, starts, n, ctx);
    std::reverse(starts.begin(), starts.end());
    probes += internal_probes(j, starts, n, ctx);
    agenda.observe_committed_v_column(j, starts, n, ctx);
  }
  std::size_t internals = 0, hairpins = 0, vm = 0;
  bool left_limit = false, right_limit = false, combined_limit = false;
  std::set<Pair> observed_wings;
  for (cand_pos_t j = 1; j <= n; ++j) {
    const std::vector<cand_pos_t> vm_starts{1, 1, 2};
    agenda.build_column(j, vm_starts, std::numeric_limits<std::size_t>::max(), ctx);
    std::vector<cand_pos_t> expected_starts;
    for (cand_pos_t i = j; i >= 1; --i) {
      const Key parent{grammar::NonTerminal::V, i, j, -1, -1};
      auto expected = deductions::generate_exact_basic_deductions(parent, Rule::V_INTERNAL, ctx);
      expected.erase(std::remove_if(expected.begin(), expected.end(), [&](const Source &source) {
        return !arrival.count({source.split.k, source.split.l});
      }), expected.end());
      std::sort(expected.begin(), expected.end(), [&](const Source &a, const Source &b) {
        return arrival.at({a.split.k, a.split.l}) < arrival.at({b.split.k, b.split.l});
      });
      for (const auto &source : expected) {
        const auto left = source.split.k - i - 1, right = j - source.split.l - 1;
        require(left + right <= MAXLOOP, "full generator exceeded V internal loop cap");
        left_limit |= left == MAXLOOP;
        right_limit |= right == MAXLOOP;
        combined_limit |= left > 0 && right > 0 && left + right == MAXLOOP;
      }
      internals += expected.size();
      auto hp = deductions::generate_exact_basic_deductions(parent, Rule::V_HAIRPIN, ctx);
      hairpins += hp.size();
      expected.insert(expected.end(), hp.begin(), hp.end());
      if (i == 1 || i == 2) {
        auto mv = deductions::generate_exact_basic_deductions(parent, Rule::V_VM, ctx);
        vm += mv.size();
        expected.insert(expected.end(), mv.begin(), mv.end());
      }
      if (!expected.empty()) expected_starts.push_back(i);
      absent_aliases(agenda, parent);
      const auto actual = agenda.take_for(parent);
      require(actual.size() == expected.size(), "V agenda source support changed under repeated arrivals");
      for (std::size_t k = 0; k < actual.size(); ++k) {
        const auto &source = actual[k];
        require(same_source(source, expected[k]),
                "V agenda changed complete source fields or first-arrival order");
        if (source.rule == Rule::V_INTERNAL) {
          observed_wings.emplace(source.split.k - i - 1, j - source.split.l - 1);
        } else {
          require(source.split.k == -1 && source.split.l == -1 &&
                      source.split.p == -1 && source.split.q == -1,
                  "V hairpin/VM acquired inactive split coordinates");
          for (std::size_t c = source.child_count; c < source.children.size(); ++c)
            require(deductions::ItemKeyEq{}(source.children[c].item, Key{}),
                    "V hairpin/VM acquired an inactive child");
        }
      }
      require(agenda.take_for(parent).empty(), "V source buffer was not consumed once");
    }
    require(agenda.parent_starts() == expected_starts, "V agenda parent order changed");
  }
  require(left_limit && right_limit && combined_limit && hairpins > 0 && vm > 0,
          "V loop cap corners or noninternal source fields were not exercised");
  // Cover every legal pair of wings, including all row endpoints and the
  // diagonal. The generator comparison above verifies support/order without
  // reproducing the queue's triangular indexing formula.
  for (cand_pos_t left = 0; left <= MAXLOOP; ++left)
    for (cand_pos_t right = 0; right <= MAXLOOP - left; ++right)
      require(observed_wings.count({left, right}) == 1,
              "V full-generator witness missed a legal wing pair");
  require(agenda.stats().internal_sources == internals &&
          agenda.stats().hairpin_sources == hairpins && agenda.stats().vm_sources == vm,
          "V agenda source counts changed under repeated arrivals");
  require(agenda.stats().internal_parent_probes == probes &&
              agenda.stats().hairpin_parent_probes == std::size_t(n) * (n + 1) / 2 &&
              agenda.stats().vm_parent_probes == 3 * std::size_t(n),
          "V agenda changed the event-probe protocol");
  std::cout << "repeated_v_events=" << internals << " loop_cap_corners=ok\n";
}

// Absolute split coordinates must survive beyond the byte range even though
// the retained queue can store bounded wing lengths in bytes.
void translated_internal_sources() {
  const cand_pos_t k = 300 + MAXLOOP + 2;
  const std::string sequence = std::string(k, 'G') + std::string(MAXLOOP + 16, 'C');
  weights::RuntimeSnapshot owner(sequence, std::string(sequence.size(), '.'));
  replay::WFinalLogWeightContext ctx(owner);
  const cand_pos_t l = k + ctx.turn() + 2;
  const Key child{grammar::NonTerminal::V, k, l, -1, -1};
  beam::VRetainedAgenda<replay::WFinalLogWeightContext> agenda;
  agenda.observe_committed_v_column(l, {k, k}, sequence.size(), ctx);
  using Pair = std::pair<cand_pos_t, cand_pos_t>;
  for (const auto &[left, right] : std::vector<Pair>{
           {0, 0}, {MAXLOOP, 0}, {0, MAXLOOP}, {1, MAXLOOP - 1},
           {MAXLOOP - 1, 1}, {MAXLOOP / 2, MAXLOOP - MAXLOOP / 2}}) {
    const Key parent{grammar::NonTerminal::V, k - left - 1, l + right + 1, -1, -1};
    auto expected = deductions::generate_exact_basic_deductions(parent, Rule::V_INTERNAL, ctx);
    expected.erase(std::remove_if(expected.begin(), expected.end(), [&](const Source &source) {
      return source.child_count != 1 ||
          !deductions::ItemKeyEq{}(source.children[0].item, child);
    }), expected.end());
    const auto actual = agenda.take_for(parent);
    require(expected.size() == 1 && actual.size() == 1 &&
                same_source(actual.front(), expected.front()) &&
                actual.front().split.k == k && actual.front().split.l == l,
            "translated V wings changed absolute split/child fields");
    require(agenda.take_for(parent).empty(), "translated V duplicate survived consumption");
  }
}

void skipped_parent_lifetime() {
  std::string sequence;
  for (int k = 0; k < 30; ++k) sequence += "AGCU";
  std::string scaffold(sequence.size(), '.');
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalLogWeightContext ctx(owner);
  beam::VRetainedAgenda<replay::WFinalLogWeightContext> history, stream;
  const cand_pos_t n = sequence.size();
  std::size_t skipped = 0, compared = 0;
  for (cand_pos_t j = 1; j <= n; ++j) {
    history.build_column(j, {}, 30, ctx);
    stream.build_column(j, {}, 30, ctx);
    require(history.parent_starts() == stream.parent_starts(),
            "V retirement changed target parent order");
    for (const auto i : history.parent_starts()) {
      const Key parent{grammar::NonTerminal::V, i, j, -1, -1};
      if (!deductions::exact_basic_item_parent_active(parent, ctx)) {
        ++skipped;
        continue; // Production also omits take_for() for these parents.
      }
      const auto a = history.take_for(parent), b = stream.take_for(parent);
      require(a.size() == b.size(), "V retirement changed future source support");
      for (std::size_t k = 0; k < a.size(); ++k) {
        require(same_source(a[k], b[k]), "V retirement changed complete future source fields/order");
        ++compared;
      }
    }
    std::vector<cand_pos_t> starts;
    for (cand_pos_t i = 1; i < j - ctx.turn() && starts.size() < 4; ++i)
      if (ctx.can_pair(i, j)) starts.push_back(i);
    history.observe_committed_v_column(j, starts, n, ctx);
    stream.observe_committed_v_column(j, starts, n, ctx);
    stream.release_column();
    stream.release_column();
    require(stream.pending_column_count() <= MAXLOOP + 1,
            "V queue exceeded its fixed future-target window");
    require(stream.stats().internal_sources == history.stats().internal_sources &&
                stream.stats().internal_parent_probes == history.stats().internal_parent_probes,
            "V retirement changed cumulative work");
  }
  require(skipped > 0 && compared > 0 && history.pending_parent_count() > 0 &&
              stream.pending_parent_count() == 0,
          "V skipped-parent lifetime witness did not exercise expired candidates");
}

void nonempty_scaffold_streams(bool numeric_pending = false) {
  const std::string sequence = "GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU";
  const std::string branching = "....(...)((....)..)((.(.((...))))(.(....)((...))))";
  auto forbidden = branching;
  forbidden[0] = 'x';
  forbidden[forbidden.find('.', forbidden.size() / 2)] = 'x';
  forbidden[forbidden.rfind('.')] = 'x';
  std::size_t compared = 0, internals = 0, cap_exceptions = 0;
  for (const auto &scaffold : {branching, forbidden}) {
    weights::RuntimeSnapshot owner(sequence, scaffold);
    for (const auto policy : {WmbClosurePolicy::RequireWeakClosure,
                              WmbClosurePolicy::LegacyCompatible}) {
      replay::WFinalLogWeightContext ctx(owner, policy);
      std::vector<Source> owned, owned_expected;
      {
        beam::VRetainedAgenda<replay::WFinalLogWeightContext> agenda;
        using Pair = std::pair<cand_pos_t, cand_pos_t>;
        std::map<Pair, std::size_t> arrival;
        const cand_pos_t n = sequence.size();
        if (numeric_pending)
          require(agenda.enable_numeric_pending(n), "nonempty-scaffold ring did not enable");
        const auto metadata_size = agenda.memory_snapshot("v").vector_size_bytes;
        constexpr cand_pos_t cap = 3;
        std::size_t probes = 0, hairpin_probes = 0, expected_internals = 0;
        std::size_t expected_hairpins = 0, expected_vm = 0;
        for (cand_pos_t j = 1; j <= n; ++j) {
          const std::vector<cand_pos_t> vm_starts{1, 1, 2};
          agenda.build_column(j, vm_starts, cap, ctx);
          const auto min_hairpin = std::max<cand_pos_t>(1, j - cap - 1);
          hairpin_probes += j - min_hairpin + 1;
          const auto scaffold_i = ctx.pair_partner(j);
          if (scaffold_i > 0 && scaffold_i < min_hairpin) ++hairpin_probes;
          std::vector<cand_pos_t> expected_starts, retained;
          for (cand_pos_t i = j; i >= 1; --i) {
            const Key parent{grammar::NonTerminal::V, i, j, -1, -1};
            auto expected = deductions::generate_exact_basic_deductions(parent, Rule::V_INTERNAL, ctx);
            expected.erase(std::remove_if(expected.begin(), expected.end(), [&](const Source &source) {
              return !arrival.count({source.split.k, source.split.l});
            }), expected.end());
            std::sort(expected.begin(), expected.end(), [&](const Source &a, const Source &b) {
              return arrival.at({a.split.k, a.split.l}) < arrival.at({b.split.k, b.split.l});
            });
            expected_internals += expected.size();
            if (i >= min_hairpin || i == scaffold_i) {
              auto hp = deductions::generate_exact_basic_deductions(parent, Rule::V_HAIRPIN, ctx);
              expected_hairpins += hp.size();
              cap_exceptions += !hp.empty() && i < min_hairpin;
              expected.insert(expected.end(), hp.begin(), hp.end());
            }
            if (i == 1 || i == 2) {
              auto vm = deductions::generate_exact_basic_deductions(parent, Rule::V_VM, ctx);
              expected_vm += vm.size();
              expected.insert(expected.end(), vm.begin(), vm.end());
            }
            if (!expected.empty()) {
              expected_starts.push_back(i);
              absent_aliases(agenda, parent);
              if (j - i - 1 >= ctx.turn() && ctx.can_pair(i, j) && retained.size() < 4)
                retained.push_back(i);
            }
            const auto actual = agenda.take_for(parent);
            require(actual.size() == expected.size(), "nonempty V source support changed");
            for (std::size_t c = 0; c < actual.size(); ++c) {
              require(same_source(actual[c], expected[c]),
                      "nonempty V complete source fields/order differ from full generator");
              ++compared;
            }
            if (owned.empty() && !actual.empty()) { owned = actual; owned_expected = expected; }
          }
          require(agenda.parent_starts() == expected_starts,
                  "nonempty V parent frontier order changed");
          if (j % 2) std::reverse(retained.begin(), retained.end());
          if (!retained.empty()) retained.push_back(retained.front());
          for (const auto i : retained) arrival.emplace(Pair{i, j}, arrival.size());
          probes += internal_probes(j, retained, n, ctx);
          agenda.observe_committed_v_column(j, retained, n, ctx);
          std::reverse(retained.begin(), retained.end());
          probes += internal_probes(j, retained, n, ctx);
          agenda.observe_committed_v_column(j, retained, n, ctx);
          agenda.release_column();
          agenda.release_column();
          require(agenda.pending_column_count() <= MAXLOOP + 1,
                  "nonempty compact V queue exceeded its future-column window");
        }
        require(agenda.stats().internal_sources == expected_internals &&
                    agenda.stats().hairpin_sources == expected_hairpins &&
                    agenda.stats().vm_sources == expected_vm &&
                    agenda.stats().internal_parent_probes == probes &&
                    agenda.stats().hairpin_parent_probes == hairpin_probes &&
                    agenda.stats().vm_parent_probes == 3 * std::size_t(n),
                "nonempty compact V work counts changed");
        require(agenda.pending_parent_count() == 0 && agenda.pending_column_count() == 0 &&
                    agenda.memory_snapshot("v").vector_size_bytes == metadata_size,
                "nonempty compact V sources survived release");
        internals += expected_internals;
      }
      require(!owned.empty() && owned.size() == owned_expected.size(),
              "V full-buffer ownership witness is vacuous");
      for (std::size_t c = 0; c < owned.size(); ++c)
        require(same_source(owned[c], owned_expected[c]),
                "returned V source buffer depended on the destroyed agenda");
    }
  }
  require(compared > 0 && internals > 0 && cap_exceptions > 0,
          "nonempty V source/internal/scaffold-cap coverage is vacuous");
  std::cout << "nonempty_v_sources=" << compared << " internals=" << internals << '\n';
}

using NumericVAgenda = beam::VRetainedAgenda<replay::WFinalLogWeightContext>;
void same_pending_state(const NumericVAgenda &a, const NumericVAgenda &b) {
  const auto &x = a.stats();
  const auto &y = b.stats();
  require(std::tie(x.hairpin_parent_probes, x.hairpin_sources, x.internal_parent_probes,
                   x.internal_sources, x.vm_parent_probes, x.vm_sources) ==
              std::tie(y.hairpin_parent_probes, y.hairpin_sources, y.internal_parent_probes,
                       y.internal_sources, y.vm_parent_probes, y.vm_sources) &&
              a.column() == b.column() && a.parent_starts() == b.parent_starts() &&
              a.pending_parent_count() == b.pending_parent_count() &&
              a.pending_column_count() == b.pending_column_count(),
          "numeric V queue changed counters, parent order or pending lifetime");
}
void drain_pending_column(NumericVAgenda &a, NumericVAgenda &b, cand_pos_t j) {
  for (cand_pos_t i = j; i >= 1; --i) {
    const Key key{grammar::NonTerminal::V, i, j, -1, -1};
    const auto actual = a.take_for(key), expected = b.take_for(key);
    require(actual.size() == expected.size(), "numeric V queue changed source support");
    for (std::size_t k = 0; k < actual.size(); ++k)
      require(same_source(actual[k], expected[k]),
              "numeric V queue changed complete source fields or arrival order");
    same_pending_state(a, b);
  }
}

// A custom context may observe queries even when it delegates numeric answers
// to the standard owner. Numeric pending storage must not opt it into caching.
struct GenericVQueries final : deductions::PairabilityContext {
  const deductions::PairabilityContext &base;
  mutable std::vector<std::tuple<char, cand_pos_t, cand_pos_t>> trace;
  explicit GenericVQueries(const deductions::PairabilityContext &ctx) : base(ctx) {}
  cand_pos_t turn() const override { trace.emplace_back('t', 0, 0); return base.turn(); }
  cand_pos_t sequence_length() const override { trace.emplace_back('n', 0, 0); return base.sequence_length(); }
  cand_pos_t pair_partner(cand_pos_t p) const override { trace.emplace_back('p', p, 0); return base.pair_partner(p); }
  cand_pos_t unpaired_prefix(cand_pos_t p) const override { trace.emplace_back('u', p, 0); return base.unpaired_prefix(p); }
  bool can_pair(cand_pos_t i, cand_pos_t j) const override {
    trace.emplace_back('c', i, j); return base.can_pair(i, j);
  }
};
void generic_pending_query_protocol(const deductions::PairabilityContext &base) {
  GenericVQueries a(base), b(base);
  beam::VRetainedAgenda<GenericVQueries> plain, ring;
  const auto n = base.sequence_length();
  require(ring.enable_numeric_pending(n), "generic numeric V fixture did not enable");
  plain.observe_committed_v_column(8, {3, 3}, n, a);
  ring.observe_committed_v_column(8, {3, 3}, n, b);
  require(a.trace == b.trace && plain.pending_parent_count() > 0 &&
              plain.pending_parent_count() == ring.pending_parent_count() &&
              plain.stats().internal_parent_probes == ring.stats().internal_parent_probes &&
              plain.stats().internal_sources == ring.stats().internal_sources,
          "numeric V storage changed generic-context queries or sources");
  for (cand_pos_t j = 1; j <= n; ++j) {
    for (cand_pos_t i = 1; i <= j; ++i) {
      const Key parent{grammar::NonTerminal::V, i, j, -1, -1};
      const auto expected = plain.take_for(parent), actual = ring.take_for(parent);
      require(actual.size() == expected.size(), "generic numeric V support differs");
      for (std::size_t c = 0; c < actual.size(); ++c)
        require(same_source(actual[c], expected[c]), "generic numeric V fields/order differ");
    }
  }
  require(!plain.pending_parent_count() && !ring.pending_parent_count(),
          "generic numeric V drain retained an owner");
}

// The ordinary packed-map path is a second owner for the same shared sources.
// Keep live colliding columns, mix out-of-window events back into the current
// window, then compare every public source field and lifetime observation.
void numeric_pending_protocol() {
  using Context = replay::WFinalLogWeightContext;
  using Agenda = beam::VRetainedAgenda<Context>;
  using NT = grammar::NonTerminal;
  std::string sequence;
  for (int i = 0; i < 128; ++i) sequence += i % 2 ? 'C' : 'G';
  weights::RuntimeSnapshot owner(sequence, std::string(sequence.size(), '.'));
  Context ctx(owner);
  generic_pending_query_protocol(ctx);
  Agenda ring, packed;
  require(!ring.enable_numeric_pending(0) && !ring.enable_numeric_pending(-1),
          "numeric V queue accepted an invalid sequence domain");
  require(ring.enable_numeric_pending(sequence.size()), "numeric V queue did not enable");
  const auto metadata_capacity = ring.memory_snapshot("v").vector_capacity_bytes;
  const auto state = same_pending_state;
  const auto drain = drain_pending_column;
  const auto build = [&](cand_pos_t j) {
    ring.build_column(j, {1, 2}, 30, ctx);
    packed.build_column(j, {1, 2}, 30, ctx);
    state(ring, packed);
  };
  const auto release = [&] {
    ring.release_column(); packed.release_column();
    ring.release_column(); packed.release_column();
    state(ring, packed);
  };
  build(8);
  require(ring.pending_parent_count() > 0, "numeric V collision witness is vacuous");
  require(!ring.enable_numeric_pending(sequence.size()), "live numeric V queue was reset");
  // End 40 occupies end 8's ring slot while those earlier owners remain live.
  build(40);
  drain(ring, packed, 8);
  drain(ring, packed, 40);
  release();
  build(72); // Reuse the now-empty slot and its touched-parent marks.
  absent_aliases(ring, {NT::V, 1, 72, -1, -1});
  absent_aliases(packed, {NT::V, 1, 72, -1, -1});
  drain(ring, packed, 72);
  release();
  build(8);
  ring.observe_committed_v_column(100, {41, 41}, 128, ctx);
  packed.observe_committed_v_column(100, {41, 41}, 128, ctx);
  state(ring, packed);
  build(101); // Existing fallback parents must not acquire a second ring owner.
  for (cand_pos_t j = 1; j <= 128; ++j) drain(ring, packed, j);
  release();
  require(ring.pending_parent_count() == 0, "mixed numeric/map V owner survived consumption");
  require(ring.memory_snapshot("v").vector_capacity_bytes -
                  packed.memory_snapshot("v").vector_capacity_bytes == metadata_capacity,
          "numeric V queue retained historical source-vector capacity");

  // Copying must preserve independently consumable owners; normal moves must
  // transfer them without leaving stale live counters in their source.
  build(48);
  auto copied = ring;
  auto expected = packed;
  drain(copied, expected, 48);
  state(ring, packed);
  copied = ring;
  expected = packed;
  copied = copied;
  auto moved = std::move(copied);
  require(copied.pending_parent_count() == 0, "moved V agenda retained pending owners");
  copied = std::move(moved);
  require(moved.pending_parent_count() == 0, "move-assigned V agenda retained pending owners");
  drain(copied, expected, 48);
  copied.release_column(); expected.release_column();
  state(copied, expected);
  // Copy assignment may also change which of the two representations owns it.
  Agenda assigned;
  assigned = ring;
  expected = packed;
  drain(assigned, expected, 48);
  copied = packed;
  expected = packed;
  drain(copied, expected, 48);
  drain(ring, packed, 48);
  release();

  // A caller-supplied smaller optional domain retains packed-map fallback.
  Agenda short_ring, short_packed;
  require(short_ring.enable_numeric_pending(16), "short numeric V domain did not enable");
  short_ring.build_column(48, {1, 2}, 30, ctx);
  short_packed.build_column(48, {1, 2}, 30, ctx);
  drain(short_ring, short_packed, 48);
  short_ring.release_column(); short_packed.release_column();
  state(short_ring, short_packed);
}

void numeric_pending_allocation_failure() {
  std::string sequence;
  for (int i = 0; i < 128; ++i) sequence += i % 2 ? 'C' : 'G';
  weights::RuntimeSnapshot owner(sequence, std::string(sequence.size(), '.'));
  replay::WFinalLogWeightContext ctx(owner);
  const auto fill = [&](NumericVAgenda &agenda, cand_pos_t j) {
    agenda.build_column(j, {1, 2}, 30, ctx);
    agenda.observe_committed_v_column(j, {3, 3}, 128, ctx);
  };
  // Packed references are complete independent owners constructed before any
  // faults; metadata, source vectors, and future parent lists must agree.
  NumericVAgenda source, destination, source_map, destination_map;
  require(source.enable_numeric_pending(128) && destination.enable_numeric_pending(64),
          "copy-fault numeric V fixture did not enable");
  fill(source, 48); fill(source_map, 48);
  fill(destination, 8); fill(destination_map, 8);
  require(source.pending_parent_count() && destination.pending_parent_count() &&
              source.pending_column_count() && destination.pending_column_count(),
          "V copy-fault current/future owners were not exercised");
  const auto compare = [](NumericVAgenda actual, NumericVAgenda expected) {
    same_pending_state(actual, expected);
    for (cand_pos_t j = 1; j <= 128; ++j) drain_pending_column(actual, expected, j);
    require(actual.pending_parent_count() == 0 && expected.pending_parent_count() == 0,
            "copy-fault numeric V drain left an owner");
  };
  compare(source, source_map); compare(destination, destination_map);
  bool success = false;
  long failures = 0;
  for (long site = 0; !success && site < 2048; ++site) {
    auto attempted = destination;
    bool threw = false;
    allocation_fault::arm(site);
    try { attempted = source; } catch (const std::bad_alloc &) { threw = true; }
    const bool fired = allocation_fault::fired;
    allocation_fault::disable();
    require(threw == fired, "copy-fault V assignment swallowed an allocation failure");
    compare(attempted, threw ? destination_map : source_map);
    compare(source, source_map);
    failures += threw;
    success = !fired;
  }
  require(success && failures >= 3, "V copy fault sweep did not reach all allocations");
  // Failure of any of the three optional metadata allocations must leave an
  // ordinary, usable map owner, with no partial numeric routing state.
  for (long site = 0; site < 3; ++site) {
    NumericVAgenda fallback;
    allocation_fault::arm(site);
    const bool enabled = fallback.enable_numeric_pending(128);
    const bool fired = allocation_fault::fired;
    allocation_fault::disable();
    require(fired && !enabled, "numeric V metadata failure did not select map fallback");
    fill(fallback, 48);
    compare(fallback, source_map);
  }
}

} // namespace

int main() {
  using namespace cparty::scfg2;
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  skipped_parent_lifetime();
  nonempty_scaffold_streams();
  nonempty_scaffold_streams(true);
  numeric_pending_protocol();
  numeric_pending_allocation_failure();
  std::string sequence;
  for (int pos = 0; pos < 64; ++pos) {
    sequence.push_back((pos % 2) == 0 ? 'G' : 'C');
  }
  std::string scaffold(sequence.size(), '.');
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const cand_pos_t j = static_cast<cand_pos_t>(sequence.size());
  const std::vector<cand_pos_t> no_vm_starts;

  beam::VRetainedAgenda<replay::WFinalExactInsideContext> exact;
  exact.build_column(
      j,
      no_vm_starts,
      std::numeric_limits<std::size_t>::max(),
      ctx);
  assert(exact.stats().hairpin_parent_probes == sequence.size());

  beam::VRetainedAgenda<replay::WFinalExactInsideContext> capped;
  capped.build_column(j, no_vm_starts, 30, ctx);
  assert(capped.stats().hairpin_parent_probes == 32);
  assert(exact.parent_starts().size() > capped.parent_starts().size());
  for (const cand_pos_t i : capped.parent_starts()) {
    assert(j - i - 1 <= 30);
    const deductions::ItemKey parent{
        grammar::NonTerminal::V, i, j, -1, -1};
    for (const auto &source : capped.take_for(parent)) {
      assert(source.rule == grammar::RuleId::V_HAIRPIN);
    }
  }

  std::cout << "exact_hairpin_probes="
            << exact.stats().hairpin_parent_probes
            << " capped_hairpin_probes="
            << capped.stats().hairpin_parent_probes << '\n';
  repeated_internal_events();
  repeated_internal_events(true);
  translated_internal_sources();
  return EXIT_SUCCESS;
}
