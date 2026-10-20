#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/beam/v_retained_agenda.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

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
void repeated_internal_events() {
  std::string sequence = std::string(32, 'G') + std::string(32, 'C');
  std::string scaffold(sequence.size(), '.');
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  replay::WFinalLogWeightContext ctx(owner);
  beam::VRetainedAgenda<replay::WFinalLogWeightContext> agenda;
  using Pair = std::pair<cand_pos_t, cand_pos_t>;
  std::map<Pair, std::size_t> arrival;
  const cand_pos_t n = sequence.size();
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
      for (std::size_t k = 0; k < actual.size(); ++k)
        require(same_source(actual[k], expected[k]),
                "V agenda changed complete source fields or first-arrival order");
      require(agenda.take_for(parent).empty(), "V source buffer was not consumed once");
    }
    require(agenda.parent_starts() == expected_starts, "V agenda parent order changed");
  }
  require(left_limit && right_limit && combined_limit, "V loop cap corners were not exercised");
  require(agenda.stats().internal_sources == internals &&
          agenda.stats().hairpin_sources == hairpins && agenda.stats().vm_sources == vm,
          "V agenda source counts changed under repeated arrivals");
  require(agenda.stats().internal_parent_probes == probes &&
              agenda.stats().hairpin_parent_probes == std::size_t(n) * (n + 1) / 2 &&
              agenda.stats().vm_parent_probes == 3 * std::size_t(n),
          "V agenda changed the event-probe protocol");
  std::cout << "repeated_v_events=" << internals << " loop_cap_corners=ok\n";
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

void nonempty_scaffold_streams() {
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
                    agenda.memory_snapshot("v").vector_size_bytes == 0,
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
} // namespace

int main() {
  using namespace cparty::scfg2;
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  skipped_parent_lifetime();
  nonempty_scaffold_streams();
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
  return EXIT_SUCCESS;
}
