#ifndef CPARTY_SCFG2_PROFILES_PROFILE_ASSEMBLY_HH
#define CPARTY_SCFG2_PROFILES_PROFILE_ASSEMBLY_HH

#include "scfg2/profiles/loop_context.hh"

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <unordered_map>
#include <vector>

namespace cparty::scfg2::profiles {

using LoopAccessibilityRow = std::array<double, unpaired_loop_context_count>;

// Structural-context profile (PK-CapR) and window accessibility (PK-Raccess)
// over the conditional ensemble of G' given the scaffold G.
struct ProfileResult {
  bool ok = false;
  std::string error;
  double log_partition_weight = -std::numeric_limits<double>::infinity();
  double ensemble_free_energy = std::numeric_limits<double>::quiet_NaN();
  // RT in kcal/mol, so that ED = -RT log P.
  double rt_kcal_mol = std::numeric_limits<double>::quiet_NaN();
  // context[pos][c], one-based pos; index zero unused.
  std::vector<std::array<double, loop_context_count>> context;
  // accessibility[w_idx][s] = P(bases s..s+w-1 all unpaired), one-based s.
  std::vector<cand_pos_t> window_lengths;
  std::vector<std::vector<double>> accessibility;
  // Optional [window][one-based start][unpaired context], in enum order
  // excluding Stem: Hairpin, Bulge, Internal, Multi, Exterior, PseudoLoop.
  std::vector<std::vector<LoopAccessibilityRow>> accessibility_by_loop;
  // Retained backend counts public/private reverse cells. It saves no terms
  // or factors; its hash-allocation byte estimate is unavailable (zero).
  // Exact/tape consumers leave these diagnostics zero.
  std::size_t reverse_records = 0;
  std::size_t reverse_terms = 0;
  std::size_t reverse_factors = 0;
  std::size_t reverse_storage_bytes = 0;
  // One per chain step examined per requested valid window length.
  std::size_t chain_window_steps = 0;
};

// Posterior mass collected from one inside/outside pass.
//
// Whole runs: a deduction (or explicit beam term) that emits whole unpaired
// runs, with its posterior.
//
// Chain steps: unpaired runs grown one base at a time along a chain keyed by
// (family, i).  Step t emits base t.  For a window [s, e] inside one chain,
//   P = entry(s) * prod_{t=s+1..e} step(t) * outside(e) / Z,
// where entry(s) is the mass of all chain derivations that reach base s,
// step(t) the weight of one more base, and outside(e)/Z = adjoint(e) /
// value(e) for the chain state after base e.  `has_prev` is false at a chain
// start or where the chain was cut (beam pruning), which ends the window.
// A mirrored chain grows leftwards (the beam's start prefixes): its steps are
// stored at -t, so the same walk applies and window starts map back.
class ProfileAccumulator {
 public:
  struct ChainStep {
    double log_entry = -std::numeric_limits<double>::infinity();
    double log_step = 0.0;
    bool has_prev = false;
    double log_state = -std::numeric_limits<double>::infinity();
    double state_adjoint = 0.0;
    // Posterior that this chain emits base t.
    double emit_probability = 0.0;
  };

  explicit ProfileAccumulator(cand_pos_t n, bool by_loop = false)
      : n_(n), by_loop_(by_loop) {}

  // Online range accumulation: keep O(n * number_of_windows) output-sized
  // differences instead of one record per regenerated deduction.
  ProfileAccumulator(cand_pos_t n, const std::vector<cand_pos_t> &windows,
                     bool by_loop = false)
      : n_(n), streaming_(true), by_loop_(by_loop), windows_(windows) {
    for (auto &d : context_diff_) d.assign(static_cast<std::size_t>(n) + 2, 0.0);
    window_diff_.assign(windows.size(), std::vector<double>(static_cast<std::size_t>(n) + 2));
    if (by_loop_)
      loop_window_diff_.assign(windows.size(),
          std::vector<LoopAccessibilityRow>(static_cast<std::size_t>(n) + 2));
  }

  void add_whole_runs(const UnpairedEmission &emission, double probability) {
    if (probability == 0.0) return;
    if (streaming_) {
      for (std::size_t r = 0; r < emission.count; ++r) {
        const auto &run = emission.runs[r];
        if (run.lo > run.hi) continue;
        auto &d = context_diff_[static_cast<std::size_t>(run.context)];
        d[run.lo] += probability;
        d[run.hi + 1] -= probability;
        for (std::size_t w = 0; w < windows_.size(); ++w) {
          if (windows_[w] < 1 || windows_[w] > run.hi - run.lo + 1) continue;
          window_diff_[w][run.lo] += probability;
          window_diff_[w][run.hi - windows_[w] + 2] -= probability;
          if (by_loop_) {
            const auto c = unpaired_index(run.context);
            loop_window_diff_[w][run.lo][c] += probability;
            loop_window_diff_[w][run.hi - windows_[w] + 2][c] -= probability;
          }
        }
      }
      return;
    }
    whole_.push_back({emission, probability});
  }

  std::size_t stored_whole_runs() const { return whole_.size(); }

  ChainStep &chain_step(int family, cand_pos_t i, cand_pos_t t,
                        LoopContext context, bool mirrored = false) {
    auto &chain = chains_[{family, i}];
    chain.context = context;
    chain.mirrored = mirrored;
    const auto pos = mirrored ? -t : t;
    // Both forward and regenerated reverse visit a chain monotonically.
    // Correct begin/end hints make these insertions amortized O(1), without
    // a log(sequence length) lookup per step. Arbitrary callers still work.
    if (chain.steps.empty() || pos < chain.steps.begin()->first)
      return chain.steps.emplace_hint(chain.steps.begin(), pos, ChainStep{})->second;
    if (pos > chain.steps.rbegin()->first)
      return chain.steps.emplace_hint(chain.steps.end(), pos, ChainStep{})->second;
    return chain.steps[pos];
  }

  ProfileResult finish(double log_partition_weight,
                       const std::vector<cand_pos_t> &window_lengths) const {
    ProfileResult out;
    if (streaming_ && window_lengths != windows_) {
      out.error = "streamed profile windows differ from initialized windows";
      return out;
    }
    out.log_partition_weight = log_partition_weight;
    const std::size_t size = static_cast<std::size_t>(n_) + 2;

    // Streaming accumulators already own these difference arrays. Read them
    // directly so finish stays repeatable without allocating scratch copies.
    std::array<std::vector<double>, loop_context_count> diff;
    if (!streaming_) for (auto &d : diff) d.assign(size, 0.0);
    for (const auto &[emission, probability] : whole_) {
      for (std::size_t r = 0; r < emission.count; ++r) {
        const auto &run = emission.runs[r];
        if (run.lo > run.hi) continue;
        auto &d = diff[static_cast<std::size_t>(run.context)];
        d[static_cast<std::size_t>(run.lo)] += probability;
        d[static_cast<std::size_t>(run.hi) + 1] -= probability;
      }
    }
    const auto &context_diff = streaming_ ? context_diff_ : diff;
    out.context.assign(static_cast<std::size_t>(n_) + 1, {});
    for (std::size_t c = 0; c < loop_context_count; ++c) {
      double running = 0.0;
      for (cand_pos_t pos = 1; pos <= n_; ++pos) {
        running += context_diff[c][static_cast<std::size_t>(pos)];
        out.context[static_cast<std::size_t>(pos)][c] = running;
      }
    }
    for (const auto &[key, chain] : chains_) {
      for (const auto &[key_t, step] : chain.steps) {
        const cand_pos_t t = chain.mirrored ? -key_t : key_t;
        out.context[static_cast<std::size_t>(t)]
                   [static_cast<std::size_t>(chain.context)] +=
            step.emit_probability;
      }
    }
    for (cand_pos_t pos = 1; pos <= n_; ++pos) {
      auto &row = out.context[static_cast<std::size_t>(pos)];
      double unpaired = 0.0;
      for (std::size_t c = 1; c < loop_context_count; ++c) unpaired += row[c];
      row[static_cast<std::size_t>(LoopContext::Stem)] = 1.0 - unpaired;
    }

    out.window_lengths = window_lengths;
    for (std::size_t wi = 0; wi < window_lengths.size(); ++wi) {
      const cand_pos_t w = window_lengths[wi];
      std::vector<double> starts;
      if (!streaming_) starts.assign(size, 0.0);
      std::vector<double> row(static_cast<std::size_t>(n_) + 1, 0.0);
      std::vector<LoopAccessibilityRow> loop_starts, loop_row;
      if (by_loop_) {
        if (!streaming_) loop_starts.resize(size);
        loop_row.resize(static_cast<std::size_t>(n_) + 1);
      }
      if (w >= 1 && w <= n_) {
        for (const auto &[emission, probability] : whole_) {
          for (std::size_t r = 0; r < emission.count; ++r) {
            const auto &run = emission.runs[r];
            if (run.hi - run.lo + 1 < w) continue;
            starts[static_cast<std::size_t>(run.lo)] += probability;
            starts[static_cast<std::size_t>(run.hi - w + 2)] -= probability;
            if (by_loop_) {
              const auto c = unpaired_index(run.context);
              loop_starts[run.lo][c] += probability;
              loop_starts[run.hi - w + 2][c] -= probability;
            }
          }
        }
        const auto &window_diff = streaming_ ? window_diff_[wi] : starts;
        const auto *loop_window_diff = by_loop_
            ? &(streaming_ ? loop_window_diff_[wi] : loop_starts) : nullptr;
        double running = 0.0;
        LoopAccessibilityRow loop_running{};
        for (cand_pos_t s = 1; s + w - 1 <= n_; ++s) {
          running += window_diff[static_cast<std::size_t>(s)];
          row[static_cast<std::size_t>(s)] = running;
          if (by_loop_)
            for (std::size_t c = 0; c < unpaired_loop_context_count; ++c) {
              loop_running[c] += (*loop_window_diff)[s][c];
              loop_row[s][c] = loop_running[c];
            }
        }
        for (const auto &[key, chain] : chains_)
          add_chain_windows(chain, w, &row, by_loop_ ? &loop_row : nullptr,
                            &out.chain_window_steps);
      }
      out.accessibility.push_back(std::move(row));
      if (by_loop_) out.accessibility_by_loop.push_back(std::move(loop_row));
    }
    out.ok = true;
    return out;
  }

 private:
  static std::size_t unpaired_index(LoopContext context) {
    const auto c = static_cast<std::size_t>(context);
    if (c == 0 || c >= loop_context_count)
      throw std::logic_error("unpaired run has no unpaired loop context");
    return c - 1;
  }

  struct Chain {
    LoopContext context = LoopContext::Exterior;
    bool mirrored = false;
    std::map<cand_pos_t, ChainStep> steps;
  };

  // Prefix log products evaluate each window in O(1) after O(E) preparation,
  // including when w grows with sequence length. Cuts and zero transitions
  // start a new stretch, so windows never bridge a missing chain link.
  static void add_chain_windows(const Chain &chain,
                                cand_pos_t w,
                                std::vector<double> *row,
                                std::vector<LoopAccessibilityRow> *loop_row,
                                std::size_t *steps_examined) {
    const auto c = loop_row ? unpaired_index(chain.context) : 0;
    struct Prefix {
      cand_pos_t position;
      const ChainStep *step;
      long double log_product;
    };
    std::vector<Prefix> prefix;
    prefix.reserve(chain.steps.size());
    std::size_t stretch_start = 0;
    long double log_product = 0.0L;
    for (const auto &[t, step] : chain.steps) {
      ++*steps_examined;
      const bool cut = prefix.empty() || t != prefix.back().position + 1 ||
                       !step.has_prev || step.log_step == -std::numeric_limits<double>::infinity();
      if (cut) {
        stretch_start = prefix.size();
        log_product = 0.0L;
      } else {
        log_product += static_cast<long double>(step.log_step);
      }
      prefix.push_back(Prefix{t, &step, log_product});
      const std::size_t end = prefix.size() - 1;
      if (end - stretch_start + 1 < static_cast<std::size_t>(w)) continue;
      const auto &start = prefix[end - static_cast<std::size_t>(w) + 1];
      const ChainStep &first = *start.step;
      if (!std::isfinite(first.log_entry)) continue;
      const ChainStep &last = step;
      if (last.state_adjoint == 0.0 || !std::isfinite(last.log_state)) continue;
      const long double log_mass = first.log_entry + log_product - start.log_product;
      const cand_pos_t window_start = chain.mirrored ? -t : start.position;
      const double probability =
          std::exp(static_cast<double>(log_mass - last.log_state)) * last.state_adjoint;
      (*row)[static_cast<std::size_t>(window_start)] += probability;
      if (loop_row) (*loop_row)[static_cast<std::size_t>(window_start)][c] += probability;
    }
  }

  cand_pos_t n_;
  bool streaming_ = false;
  bool by_loop_ = false;
  std::vector<cand_pos_t> windows_;
  std::array<std::vector<double>, loop_context_count> context_diff_;
  std::vector<std::vector<double>> window_diff_;
  std::vector<std::vector<LoopAccessibilityRow>> loop_window_diff_;
  std::vector<std::pair<UnpairedEmission, double>> whole_;
  struct ChainKeyHash {
    std::size_t operator()(const std::pair<int, cand_pos_t> &key) const {
      return std::hash<int>{}(key.first) * 1315423911u ^ std::hash<cand_pos_t>{}(key.second);
    }
  };
  std::unordered_map<std::pair<int, cand_pos_t>, Chain, ChainKeyHash> chains_;
};

}  // namespace cparty::scfg2::profiles

#endif  // CPARTY_SCFG2_PROFILES_PROFILE_ASSEMBLY_HH
