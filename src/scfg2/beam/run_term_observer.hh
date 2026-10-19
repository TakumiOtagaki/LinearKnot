#ifndef CPARTY_SCFG2_BEAM_RUN_TERM_OBSERVER_HH
#define CPARTY_SCFG2_BEAM_RUN_TERM_OBSERVER_HH

#include "scfg2/types.hh"
#include "scfg2/grammar/rules.hh"
#include "scfg2/semirings/value_boundary.hh"

namespace cparty::scfg2::beam {

// Observation-only hook for consumers of deduction posteriors (PK-CapR,
// PK-Raccess).  The beam folds several run-emitting rules into private
// aggregates, so the generic candidate observer sees one aggregate per
// parent instead of one contribution per run.  These callbacks expose the
// pieces the aggregates are built from; they never change a value.
template <class Value>
class RunTermObserver {
 public:
  virtual ~RunTermObserver() = default;

  // One explicit term of a selected start-prefix cell.  Its rule emits the
  // unpaired run [lo, hi] (empty when lo > hi).
  virtual void on_run_term(grammar::RuleId rule,
                           cand_pos_t lo,
                           cand_pos_t hi,
                           const Value &term) = 0;

  // VPR right-run aggregate C(i,t) = sum_k VP(i,k) * cp^(t-k) over the
  // retained seeds k < t.  `carried` is true when C(i,t) includes
  // C(i,t-1) * cp.
  virtual void on_vpr_run_cell(cand_pos_t i,
                               cand_pos_t t,
                               const Value &total,
                               bool carried) = 0;

  // KeepAll start-prefix columns grow runs leftwards one base at a time.
  // `unit` is the mass of all runs of `rule`'s family at column j that cover
  // `base` and extend to its right; `alias` is a second node with the same
  // value that also carries this mass (or nullptr).  `linked` is true when
  // `unit` includes the unit of base + 1.
  virtual void on_left_run_unit(grammar::RuleId rule,
                                cand_pos_t j,
                                cand_pos_t base,
                                const Value &unit,
                                const Value *alias,
                                bool linked) = 0;
};

template <class Value>
RunTermObserver<Value> *&run_term_observer();

template <class Value>
inline void notify_left_run_unit(grammar::RuleId rule,
                                 cand_pos_t j,
                                 cand_pos_t base,
                                 const Value &unit,
                                 const Value *alias,
                                 bool linked) {
  if (auto *observer = run_term_observer<Value>()) {
    observer->on_left_run_unit(rule, j, base, unit, alias, linked);
  }
}

template <class Value>
RunTermObserver<Value> *&run_term_observer() {
  thread_local RunTermObserver<Value> *observer = nullptr;
  return observer;
}

template <class Value>
inline void notify_run_term(grammar::RuleId rule,
                            cand_pos_t lo,
                            cand_pos_t hi,
                            const Value &term) {
  if (auto *observer = run_term_observer<Value>()) {
    observer->on_run_term(rule, lo, hi, term);
  }
}

// Installs an observer for the lifetime of the scope.
template <class Value>
class ScopedRunTermObserver {
 public:
  explicit ScopedRunTermObserver(RunTermObserver<Value> *observer)
      : previous_(run_term_observer<Value>()) {
    run_term_observer<Value>() = observer;
  }
  ~ScopedRunTermObserver() { run_term_observer<Value>() = previous_; }
  ScopedRunTermObserver(const ScopedRunTermObserver &) = delete;
  ScopedRunTermObserver &operator=(const ScopedRunTermObserver &) = delete;

 private:
  RunTermObserver<Value> *previous_;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_RUN_TERM_OBSERVER_HH
