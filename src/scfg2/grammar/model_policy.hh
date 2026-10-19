#ifndef CPARTY_SCFG2_GRAMMAR_MODEL_POLICY_HH
#define CPARTY_SCFG2_GRAMMAR_MODEL_POLICY_HH

namespace cparty::scfg2 {

// Only the domain of completed WMB changes. WMBP/WMBW remain partial
// band-chain carriers; local weights (including Psm ownership) are unchanged.
enum class WmbClosurePolicy {
  RequireWeakClosure,
  LegacyCompatible,
};

inline constexpr const char *wmb_closure_policy_name(WmbClosurePolicy policy) {
  return policy == WmbClosurePolicy::LegacyCompatible
             ? "legacy-compatible-wmb"
             : "weakly-closed-wmb";
}

}  // namespace cparty::scfg2

#endif  // CPARTY_SCFG2_GRAMMAR_MODEL_POLICY_HH
