#ifndef CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_PROFILE_HH
#define CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_PROFILE_HH

#include "scfg2/deductions/deduction.hh"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <map>
#include <set>
#include <tuple>
#include <vector>

namespace cparty::scfg2::beam {

enum class WmbpBeVpPrefixMode : unsigned char {
  Wmbp,
  Wmbw,
  WiEmpty,
  Wi,
};

inline constexpr bool is_wmbp_be_vp_rule(grammar::RuleId rule) {
  return rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP ||
         rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP ||
         rule == grammar::RuleId::WMBP_SPLIT_BE_WI_VP;
}

struct WmbpBeVpPrefixKey {
  WmbpBeVpPrefixMode mode = WmbpBeVpPrefixMode::Wmbp;
  deductions::ItemKey be{};
  deductions::ItemKey left{};
  bool has_left = false;
};

struct WmbpBeVpPrefixKeyHash {
  std::size_t operator()(const WmbpBeVpPrefixKey &key) const noexcept {
    std::size_t out = static_cast<std::size_t>(key.mode);
    out = (out * 1315423911u) ^ deductions::ItemKeyHash{}(key.be);
    out = (out * 1315423911u) ^ static_cast<std::size_t>(key.has_left);
    if (key.has_left) {
      out = (out * 1315423911u) ^ deductions::ItemKeyHash{}(key.left);
    }
    return out;
  }
};

inline auto item_tie_key(const deductions::ItemKey &item) {
  return std::make_tuple(static_cast<unsigned char>(item.nonterminal),
                         item.i,
                         item.j,
                         item.ip,
                         item.jp);
}

inline auto wmbp_be_vp_prefix_tie_key(const WmbpBeVpPrefixKey &key) {
  return std::make_tuple(static_cast<unsigned char>(key.mode),
                         item_tie_key(key.be),
                         key.has_left,
                         item_tie_key(key.left));
}

inline bool operator==(const WmbpBeVpPrefixKey &lhs,
                       const WmbpBeVpPrefixKey &rhs) {
  return lhs.mode == rhs.mode &&
         deductions::ItemKeyEq{}(lhs.be, rhs.be) &&
         lhs.has_left == rhs.has_left &&
         (!lhs.has_left || deductions::ItemKeyEq{}(lhs.left, rhs.left));
}

inline WmbpBeVpPrefixKey wmbp_be_vp_prefix_key(
    const deductions::Deduction &deduction) {
  WmbpBeVpPrefixKey out;
  out.be = deduction.children[0].item;
  if (deduction.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBP_VP) {
    out.mode = WmbpBeVpPrefixMode::Wmbp;
    out.left = deduction.children[1].item;
    out.has_left = true;
  } else if (deduction.rule == grammar::RuleId::WMBP_SPLIT_BE_WMBW_VP) {
    out.mode = WmbpBeVpPrefixMode::Wmbw;
    out.left = deduction.children[1].item;
    out.has_left = true;
  } else if (deduction.child_count == 2) {
    out.mode = WmbpBeVpPrefixMode::WiEmpty;
  } else {
    out.mode = WmbpBeVpPrefixMode::Wi;
    out.left = deduction.children[1].item;
    out.has_left = true;
  }
  return out;
}

inline const deductions::ItemKey &wmbp_be_vp_right_child(
    const deductions::Deduction &deduction) {
  return deduction.children[deduction.child_count - 1].item;
}

struct WmbpBeVpPrefixProfile {
  std::size_t source_count = 0;
  std::size_t wmbp_count = 0;
  std::size_t wmbw_count = 0;
  std::size_t wi_empty_count = 0;
  std::size_t wi_count = 0;
  std::size_t prefix_count = 0;
  std::size_t reused_prefix_count = 0;
  std::size_t right_vp_count = 0;
  std::size_t reused_right_vp_count = 0;
  std::size_t multi_source_prefix_count = 0;
  std::size_t maximum_sources_per_prefix = 0;
  std::size_t local_variant_prefix_count = 0;
  std::size_t maximum_locals_per_prefix = 0;
  std::size_t maximum_prefixes_per_right_vp = 0;
  std::size_t column_count = 0;
  std::size_t maximum_sources_per_column = 0;
};

class WmbpBeVpPrefixProfiler {
 private:
  struct PrefixLess {
    bool operator()(const WmbpBeVpPrefixKey &lhs,
                    const WmbpBeVpPrefixKey &rhs) const {
      return wmbp_be_vp_prefix_tie_key(lhs) <
             wmbp_be_vp_prefix_tie_key(rhs);
    }
  };
  struct ItemLess {
    bool operator()(const deductions::ItemKey &lhs,
                    const deductions::ItemKey &rhs) const {
      return item_tie_key(lhs) < item_tie_key(rhs);
    }
  };

 public:
  template <class Chart>
  void observe(const deductions::Deduction &deduction,
               const Chart &chart,
               double local_weight,
               double direct_contribution) {
    if (!is_wmbp_be_vp_rule(deduction.rule)) return;

    const auto prefix = wmbp_be_vp_prefix_key(deduction);
    const auto &right_vp = wmbp_be_vp_right_child(deduction);
    double prefix_factor = local_weight * chart.get(prefix.be);
    if (prefix.has_left) prefix_factor *= chart.get(prefix.left);
    const double reconstructed = prefix_factor * chart.get(right_vp);
    const double tolerance =
        std::max(1e-14,
                 std::max(std::fabs(direct_contribution),
                          std::fabs(reconstructed)) *
                     1e-12);
    if (std::fabs(direct_contribution - reconstructed) <= tolerance) {
      ++factor_identity_count_;
    }

    switch (prefix.mode) {
      case WmbpBeVpPrefixMode::Wmbp: ++wmbp_count_; break;
      case WmbpBeVpPrefixMode::Wmbw: ++wmbw_count_; break;
      case WmbpBeVpPrefixMode::WiEmpty: ++wi_empty_count_; break;
      case WmbpBeVpPrefixMode::Wi: ++wi_count_; break;
    }
    ++sources_per_prefix_[prefix];
    auto &locals = locals_per_prefix_[prefix];
    const auto same_local = [&](double value) {
      const double local_tolerance =
          std::max(1e-14,
                   std::max(std::fabs(value), std::fabs(local_weight)) *
                       1e-12);
      return std::fabs(value - local_weight) <= local_tolerance;
    };
    if (std::none_of(locals.begin(), locals.end(), same_local)) {
      locals.push_back(local_weight);
    }
    prefixes_per_right_vp_[right_vp].insert(prefix);
    ++sources_per_column_[deduction.parent.j];
    ++source_count_;
  }

  WmbpBeVpPrefixProfile profile() const {
    WmbpBeVpPrefixProfile out;
    out.source_count = source_count_;
    out.wmbp_count = wmbp_count_;
    out.wmbw_count = wmbw_count_;
    out.wi_empty_count = wi_empty_count_;
    out.wi_count = wi_count_;
    out.prefix_count = sources_per_prefix_.size();
    out.reused_prefix_count = source_count_ - sources_per_prefix_.size();
    out.right_vp_count = prefixes_per_right_vp_.size();
    out.reused_right_vp_count = source_count_ - prefixes_per_right_vp_.size();
    for (const auto &[prefix, count] : sources_per_prefix_) {
      (void)prefix;
      if (count > 1) ++out.multi_source_prefix_count;
      out.maximum_sources_per_prefix =
          std::max(out.maximum_sources_per_prefix, count);
    }
    for (const auto &[prefix, locals] : locals_per_prefix_) {
      (void)prefix;
      if (locals.size() > 1) ++out.local_variant_prefix_count;
      out.maximum_locals_per_prefix =
          std::max(out.maximum_locals_per_prefix, locals.size());
    }
    for (const auto &[vp, prefixes] : prefixes_per_right_vp_) {
      (void)vp;
      out.maximum_prefixes_per_right_vp =
          std::max(out.maximum_prefixes_per_right_vp, prefixes.size());
    }
    out.column_count = sources_per_column_.size();
    for (const auto &[j, count] : sources_per_column_) {
      (void)j;
      out.maximum_sources_per_column =
          std::max(out.maximum_sources_per_column, count);
    }
    return out;
  }

  bool all_factor_identities_hold() const {
    return factor_identity_count_ == source_count_;
  }

 private:
  std::map<WmbpBeVpPrefixKey, std::size_t, PrefixLess> sources_per_prefix_;
  std::map<WmbpBeVpPrefixKey, std::vector<double>, PrefixLess>
      locals_per_prefix_;
  std::map<deductions::ItemKey,
           std::set<WmbpBeVpPrefixKey, PrefixLess>,
           ItemLess>
      prefixes_per_right_vp_;
  std::map<cand_pos_t, std::size_t> sources_per_column_;
  std::size_t source_count_ = 0;
  std::size_t wmbp_count_ = 0;
  std::size_t wmbw_count_ = 0;
  std::size_t wi_empty_count_ = 0;
  std::size_t wi_count_ = 0;
  std::size_t factor_identity_count_ = 0;
};

}  // namespace cparty::scfg2::beam

#endif  // CPARTY_SCFG2_BEAM_WMBP_BE_VP_PREFIX_PROFILE_HH
