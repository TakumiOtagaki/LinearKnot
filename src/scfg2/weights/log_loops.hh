#ifndef CPARTY_SCFG2_WEIGHTS_LOG_LOOPS_HH
#define CPARTY_SCFG2_WEIGHTS_LOG_LOOPS_HH

#include "scfg2/weights/local_weight.hh"
#include "scfg2/weights/internal_loop.hh"
#include <cstring>

namespace cparty::scfg2::weights {

// Log-space translation of the bundled ViennaRNA PF helpers in
// loops/hairpin.h, loops/internal.h, loops/multibranch.h and loops/external.h.
// Origin and terms: THIRD_PARTY_NOTICES.md and src/ViennaRNA/COPYING.
// P is built by cparty_log_pf_parameters: exp* fields store log weights.
// Preserve the bundled PF helper branches, not integer MFE energy semantics.
inline double log_hairpin(int u, int type, short si, short sj,
                          const char *sequence, const vrna_exp_param_t &p) {
  double q = p.exphairpin[std::min(u, 30)];
  if (u > 30) q -= p.lxc * std::log(u / 30.) * 10. / p.kT;
  if (u < 3) return q;
  if (sequence && p.model_details.special_hp) {
    if (u == 4) {
      char loop[7]{}; std::memcpy(loop, sequence, 6);
      if (const char *s = std::strstr(p.Tetraloops, loop)) {
        const double tetra = p.exptetra[(s - p.Tetraloops) / 7];
        if (type != 7) return tetra;
        q += tetra;
      }
    } else if (u == 6) {
      char loop[9]{}; std::memcpy(loop, sequence, 8);
      if (const char *s = std::strstr(p.Hexaloops, loop))
        return p.exphex[(s - p.Hexaloops) / 9];
    } else if (u == 3) {
      char loop[6]{}; std::memcpy(loop, sequence, 5);
      if (const char *s = std::strstr(p.Triloops, loop))
        return p.exptri[(s - p.Triloops) / 6];
      return q + (type > 2 ? p.expTermAU : 0.);
    }
  }
  return q + p.expmismatchH[type][si][sj];
}

inline double log_internal_loop(int u1, int u2, int type, int type2,
    short si, short sj, short sp, short sq, const vrna_exp_param_t &p,
    const InternalLoopExtrapolation &extended) {
  const int ul = std::max(u1, u2), us = std::min(u1, u2), size = u1 + u2;
  if (ul == 0) return p.expstack[type][type2];
  if (p.model_details.noGUclosure &&
      (type == 3 || type == 4 || type2 == 3 || type2 == 4)) return log_zero;
  const auto round = [&](double value) {
    return p.model_details.pf_smooth ? value : std::trunc(value);
  };
  if (size > MAXLOOP) {
    const double initiation = round((us == 0 ? extended.bulge_energy_30 :
        extended.internal_energy_30) + p.lxc * std::log(size / 30.));
    double q = -initiation * 10. / p.kT;
    if (us == 0) return q + (type > 2 ? p.expTermAU : 0.) +
                            (type2 > 2 ? p.expTermAU : 0.);
    q -= std::min(extended.max_ninio, (ul-us)*round(extended.ninio_per_base)) * 10. / p.kT;
    if (us == 1) return q + p.expmismatch1nI[type][si][sj] + p.expmismatch1nI[type2][sq][sp];
    return q + p.expmismatchI[type][si][sj] + p.expmismatchI[type2][sq][sp];
  }
  if (us == 0) {
    double q = p.expbulge[ul];
    if (ul == 1) return q + p.expstack[type][type2];
    return q + (type > 2 ? p.expTermAU : 0.) + (type2 > 2 ? p.expTermAU : 0.);
  }
  if (us == 1) {
    if (ul == 1) return p.expint11[type][type2][si][sj];
    if (ul == 2) return u1 == 1 ? p.expint21[type][type2][si][sq][sj]
                              : p.expint21[type2][type][sq][si][sp];
    return p.expinternal[size] + p.expmismatch1nI[type][si][sj] +
           p.expmismatch1nI[type2][sq][sp] + p.expninio[2][ul-us];
  }
  if (us == 2 && ul == 2) return p.expint22[type][type2][si][sp][sq][sj];
  if (us == 2 && ul == 3) return p.expinternal[5] + p.expmismatch23I[type][si][sj] +
                               p.expmismatch23I[type2][sq][sp] + p.expninio[2][1];
  return p.expinternal[size] + p.expmismatchI[type][si][sj] +
         p.expmismatchI[type2][sq][sp] + p.expninio[2][ul-us];
}

inline double log_stem(int type, int si, int sj, bool multiloop,
                       const vrna_exp_param_t &p) {
  double q = type > 2 ? p.expTermAU : 0.;
  if (multiloop) q += p.expMLintern[type];
  if (si >= 0 && sj >= 0) q += multiloop ? p.expmismatchM[type][si][sj]
                                                       : p.expmismatchExt[type][si][sj];
  else if (si >= 0) q += p.expdangle5[type][si];
  else if (sj >= 0) q += p.expdangle3[type][sj];
  return q;
}

}  // namespace cparty::scfg2::weights
#endif
