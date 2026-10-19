#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/weights/log_loops.hh"
#include "rna_pair_table.hh"
#include "scfg2/weights/internal_loop.hh"
#include "scfg2/replay/w_final_exact_inside.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

extern "C" {
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/params/io.h"
}

namespace {
using cparty::scfg2::weights::InternalLoopExtrapolation;
using cparty::scfg2::weights::internal_loop_boltzmann_weight;

bool close(double a, double b) {
  return std::abs(a - b) <= 2e-12 * std::max(std::abs(a), std::abs(b));
}

InternalLoopExtrapolation extrapolation(const vrna_exp_param_t &p) {
  const double ratio = (p.temperature + K0) / Tmeasure;
  const auto rescale = [ratio](double g, double h) {
    return h - (h - g) * ratio;
  };
  return {rescale(bulge37[30], bulgedH[30]),
          rescale(internal_loop37[30], internal_loopdH[30]),
          rescale(ninio37, niniodH), static_cast<double>(MAX_NINIO)};
}

void check_parameters(double temperature, bool smooth) {
  vrna_md_t md;
  vrna_md_set_default(&md);
  md.temperature = temperature;
  md.pf_smooth = smooth;
  auto *p = vrna_exp_params(&md);
  const auto extended = extrapolation(*p);
  for (int u1 = 0; u1 <= MAXLOOP; ++u1)
    for (int u2 = 0; u1 + u2 <= MAXLOOP; ++u2)
      for (int type = 1; type <= 6; ++type)
        for (int type2 = 1; type2 <= 6; ++type2)
          for (short base = 1; base <= 4; ++base) {
            const double before = exp_E_IntLoop(u1, u2, type, type2,
                                                base, 1, 2, 3, p);
            const double after = internal_loop_boltzmann_weight(
                u1, u2, type, type2, base, 1, 2, 3, p, extended);
            assert(before == after);
          }

  // At 37 C, the independent integer-energy helper supplies an oracle for
  // all loop branches. Add back only the known fractional initiation term
  // for smooth PF parameters; small-loop weights above remain bit-identical.
  if (temperature == 37.0) {
    auto *energy_params = vrna_params(&md);
    for (const auto &[u1, u2] : std::vector<std::pair<int, int>>{
             {16,15}, {15,16}, {31,0}, {0,31}, {1,30}, {30,1},
             {2,40}, {40,2}, {500,500}, {1,1000}, {1000,1}})
      for (int type = 1; type <= 6; ++type)
        for (int type2 = 1; type2 <= 6; ++type2)
          for (short base = 1; base <= 4; ++base) {
            double expected_energy = E_IntLoop(u1, u2, type, type2,
                                               base, 1, 2, 3, energy_params);
            const double increment = p->lxc * std::log((u1 + u2) / 30.0);
            if (smooth) expected_energy += increment - std::trunc(increment);
            const double expected = std::exp(-expected_energy * 10.0 / p->kT);
            const double actual = internal_loop_boltzmann_weight(
                u1, u2, type, type2, base, 1, 2, 3, p, extended);
            assert(close(actual, expected));
          }
    free(energy_params);
  }

  p->model_details.noGUclosure = 1;
  assert(internal_loop_boltzmann_weight(31, 0, 3, 1, 1, 1, 1, 1,
                                       p, extended) == 0.0);
  assert(internal_loop_boltzmann_weight(16, 15, 1, 4, 1, 1, 1, 1,
                                       p, extended) == 0.0);
  assert(internal_loop_boltzmann_weight(16, 15, 2, 1, 1, 1, 1, 1,
                                       p, extended) > 0.0);
  free(p);
}

void check_fractional_rounding_and_asymmetry() {
  vrna_md_t md;
  vrna_md_set_default(&md);
  auto *p = vrna_exp_params(&md);
  // A custom small NINIO slope must continue beyond index 30, before the
  // energy cap is reached. Clamping the asymmetry index to 30 is incorrect.
  InternalLoopExtrapolation e{610.75, 370.25, 1.75, 300.0};
  for (const bool smooth : {true, false}) {
    p->model_details.pf_smooth = smooth;
    const double weight = internal_loop_boltzmann_weight(
        101, 1, 2, 1, 1, 1, 1, 1, p, e);
    const double mismatch = p->expmismatch1nI[2][1][1] *
                            p->expmismatch1nI[1][1][1];
    const double recovered = -std::log(weight / mismatch) * p->kT / 10.0;
    const double initiation = 370.25 + p->lxc * std::log(102.0 / 30.0);
    const double expected = smooth ? initiation + 175.0
                                   : std::trunc(initiation) + 100.0;
    assert(std::abs(recovered - expected) < 1e-10);
    e.ninio_per_base = 0.0;
    const double zero_asymmetry = internal_loop_boltzmann_weight(
        101, 1, 2, 1, 1, 1, 1, 1, p, e);
    assert(close(weight / zero_asymmetry,
                 std::exp(-(smooth ? 175.0 : 100.0) * 10.0 / p->kT)));
    e.ninio_per_base = 1.75;
  }
  free(p);
}

double raw_internal_loop_factor(const cparty::scfg2::weights::RuntimeSnapshot &owner,
                                cand_pos_t i, cand_pos_t j, cand_pos_t k, cand_pos_t l) {
  return std::exp(cparty::scfg2::weights::log_internal_loop(k-i-1, j-l-1,
      cparty::rna::pair_type_for(owner.encoded_base(i), owner.encoded_base(j)),
      cparty::rna::reverse_pair_type(cparty::rna::pair_type_for(owner.encoded_base(k), owner.encoded_base(l))),
      owner.adjacent_base(i+1), owner.adjacent_base(j-1), owner.adjacent_base(k-1), owner.adjacent_base(l+1),
      owner.log_pf_parameters(), owner.loop_extrapolation()));
}

void check_runtime_be(int left, int right, double expected_kcal) {
  using namespace cparty::scfg2;
  const int n = left + right + 11, k = left + 2, l = left + 10;
  std::string sequence(n, 'A');
  for (int i : {1,k}) sequence[i-1] = 'G';
  for (int i : {l,n}) sequence[i-1] = 'C';
  std::string scaffold = "(" + std::string(left, '.') + "(" +
      std::string(7, '.') + ")" + std::string(right, '.') + ")";
  weights::RuntimeSnapshot owner(sequence, scaffold);
  replay::WFinalExactInsideContext ctx(owner);
  const deductions::ItemKey parent{grammar::NonTerminal::BE,1,n,k,l};
  assert(deductions::exact_be_split_allows(
      parent, grammar::RuleId::BE_INTERNAL_LOOP, k, l, ctx));
  const auto sources = deductions::generate_exact_basic_deductions(
      parent, grammar::RuleId::BE_INTERNAL_LOOP, ctx);
  assert(sources.size() == 1);
  const double raw = raw_internal_loop_factor(owner,1,n,k,l);
  const double energy = -std::log(raw) * owner.parameters().kT / 1000.0;
  assert(std::abs(energy - expected_kcal) < 1e-12);
  const double expected = std::pow(raw, owner.pseudoknot_model().internal_multiplier);
  assert(close(ctx.be_internal_loop_weight(1,n,k,l), expected));
  // Parameter ownership: a new global parameter set must not change the
  // existing owner's long-loop weights while its short-loop PF tables stay
  // unchanged. Poison the raw base/slope fields and restore them afterward.
  const int old_internal = internal_loop37[30], old_bulge = bulge37[30];
  const int old_ninio = ninio37;
  internal_loop37[30] += 123;
  bulge37[30] += 234;
  ninio37 += 345;
  assert(raw_internal_loop_factor(owner,1,n,k,l) == raw);
  internal_loop37[30] = old_internal;
  bulge37[30] = old_bulge;
  ninio37 = old_ninio;
}
}  // namespace

int main() {
  read_parameter_file("params/rna_Turner04.par");
  for (double temperature : {10.0, 37.0, 55.0})
    for (bool smooth : {false, true}) check_parameters(temperature, smooth);
  check_fractional_rounding_and_asymmetry();
  check_runtime_be(16,14,4.9);
  check_runtime_be(16,15,4.335365791303965);
  check_runtime_be(31,0,6.135365791303965);
  std::cout << "scfg2_internal_loop_weight=ok\n";
}
