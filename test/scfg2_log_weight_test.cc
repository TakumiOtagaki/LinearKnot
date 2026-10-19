#include "scfg2/weights/runtime_snapshot.hh"
#include "rna_pair_table.hh"

extern "C" {
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/default.h"
}
#include "scfg2/api/beam.hh"
#include "scfg2/replay/w_final_exact_inside.hh"
#include "scfg2/semirings/logsumexp_tape.hh"
#include "scfg2/semirings/factorized_logsumexp.hh"

#include <cstdlib>
#include <array>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <type_traits>
#include <tuple>

namespace s = cparty::scfg2;
using Params = std::unique_ptr<vrna_exp_param_t, decltype(&std::free)>;

static void require(bool condition, const char *message) {
  if (!condition) throw std::runtime_error(message);
}
static void close(double a, double b, double tolerance = 1e-10) {
  require(std::isfinite(a) && std::isfinite(b) && std::abs(a-b) <= tolerance,
          "log weight / PF reference mismatch");
}
template <class T>
static void compare_table(const T &raw, const T &log) {
  if constexpr (std::is_array_v<T>) {
    for (std::size_t i=0; i<std::extent_v<T>; ++i) compare_table(raw[i],log[i]);
  } else if (raw > 0. && std::isfinite(raw)) close(std::log(raw),log);
}

static void pf_parameters_and_loops() {
  for (double temperature : {10.,37.,60.}) for (int smooth : {0,1})
    for (int dangles : {0,2}) {
      vrna_md_t md; vrna_md_set_default(&md);
      md.temperature=temperature; md.pf_smooth=smooth; md.dangles=dangles;
      Params raw(vrna_exp_params(&md), &std::free);
      Params log(cparty_log_pf_parameters(&md), &std::free);
      const double ratio=(temperature+K0)/Tmeasure;
      const auto rescale=[&](double g,double h) { return h-(h-g)*ratio; };
      const s::weights::InternalLoopExtrapolation extended{
          rescale(bulge37[30],bulgedH[30]),
          rescale(internal_loop37[30],internal_loopdH[30]),
          rescale(ninio37,niniodH),static_cast<double>(MAX_NINIO)};
#define COMPARE(FIELD) compare_table(raw->FIELD,log->FIELD)
      COMPARE(expMLbase); COMPARE(expMLclosing); COMPARE(expTermAU);
      COMPARE(expMLintern); COMPARE(expstack); COMPARE(expbulge);
      COMPARE(expinternal); COMPARE(exphairpin); COMPARE(expninio);
      COMPARE(expdangle5); COMPARE(expdangle3); COMPARE(expmismatchH);
      COMPARE(expmismatchI); COMPARE(expmismatch1nI); COMPARE(expmismatch23I);
      COMPARE(expmismatchM); COMPARE(expmismatchExt);
      COMPARE(expint11); COMPARE(expint21); COMPARE(expint22);
#undef COMPARE
      for (int type=1; type<=7; ++type) {
        for (int u : {3,4,6,10,30,31,100,10000}) {
          std::string loop(u+2,'A');
          const double q=exp_E_Hairpin(u,type,1,1,loop.c_str(),raw.get());
          if (q>0. && std::isfinite(q)) close(std::log(q),
              s::weights::log_hairpin(u,type,1,1,loop.c_str(),*log));
        }
        // Every loaded special loop: lookup and type-7 tetra behavior.
        for (auto [catalogue,u,stride] : {
            std::tuple{raw->Triloops,3,6}, std::tuple{raw->Tetraloops,4,7},
            std::tuple{raw->Hexaloops,6,9}}) {
          for (std::size_t pos=0;pos+u+2<=std::strlen(catalogue);pos+=stride) {
            std::string loop(catalogue+pos,u+2);
            const double q=exp_E_Hairpin(u,type,1,1,loop.c_str(),raw.get());
            if (q>0. && std::isfinite(q)) close(std::log(q),
                s::weights::log_hairpin(u,type,1,1,loop.c_str(),*log));
          }
        }
        for (int si=-1;si<=4;++si) for(int sj=-1;sj<=4;++sj) {
          close(std::log(exp_E_ExtLoop(type,si,sj,raw.get())),
                s::weights::log_stem(type,si,sj,false,*log));
          close(std::log(exp_E_MLstem(type,si,sj,raw.get())),
                s::weights::log_stem(type,si,sj,true,*log));
        }
        for (int inner=1;inner<=7;++inner) for (int u1=0;u1<=15;++u1)
          for(int u2=0;u2<=15;++u2) for(int base=1;base<=4;++base) {
            const double q=exp_E_IntLoop(u1,u2,type,inner,base,2,3,4,raw.get());
            if(q>0. && std::isfinite(q)) close(std::log(q),
                s::weights::log_internal_loop(u1,u2,type,inner,base,2,3,4,*log,{}));
          }
        for (int inner=1;inner<=7;++inner)
          for (auto [u1,u2] : {std::pair{0,31},std::pair{31,0},
              std::pair{1,30},std::pair{30,1},std::pair{15,16},
              std::pair{4,100},std::pair{0,10000}}) {
            const double q=s::weights::internal_loop_boltzmann_weight(
                u1,u2,type,inner,1,2,3,4,raw.get(),extended);
            close(std::log(q),s::weights::log_internal_loop(
                u1,u2,type,inner,1,2,3,4,*log,extended));
          }
      }
      raw->model_details.noGUclosure=log->model_details.noGUclosure=1;
      require(s::weights::log_internal_loop(1,30,3,1,1,2,3,4,*log,extended)
                  == s::weights::log_zero, "GU closure became reachable");
      close(std::log(exp_E_IntLoop(0,0,3,1,1,2,3,4,raw.get())),
            s::weights::log_internal_loop(0,0,3,1,1,2,3,4,*log,extended));
      // Vienna's INF sentinel is impossible, not a tiny positive branch.
      require(log->exphairpin[0] == s::weights::log_zero, "INF sentinel became reachable");
    }
}

static void typed_boundary_and_owner() {
  {
    const int saved_g=hairpin37[3], saved_h=hairpindH[3];
    hairpin37[3]=hairpindH[3]=200000;
    vrna_md_t md; vrna_md_set_default(&md);
    Params raw(vrna_exp_params(&md),&std::free), log(cparty_log_pf_parameters(&md),&std::free);
    require(raw->exphairpin[3] == 0., "finite-energy underflow probe ineffective");
    close(log->exphairpin[3], -200000.*10./log->kT,1e-8);
    hairpin37[3]=saved_g; hairpindH[3]=saved_h;
  }
  static_assert(!std::is_convertible_v<s::weights::LogWeight,double>);
  using S=s::semirings::LogSumExpSemiring;
  require(S::from_weight(s::weights::LogWeight{0.}) == S::one(), "neutral log weight lost");
  require(S::is_zero(S::from_weight(s::weights::LogWeight{s::weights::log_zero})), "impossible log weight admitted");
  s::semirings::LogSumExpTape tape;
  { s::semirings::TapeLogSumExpSemiring::Scope scope(tape);
    auto x=s::semirings::TapeLogSumExpSemiring::from_weight(s::weights::LogWeight{-800.});
    close(x.log_value,-800.);
  }
  s::semirings::FactorizedInsideGraph graph;
  { s::semirings::FactorizedLogSumExpSemiring::Scope scope(graph);
    auto x=s::semirings::FactorizedLogSumExpSemiring::from_weight(s::weights::LogWeight{-800.});
    close(x.log_value,-800.);
  }
  std::string sequence(16386,'A'), scaffold(sequence.size(),'.');
  cparty::scfg2::weights::RuntimeSnapshot owner(
      sequence, scaffold);
  s::replay::WFinalLogWeightContext ctx(owner);
  const double cp=ctx.vpl_split_vp_weight(4096);
  require(std::isfinite(cp) && cp < -745., "long factor underflowed before log injection");
  close(cp,4096*owner.log_pk_parameters().cp);
  s::engine::AllLocalWeightModel<s::replay::WFinalLogWeightContext> model(&ctx), absent(nullptr);
  const auto unpaired=s::deductions::expand_basic_rule(
      {s::grammar::NonTerminal::W,1,1,-1,-1},s::grammar::RuleId::W_EXTEND_UNPAIRED,{});
  close(model.log_weight_for(unpaired),0.);
  require(s::weights::is_zero(absent.local_weight_for(unpaired)), "null log context admitted");
  const double before=ctx.w_split_wmb_weight();
  const auto saved=s::weights::pseudoknot_energy_defaults().exterior;
  s::weights::pseudoknot_energy_defaults().exterior += 100;
  std::string small="AAAA", g="....";
  cparty::scfg2::weights::RuntimeSnapshot later(
      small, g);
  close(ctx.w_split_wmb_weight(),before);
  require(std::abs(later.log_pk_parameters().ps-before)>0.1,"owner isolation probe ineffective");
  s::weights::pseudoknot_energy_defaults().exterior=saved;
}

static void multiloop_closing_boundaries() {
  const std::string core="GAGGGAAAACCCAGGGAAAACCCAC";
  std::size_t cases=0;
  for (const auto &pair : {"GC","CG","AU","UA","GU","UG"})
    for (int left : {0,1}) for (int right : {0,1}) for (int dangles : {0,2}) {
      auto sequence=std::string(left,'A')+core+std::string(right,'A');
      const cand_pos_t i=left+1, j=left+core.size();
      sequence[i-1]=pair[0]; sequence[j-1]=pair[1];
      s::weights::RuntimeSnapshot owner(sequence,std::string(sequence.size(),'.'),
          s::weights::RuntimeOptions{dangles});
      s::replay::WFinalLogWeightContext ctx(owner);
      s::replay::WFinalExactInsideContext raw_ctx(owner);
      auto md=owner.log_pf_parameters().model_details;
      Params raw(vrna_exp_params(&md),&std::free);
      const int type=cparty::rna::pair_type_for(owner.encoded_base(j),owner.encoded_base(i));
      // Independent Vienna PF helper: closing neighbors are inside the pair,
      // including when the pair touches either or both sequence endpoints.
      const int si=dangles ? owner.encoded_base(j-1) : -1;
      const int sj=dangles ? owner.encoded_base(i+1) : -1;
      const double expected=std::log(raw->expMLclosing * exp_E_MLstem(type,si,sj,raw.get()));
      close(ctx.split_wm_wmv_weight(i,j),expected);
      close(ctx.split_wm_wmp_weight(i,j),expected);
      close(ctx.split_wmp_base_weight(i,j,2),expected+2*std::log(raw->expMLbase));
      close(std::log(raw_ctx.split_wm_wmv_weight(i,j)),expected);
      // Branch/exterior stems retain their outer-neighbor endpoint guards.
      const int branch_type=cparty::rna::pair_type_for(owner.encoded_base(i),owner.encoded_base(j));
      const int outer5=dangles && i>1 ? owner.encoded_base(i-1) : -1;
      const int outer3=dangles && j<owner.length() ? owner.encoded_base(j+1) : -1;
      close(ctx.wmv_stem_weight(i,j),std::log(exp_E_MLstem(branch_type,outer5,outer3,raw.get())));
      close(ctx.exp_extloop(i,j),std::log(exp_E_ExtLoop(branch_type,outer5,outer3,raw.get())));
      ++cases;
    }
  std::cout << "multiloop_closing_boundaries=ok cases=" << cases << '\n';
}

static void multiloop_fixed_structure_witnesses() {
  const std::string core="GAGGGAAAACCCAGGGAAAACCCAC";
  const std::string target="(.(((....))).(((....))).)";
  // Turner04, 37 C, dangles=2. Independent single-structure Vienna PF scores;
  // these are fixed targets, not a claim about their MFE optimality.
  const std::array<double,4> expected{{-1.7,-3.4,-1.9,-3.2}};
  std::size_t cases=0;
  for (int left : {0,1}) for (int right : {0,1})
    for (auto policy : {s::WmbClosurePolicy::RequireWeakClosure,
                        s::WmbClosurePolicy::LegacyCompatible}) {
      const auto sequence=std::string(left,'A')+core+std::string(right,'A');
      const auto full=std::string(left,'.')+target+std::string(right,'.');
      for (int scaffold_mode : {0,1,2}) {
        std::string scaffold(sequence.size(),'.'), generated=full;
        if (scaffold_mode==1) {
          scaffold[left]='('; scaffold[left+core.size()-1]=')';
          generated[left]=generated[left+core.size()-1]='.';
        } else if (scaffold_mode==2) {
          scaffold=full; generated.assign(sequence.size(),'.');
        }
        s::api::ExactSession session(sequence,scaffold,policy);
        const auto evaluated=session.evaluate_energy(generated);
        require(evaluated.ok && evaluated.reachable,"multiloop target unreachable");
        close(evaluated.energy,expected[left*2+right]);
        const auto exact=session.outside();
        require(exact.ok && exact.reachable,"multiloop exact outside failed");
        for (auto ranking : {s::api::BeamRanking::Baseline,s::api::BeamRanking::Outer}) {
          s::api::BeamConfig config{100}; config.ranking=ranking;
          config.wmb_closure_policy=policy;
          const auto beam=s::api::beam_outside(sequence,scaffold,config);
          require(beam.ok && beam.reachable,"multiloop beam outside failed");
          close(beam.log_partition_weight,exact.log_partition_weight);
          require(beam.pair_probabilities.size()==exact.pair_probabilities.size(),
                  "multiloop pair output length mismatch");
          for (std::size_t k=0;k<exact.pair_probabilities.size();++k) {
            const auto &a=exact.pair_probabilities[k], &b=beam.pair_probabilities[k];
            require(a.i==b.i && a.j==b.j,"multiloop pair coordinates differ");
            close(a.round_probability,b.round_probability);
            close(a.square_probability,b.square_probability);
            close(a.union_probability,b.union_probability);
          }
          require(beam.unpaired_probabilities.size()==exact.unpaired_probabilities.size(),
                  "multiloop unpaired output length mismatch");
          for (std::size_t k=0;k<exact.unpaired_probabilities.size();++k)
            close(beam.unpaired_probabilities[k],exact.unpaired_probabilities[k]);
        }
        ++cases;
      }
    }
  std::cout << "multiloop_fixed_structures=ok cases=" << cases << '\n';
}

static void public_numeric_limits() {
  const std::string sequence=std::string(150,'G')+"AAA"+std::string(150,'C');
  const std::string scaffold=std::string(150,'(')+"..."+std::string(150,')');
  s::api::ExactSession session(sequence,scaffold);
  const auto inside=session.inside();
  const auto best=session.viterbi();
  const auto probability=session.probability(std::vector<s::api::BasePair>{});
  require(inside.ok && inside.reachable && best.ok && best.reachable &&
          probability.ok && probability.reachable, "303-nt public exact log result failed");
  require(std::isinf(inside.partition_weight) && std::isinf(best.weight), "overflow witness no longer exercises projection");
  close(inside.ensemble_free_energy,-486.3,1e-8);
  close(best.energy,-486.3,1e-8);
  close(probability.log_probability,0.,1e-9);
  const auto beam=s::api::beam_inside(sequence,scaffold,s::api::BeamConfig{20});
  require(beam.ok && beam.reachable,"long beam log result failed");
  close(beam.log_partition_weight,inside.log_partition_weight,1e-9);
  const auto neutral=s::api::beam_inside("AAAAAAAAAA","..........",s::api::BeamConfig{20});
  require(neutral.ok && neutral.reachable,"zero log weight was pruned");
  close(neutral.log_partition_weight,0.);
  setenv("CPARTY_PF_OVERRIDE_SCALE","0.5",1);
  const auto forced=s::api::inside("GGGAAACCC","(((...)))");
  unsetenv("CPARTY_PF_OVERRIDE_SCALE");
  require(forced.ok,"legacy scale override affected public log API");
  close(forced.ensemble_free_energy,-1.2,1e-9);
}

int main() {
  try {
    setenv("CPARTY_FIXED_ENERGY_PARAM_FILE","params/rna_Turner04.par",1);
    unsetenv("CPARTY_PF_OVERRIDE_SCALE");
    pf_parameters_and_loops();
    multiloop_closing_boundaries();
    multiloop_fixed_structure_witnesses();
    typed_boundary_and_owner();
    public_numeric_limits();
    std::cout << "scfg2_log_weight=ok\n";
  } catch(const std::exception &e) {
    std::cerr << e.what() << '\n'; return 1;
  }
}
