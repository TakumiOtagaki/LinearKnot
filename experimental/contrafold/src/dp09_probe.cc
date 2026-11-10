// Exact BPP and arbitrary pair-adjoint VJP for the twelve consumed DP09
// coordinates, on the unchanged shared generator, executor and local weights.
//
// With nearest-neighbor tables fixed, every local log weight is affine in the
// DP09 scalars and in the two band multipliers. Weights at a displacement
// delta are read from one snapshot built at the displaced (continuous) DP09
// values. The feature of coordinate k is the local energy change per unit
// step (1 kcal/mol, or +1 for a multiplier), read from unit-perturbed
// snapshots. For a BPP objective with fixed pair adjoints a,
// dL/d delta_k = -Cov(A, F_k) / RT.
//
// Fast feature path: a deduction whose log weight is bitwise unchanged under
// one snapshot that moves all twelve coordinates at once (incommensurate
// steps) reads no DP09 coordinate, and its features are zero; only the other
// deductions evaluate the twelve unit snapshots. controls.py checks the VJP
// against finite differences.
//
// Protocol, one request per line, one JSON response per line:
//   R id sequence scaffold          register a case
//   E id d1..d12                    BPP and expected features at displacement d
//   B id d1..d12                    BPP only (fitting path)
//   V id version count (i j a)*     VJP at the evaluated displacement
//   T id k step                     log Z from a truly perturbed snapshot
//   D id                            release a registered case
#include "learning.hh"
#include "scfg2/engine/outside_marginals.hh"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

namespace sc=cparty::scfg2;
namespace ex=linearknot_experiment;
using Item=sc::deductions::ItemKey;
using Deduction=sc::deductions::Deduction;
using Sum=sc::semirings::LogSumExpSemiring;
constexpr int K=12;
using Vec=std::array<double,K>;
static const char *NAMES[K]={"exterior","in_multiloop","in_pseudoloop","band","unpaired","nested_region",
  "multiloop_branch","band_multiloop_closing","band_multiloop_branch","band_multiloop_unpaired",
  "stack_multiplier","internal_multiplier"};

void require(bool value,const char *message){if(!value)throw std::runtime_error(message);}
template<class T>void vector_json(const T &v){std::cout<<'[';bool first=true;for(double x:v){if(!first)std::cout<<',';first=false;std::cout<<x;}std::cout<<']';}

// Step applied to coordinate k of model m; scalar fields use 10 cal/mol units.
void step_into(ex::Model &m,int k,double step){
  const double units=step*100.;
  switch(k){
    case 0:m.exterior+=units;break;case 1:m.in_multiloop+=units;break;case 2:m.in_pseudoloop+=units;break;
    case 3:m.band+=units;break;case 4:m.unpaired+=units;break;case 5:m.nested_region+=units;break;
    case 6:m.multiloop_branch+=units;break;case 7:m.band_multiloop_closing+=units;break;
    case 8:m.band_multiloop_branch+=units;break;case 9:m.band_multiloop_unpaired+=units;break;
    case 10:m.stack_multiplier+=step;break;case 11:m.internal_multiplier+=step;break;
    default:throw std::invalid_argument("coordinate index");
  }
}
ex::Model stepped(int k,double step){ex::Model m;step_into(m,k,step);return m;}
ex::Model displaced(const std::array<double,12> &d){ex::Model m;for(int k=0;k<12;++k)step_into(m,k,d[k]);return m;}
// Incommensurate steps for the all-coordinate zero test.
ex::Model combined_probe(){std::array<double,12> d;for(int k=0;k<12;++k)d[k]=0.731+0.1*std::sqrt(2.+k);return displaced(d);}

// Scalar observable moment: log mass and the mass-normalized expectation of A.
struct Local {double log_weight=-INFINITY;double a=0.;};
struct ScalarMoment {
  using value_type=Local;
  static Local zero(){return {};}
  static Local one(){return {0.,0.};}
  static bool is_zero(const Local &v){return v.log_weight==-INFINITY;}
  static Local from_weight(const Local &v){return v;}
  static Local plus(const Local &x,const Local &y){
    if(is_zero(x))return y;if(is_zero(y))return x;
    const double hi=std::max(x.log_weight,y.log_weight);
    const double log_weight=hi+std::log1p(std::exp(std::min(x.log_weight,y.log_weight)-hi));
    return {log_weight,std::exp(x.log_weight-log_weight)*x.a+std::exp(y.log_weight-log_weight)*y.a};
  }
  static Local times(const Local &x,const Local &y){
    if(is_zero(x)||is_zero(y))return zero();return {x.log_weight+y.log_weight,x.a+y.a};
  }
};

struct Case {
  std::string sequence,g;int n;ex::Snapshot baseline,combined;std::vector<std::unique_ptr<ex::Snapshot>> unit;
  std::unique_ptr<ex::Snapshot> current;double rt;
  std::vector<Item> schedule;Item root;std::vector<double> probability,adjoint;
  Vec delta{},expected{};bool has_expected=false;double log_z=-INFINITY;std::size_t version=0;
  Case(std::string seq,std::string scaffold):sequence(std::move(seq)),g(std::move(scaffold)),n(sequence.size()),
      baseline(sequence,g,ex::Model{}),combined(sequence,g,combined_probe()),root{sc::grammar::NonTerminal::W,1,n,-1,-1},
      probability((n+1)*(n+1)),adjoint((n+1)*(n+1)){
    std::string error;require(sc::api::detail::validate_session_inputs(sequence,g,&error),"invalid sequence/G");
    rt=baseline.owner->partition_free_energy_from_log_weight(-1.);
    for(int k=0;k<K;++k)unit.push_back(std::make_unique<ex::Snapshot>(sequence,g,stepped(k,1.)));
    schedule=sc::engine::build_pkfree_basic_exact_schedule(n,*baseline.context);
  }
  auto generate(const Item &item){
    auto all=sc::deductions::generate_exact_basic_deductions_for_item(item,*baseline.context);
    for(const auto &d:all)sc::engine::require_legacy_span_schedule_support(d);return all;
  }
  double pair_score(const Deduction &d) const {
    std::vector<sc::engine::TracebackPair> emitted;
    sc::engine::append_deduction_generated_pairs(d,*baseline.context,&emitted);
    double value=0.; for(auto p:emitted) value+=pkfit::pair(sequence,p.i,p.j); return value;
  }
  double weight_at(const Deduction &d)const{return current->local->log_weight_for(d)+pair_score(d);}
  // Features of d; zero without unit evaluations when d reads no DP09 coordinate.
  void features(const Deduction &d,Vec &f)const{
    f={};const double base=baseline.local->log_weight_for(d);
    if(!std::isfinite(base)||combined.local->log_weight_for(d)==base)return;
    for(int k=0;k<K;++k){
      const double changed=unit[k]->local->log_weight_for(d);
      require(std::isfinite(changed),"unit step changes admission");f[k]=-rt*(changed-base);
    }
  }
  Local local(const Deduction &d){
    Local out;out.log_weight=weight_at(d);if(!std::isfinite(out.log_weight))return {};
    thread_local std::vector<sc::engine::TracebackPair> emitted;emitted.clear();
    sc::engine::append_deduction_generated_pairs(d,*baseline.context,&emitted);
    for(auto p:emitted)out.a+=adjoint[p.i*(n+1)+p.j];
    return out;
  }
  void evaluate(const std::string &id,const Vec &d,bool with_features){
    const auto start=std::chrono::steady_clock::now();delta=d;current=std::make_unique<ex::Snapshot>(sequence,g,displaced(d));
    const auto weight=[&](const Deduction &x){return sc::weights::LogWeight{weight_at(x)};};
    const auto inside=sc::engine::run_schedule<Sum>(schedule,[&](const Item &i){return generate(i);},weight);
    log_z=inside.get(root);require(std::isfinite(log_z),"empty support");
    expected={};has_expected=with_features;
    const auto observe=[&](const Deduction &x,double beta){
      if(!with_features)return;
      double mass=beta+weight_at(x);for(std::size_t c=0;c<x.child_count;++c)mass+=inside.get(sc::engine::RuntimeDeductionChildItem{}(x,c));
      if(!std::isfinite(mass))return;const double posterior=std::exp(mass-log_z);
      Vec f;features(x,f);for(int k=0;k<K;++k)expected[k]+=posterior*f[k];
    };
    const auto outside=sc::engine::run_regenerated_outside<Sum>(schedule,inside,root,[&](const Item &i){return generate(i);},
      weight,sc::engine::RuntimeDeductionChildItem{},observe,n);
    const auto bpp=sc::engine::collect_outside_marginals_by_item(inside,[&](const Item &item,double value)->std::optional<double>{
      const double mass=outside.get(item)+value-log_z;if(!std::isfinite(mass))return std::nullopt;return std::exp(mass);},n,*baseline.context);
    require(bpp.ok,"invalid BPP");std::fill(probability.begin(),probability.end(),0.);
    for(auto p:bpp.pairs)probability[p.i*(n+1)+p.j]=p.union_probability;
    ++version;
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"{\"id\":"<<std::quoted(id)<<",\"version\":"<<version<<",\"logZ\":"<<log_z<<",\"RT\":"<<rt<<",\"seconds\":"<<seconds;
    if(with_features){std::cout<<",\"expected_feature\":";vector_json(expected);}
    std::cout<<",\"bpp\":[";bool first=true;
    for(int i=1;i<=n;++i)for(int j=i+1;j<=n;++j)if(double p=probability[i*(n+1)+j];p>0.){
      if(!first)std::cout<<',';first=false;std::cout<<'['<<i<<','<<j<<','<<p<<']';}std::cout<<"]}";
  }
  void vjp(const std::string &id,std::size_t expected_version){
    require(version&&version==expected_version,"stale/missing BPP evaluation");
    const auto start=std::chrono::steady_clock::now();double expected_a=0.;
    // Fixed scaffold pairs are constants (probability one) and are never
    // generated by a deduction; their adjoints cannot change the gradient.
    const auto &ctx=*baseline.context;
    for(int i=1;i<=n;++i)for(int j=i+1;j<=n;++j){
      auto &a=adjoint[i*(n+1)+j];if(!ctx.is_unpaired(i)||!ctx.is_unpaired(j))a=0.;
      expected_a+=probability[i*(n+1)+j]*a;
    }
    // Cov(A,F_k) = sum_d f_dk P(d) (E[A | d] - E[A]): scalar-moment inside and
    // outside, with features read only where a deduction carries DP09 terms.
    const auto gen=[&](const Item &i){return generate(i);};
    const auto loc=[&](const Deduction &d){return local(d);};
    const auto inside=sc::engine::run_schedule<ScalarMoment>(schedule,gen,loc);
    const auto value=inside.get(root);require(std::abs(value.log_weight-log_z)<1e-8,"moment partition mismatch");
    require(std::abs(value.a-expected_a)<1e-8*std::max(1.,std::abs(expected_a)),"observable ownership mismatch");
    Vec covariance{},gradient{},mean_feature{};
    const auto observe=[&](const Deduction &d,const Local &beta){
      Vec f;features(d,f);bool any=false;for(double x:f)any|=x!=0.;if(!any)return;
      Local total=ScalarMoment::times(beta,local(d));
      for(std::size_t c=0;c<d.child_count;++c)total=ScalarMoment::times(total,inside.get(sc::engine::RuntimeDeductionChildItem{}(d,c)));
      if(ScalarMoment::is_zero(total))return;
      const double posterior=std::exp(total.log_weight-log_z),centered=total.a-value.a;
      for(int k=0;k<K;++k){covariance[k]+=posterior*centered*f[k];mean_feature[k]+=posterior*f[k];}
    };
    sc::engine::run_regenerated_outside<ScalarMoment>(schedule,inside,root,gen,loc,sc::engine::RuntimeDeductionChildItem{},observe,n);
    for(int k=0;k<K;++k){
      if(has_expected)require(std::abs(mean_feature[k]-expected[k])<1e-7*std::max(1.,std::abs(expected[k])),"feature expectation mismatch");
      gradient[k]=-covariance[k]/rt;
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    std::cout<<"{\"id\":"<<std::quoted(id)<<",\"version\":"<<version<<",\"seconds\":"<<seconds<<",\"covariance\":";vector_json(covariance);
    std::cout<<",\"gradient\":";vector_json(gradient);std::cout<<'}';
  }
  void target(const std::string &id, const std::string &structure) {
    require(version>0,"evaluate first");
    std::string error;std::vector<sc::api::BasePair> pairs;
    require(structure.size()==sequence.size() && sc::api::detail::parse_dot_bracket_pairs(structure,&pairs,&error),"invalid target");
    bool compatible=true;
    for(int i=1;i<=n;++i){int j=baseline.context->pair_partner(i);
      if(j>i && std::none_of(pairs.begin(),pairs.end(),[&](const auto&p){return p.i==i&&p.j==j;}))compatible=false;}
    double target_z=-INFINITY;Vec target_f{},grad{};
    if(compatible){
      const sc::engine::GeneratedTargetIndex target(sc::engine::generated_traceback_pairs(sc::api::detail::to_traceback_pairs(pairs),*baseline.context));
      const auto gen=[&](const Item &i){auto all=generate(i);std::vector<Deduction> out;
        for(const auto&d:all)if(target.exactly_partitions(d,*baseline.context))out.push_back(d);return out;};
      const auto weight=[&](const Deduction&d){return sc::weights::LogWeight{weight_at(d)};};
      const auto inside=sc::engine::run_schedule<Sum>(schedule,gen,weight);target_z=inside.get(root);
      if(std::isfinite(target_z)){
        const auto observe=[&](const Deduction&d,double beta){double mass=beta+weight_at(d);
          for(std::size_t c=0;c<d.child_count;++c)mass+=inside.get(sc::engine::RuntimeDeductionChildItem{}(d,c));
          if(!std::isfinite(mass))return;Vec f;features(d,f);double p=std::exp(mass-target_z);
          for(int k=0;k<K;++k)target_f[k]+=p*f[k];};
        sc::engine::run_regenerated_outside<Sum>(schedule,inside,root,gen,weight,sc::engine::RuntimeDeductionChildItem{},observe,n);
        require(has_expected,"target likelihood requires E evaluation");
        for(int k=0;k<K;++k)grad[k]=(expected[k]-target_f[k])/rt;
      }
    }
    std::cout<<"{\"id\":"<<std::quoted(id)<<",\"reachable\":"<<(std::isfinite(target_z)?"true":"false")<<",\"log_probability\":";
    if(std::isfinite(target_z))std::cout<<target_z-log_z;else std::cout<<"null";
    std::cout<<",\"gradient_log_probability\":";vector_json(grad);std::cout<<'}';
  }
  void truth(const std::string &id,int k,double step){
    ex::Snapshot changed(sequence,g,stepped(k,step));
    const auto inside=sc::engine::run_schedule<Sum>(schedule,[&](const Item &i){return generate(i);},
      [&](const Deduction &d){return sc::weights::LogWeight{changed.local->log_weight_for(d)+pair_score(d)};});
    std::cout<<"{\"id\":"<<std::quoted(id)<<",\"coordinate\":"<<std::quoted(NAMES[k])<<",\"step\":"<<step<<",\"logZ\":"<<inside.get(root)<<'}';
  }
};

int main(){
  std::cout<<std::setprecision(17);std::map<std::string,std::unique_ptr<Case>> cases;std::string line;
  while(std::getline(std::cin,line)){if(line.empty())continue;std::istringstream in(line);std::string command,id,extra;
    std::ostringstream buffer;auto *saved=std::cout.rdbuf(buffer.rdbuf());
    try{
      require(bool(in>>command>>id),"missing command/id");
      if(command=="M"){
        require(!(in>>extra),"extra metadata fields");
        sc::weights::RuntimeSnapshot owner("AAAAA",".....");
        const double rt=owner.partition_free_energy_from_log_weight(-1.);
#ifdef PK_FIXED_CF_BRANCH
        const char *backend="contrafold-pk-v2";
        const double branch=pkfit::multi_paired_score();
#else
        const char *backend="contrafold-pk-legacy-branch";
        const double branch=-0.03/rt;
#endif
        std::cout<<"{\"id\":"<<std::quoted(id)<<",\"backend\":"<<std::quoted(backend)
          <<",\"RT\":"<<rt<<",\"ordinary_branch_log_score\":"<<branch
          <<",\"ordinary_branch_energy_equivalent\":"<<-rt*branch<<'}';
      }else if(command=="R"){
        std::string seq,g;require(bool(in>>seq>>g)&&!cases.count(id),"invalid registration");require(!(in>>extra),"extra fields");
        cases[id]=std::make_unique<Case>(seq,g);std::cout<<"{\"id\":"<<std::quoted(id)<<",\"registered\":true}";
      }else if(command=="E"||command=="B"){
        require(cases.count(id),"unknown case");Vec d;for(double &x:d)require(bool(in>>x)&&std::isfinite(x),"invalid displacement");
        require(!(in>>extra),"extra fields");cases[id]->evaluate(id,d,command=="E");
      }else if(command=="V"){
        require(cases.count(id),"unknown case");auto &c=*cases[id];std::size_t version;int count;
        require(bool(in>>version>>count)&&count>=0,"invalid VJP header");std::fill(c.adjoint.begin(),c.adjoint.end(),0.);
        for(int k=0;k<count;++k){int i,j;double a;require(bool(in>>i>>j>>a)&&1<=i&&i<j&&j<=c.n&&std::isfinite(a),"invalid adjoint");
          auto &slot=c.adjoint[i*(c.n+1)+j];require(slot==0.,"duplicate adjoint");slot=a;}
        require(!(in>>extra),"extra fields");c.vjp(id,version);
      }else if(command=="L"){
        require(cases.count(id),"unknown case");std::string structure;require(bool(in>>structure)&&!(in>>extra),"invalid target request");cases[id]->target(id,structure);
      }else if(command=="D"){
        require(cases.erase(id)==1,"unknown case");std::cout<<"{\"id\":"<<std::quoted(id)<<",\"released\":true}";
      }else{
        require(command=="T"&&cases.count(id),"unknown command/case");int k;double step;
        require(bool(in>>k>>step)&&0<=k&&k<K,"invalid truth request");cases[id]->truth(id,k,step);
      }
      std::cout.rdbuf(saved);std::cout<<buffer.str()<<std::endl;
    }catch(const std::exception &error){std::cout.rdbuf(saved);std::cout<<"{\"id\":"<<std::quoted(id)<<",\"error\":"<<std::quoted(error.what())<<'}'<<std::endl;}
  }
}
