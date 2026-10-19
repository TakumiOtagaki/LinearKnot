#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/api/beam_detail.hh"
#include "scfg2/api/detail.hh"
#include "scfg2/engine/chosen_pair_output.hh"

#include <cstdlib>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cparty::scfg2;
void require(bool ok, const std::string &message) {
  if (!ok) { std::cerr << message << '\n'; std::exit(EXIT_FAILURE); }
}
bool same(const std::vector<api::BasePair> &a, const std::vector<api::BasePair> &b) {
  if (a.size() != b.size()) return false;
  for (std::size_t i = 0; i < a.size(); ++i)
    if (a[i].i != b[i].i || a[i].j != b[i].j) return false;
  return true;
}
int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", SCFG2_TEST_PARAM_FILE, 1);
  std::mt19937 rng(20261014);
  std::size_t ordering_checks = 0, api_checks = 0;
  for (int n : {0,1,10,100,100000}) {
    for (int count : {0,1,2000}) {
      if (!n && count) continue;
      std::vector<engine::TracebackPair> pairs;
      for (int i = 0; i < count; ++i)
        pairs.push_back({int(rng() % n) + 1, int(rng() % n) + 1,
            i % 2 ? engine::TracebackPairLane::Round : engine::TracebackPairLane::Square});
      require(same(api::detail::to_base_pairs(pairs), api::detail::to_base_pairs_linear(pairs,n)),
              "coordinate pair ordering differs from independent comparison sort");
      ++ordering_checks;
    }
  }
  bool invalid = false;
  try { api::detail::to_base_pairs_linear({{0,2,engine::TracebackPairLane::Round}},3); }
  catch (const std::invalid_argument &) { invalid = true; }
  require(invalid, "out-of-domain pair coordinate accepted by linear output contract");
  const std::vector<std::pair<std::string,std::string>> inputs{
      {"GCGCGCGCGCGCGC", ".............."},
      {"GCGCGCGCGCGCGC", "(.(......)...)"},
      {"GGUUUUCCGGGUACCUGGUGGCAGGUCACGUUUAUGCCGGCUGCUUUAUU",
       "....(...)((....)..)((.(.((...))))(.(....)((...))))"}};
  for (auto [seq,g] : inputs) for (std::size_t b : {3,20}) {
    if (b > seq.size()) continue;
    for (auto ranking : {api::BeamRanking::Baseline,api::BeamRanking::Outer})
      for (auto exponent : {api::OuterPilotExponent::Half,api::OuterPilotExponent::TwoThirds})
      for (auto policy : {WmbClosurePolicy::RequireWeakClosure,WmbClosurePolicy::LegacyCompatible}) {
        if (ranking == api::BeamRanking::Baseline && exponent == api::OuterPilotExponent::Half) continue;
        api::BeamConfig config;
        config.beam_size = b; config.ranking = ranking; config.wmb_closure_policy = policy;
        config.outer_pilot_exponent = exponent;
        weights::RuntimeSnapshot owner(seq, g);
        replay::WFinalLogWeightContext ctx(owner,policy);
        const auto prepared = api::detail::prepare_beam_ranking(ctx,seq.size(),config);
        auto old = replay::run_w_final_vp_wmbp_online_viterbi_with_semiring<semirings::LogMaxSemiring>(
            ctx,seq.size(),prepared.widths,prepared.ranking);
        const deductions::ItemKey root{grammar::NonTerminal::W,1,int(seq.size()),-1,-1};
        const auto expected_weight = old.dp.chart.get(root);
        auto tree = engine::build_viterbi_traceback_result(root,
            engine::ViterbiRunResult<double>{std::move(old.dp.chart),std::move(old.choices)},ctx);
        const auto result = api::beam_viterbi(seq,g,config);
        require(result.ok && result.reachable,"public MFE failed: " + result.error);
        require(result.log_weight == expected_weight &&
                same(result.generated_pairs,api::detail::to_base_pairs(tree.generated_pairs)) &&
                result.generated_structure == api::detail::render_generated_structure(seq.size(),tree.pairs,ctx),
                "choice-DAG public MFE differs from old tree output");
        ++api_checks;
      }
  }
  std::cout << "scfg2_mfe_output=ok ordering=" << ordering_checks << " api=" << api_checks << '\n';
}
