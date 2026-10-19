#include "scfg2/weights/runtime_snapshot.hh"
#include "scfg2/replay/w_final_vp_wmbp_online.hh"

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <string>

namespace {

constexpr const char *kSequence =
    "CGGCAUACUCGAGGAUUUUAAUCGAUGGGUGUCGCUCUUUCUGGGGCCCUUCGCCUCCUUGGGGGGUCGG"
    "AUAUAAAGGAUUUGUCGUUAUGGACCCGGCUUCUGUAUGAUUGAAAUGACCUCCGAGUGUGAUUCCUCGU"
    "GUACACCGCCCCUUUUAUUGGUAGGUCUUAGGACGCUCUGUGGGCUUGGUAUCUGGGUUCAUGCAAAUCG"
    "GGCUCAGAUGUUCGCCCGCGUGGUUAUCGUGAAAGAUAGCCUUUCGGAUUUGAUCGGACCUCCCACCCGG"
    "GGAGUCAUGCUGUUAAACGCCGUGGGCGACCCUCGAGACGAUUUCCAGGGUGUGUAGGAUUCGUCUAUAA"
    "GUGAAAACACCCAGAAUUAUCUGGCCCAUCGGGUAAGUUUAUGAGGUCCG";

constexpr const char *kScaffold =
    ".((.....))(....)(...(....)(.((...))((((...)))((.....))((......)).))((."
    ".(...(.(...))))(..(...(...))(....))((...)((...)(....).))(.(...)(....(."
    "(....))))(....(((....)))(....(...((...))((....))(..(...).)))((.......)"
    ".(.(.....))(.(....)))).((.....)(....)))(..(..(...).(...(.....))(.(...)"
    ")(...((...))(...))(..(....)((..(..(...)).))((......).))))(..(.....(..."
    ".)(....)(....(.(...))))(.(....)(.(...(.(...)))))))";

bool close(double lhs, double rhs, double tolerance = 1e-12) {
  return std::fabs(lhs - rhs) <= tolerance;
}

}  // namespace

int main() {
  setenv("CPARTY_FIXED_ENERGY_PARAM_FILE", "params/rna_Turner04.par", 1);
  const std::string sequence = kSequence;
  const std::string scaffold = kScaffold;
  assert(sequence.size() == 400);
  assert(scaffold.size() == sequence.size());

  std::string owner_sequence = sequence;
  std::string owner_scaffold = scaffold;
  cparty::scfg2::weights::RuntimeSnapshot owner(
      owner_sequence, owner_scaffold);
  cparty::scfg2::replay::WFinalExactInsideContext ctx(owner);
  const auto n = static_cast<cand_pos_t>(sequence.size());
  const cparty::scfg2::deductions::ItemKey root{
      cparty::scfg2::grammar::NonTerminal::W, 1, n, -1, -1};
  const auto widths =
      cparty::scfg2::replay::WFinalVpWmbpOnlineWidths::uniform(80);

  const auto independent_inside =
      cparty::scfg2::replay::run_w_final_vp_wmbp_online<
          cparty::scfg2::semirings::LogSumExpSemiring>(ctx, n, widths);
  const auto independent_viterbi =
      cparty::scfg2::replay::run_w_final_vp_wmbp_online<
          cparty::scfg2::semirings::LogMaxSemiring>(ctx, n, widths);
  const auto coupled = cparty::scfg2::replay::
      run_w_final_vp_wmbp_online_coupled_viterbi_traceback(ctx, n, widths);

  const double independent_log_z = independent_inside.chart.get(root);
  const double independent_log_v = independent_viterbi.chart.get(root);
  const auto coupled_root = coupled.run.chart.get(root);

  // This historical input no longer violates the independent max/sum bound
  // at baseline 7614910d. Keep the long input and frozen baseline scores as
  // secondary regression controls, not as an energy-model semantic oracle.
  // Independent runs still have no general cross-DAG max <= sum guarantee.
  assert(std::isfinite(independent_log_z));
  assert(std::isfinite(independent_log_v));
  assert(std::isfinite(coupled_root.log_sum));
  assert(std::isfinite(coupled_root.log_max));
  assert(coupled_root.log_max <= coupled_root.log_sum + 1e-12);
  assert(close(coupled_root.log_sum, independent_log_z));
  assert(cparty::scfg2::deductions::ItemKeyEq{}(coupled.root, root));
  assert(coupled.trace.chosen_deduction.has_value());
  assert(close(coupled.trace.chosen_contribution.log_max, coupled_root.log_max));
  assert(close(independent_log_z, -481.5498157235681, 1e-9));
  assert(close(independent_log_v, -484.95038927987594, 1e-9));
  assert(close(coupled_root.log_max, -484.95038927987594, 1e-9));

  std::cout << "scfg2_coupled_log_probability_regression=ok\n";
}
