// Research observations use the same pilot/ranking implementation as production.
#ifndef SCFG2_PILOT_GUIDED_INSIDE_HH
#define SCFG2_PILOT_GUIDED_INSIDE_HH
#include "scfg2_bounded_outer.hh"
#include "scfg2/beam/pilot_outer.hh"
using sc::beam::PilotExponent;
using sc::beam::outer_detail::sqrt_pilot_width;
using sc::beam::outer_detail::pilot_exponent_name;
using sc::beam::outer_detail::pilot_width_for;
using sc::beam::outer_detail::PilotRankingCounts;
using sc::beam::outer_detail::PilotOuterRanking;
using sc::beam::outer_detail::PilotGuidedInsideRun;
using sc::beam::outer_detail::run_pilot_guided_inside;
#endif
