# deductions

Own the shared deduction generator:

- item/key definitions
- candidate expansion from a parent state/span/context
- child extraction for each applicable rule

This layer should be independent of semiring choice.

## VP/WMB geometry contract

The scientific source is CParty's [supplementary information](https://doi.org/10.1093/bioinformatics/btae748),
§§1.2--1.4, eqs. (ii)--(v). The implementation also preserves explicit
patched-legacy support conditions. The table distinguishes those conditions
from the printed recurrence; reorganizing them is not a clean-room claim.

Coordinates are inclusive and 1-based. Write `c_i=parent_index(i)` and
`c_j=parent_index(j)`. Missing cover/border values are negative; root cover
is zero. A normal free endpoint has partner `-2`, while `-1` forbids pairing.
The published `b'` is the API's `border_bp`, and `B'` is `border_Bp`.

| Rule/view | Geometry and current admission | Expansion |
| --- | --- | --- |
| VP WI Case 1 | `c_i>0`, `0<=c_j<c_i`, `Bp,B>=0`, `bp<0` | WI intervals `[i+1,Bp-1]`, `[B+1,j-1]` |
| VP WI Case 2 | `c_j>0`, `0<=c_i<c_j`, `b,bp>=0`, `Bp<0` | WI intervals `[i+1,b-1]`, `[bp+1,j-1]` |
| VP WI Case 3 | `c_i,c_j>0`, all four borders nonnegative | WI intervals `[i+1,Bp-1]`, `[B+1,b-1]`, `[bp+1,j-1]` |
| VP internal loop, SI (v.5) | `i<k<l<j`; the two flanks are empty in G; `(k,l)` is pairable; `l>max(bp,B)`; exclude the zero-gap stack | One VP child `(k,l)` |
| Partial WMBP/WMBW region | `j-i-1>TURN`; each endpoint is ordinary free or has its fixed partner inside `[i,j]` | Interior G pairs may still leave the interval |
| Completed WMB | Same endpoint/span gate; additionally weak closure under `RequireWeakClosure` | Open partial prefixes remain available; `LegacyCompatible` omits the completed-region closure check |
| WMBW join | `l` is non-fixed; `parent_index(l)=parent_index(j)>=0`; `pair_partner(j)<j` | WMBP `(i,l)` and WI `(l+1,j)` |
| Shared BE/VP WMBP join | Existing border/parent-anchor checks **and** `BE.outer_open >= WMBP.i` | The last condition preserves oracle fill availability and prevents overlapping V ownership; it is not inferred solely from SI geometry |

`vp_geometry.hh` observes the WI border partition without weights or chart
ownership. The printed SI describes the first three cases in terms of covers;
the more specific parent/border tests above are the retained implementation
contract. Case 3 gains no additional parent-inequality guard in this cleanup.

For VP internal loops, define the unpaired gaps `u_L=k-i-1` and `u_R=j-l-1`.
The existing implementation search budget is **`u_L+u_R<MAXLOOP`**, not an
inclusive limit; the SI itself does not impose this implementation cap.
The interval view derives the lowest `l` from
`u_R <= MAXLOOP-1-u_L` and the strict right-border inequality. For the left
edge it preserves the current observation-context convention: apply
`k<min(Bp,b)` when both borders are positive; otherwise the loop budget and
span bound apply. The outer VP admission/generator remains the sole support
surface. This does not redefine absent borders or broaden candidate support.

Candidates retain ascending `k`, descending `l` order. This preserves sum
order and historical traceback ties. The old final `k` iteration with no
feasible `l` is omitted by the gap budget. No local weights change.

Verification includes 675 rectangle-enumeration comparisons with missing,
zero and positive borders and blocked flanks in
`scfg2_vp_internal_border_regression_test`, both fixed VP boundary targets,
and the independently measured long-stem BE-availability/ownership witnesses
in `scfg2_wmbp_be_availability_test`. These are finite regression checks,
not a global proof of the published grammar.
