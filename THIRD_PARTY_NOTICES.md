# Third-party notices

## CParty

LinearKnot derives from [CParty](https://github.com/TheCOBRALab/CParty).
The existing GPL version 3 text is preserved in [LICENSE](LICENSE).
The original CParty authors and contributors retain their rights in their
contributions. This notice does not replace or amend their license.

Scientific reference: Mateo Gray, Luke Trinity, Ulrike Stege, Yann Ponty,
Sebastian Will, and Hosna Jabbari. CParty: hierarchically constrained
partition function of RNA pseudoknots. *Bioinformatics* 41(1), btae748 (2025).
<https://doi.org/10.1093/bioinformatics/btae748>

## ViennaRNA

`src/ViennaRNA/` contains ViennaRNA code and thermodynamic parameter tables
inherited through CParty. Credit belongs to the ViennaRNA authors and the
Institute for Theoretical Chemistry, University of Vienna.

The upstream [COPYING](src/ViennaRNA/COPYING) and
[AUTHORS](src/ViennaRNA/AUTHORS) files are included verbatim from
ViennaRNA v2.4.17, commit `5d0cf4fb5e794a9859c893fad2a7f34ac64cc403`.
Representative bundled files (`eval.c`, `model.c`, `params/default.c`, and
`utils/utils.c`) are byte-identical to that release. This comparison does
not establish an exact version for every file in the inherited subtree.

ViennaRNA has its own terms, including a restriction on redistribution for
a fee beyond media costs and a request to contact the authors for inclusion
in a commercial product. Read the full COPYING text. The repository's GPL
text does not relicense the bundled ViennaRNA files.

Scientific reference: Ronny Lorenz et al. ViennaRNA Package 2.0.
*Algorithms for Molecular Biology* 6, 26 (2011).
<https://doi.org/10.1186/1748-7188-6-26>

## Component origin map

| Component | Recorded origin and local changes | Terms / scope |
| --- | --- | --- |
| SCFG2 grammar and geometry | Source-informed re-expression of CParty's model and deduction system using a shared executor, explicit deductions, semirings and static scaffold geometry. VP/WMB support conditions are specified in `src/scfg2/deductions/README.md`, including retained compatibility conditions. | The current root LICENSE remains in place. This description does not certify independent authorship or grant an MIT license. |
| `src/ViennaRNA/` | ViennaRNA inherited through CParty, including local changes to PF raw/log parameter construction and parameter-file envelope validation. | ViennaRNA's COPYING and AUTHORS above; the representative v2.4.17 match does not certify every file. |
| `src/scfg2/weights/log_loops.hh` | Log-space translation of the bundled ViennaRNA PF helpers in `loops/hairpin.h`, `loops/internal.h`, `loops/multibranch.h` and `loops/external.h`. | Vienna-derived implementation outside the vendor directory; its relocation/translation does not establish a new license grant. See the ViennaRNA terms above. |
| `src/scfg2/weights/internal_loop.hh`, runtime and pair-table adapters | Local long-loop extrapolation and backend adapters using ViennaRNA's model tables/API. The small-loop reference delegates to Vienna's helper; the log consumer uses the translated helper above. | Keep the Vienna/model provenance in scope. These wrappers are not a replacement thermodynamic backend or an independently cleared MIT component. |
| `cmake/BundledViennaRNA.cmake` | Locally organized build module preserving the inherited Vienna source inventory and exported static-library contract. | Build refactoring does not change the terms of the compiled sources or remove their attribution. |
| `src/scfg2/weights/pseudoknot_parameters.hh` and `params/*.par` | Retained CParty/DP09 parameterization, band multipliers and inherited parameter files. | Model/data origin is separate from authorship of the executor; the table provenance is recorded below. |
| `test/oracle/legacy_cparty/patches/` | Local corrections with small original CParty context hunks. The checked-out oracle is external to production. | Retain upstream context attribution; these patches are excluded from the Python sdist/wheel by the current explicit packaging rules. |

## Parameter-file provenance

The six top-level parameter files were inherited byte-for-byte from
[CParty at `800e95a`](https://github.com/TheCOBRALab/CParty/tree/800e95af9d1a97555233e2c7b37d73dbaa1a92b4/params).
The SHA-256 values identify the exact inherited data. This records the
immediate source, not a claim that CParty authored every table or that each
original dataset's redistribution terms have been independently verified.

| File under `params/` | SHA-256 |
| --- | --- |
| `dna_Matthews04.par` | `71e370178328e563ce8845c1477a369accf4a4bb5af96936978f25e1aa6c6ffb` |
| `rna_CaoChen06.par` | `b3d69be194674246e518a0ccf948ee427c8b3e7947d773437b7adda0bceef404` |
| `rna_CaoChen09.par` | `5253c6f76df3d287831b502b038eda4aa2545c323112b2bc101fae8c7e638eab` |
| `rna_DirksPierce03.par` | `b3d69be194674246e518a0ccf948ee427c8b3e7947d773437b7adda0bceef404` |
| `rna_DirksPierce09.par` | `46e750b6b9f09ebd59cf2f65b388668f4cc2cce6726f2f151c65b50e06a4f6d6` |
| `rna_Turner04.par` | `0c671a7455f5686bc160fc655194b11888594ab7d506626cdad1470bcadd17d9` |

Top-level `params/` is not included in the Python sdist/wheel. Native commands
still incorporate the bundled ViennaRNA default tables; absence of those
standalone files does not remove the thermodynamic data/backend boundary.

## Distribution review status

Adding these notices restores attribution and includes the upstream terms;
it does not establish compatibility between the CParty GPL and ViennaRNA
terms. Additional permissions or exceptions have not been verified. The
inherited `GPL-3.0-only` package metadata must be reviewed for the combined
distribution before public release. Public publication is deferred.
