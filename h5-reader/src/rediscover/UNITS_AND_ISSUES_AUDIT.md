# Units And Issues Audit

Status: investigation/report only. No code fixes, no commit, no merge. Sources are current
docs/code/comments only; no git archaeology. The IUPAC/topology episode is not investigated
beyond current-file boundary notes, per `/tmp/codex-units-issues-investigation-brief.md:42-46`.

## Executive Verdict

The ring-current and McConnell "literature-fixed" de-circularisation rows are not
comparing matched ppm predictions to DFT. They are comparing producer bare kernels
against ORCA shielding targets.

- DFT target: ORCA `CHEMICAL SHIELDINGS` total shielding tensor, in ppm, stored as
  `DftAtomShielding`/`dft_sigma_iso`; this is sigma, not chemical shift delta
  (`src/io/OrcaShieldingParser.cpp:72-108`, `src/model/DftShielding.h:20-35`,
  `src/rediscover/ExtractionSupport.cpp:45-67`, `src/rediscover/RecordSink.cpp:120-128`).
- Ring current: producer `bs_shielding` is a unit-current Biot-Savart geometric
  kernel with `PPM_FACTOR` and shielding sign baked in, units `ppm_T_per_nA`, not
  literature-scaled ppm. A ring-current coefficient/intensity must multiply it to
  produce ppm (`../GEOMETRIC_KERNEL_CATALOGUE.md:18-31`,
  `../GEOMETRIC_KERNEL_CATALOGUE.md:40-44`, `../doc/calculators/biotsavart.tex:24-30`,
  `../doc/calculators/biotsavart.tex:186-190`, `../src/BiotSavartResult.cpp:228-239`,
  `../src/BiotSavartResult.cpp:289-332`).
- McConnell: producer `mc_shielding` is the full asymmetric `M/r^3` kernel in
  `Angstrom^-3`, not ppm and not multiplied by `Delta_chi`, `mu0/4pi`, molar-unit
  conversion, or any literature coefficient (`../GEOMETRIC_KERNEL_CATALOGUE.md:373-385`,
  `../doc/calculators/mcconnell.tex:23-28`, `../doc/calculators/mcconnell.tex:215-235`,
  `../doc/calculators/mcconnell.tex:261-273`, `../src/McConnellResult.cpp:26-38`,
  `../src/McConnellResult.cpp:84-97`, `../src/McConnellResult.cpp:196-284`,
  `../src/McConnellShieldingTimeSeriesTrajectoryResult.cpp:101-106`).

Therefore:

- Ring T0 gamma = -11.3 is not evidence that a fixed-literature ppm kernel failed.
  It is intensity-like, with units roughly `nA/T`, and its negative sign is consistent
  with the documented diamagnetic ring-current sign convention (`../GEOMETRIC_KERNEL_CATALOGUE.md:158-171`,
  `../spec/MATHS_GOALS.md:155-162`, `../spec/CONSTITUTION.md:243-254`).
- McConnell T0 gamma = -4.75 is not dimensionless. It is a fitted conversion from
  `Angstrom^-3` geometry to ppm shielding, absorbing susceptibility and unit/sign
  prefactors. The sign belongs to that omitted prefactor/convention, not to a
  sigma-vs-delta flip (`../GEOMETRIC_KERNEL_CATALOGUE.md:373-385`,
  `../GEOMETRIC_KERNEL_CATALOGUE.md:914-966`).
- A sigma-vs-delta sign convention does not explain the negative gammas in the
  current pipeline. The parsed target is shielding sigma, named/stored as such, and
  no reference subtraction or chemical-shift delta transform is present
  (`src/io/OrcaShieldingParser.cpp:72-108`, `src/model/DftShielding.h:1-10`,
  `src/rediscover/ExtractionSupport.cpp:56-58`, `src/rediscover/RecordSink.cpp:123`).

## Units Map

| Item | Actual unit/convention | Consumer metadata/status |
|---|---|---|
| ORCA target total/dia/para | ppm shielding sigma; total shielding tensor is parsed from ORCA and decomposed | Correct in target sidecars/CSV (`../python/nmr_extract/_catalog.py:199-222`, `src/io/QtFieldCatalog.gen.h:254-256`, `src/io/OrcaShieldingParser.cpp:94-108`) |
| `bs_shielding` producer kernel | Unit-current BS kernel, `ppm_T_per_nA`, sign `sigma_ab = -dB_sec/dB0`; multiply by ring-current coefficient to get ppm | Correct in `_catalog.py` and generated field catalog (`../python/nmr_extract/_catalog.py:148-156`, `src/io/QtFieldCatalog.gen.h:198-201`), correct reader comment (`src/model/QtBiotSavartGroup.h:6-10`), but producer H5 group attr incorrectly writes `units="ppm"` (`../src/BsShieldingTimeSeriesTrajectoryResult.cpp:138-143`) |
| Rediscover `kernel_bs` catalog | Should be `ppm_T_per_nA` | Incorrectly labeled `ppm` in local rediscover catalog (`src/rediscover/Catalog.cpp:67-72`) |
| Ring source `jb_unit_*` | Unit-current local Johnson-Bovey kernel, `ppm_T_per_nA` | Correctly labeled in CSV schema (`src/rediscover/RingCurrentNeighborhood.cpp:80-86`) |
| Ring source `jb_T*` | Literature-scaled ppm because `ring.LiteratureIntensity()` is multiplied | Correctly labeled ppm in CSV schema and code (`src/rediscover/RingCurrentNeighborhood.cpp:87-92`, `src/rediscover/RingCurrentNeighborhood.cpp:260-273`) |
| Rediscover ring `bare_T0/T1/T2` and bare sidecar | The H5 `bsShielding()` bare kernel, still `ppm_T_per_nA` | `_catalog.py` sidecar units are correct (`../python/nmr_extract/_catalog.py:203-210`), but generic CSV `BareKernelColumns()` labels all bare components as ppm (`src/rediscover/ExtractionSupport.cpp:91-98`), and `literature_fixed_decirc.py` reads these as fixed-literature (`src/rediscover/analysis/literature_fixed_decirc.py:32-39`, `src/rediscover/analysis/literature_fixed_decirc.py:381-387`) |
| `mc_shielding` producer kernel | Full `M/r^3`, `Angstrom^-3`, T0+T1+T2, no chi/prefactor | Correct in producer writer, generated catalog, and reader comment (`../src/McConnellShieldingTimeSeriesTrajectoryResult.cpp:101-106`, `src/io/QtFieldCatalog.gen.h:213-215`, `src/model/QtMcConnellGroup.h:1-8`) |
| Rediscover `kernel_mc` catalog | Should be `Angstrom^-3` | Incorrectly labeled `ppm` in local rediscover catalog (`src/rediscover/Catalog.cpp:70-72`) |
| Rediscover McConnell `bare_T0/T1/T2` and bare sidecar | The H5 `mcShielding()` bare kernel, `Angstrom^-3` | `_catalog.py` sidecar units are correct (`../python/nmr_extract/_catalog.py:215-222`), but generic CSV bare columns label ppm (`src/rediscover/ExtractionSupport.cpp:91-98`), and `literature_fixed_decirc.py` calls it "fixed literature-chi emitted kernel" (`src/rediscover/analysis/literature_fixed_decirc.py:40-46`) |
| HM/ringchi | HM `Angstrom^-1`; ring susceptibility `Angstrom^-3`; both are geometric kernels until calibrated | Consistent in `_catalog.py` and calculator docs (`../python/nmr_extract/_catalog.py:160-186`, `../doc/calculators/haighmallion.tex:25-29`) |
| APBS E/EFG | APBS reader/catalog pin E to `V/A`, EFG to `V/A^2` | Correct current convention (`src/model/QtApbsGroup.h:1-10`, `src/io/QtFieldCatalog.gen.h:252-253`, `../python/nmr_extract/_catalog.py:360-364`) |
| Substrate charge multipole field/gradient | Conventions doc describes raw charge-derived substrate fields as `e/Angstrom^2` and gradients `e/Angstrom^3` | Historical Q6 flagged this against APBS V-units (`spec/substrate_conventions_2026-05-30.md:46-48`, `src/rediscover/SURFACE_DESIGN.md:186-189`). Current state distinguishes APBS V-units from FF14SB broad local field `e/Angstrom^2` (`../python/nmr_extract/_catalog.py:230-235`) |
| EFG feature | APBS EFG T2 in `V/A^2`; DFT target T2 in ppm | Correct in `_catalog.py` (`../python/nmr_extract/_catalog.py:237-245`) |
| Buckingham E-field feature/T1 | APBS E-field in `V/A`; T1 target emitted in ppm but convention-unverified | Correctly flagged in `_catalog.py` (`../python/nmr_extract/_catalog.py:247-256`) |

## Gamma Trace

`literature_fixed_decirc.py` claims it reads "fixed-literature" kernels and reports
dimensionless gamma (`src/rediscover/analysis/literature_fixed_decirc.py:1-15`). It
then selects:

- ring T0: `bare_T0` from `ring_current_aggregated.csv`
  (`src/rediscover/analysis/literature_fixed_decirc.py:32-39`,
  `src/rediscover/analysis/literature_fixed_decirc.py:381-387`);
- McConnell T0: `bare_T0` from `mcconnell_aggregated.csv`
  (`src/rediscover/analysis/literature_fixed_decirc.py:40-46`,
  `src/rediscover/analysis/literature_fixed_decirc.py:381-387`);
- T2: `rediscover_*_aggregated_bare_kernel_T2.npy`
  (`src/rediscover/analysis/literature_fixed_decirc.py:402-403`).

Those are producer bare kernels. The generated sidecar metadata even says ring bare
T2 is `ppm_T_per_nA` and McConnell bare T2 is `Angstrom^-3`
(`../python/nmr_extract/_catalog.py:203-222`). The report line saying "matched ppm
units" is therefore false for the columns used
(`src/rediscover/analysis/LITERATURE_FIXED_DECIRC.md:7`).

The reported results match this diagnosis:

- ring T0 gamma `-11.3509` / `-11.3042` is intensity-like, not dimensionless
  (`src/rediscover/analysis/LITERATURE_FIXED_DECIRC.md:17-18`);
- McConnell T0 gamma `-4.7778` / `-4.7514` is a ppm-per-`Angstrom^-3` coefficient,
  not dimensionless (`src/rediscover/analysis/LITERATURE_FIXED_DECIRC.md:21-22`);
- ring T2 gamma near 1 has negligible `gamma_scaled_R2` and does not prove ppm
  units; the code still feeds the bare sidecar (`src/rediscover/analysis/LITERATURE_FIXED_DECIRC.md:19-20`,
  `src/rediscover/analysis/literature_fixed_decirc.py:402-403`);
- McConnell T2 gamma near 0 likewise has near-null tensor score and carries no unit
  verdict (`src/rediscover/analysis/LITERATURE_FIXED_DECIRC.md:23-24`).

The broad-backbone "fixed literature" T2 path has the same naming problem. The sink
writes `*_literature_kernel_T2` sidecars (`src/rediscover/BroadBackboneSink.cpp:139-147`),
but the reducer fills ring/bond from H5 bare kernels or geometric fallback, not from
fixed ppm literature predictions (`src/rediscover/BroadBackbone.cpp:65-77`,
`src/rediscover/BroadBackbone.cpp:79-105`, `src/rediscover/BroadBackbone.cpp:107-127`,
`src/rediscover/BroadBackbone.cpp:468-477`). The calibration script then labels ring
and McConnell sidecars as dimensionless fixed-literature multipliers with coefficient
1 (`src/rediscover/analysis/static_environment_calibration.py:493-523`), which the
current emit does not support.

## Issue History From Docs/Comments

- The original rediscover design intentionally used DFT sigma as the target and
  producer bare `bs`/`mc` kernels as cross-checks, not as already-calibrated
  literature predictions (`src/rediscover/DESIGN.md:57-61`,
  `src/rediscover/DESIGN.md:97-108`).
- The consumer discipline repeatedly says Python reads emitted substrate and does
  not recompute physics; if a fixed-literature quantity is missing, C++ emit must
  be extended (`src/rediscover/GUIDANCE.md:22-34`,
  `src/rediscover/analysis/PATTERNS.md:9-19`,
  `src/rediscover/analysis/PATTERNS.md:42-56`).
- The substrate conventions already warn that kernels are not shielding and
  residuals need explicit baselines (`spec/substrate_conventions_2026-05-30.md:376-392`,
  `spec/substrate_conventions_2026-05-30.md:719-729`).
- Early scalar docs knew the "producer kernel vs DFT" sign was convention-bearing
  and that "coefficient to literature" was still open (`src/rediscover/analysis/FINDINGS.md:46-48`,
  `src/rediscover/analysis/FINDINGS.md:224-237`).
- Later de-circularisation docs/scripts drifted from "bare kernel cross-check" into
  "fixed-literature ppm" language (`src/rediscover/analysis/literature_fixed_decirc.py:32-46`,
  `src/rediscover/analysis/LITERATURE_FIXED_DECIRC.md:7`,
  `src/rediscover/analysis/static_environment_calibration.py:493-523`).

IUPAC/topology note: current docs/code mention typed selectors and "IUPAC trap"
avoidance (`src/rediscover/SURFACE_DESIGN.md:148-157`,
`spec/substrate_conventions_2026-05-30.md:422-425`,
`spec/substrate_conventions_2026-05-30.md:728`). Per brief, this was not
investigated further; no archives were extracted and no archaeology was done.

## Consolidated Open/Resolved Issues

| Issue | Location/history | Status | Impact / recommended resolution |
|---|---|---|---|
| Gamma-units thread | Triggered by ring T0 `-11.3`, McConnell T0 `-4.75` (`/tmp/codex-units-issues-investigation-brief.md:8-14`); actual code reads bare kernels (`src/rediscover/analysis/literature_fixed_decirc.py:381-403`) | Open metadata/analysis bug, physics verdict resolved | Do not call these dimensionless fixed-literature tests. For ring, either use `jb_T*` ppm source columns or multiply H5 BS by `LiteratureIntensity()`. For McConnell, either emit/multiply a named `Delta_chi`+unit prefactor or report gamma with units. |
| BS H5 writer unit attr | Producer H5 attr says `ppm` while producer catalogue/reader say `ppm_T_per_nA` (`../src/BsShieldingTimeSeriesTrajectoryResult.cpp:138-143`, `../python/nmr_extract/_catalog.py:148-156`, `src/model/QtBiotSavartGroup.h:6-10`) | Open | H5 metadata can mislead any non-catalog consumer. Fix writer attr when code changes are allowed. |
| Rediscover `Catalog.cpp` kernel units | `KernelBs` and `KernelMc` both labeled `ppm` (`src/rediscover/Catalog.cpp:67-72`) | Open | Propagates false assumptions into broad/reducer metadata. Should be `ppm_T_per_nA` and `Angstrom^-3`. |
| Generic `BareKernelColumns()` ppm labels | Applies ppm to ring and McConnell bare columns (`src/rediscover/ExtractionSupport.cpp:91-98`) | Open | CSV headers lie for both central mechanisms. Make bare-kernel units mechanism-specific or carry a unit field. |
| Broad "literature_kernel" sidecars | Ring/bond sidecars are H5 bare/fallback geometric kernels (`src/rediscover/BroadBackbone.cpp:468-477`); script says fixed literature coefficient 1 (`src/rediscover/analysis/static_environment_calibration.py:493-523`) | Open naming/metadata bug | Rename to `bare/geometric_kernel_T2` or actually emit fixed ppm literature predictions. |
| Sigma-vs-delta sign | ORCA parser reads shielding tensors; target columns are `dft_sigma_iso` (`src/io/OrcaShieldingParser.cpp:72-108`, `src/rediscover/RecordSink.cpp:123`) | Resolved: not the cause | Do not flip DFT target sign. Chemical-shift delta would be a new, explicit target transform. |
| Ring-current sign convention | `G=-n*B*PPM_FACTOR`; negative diamagnetic intensity gives shielding above ring (`../GEOMETRIC_KERNEL_CATALOGUE.md:158-171`, `../spec/MATHS_GOALS.md:486-490`) | Resolved | Negative ring gamma is expected if fitting the omitted intensity. Keep sign in reports. |
| McConnell tensor formula provenance | Initial applied-maths audit called broad bond tensor home-rolled (`src/rediscover/APPLIED_MATHS_AUDIT.md:24-40`); independent audit resolves it in favor of project-canonical `K * chi_aniso` (`src/rediscover/FIXES_AUDIT_opus.md:18-32`, `src/rediscover/FIXES_AUDIT_opus.md:36-90`) | Resolved formula; coefficient units still open | The tensor expression is not the bug. Unit scaling/name is the bug. |
| McConnell producer-kernel reconstruction gap | FINDINGS reports producer MC kernel only approx 55% reconstructable from scalar `(r, cos theta, category)` (`src/rediscover/analysis/FINDINGS.md:68-80`); STATE still lists gap (`src/rediscover/STATE.md:352-356`) | Partly open | Scalar source summaries do not fully reconstruct full producer tensor. Use full axis/tensor emissions when making tensor claims. |
| T1 convention | H5 time-series attr labels T1 as m-basis while NPY/group views carry Cartesian antisymmetric pseudovector; parity mismatch noted (`src/model/Types.h:512-526`); T1 flagged unverified in conventions/docs (`src/rediscover/SURFACE_DESIGN.md:167-172`, `src/rediscover/STUB_LANGUAGE.md:941-953`) | Open | Do not fit/report T1 physics as verified. Keep T1 emitted for audit only until basis/parity is pinned. |
| T2 basis | Library T2 order `[xy,yz,zz,xz,xx-yy]` pinned in conventions and code (`spec/substrate_conventions_2026-05-30.md:64-125`, `src/rediscover/SphericalBasis.cpp:7-37`, `src/rediscover/analysis/change_of_basis.py:40-57`) | Resolved for rediscover | Continue using frozen `change_of_basis.get_C()`; do not rederive in Python. |
| T2 Cartesian frame | Original caveat said ORCA/H5 frame comparison was unverified (`src/rediscover/DESIGN.md:80-83`, `src/rediscover/GUIDANCE.md:98-103`); STATE records Kabsch check mean `8.9e-5` deg, max `2.4e-4` deg (`src/rediscover/STATE.md:342-350`) | Resolved for ORCA/H5 tensor frame | T2 components are comparable as emitted for that substrate. |
| EFG lab-frame confound | Applied audit flagged EFG lab-frame de-meaning (`src/rediscover/APPLIED_MATHS_AUDIT_codex.md:60-80`); corrected evidence says local-frame emitter removes confound and old O/C artifacts (`src/rediscover/analysis/EFG_ARC_EVIDENCE.md:22-48`, `src/rediscover/analysis/EFG_ARC_EVIDENCE.md:78-104`) | Resolved in corrected EFG rerun; keep caveat for old outputs | Only cite corrected local-frame EFG outputs as evidence. |
| APBS units thread | Substrate conventions say fields/gradients in raw charge units (`spec/substrate_conventions_2026-05-30.md:46-48`); SURFACE_DESIGN Q6 says pin APBS from producer (`src/rediscover/SURFACE_DESIGN.md:186-189`); current APBS is V/A and V/A^2 (`src/model/QtApbsGroup.h:1-10`) | Mostly resolved, with FF14SB/APBS scale caveat | Do not compare APBS and FF14SB fitted coefficients without unit conversion/labels. |
| FF14SB broad field units | Broad local field catalog says `e/Angstrom^2` (`../python/nmr_extract/_catalog.py:230-235`); applied audit says FF14SB field lacks Coulomb prefactor while charge T2 uses it (`src/rediscover/APPLIED_MATHS_AUDIT_codex.md:46-58`) | Open label/scale caveat | Correlations scale-invariant; fitted coefficients and APBS-vs-FF14SB comparisons are not. |
| CPCM solvation caveat | DFT is r2SCAN/def2-SVP CPCM(water), protein-only; fixed-charge/APBS/water predictors are not the same electrostatic treatment (`src/rediscover/analysis/VARIANCE_DECOMPOSITION_METHOD.md:240-266`, `src/rediscover/analysis/VARIANCE_DECOMPOSITION_METHOD.md:270-306`, `src/rediscover/analysis/EFG_ARC_EVIDENCE.md:112-115`) | Open interpretive caveat | Magnetic ring/MC claims are not compromised; electrostatic results must be labeled mismatched-treatment correlations. |
| Thin N / effective N | Docs require atom-count and autocorrelation-deflated N, with thin flags for small strata (`src/rediscover/analysis/VARIANCE_DECOMPOSITION_METHOD.md:127-149`, `src/rediscover/analysis/PATTERNS.md:57-63`) | Open disclosure requirement | Report per-stratum atom counts, N_eff, lag-1, and thin flags; do not bootstrap-inflate 4-atom strata. |
| Variance axis ambiguity | Docs find fitters mixing within-only and pooled/total axes (`src/rediscover/analysis/VARIANCE_DECOMPOSITION_METHOD.md:20-65`); codex audit gives full between/within design (`src/rediscover/APPLIED_MATHS_AUDIT_codex.md:198-256`) | Open/reporting discipline | Every result should state between, within, or total axis and use train-only centering for split metrics. |
| DFT validation | STATE said element equality and raw total ~= dia+para deferred (`src/rediscover/STATE.md:352-356`); current loader checks atom count, parsed element known, and T0 total==dia+para (`src/io/DftShieldingLoader.cpp:78-104`) | Partially resolved | Total identity is now checked; explicit parsed-element vs topology-element equality still not visible in current loader. |
| Unwired/fail-loud stubs | STATE says seven relationships fail loud (`src/rediscover/STATE.md:299-306`); REDISCOVERY_MAP still lists several as stubs (`src/rediscover/REDISCOVERY_MAP.md:33-38`, `src/rediscover/REDISCOVERY_MAP.md:65-73`); later docs show some progressed (EFG corrected evidence) | Mixed/stale docs | Treat scenario availability as runtime `ValidateScenario`/catalog fact, not doc grid fact. Update status docs when code changes are allowed. |
| Missing producer bare-kernel TS for many relationships | STUB_LANGUAGE says only ring, McConnell, and APBS EFG have real cross-check TS; Buckingham/charge/CRG/embedding do not (`src/rediscover/STUB_LANGUAGE.md:954-964`) | Open | Cross-check kernel must remain optional; no placeholder `KernelBs` for unrelated mechanisms. |
| MOPAC per-frame charge source | STUB_LANGUAGE says MOPAC charge per-atom TS absent, should validate-fail rather than emit zeros (`src/rediscover/STUB_LANGUAGE.md:828-839`); local catalog still marks `MopacCharge` absent (`src/rediscover/Catalog.cpp:95-97`) | Open until data lands | Keep fail-loud; no fallback chain. |
| AIMNet2 embedding carrier awkwardness | STUB_LANGUAGE notes embedding is a feature, not a source-slot value (`src/rediscover/STUB_LANGUAGE.md:928-931`) | Open design issue | Needs a non-source carrier or aggregated-row feature layout before fitting. |
| PBC substrate rule | Conventions corrected prior minimum-image guidance; extractor makes protein whole upstream and rediscover uses `PbcMode=None` with sanity assert (`spec/substrate_conventions_2026-05-30.md:660-684`) | Resolved, assert still needed | Do not reopen reader-side PBC unless solvent/cross-box queries enter scope. |
| Differencing system-ID | WORK_CATALOG parks differencing because every-other-frame DFT spacing fails smoothness gate (`src/rediscover/WORK_CATALOG.md:42-58`) | Parked | Needs dense consecutive-frame DFT burst plus more contributors; full 750 with same spacing does not fix it. |
| Documentation staleness | BACKBONE/STATIC reports contain stale "fixed literature" and emitted-sidecar language (`src/rediscover/analysis/BACKBONE_LAW_EVIDENCE.md:115-117`, `src/rediscover/analysis/STATIC_ENVIRONMENT_CALIBRATION.md:35-52`, `src/rediscover/analysis/STATIC_ENVIRONMENT_CALIBRATION.md:65-68`) | Open | Mark affected reports as pre-units-audit or revise wording once fixes are allowed. |

## Minimal Next Fix Set (Not Performed Here)

1. Correct unit metadata: BS H5 attr, rediscover `Catalog.cpp` `KernelBs/KernelMc`,
   and mechanism-specific `BareKernelColumns()`.
2. Rename current `literature_fixed_decirc` inputs/results to bare-kernel scale
   diagnostics, or feed actual literature-scaled ppm arrays.
3. For ring fixed-literature, consume `jb_T*` source columns or emit an aggregated
   literature-scaled ppm kernel.
4. For McConnell fixed-literature, choose and document the `Delta_chi`/unit prefactor
   convention, then emit a ppm prediction or keep gamma explicitly unitful.
5. Update broad-backbone `*_literature_kernel_T2` naming/metadata to distinguish
   bare BS, bare McConnell, fallback geometric McConnell, and charge EFG.
