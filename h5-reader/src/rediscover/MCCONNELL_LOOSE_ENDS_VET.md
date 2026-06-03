# McConnell Loose Ends Vet

Scope: read-only second-opinion pass against code. No fixes, no re-emit, no ORCA, no `trajectory.h5` read. Final scientific call stays with the lead.

C1 note: I agree the structural ceiling is that McConnell-alone is a minority contributor to the full DFT T2; the fair test is joint/residualized, not standalone McConnell-only.

## C6 - QtBond A/B Endpoint Ordering

Verdict: DISPUTED as a consumed-T2 sign bug.

Evidence:

- Producer endpoint creation is index-canonical, not chemically oriented: OpenBabel bonds are stored as `{min(a,b), max(a,b)}` in `/shared/2026Thesis/nmr-shielding/src/CovalentTopology.cpp:45-52`; the fallback detector loops `j = i + 1` and stores `{i,j}` in `/shared/2026Thesis/nmr-shielding/src/CovalentTopology.cpp:114-120`; classification copies those endpoints to `Bond.atom_index_a/b` in `/shared/2026Thesis/nmr-shielding/src/CovalentTopology.cpp:128-133`. Force-added disulfides also use min/max in `/shared/2026Thesis/nmr-shielding/src/CovalentTopology.cpp:344-350`.
- The sidecar preserves those endpoints verbatim: `bonds.npy` writes `atom_index_a/b` from `Bond.atom_index_a/b` in `/shared/2026Thesis/nmr-shielding/src/TopologySidecar.cpp:265-279`; the reader row layout has the same fields in `src/io/QtNpyRecords.h:62-66`; `DecodeBond` copies them to `QtBond.atomIndexA/B` in `src/io/QtTopologySidecar.cpp:254-258`.
- Rediscover builds the bond axis as `posB - posA`: standalone McConnell in `src/rediscover/McConnellNeighborhood.cpp:163-173` and broad backbone in `src/rediscover/BroadBackbone.cpp:320-334`; both emit/store `bond_atom_a/b` and `bond_axis_local` from that axis in `src/rediscover/McConnellNeighborhood.cpp:188-198` and `src/rediscover/BroadBackbone.cpp:349-356`.
- The consumed McConnell tensor is even under endpoint reversal. `McConnellSourceLiteratureKernelLocal` computes `cosTheta = dHat.dot(bHat)` and then `(9*cosTheta*dHat*bHat^T - 3*bHat*bHat^T - (3*dHat*dHat^T - I))/r3` in `src/rediscover/McConnellLiteratureKernel.cpp:65-75`; the broad fallback uses the same structure in `src/rediscover/BroadBackbone.cpp:90-97`. If endpoints swap, `bHat -> -bHat` and `cosTheta -> -cosTheta`, so `cosTheta*dHat*bHat^T` is unchanged, and `bHat*bHat^T` is also unchanged.

Assessment: A/B is stable and index-canonical through the loader. It is not chemically canonical C->O/C->N, so the raw `bond_axis_local_*` vector is arbitrary as a polar vector. But the emitted and consumed `mc_lit_T2_local_*` tensors and `bondKernelT2FromSources` do not flip sign when A/B flips. The audit's term-1 sign-flip claim is over-called for this path.

Warrants a fix? No for the McConnell T2 sign claim. Optional documentation/future-proofing: label `bond_axis_local_*` as index-oriented, or add a chemically oriented axis only if a future odd-vector consumer needs a polar bond direction.

## C2 - Rediscover 8 A Cutoff vs Producer 10 A

Verdict: AGREED, with the 8-10 A magnitude not numerically quantified in this read-only pass.

Evidence:

- Standalone rediscover McConnell defaults to 8.0 A in `src/rediscover/McConnellNeighborhood.h:26-30` and the CLI `--mc-cutoff` default in `src/rediscover/main_extract.cpp:165-170`; `main_extract` passes that value into the composed/procedural McConnell paths in `src/rediscover/main_extract.cpp:504-515`.
- Broad backbone also defaults the anisotropic-bond cutoff to 8.0 A in `src/rediscover/main_extract.cpp:177-185` and passes it to `MakeBroadBackboneRelationship` in `src/rediscover/main_extract.cpp:430-433`, which uses it for `CloudKind::BondMidpoints` in `src/rediscover/BroadBackbone.cpp:447-449`.
- The producer default is 10.0 A: `/shared/2026Thesis/nmr-shielding/src/CalculatorConfig.cpp:49-51`, consumed by `McConnellResult` through `BondsWithinRadius(... mcconnell_bond_anisotropy_cutoff)` in `/shared/2026Thesis/nmr-shielding/src/McConnellResult.cpp:139-140`.
- The producer also applies self/near-field filters in `/shared/2026Thesis/nmr-shielding/src/McConnellResult.cpp:120-125`; the near-field criterion is `distance > near_field_exclusion_ratio * source_extent` in `/shared/2026Thesis/nmr-shielding/src/KernelEvaluationFilter.h:174-179`. Broad backbone has a producer-valid analogue via `mc_source_is_self_or_bonded` in `src/rediscover/BroadBackbone.cpp:357-360` and valid sums in `src/rediscover/BroadBackboneSink.h:111-118`; standalone McConnell still declares valid == all in `src/rediscover/McConnellNeighborhood.cpp:102-104`.
- The tensor/source law is explicitly `1/r^3`: `McConnellLiteratureKernel.cpp:68-75` and `BroadBackbone.cpp:93-97`.

Assessment: The mismatch is real. I did not measure the actual 8-10 A shell, but calling it negligible is not supported by code: it is inside the producer's canonical radius, is beyond the near-field/self-source region, and a summed `1/r^3` field over a growing shell can contribute materially after angular cancellation. The concern is especially valid for comparisons to producer `mc_shielding`, which includes 10 A sources.

Warrants a fix? Yes. Align rediscover McConnell/broad `bond-cutoff` default to 10.0 A, or bind it to the producer config value, keep recording the cutoff, and rerun/report 8 vs 10 as a sensitivity check.

## C5 - Within-Axis De-Mean vs Between-Axis

Verdict: AGREED as a reporting/weighting problem, not as "between cannot be computed."

Evidence:

- The older fixed-literature consumer is within-only by design: it says scalar and tensor comparisons use within-atom centering in `src/rediscover/analysis/literature_fixed_decirc.py:12-14`, centers per atom in `src/rediscover/analysis/literature_fixed_decirc.py:120-136`, and both its frame-split and leave-atoms-out protocols operate on centered arrays in `src/rediscover/analysis/literature_fixed_decirc.py:449-479`.
- The current `mc_lit_*` decirc consumer reads the valid columns by default in `src/rediscover/analysis/mcconnell_literature_decirc.py:58-61` and `src/rediscover/analysis/mcconnell_literature_decirc.py:609-614`. It computes a within row by per-atom de-meaning in `src/rediscover/analysis/mcconnell_literature_decirc.py:118-134`, no-intercept gamma in `src/rediscover/analysis/mcconnell_literature_decirc.py:160-169`, and row append in `src/rediscover/analysis/mcconnell_literature_decirc.py:384-418`.
- The same script does compute between rows: atom means in `src/rediscover/analysis/mcconnell_literature_decirc.py:137-157`, intercept fit in `src/rediscover/analysis/mcconnell_literature_decirc.py:172-189`, and row append in `src/rediscover/analysis/mcconnell_literature_decirc.py:422-458`. The report prints both axes in the T2 lead table in `src/rediscover/analysis/mcconnell_literature_decirc.py:551-571`; the current generated report shows both axes for each stratum in `src/rediscover/analysis/MCCONNELL_LITERATURE_DECIRC.md:17-28`.
- The dchi calibration likewise writes both axes in its lead table in `src/rediscover/analysis/MCCONNELL_DCHI_CALIBRATION.md:17-28` and diagnostics in `src/rediscover/analysis/MCCONNELL_DCHI_CALIBRATION.md:36-47`; it explicitly calls out N-between and O-between in `src/rediscover/analysis/MCCONNELL_DCHI_CALIBRATION.md:68`.
- The method document defines between as the static-environment axis and within as frame-to-frame motion in `src/rediscover/analysis/VARIANCE_DECOMPOSITION_METHOD.md:31-36`, and says bond anisotropy on carbonyl C is plausibly between-led in `src/rediscover/analysis/VARIANCE_DECOMPOSITION_METHOD.md:349-353`.

Assessment: The audit is right that a near-static McConnell mechanism should not be led by a uniform per-atom de-meaned read. The current `mc_lit_*` scripts are better than the old within-only consumer because they compute both axes, so "between is absent" would be over-called. The remaining issue is presentation/default weighting: within is produced first, uses many frame rows, and inherited thin/confidence flags are based on within-style atom modulation even for between rows (`mcconnell_literature_decirc.py:380-408` and `:442-447`; `mcconnell_dchi_calibration.py:541-565` and `:591-596`). That can understate the between read instead of making it the primary static-mechanism diagnostic.

Warrants a fix? Yes. Make McConnell reports axis-explicit with between as the lead/default static read, use `atom_mean_signal_neff` or atom count for between thin flags, and keep within rows as dynamic diagnostics rather than the headline.

## SDK Gap - `mc_lit_*` Columns Missing From `_catalog.py`

Verdict: AGREED.

Evidence:

- `_catalog.py` describes `ArraySpec` as metadata for one NPY file in `/shared/2026Thesis/nmr-shielding/python/nmr_extract/_catalog.py:44-46`.
- The rediscover catalog lists old McConnell standalone sidecars in `/shared/2026Thesis/nmr-shielding/python/nmr_extract/_catalog.py:211-222` and broad-backbone NPYs only for target T2 and field-local in `/shared/2026Thesis/nmr-shielding/python/nmr_extract/_catalog.py:224-235`. There are no `mc_lit_*` specs.
- Broad-backbone explicitly says only the three NPY payloads have ArraySpec entries in `src/rediscover/BroadBackboneSink.h:21-24`.
- The consumed columns are CSV-only today: source `mc_lit_kernel_present`, `mc_lit_T0`, `mc_lit_T2_local_0..4` in `src/rediscover/BroadBackboneSink.cpp:139-146`; aggregate all-source and valid `mc_lit_*` columns in `src/rediscover/BroadBackboneSink.cpp:149-187`.
- The analyses require and consume those CSV columns: `mcconnell_literature_decirc.py` defines `mc_lit_T2_local_*` and valid variants in `src/rediscover/analysis/mcconnell_literature_decirc.py:28-29`, selects valid/all in `src/rediscover/analysis/mcconnell_literature_decirc.py:609-614`, and requires them in `src/rediscover/analysis/mcconnell_literature_decirc.py:621-631`. The dchi calibration source-row consumer similarly defines all/valid T2 column names in `src/rediscover/analysis/mcconnell_dchi_calibration.py:29-31`.

Assessment: The columns are genuinely absent from the SDK catalog. Because `_catalog.py` is currently NPY-oriented, this is either a missing CSV-column schema or a signal that the consumed `mc_lit_*` tensors should be emitted as NPY sidecars. Either way, downstream consumers have no SDK-level units/axis/feature contract for the McConnell literature tensors they now use.

Warrants a fix? Yes. Add a rediscover CSV-column schema for `mc_lit_T0`, `mc_lit_T2_local_0..4`, `mc_lit_valid_kernel_present`, `mc_lit_T0_valid`, and `mc_lit_T2_local_valid_0..4`, or emit aggregate/source `mc_lit` tensors as NPY sidecars and add `ArraySpec` entries with axis, units `ppm`, irreps `2e`, and mechanism `bond_anisotropy`.
