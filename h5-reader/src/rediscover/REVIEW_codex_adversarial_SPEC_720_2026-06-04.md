# Codex adversarial review - SPEC_720_STATIC_INGEST_2026-06-04.md

> **Historical review record — not current truth (trued 2026-06-04).** Session
> provenance only; current truth is the relevant `SPEC_*`, `NOW.md`, and
> corrected `STATE.md`.

Scope: thoughtful-scientist review of the described 720 static-ingest design. I read the spec, prior review, guidance/contract docs, target/basis/reducer code, and the C++/Python producer catalogs. I did not run 720/build/fit, and I did not touch code or git.

## Verdict

**Indicative maths: yes, OK for statistical position, not proof. Basic maths errors: no hard error found in the target decomposition, T2 basis, dia/para handling, additive source reductions, or C++ aggregation arithmetic.**

That verdict is conditional in the normal scientific sense: the 720 read must be reported as indicative between-axis position with protein/block and stratum support, not as row-iid proof. The spec is honest about this: it keeps the between-axis estimator open, gives the accumulator hooks for centering/null/strata/support, and keeps Python to lean fitted substrate only (`SPEC_720_STATIC_INGEST_2026-06-04.md:162`, `:177`, `:187`, `:493`).

The one place I would not let implementation drift is **T2 orientation across proteins**. T0 is rotation-invariant. T2 components are not. A plain pooled scalar fit over lab-frame T2 components from arbitrarily oriented static proteins would be a basic tensor-statistics error. I do not call the design wrong because it leaves room for local-frame/equivariant handling, but the build spec should explicitly require the C++ emit to provide either stable local-frame T2 where valid, equivariant/node-store consumption, or invariant T2 summaries before any between-protein component fit.

## Maths Check

The T0/T2 basis is sound. `DecomposeLibrary` uses `trace/3` for T0, the antisymmetric dual for T1, and the library T2 order `[xy, yz, zz, xz, xx-yy]` with the expected isometric sqrt factors (`SphericalBasis.cpp:7`, `SphericalBasis.cpp:27`). That matches the stated frozen basis bridge, not the viewer-only order (`SphericalBasis.h:1`, `SphericalBasis.h:11`).

The target build is also sound for the stated purpose. `BuildTarget` reads raw ORCA total/dia/para tensors, decomposes all three through `DecomposeLibrary`, and keeps total-local only when a valid local frame exists (`ExtractionSupport.cpp:45`, `ExtractionSupport.cpp:58`, `ExtractionSupport.cpp:63`). The per-atom writer emits total T0/T1/T2 plus dia/para T0/T1/T2 sidecars (`PerAtomSubstrate.cpp:2657`, `PerAtomSubstrate.cpp:2799`). The old logical gap that the loader validated dia+para using T0 only is real (`DftShieldingLoader.cpp:95`), but the documented follow-up found componentwise dia+para=total to ORCA print rounding, including T2 max `1.632993e-3` ppm (`POSTMORTEM_DIAPARA_CHECK_2026-06-04.md:6`, `:14`). So total-T2 is clean, and dia/para is acceptable as diagnostic split.

The C++ reductions are additive where shielding physics is additive. Charge EFG uses `q * (3 rr/r^5 - I/r^3)`, projects traceless, and decomposes in the same basis (`PerAtomSubstrate.cpp:641`, `:646`, `:649`). MOPAC field uses the correct target-minus-source sign via `(-q/r^3) * disp` where `disp = source - target` (`PerAtomSubstrate.cpp:686`, `:691`). McConnell uses bond midpoint/bond axis, rejects endpoint/near-field self sources, and only contributes valid literature kernels (`PerAtomSubstrate.cpp:831`, `:850`, `:852`, `:885`). The geometric H-bond shadow uses the expected dipolar form times `(h h^T - I/3)` before decomposition (`PerAtomSubstrate.cpp:803`, `:808`, `:811`). Dominance fractions are `max(abs contribution)/sum(abs contribution)` and are isolation descriptors, not signed physics claims (`PerAtomSubstrate.cpp:118`).

The row aggregation is mathematically honest as a lean substrate: one row per target atom, sums in C++, support/presence flags in the row/manifest, no default source table (`PerAtomSubstrate.h:1`, `:25`; `PerAtomSubstrate.cpp:2757`, `:2796`; `PerAtomSubstrate.cpp:3388`, `:3421`). For 720 statics, this is a valid indicative substrate if later statistics respect protein clustering and stratum support.

## Relationships At Risk

These are not reasons to abandon the design. They are the minimal C++ emissions I would require so the lean design does not quietly miss relationships the 720 sample could actually show.

1. **Cross-protein T2 orientation.** Current all-atom substrate uses a lab-frame function by default (`PerAtomSubstrate.cpp:246`, `:3837`). Across many static proteins, global lab axes are not a scientific stratum. A real tensor relationship can be washed out by arbitrary protein orientation.

Minimal fix: C++ emit `static_target_local_T2`, `static_feature_local_T2` for atom classes with stable local frames, plus `frame_valid/frame_variant/frame_anchor_atom_index`; for atoms without a stable frame, emit `|T2|`, tensor inner products/norms, or route that family only to an equivariant/node-store consumer. Do not let Python rotate tensors or rebuild frames.

2. **Ring self/through-space split.** The legacy ring extractor explicitly kept both all-ring and producer-valid through-space sums (`RingCurrentNeighborhood.cpp:211`, `:278`). The all-atom aggregate currently sums ring T2 including self/bonded rings while only counting the self/bonded rows (`PerAtomSubstrate.cpp:969`, `:978`), and the classical columns expose only `ring_jb_T0/T2`, not a valid/self split (`PerAtomSubstrate.cpp:2192`, `:2196`). The 720 sample could distinguish own-ring aromatic shielding from through-space ring current.

Minimal fix: C++ emit `ring_jb_valid_T0/T2` and `ring_jb_self_or_bonded_T0/T2`, with counts, using the existing `ringSourceIsSelfOrBonded` rule. This is a scalar/sidecar split, not a pair dump.

3. **Per-type/per-category mechanism structure.** Build history already says pooled-then-sliced was misleading and per-type rescued H/heavy strata (`NOW.md:20`, `:23`). The producer/catalog and current substrate can carry ring per-type T2 and Mc/MOPAC-Mc category T2 (`PerAtomSubstrate.cpp:1555`, `:1589`; `PerAtomSubstrate.cpp:2116`, `:2150`). If the final lean menu collapses these to one total per mechanism, it can miss a relationship that 720 has enough atom rows to show.

Minimal fix: require compact C++ grouped sidecars, at least `ring_{bs,hm,pq,disp}_per_type_T2` where enabled and `mc_category_T2` / `mopac_mc_category_T2`, plus category counts/support. Keep them grouped, not source-expanded.

4. **Electrostatic source distinctions.** The catalog separates FF14SB/MOPAC/EEQ/AIMNet2 charges, APBS E/EFG, MOPAC Coulomb E/EFG/shielding, AIMNet2 EFG/backbone/aromatic, and water E/EFG (`QtFieldCatalog.gen.h:240`, `:253`; `QtFieldCatalog.gen.h:268`, `:274`). These are not interchangeable: solvent PB, semiempirical charge, learned charge, and explicit water can differ in sign/convention/support. Collapsing them to a single "field" magnitude would leave source-comparison relationships on the floor.

Minimal fix: C++ emit separate compact electrostatic slabs: APBS `E` and `EFG_T2`, MOPAC Coulomb `E`, raw EFG backbone/aromatic and shielding T2, AIMNet2 charge/CRG/EFG splits, FF14SB/MOPAC/EEQ charge scalars, and per-source presence/support flags. The current direct-feature path already demonstrates this shape (`PerAtomSubstrate.cpp:1435`, `:1474`; `PerAtomSubstrate.cpp:1598`, `:1619`).

5. **H-bond/solvation/exposure.** The producer exposes geometric H-bond, Larsen per-class H-bond tensors, DSSP energies/SS/chi, water fields/EFGs, hydration shell, SASA and surface normals (`QtFieldCatalog.gen.h:221`, `:235`; `QtFieldCatalog.gen.h:289`, `:296`). These are static-relevant and high-support enough across 720 to show relationships. Treating them as optional is fine for disk, but omitting the whole family would blind the review to solvent/exposure and H-bond signal.

Minimal fix: one C++ compact sidecar for `hbond_T2`, Larsen per-class T2 plus water term, DSSP hbond energy/count/SS flags, SASA, SASA normal, water E/EFG total/first-shell, and hydration shell scalars. The current hbond/conditioning copy path is already C++ (`PerAtomSubstrate.cpp:1621`, `:1665`; `PerAtomSubstrate.cpp:2163`, `:2189`).

6. **AIMNet2 electronic environment.** The 256-d embedding is not necessary for the classical law read, but it is a real electronic-environment family and remains under the stated disk budget even for high atom counts (`SPEC_720_STATIC_INGEST_2026-06-04.md:326`, `:362`). If disabled, the design cannot see whether learned electronic state carries between-axis position.

Minimal fix: keep `static_optional_aimnet2_embedding.npy` as a separable f32 sidecar with `embedding_present` and support flags, exactly as the spec sketches (`SPEC_720_STATIC_INGEST_2026-06-04.md:227`, `:269`). Do not derive embedding features in Python.

I do **not** recommend adding mutation-delta fields as a default fix. The design is explicitly WT absolute sigma, not a mutation-pair delta model (`SPEC_720_STATIC_INGEST_2026-06-04.md:13`, `:133`). Delta surfaces in the catalog are real, but pulling them into this substrate would change the question rather than rescue a relationship inside the chosen static-ingest design.

## Rules Check

The foreclosure/discipline is respected. The spec states no terabyte Python, no second Python protein model, C++ owns geometry/reductions/source association, extractor output is consumed rather than reimplemented, exact-path loading/log-and-stop is required, memory strategy remains open, and Python only consumes lean row-aligned substrate (`SPEC_720_STATIC_INGEST_2026-06-04.md:21`, `:23`, `:91`, `:197`, `:477`). This matches `GUIDANCE.md:17`, `GUIDANCE.md:22`, `PARTITION_FILTER_DESIGN.md:16`, and `NODE_STORE_CONTRACT_2026-06-02.md:24`. The rule check is secondary to the maths, but it passes.

## Summary

**(a)** Indicative-OK + error-free: **yes**, as indicative statistical position and not proof; **no basic maths error found**. Required guard: do not pool arbitrary lab-frame T2 components across proteins without local-frame/equivariant/invariant handling.

**(b)** Real relationships left on the floor + minimal fixes: at risk are cross-protein T2 orientation, ring self vs through-space, per-type/per-category mechanism structure, electrostatic source distinctions, H-bond/solvation/exposure, and optional AIMNet2 electronic environment. Minimal fixes are all C++ emits: local/invariant T2 sidecars; ring valid/self split; grouped ring/mc category T2 sidecars; separated electrostatic slabs; compact hbond/solvation sidecar; separable AIMNet2 embedding sidecar.

**(c)** Rules check: foreclosure holds. The design keeps model/spatial work and reductions in C++, consumes exact producer outputs, prevents Python protein reconstruction and pair dumps, keeps memory strategy/statistics open, and emits a lean row-aligned substrate under a fail-loud disk budget.
