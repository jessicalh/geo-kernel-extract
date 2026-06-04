# Maths Walk Postmortem - 2026-06-04

## SUPPORTS

S1 nmr_extract producer kernels:
- Coulomb/charge D_ab uses ke*q*(3rr/r5 - I/r3), with the stated E/EFG convention and traceless projection: ../src/CoulombResult.cpp:30, ../src/CoulombResult.cpp:156, ../src/CoulombResult.cpp:165.
- MOPAC-Coulomb mirrors the same EFG formula with MOPAC Mulliken charges and traceless projection, so the aimnet2/mopac-frame shadow is an EFG, not a ppm tensor: ../src/MopacCoulombResult.cpp:29, ../src/MopacCoulombResult.cpp:145, ../src/MopacCoulombResult.cpp:177.
- Ring current sign is explicit as sigma_ab = -n_b B_a, and the JB double-loop is unit-converted at the SI boundary: ../src/BiotSavartResult.cpp:90, ../src/BiotSavartResult.cpp:122, ../src/BiotSavartResult.cpp:234.
- McConnell/bond anisotropy implements the full tensor M_ab/r3 and then derives symmetric-traceless T2 category shadows from the full-M sums: ../src/McConnellResult.cpp:29, ../src/McConnellResult.cpp:84, ../src/McConnellResult.cpp:269.
- Larsen H-bond ppm tensors are generated from geometric donor/acceptor enumeration, rotated from the canonical frame to lab, summed by class, and emitted as full 9 spherical components: ../src/LarsenHBondShieldingResult.cpp:380, ../src/LarsenHBondShieldingResult.cpp:582, ../src/LarsenHBondShieldingResult.cpp:858, ../src/LarsenHBondShieldingResult.cpp:920.
- DFT shielding parse keeps ORCA raw 3x3 dia/para/total tensors and input coordinates, so rediscover is not reconstructing target tensors from scalar summaries: src/io/OrcaShieldingParser.cpp:94, src/io/OrcaShieldingParser.cpp:105, src/io/OrcaShieldingParser.cpp:119.

S2 substrate shadows and targets:
- The emitted charge shadow is q*D_ab in the same T2 library basis used by the producers: src/rediscover/PerAtomSubstrate.cpp:641, src/rediscover/PerAtomSubstrate.cpp:646, src/rediscover/PerAtomSubstrate.cpp:649.
- The geometric H-bond shadow is the requested dipolar*(hh - I/3) tensor before library decomposition: src/rediscover/PerAtomSubstrate.cpp:803, src/rediscover/PerAtomSubstrate.cpp:808, src/rediscover/PerAtomSubstrate.cpp:811.
- The T0/T1/T2 decomposition is isometric for T2 and has a single implementation used by rediscover targets and shadows: src/rediscover/SphericalBasis.cpp:7, src/rediscover/SphericalBasis.cpp:27.
- Targets are the raw DFT dia/para/total tensors decomposed in the same library basis, then written as total/dia/para T0/T1/T2 sidecars: src/rediscover/ExtractionSupport.cpp:51, src/rediscover/ExtractionSupport.cpp:58, src/rediscover/PerAtomSubstrate.cpp:2799.
- APBS/Water/MOPAC EFG convention reconciliation is explicit: APBS writes Hessian(phi) T2, Water writes the same Hessian form, and the fit flips water into the APBS/MOPAC -Hessian convention before unified use: ../src/ApbsEfgTimeSeriesTrajectoryResult.cpp:160, ../src/WaterFieldResult.cpp:125, src/rediscover/analysis/stage2_law_fits.py:145.

S3 unified D_ab-sum fit and statistics:
- The unified design matrix is a per-term scalar linear sum over flattened 5-component T2 features after the frozen library->e3nn basis map: src/rediscover/analysis/stage2_law_fits.py:595, src/rediscover/analysis/stage2_law_fits.py:607.
- The frozen change-of-basis is checked before loading the fit inputs, with a hard failure if |C.T C - I| is not tiny: src/rediscover/analysis/stage2_law_fits.py:722.
- Ridge standardization/imputation is train-only with an unpenalized intercept, and allatom PCA helpers fit PCA on train indices before transforming eval rows: src/rediscover/analysis/stage2_law_fits.py:1375, src/rediscover/analysis/allatom_fit_common.py:401, src/rediscover/analysis/allatom_fit_common.py:736, src/rediscover/analysis/allatom_fit_common.py:776.
- Within-frame alpha selection uses blocked frame folds inside the outer training rows, and the common helper records outer test/purged rows used for alpha selection: src/rediscover/analysis/stage2_law_fits.py:1463, src/rediscover/analysis/allatom_fit_common.py:1423, src/rediscover/analysis/allatom_fit_common.py:1502.
- Stage 2.3 nulls shuffle the axis-appropriate target structure: within atom for frame modulation, across atoms for LOAO, then score through the same design functions: src/rediscover/analysis/stage2_law_fits.py:2085, src/rediscover/analysis/stage2_law_fits.py:2129, src/rediscover/analysis/stage2_law_fits.py:2171.
- Drop-one ablation recomputes the full fit without each term and reports loss versus the all-term fit, i.e. marginal-given-others attribution: src/rediscover/analysis/stage2_law_fits.py:2294, src/rediscover/analysis/stage2_law_fits.py:2308, src/rediscover/analysis/stage2_law_fits.py:2317.

S4 e3nn path:
- The library<->e3nn 2e map is frozen, orthogonal, round-tripped, and Wigner-D tested against rotated Cartesian symmetric-traceless tensors: src/rediscover/analysis/change_of_basis.py:162, src/rediscover/analysis/test_change_of_basis.py:56, src/rediscover/analysis/test_change_of_basis.py:65.
- The ring e3nn model uses e3nn spherical harmonics plus the intended 1o x 1o -> 2e TensorProduct and scatter-add pooling, so the angular map itself is equivariant: src/rediscover/analysis/equiv_t2_e3nn.py:120, src/rediscover/analysis/equiv_t2_e3nn.py:131, src/rediscover/analysis/equiv_t2_e3nn.py:139.
- The EFG e3nn path maps feature and target with the same C.T, uses a blocked/purged frame split, and normalizes the scalar magnitude from training rows only: src/rediscover/analysis/equiv_t2_efg_e3nn.py:186, src/rediscover/analysis/equiv_t2_efg_e3nn.py:250, src/rediscover/analysis/equiv_t2_efg_e3nn.py:263.

S5 double-counting checks:
- The unified D_ab-sum includes charge, MC categories, MOPAC/water EFG field terms, geometric H-bond, pi-quadrupole, and dispersion, but excludes Larsen ppm tensors and excludes ring current tensors: src/rediscover/analysis/stage2_law_fits.py:140, src/rediscover/analysis/stage2_law_fits.py:146.
- The convention ledger states the intended exclusions: ringchi is separate-convention, Larsen ppm tensors are not double-counted with geometric hbond_T2, and ring current is not in the unified symmetric D_ab sum: src/rediscover/analysis/stage2_law_fits.py:2580, src/rediscover/analysis/stage2_law_fits.py:2582, src/rediscover/analysis/stage2_law_fits.py:2583.

## OPPOSES

S1 nmr_extract producer kernels:
- No concrete formula/sign/unit error found in the producer kernels inspected.

S2 substrate shadows and targets:
- Concrete validation gap: DftShieldingLoader says T0 identity stands in for all components, but trace equality does not imply dia+para==total for T1/T2; a component-level split mismatch could pass and then be emitted as target dia/para/total sidecars: src/io/DftShieldingLoader.cpp:95.

S3 unified D_ab-sum fit and statistics:
- Concrete procedure mismatch: the stage2_law_fits LOAO path trains with one atom held out but centers the held-out feature and target rows by the held-out atom's own mean; this measures within-held-atom modulation, not LOAO atom-mean recovery: src/rediscover/analysis/stage2_law_fits.py:657, src/rediscover/analysis/stage2_law_fits.py:1678, src/rediscover/analysis/stage2_law_fits.py:1681, src/rediscover/analysis/stage2_law_fits.py:2066.

S4 e3nn path:
- Concrete e3nn frame-split leakage/suspect: ring and broad-backbone e3nn de-mean targets over all groups before the split and de-mean predictions using all groups in the forward pass; their frame gates are also random/unpurged, unlike the stricter EFG path: src/rediscover/analysis/equiv_t2_e3nn.py:185, src/rediscover/analysis/equiv_t2_e3nn.py:141, src/rediscover/analysis/equiv_t2_e3nn.py:191, src/rediscover/analysis/equiv_t2_backbone_e3nn.py:309, src/rediscover/analysis/equiv_t2_backbone_e3nn.py:261, src/rediscover/analysis/equiv_t2_backbone_e3nn.py:372.

S5 double-counting checks:
- No concrete double-counting error found in the unified D_ab-sum path; the only double-countable ppm tensors/ring tensors are excluded from the unified specification.

maths-walk verdict: 3 concrete issues found
