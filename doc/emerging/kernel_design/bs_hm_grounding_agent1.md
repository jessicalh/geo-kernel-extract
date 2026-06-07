# Biot-Savart + Haigh-Mallion grounding — first-stage (codex xhigh; read-only; web + corpus), 2026-06-06

Pre-staged for the BS/HM work-through. Brief: /tmp/bs_hm_grounding_brief.txt. Companion to mcconnell_grounding_agent1.md + efg_grounding_agent1.md.

---

## Independent corroboration — lit sweep on the ring code (2026-06-06)

A separate lit sweep (Jessica's, run independently of this grounding) landed on the same two issues —
two-path agreement on the *findings*, not just the physics:

1. **HM naming.** Our `HaighMallionResult` is a **surface-integral *tensor* variant**, not the literal
   published *scalar* Haigh–Mallion signed-area bond-sum `ΣS_ij(1/r_i³ + 1/r_j³)` (Case eq.5 / Moyna /
   Sahakyan). Same family, not provably the same formula. **In-scope move:** benchmark our tensor's
   **scalar part (`T0`)** against the published bond-sum in populated geometries — *correlate, don't
   match*, up to scale/sign. If it correlates → document it as **"Haigh–Mallion extended to the full
   tensor"** (the published HM is scalar-only; our `T2` is the principled extension — the name becomes
   earned, and it's the richer object the equivariant emit wants). If not → rename **"HM-style tensor
   surface kernel."** Honest-naming, fits *correlate-don't-match*.
2. **BS units/sign labels [SOLID — pure label/doc fix].** `bs_shielding_contribution`,
   `bs_shielding.npy`, the `"ppm"` label, and `GeometryChoice "nA"` (should be `nA/T`) all read as
   *final shielding* when the value is a **unit-current geometric kernel**. Matches this grounding's
   "BS emits unit-current kernels, name it so," and the long-standing `nA` vs `nA·T⁻¹` source-label
   bug from the calculator-exposition pass. No physics change — stop the labels lying.

## Cross-kernel: ring ↔ EFG overlap (2026-06-06)

Ring (BS/HM) and charge/EFG do **not** overlap at the *mechanism* level — ring is **magnetic**
(induced ring current → secondary magnetic field); EFG is **electric** (static field + gradient from
charges). Neither uses the other's source, so there is **no double-count to zero** (contrast
ring↔McConnell, which *is* the same magnetic physics → aromatic-zero). They **coexist** as
complementary channels.

They **do** correlate at the *geometric-descriptor* level near aromatic groups: the ring `2e` and the
aromatic part of the EFG `2e` are both low-order multipole moments of the *same* neighbourhood (the
dual-nature — a kernel is a physical field AND a geometric descriptor). That is **correlation, not
double-counting** — the law study (NEW step 1) splits shared vs unique variance by ablation and
*reports* the collinearity; it is NOT zeroed.

Real overlaps live **within** families: ring↔McConnell (magnetic — aromatic-zero), **EFG↔π-quadrupole**
(electric — partition-contested), the EFG source-fork (MOPAC/FF14SB/APBS), and
MOPAC-populations↔AIMNet2.

---

**Verdict**

Keep both calculators. The core two-path idea is well supported: Case 1995 calibrates Johnson-Bovey and Haigh-Mallion against DFT and finds no significant scalar-fit difference for protein aromatic rings; Moyna 1998 finds HM and JB within about a 5% RMS gap on protein proton shifts. What needs changing is the output contract, calibration layer, and parity metadata, not the decision to keep both.

The ring shielding emit must be full even `0e ⊕ 1e ⊕ 2e`, never `1o`. Current code decomposes full tensors correctly in [Types.cpp](/shared/2026Thesis/nmr-shielding/src/Types.cpp:25), but trajectory metadata says `0e+1o+2e` in BS/HM/ring-chi headers, for example [BsShieldingTimeSeriesTrajectoryResult.h](/shared/2026Thesis/nmr-shielding/src/BsShieldingTimeSeriesTrajectoryResult.h:24) and [HmShieldingTimeSeriesTrajectoryResult.h](/shared/2026Thesis/nmr-shielding/src/HmShieldingTimeSeriesTrajectoryResult.h:20). That is wrong for e3nn.

**Biot-Savart**

`BiotSavartResult` holds up as the production Johnson-Bovey-style distributed-source path. It uses two loops at `±lobe_offset`, integrates wire segments in SI, then builds `G_ab = -n_b B_a * PPM_FACTOR` and decomposes to `T0/T1/T2` [BiotSavartResult.cpp](/shared/2026Thesis/nmr-shielding/src/BiotSavartResult.cpp:69), [BiotSavartResult.cpp](/shared/2026Thesis/nmr-shielding/src/BiotSavartResult.cpp:228). Boyd-Skrynnikov 2002 is the right citation for the tensor lift from isotropic Johnson-Bovey to full rank-2 shielding; their eq. 3 assembles the tensor, not a separate calculator.

Changes I would make:
- Rename/metadata-contract the current emitted BS values as unit-current geometric kernels unless they are multiplied by ring intensity. The code records `ring.Intensity()` in provenance but computes with `current=1.0` and does not apply the intensity in the tensor sum [BiotSavartResult.cpp](/shared/2026Thesis/nmr-shielding/src/BiotSavartResult.cpp:190), [BiotSavartResult.cpp](/shared/2026Thesis/nmr-shielding/src/BiotSavartResult.cpp:230).
- Emit per-ring-type full 9 components, not only `T0` and `T2`; the full atom total exists, but per-type `T1` is dropped [BiotSavartResult.cpp](/shared/2026Thesis/nmr-shielding/src/BiotSavartResult.cpp:311), [BiotSavartResult.cpp](/shared/2026Thesis/nmr-shielding/src/BiotSavartResult.cpp:470).
- Revisit calibration constants. Current defaults look like older Giessner-Prettre/Pullman-style factors scaled to Phe = -12 nA/T [CalculatorConfig.cpp](/shared/2026Thesis/nmr-shielding/src/CalculatorConfig.cpp:27). Case 1995 reports larger, model-specific fitted factors, especially His and Trp-5, and explicitly separates JB from HM calibration.
- Revisit five-member-ring lobe offsets. Perkins 1977 and Case 1995 describe the Johnson-Bovey loop height as 0.64 Å from the ring plane; current defaults use 0.50-0.52 Å for His/Trp-5 [CalculatorConfig.cpp](/shared/2026Thesis/nmr-shielding/src/CalculatorConfig.cpp:38). If there is no stronger cite, treat those as hypotheses, not best-practice constants.
- Benchmark against analytic circular-loop JB/Boyd values. The polygonal loop over actual vertices is defensible, but it is a geometry generalization of the published circular-loop tables.

**Haigh-Mallion**

Keep it, but be precise. The code computes a surface integral of the dipolar kernel over a fan-triangulated ring area, contracts `H.n`, then builds the same rank-1 shielding tensor `G_ab=-n_b V_a` [HaighMallionResult.h](/shared/2026Thesis/nmr-shielding/src/HaighMallionResult.h:7), [HaighMallionResult.cpp](/shared/2026Thesis/nmr-shielding/src/HaighMallionResult.cpp:281). That is an excellent independent geometric path.

The thin point is naming/calibration: Case quotes the published HM geometric factor as a bond/triangle-area sum and then fits HM-specific intensity factors. I would either validate this quadrature implementation directly against HM/Case tabulations, or call it “HM-style surface-integral” until that check is done. It also emits Å^-1 geometric values, not calibrated ppm shielding [HmShieldingTimeSeriesTrajectoryResult.h](/shared/2026Thesis/nmr-shielding/src/HmShieldingTimeSeriesTrajectoryResult.h:26).

Changes I would make:
- Add HM-specific Case 1995 calibration, separate from BS/JB.
- Emit per-type full 9, including `T1`, just like BS should.
- Keep both `hm_H` and `hm_G`: `H` is the raw symmetric/traceless surface kernel; `G` is the shielding tensor. The sparse output already captures both [ConformationResult.cpp](/shared/2026Thesis/nmr-shielding/src/ConformationResult.cpp:50).

**Agreement**

Do not compare raw BS and HM as equal numbers. Compare after a named scale/calibration, and stratify by geometry. Expected agreement: sign, spatial lobe structure, scalar `T0`, and major `T2` orientation should correlate well outside pathological near-ring regions. Legitimate divergence occurs directly above/below the pi cloud, close to vertices, in-plane, and in hetero/fused rings. Perkins 1977 explicitly preferred JB over HM above/below aromatic amino-acid rings; Agarwal 1977 showed local anisotropy can be about 35% of the observed shift and that JB/HM differ by scale and region.

What agreement buys you: a strong sign/topology/unit sanity check. It will catch flipped ring normals, inconsistent ring walks, bad intensity application, geometry pathologies, and bogus `T2` orientation. Use diagnostics such as scaled `T0` R², full-9 cosine, `T2` cosine, norm ratios, sign agreement, and residuals binned by `rho,z,theta,ring_type`.

**Emit Contract**

Defensible and generous shape:

- Per atom, per calculator: full packed tensor `[T0, T1[3], T2[5]]`, parity `0e+1e+2e`.
- Per atom and ring type: full 8 × 9, not current 8 × `T0` plus 8 × `T2`.
- Per atom-ring sparse rows: atom/ring ids, ring type, distance, `rho`, signed `z`, `theta`, azimuth, `(3cos²θ-1)/r³`, lobe offset, intensity/coefficient set id, BS full 9, HM `H` full 9, HM `G` full 9, and agreement diagnostics.
- Scaled variants should be named by calibration set: raw unit-current/unit-geometric, `Case95-JB`, `Case95-HM`, and optionally legacy GP/Pullman for reproducibility.

**Three New Steps**

Step 1, law/correlation: feed raw and calibrated BS/HM metrics, geometry scalars, full `T0/T1/T2`, `|T2|`, and BS↔HM agreement. Run DFT R² with `T0` only, with angular `T2`, and with full tensor features. Use BS↔HM agreement as a validation diagnostic and a stratification variable, not just another feature.

Step 2, equivariant shielding-tensor model: feed full even irreps from both calculators: `1x0e + 1x1e + 1x2e`. The `2e` is the Boyd-Skrynnikov tensor lift. If a downstream observable wants symmetric shielding only, expose a secondary `0e+2e` view, but do not replace the full emit.

Step 3, experimental-shift fido: feed everything ablatable: raw tensors, calibrated tensors, per-type tensors, sparse pair rows, counts/shells, geometry, ring constants, filter flags, and BS↔HM residual/coherence. Experimental shifts are confounded by local anisotropy, electric field, solvent, hydrogen bonding, and assignment noise; fido should get the rich feature record and let ablations decide.

**Cites Used**

Core: Boyd & Skrynnikov 2002 DOI 10.1021/ja016476f; Perkins et al. 1977 DOI 10.1042/bj1650223; Case 1995; Moyna et al. 1998; Agarwal et al. 1977. Tensor/equivariant support: Plasser & Glöcklhofer 2021 DOI 10.1002/ejoc.202100352; Gershoni-Poranne & Stanger 2015 DOI 10.1039/c5cs00114e; Ben Mahmoud et al. 2024; Geiger & Smidt 2022.

Gap: I did not find held originals for Johnson-Bovey 1958 or Haigh-Mallion 1971/1980; Case/Perkins/Moyna are strong held bridges, but exact HM formula equivalence should still be checked against the originals or reproduced tables.
