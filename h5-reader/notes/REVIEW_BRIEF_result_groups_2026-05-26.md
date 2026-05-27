# Review brief — h5-reader result-group mirror (per-step)

**For:** a second set of eyes (codex / any reviewer). **Posture:** *cast a
wide net.* We verify every lead deterministically against the source before
acting, so **err toward flagging** — a false positive costs us one grep; a
missed real issue costs a silent mis-read. Uncertain hunches, smells, "this
looks off, please check X" are all welcome. Thank you for the careful read.

**Re-run cadence:** this brief is re-run after each result-group step. The
**current step** is the *ring-kernel family + the snapshot/catalog wiring*
it reads through; later steps add the bond/electrostatic, scalar/vector,
and MOPAC/Orca families (each will extend the file list below).

## What this is

The reader's per-frame **result-group views** mirror the SDK's `Protein`
groups (`../python/nmr_extract/_protein.py`, `_tensors.py`). Each view is a
thin, const handle over a `QtConformationSnapshot`, which holds a dense
per-`FieldKind` store of raw NPY columns (`NpyColumn`) filled at load via
the generated catalog. The views decode those columns into typed value
objects (`SphericalTensor`, `Vec3`, the `QtResultBlocks`). **Ground truth:**
`_catalog.py` (cols / native_axis / irreps), the SDK group/block classes
(shapes + accessor semantics), `../GEOMETRIC_KERNEL_CATALOGUE.md` (physics).

## Files in this step

- `src/model/QtConformationSnapshot.h` — the per-`FieldKind` `NpyColumn` store.
- `src/model/QtResultBlocks.h` — `PerRingTypeT0`, `PerRingTypeT2`, `RingCounts`.
- `src/model/QtRingKernelGroup.h` — base (shielding + perTypeT0 + perTypeT2).
- `src/model/QtBiotSavartGroup.h` — base + `totalB` + `ringCounts`.
- `src/model/QtHaighMallionGroup.h`, `QtPiQuadrupoleGroup.h`,
  `QtDispersionGroup.h` — trivial base subclasses (FieldKinds only).
- `src/model/QtRingSusceptibilityGroup.h` — standalone (shielding only).

## What to scrutinise (lead-maximising)

1. **Row-decode layout.** `QtRingKernelGroup::tensorAt` decodes a 9-col row
   as `T0=r[0]`, `T1=r[1..3]`, `T2=r[4..8]`. Does that match the SDK
   `SphericalTensor` packing `[T0, T1[3], T2[5]]`? `perTypeT2` reads a 40-col
   row as `byType[t][i] = r[t*5 + i]` (8 types × 5, row-major, ordered by
   `RingTypeIndex`) — correct vs `bs_per_type_T2` and `PerRingTypeT2`?
   `perTypeT0` as `byType[t] = r[t]` (8)? `RingCounts` as `{r[0..3]}` =
   3/5/8/12 Å — right order? `totalB` as `Vec3(r[0],r[1],r[2])`?
2. **FieldKind mappings.** Each group's ctor passes the right enumerators
   (`BSShielding`/`BSPerTypeT0`/`BSPerTypeT2`/`BSTotalB`/`BSRingCounts`;
   `HMShielding`/…; `PQ…`; `Disp…`; `RingChiShielding`). Cross-check the
   names + the convention against `QtFieldCatalog.gen.h`.
3. **Storage + indexing.** `NpyColumn::row(atomIdx) = data + atomIdx*cols`.
   Is `cols` always the catalog value (with `-1`→1 at load)? Any atomIdx
   bounds risk? Is the dense `std::array<NpyColumn, FieldKind::Count>` sound?
4. **"Absent, not faked."** Every accessor returns `std::nullopt` when
   `!snap.has(kind)`. Is that the right contract vs the conditional-attach
   discipline (a present 0.0 is real; absence is nullopt)?
5. **Anything peripheral** — base/subclass design, const-correctness, a
   stale comment, a unit mislabel, an off-by-one, a smell. Flag it.

## Step 2 (2026-05-26) — bond/electrostatic family

Files: `src/model/QtMcConnellGroup.h`, `QtCoulombGroup.h`, `QtHBondGroup.h`,
and the `QtResultBlocks.h` additions (`McConnellCategory` enum,
`PerBondCategoryT2`, `QtEfg`, `McConnellScalars`, `CoulombScalars`,
`HBondScalars`, `UnpackSphericalTensor`, `UnpackEfg`).

Scrutinise (wide net):
- **McConnellCategory** is a 5-value enum (backbone/sidechain/aromatic totals
  + CO_nearest + CN_nearest) — verify order vs SDK `_types.py` BondCategory,
  and that it is NOT conflated with the 8-value topology `BondCategory`.
- **PerBondCategoryT2** decode `byCategory[c][i] = r[c*5 + i]` (5×5 = 25,
  row-major) vs `mc_category_T2` (cols 25) + the writer `McConnellResult.cpp`
  + SDK `PerBondCategoryT2.for_category` / `as_block`.
- **Scalar-block field order** vs `_tensors.py`:
  `McConnellScalars{co_sum,cn_sum,sidechain_sum,aromatic_sum,nearest_CO_dist,
  nearest_CN_dist}`, `CoulombScalars{E_magnitude,E_bond_proj,E_backbone_frac,
  aromatic_E_magnitude}`, `HBondScalars{nearest_dist,inv_r3,count_3_5A,
  mcconnell_scalar}`.
- **QtEfg** is 5-component T2-only (irreps `1x2e`); `UnpackEfg` reads 5 cols.
  Verify vs SDK `EFGTensor` + catalog (`coulomb_efg_*` cols 5).
- **FieldKind mappings**: `McShielding/McCategoryT2/McScalars`,
  `CoulombShielding/CoulombE/CoulombEFGBackbone/CoulombEFGAromatic/CoulombScalars`,
  `HBondShielding/HBondScalars` — names + units vs `gen.h`.
- Coulomb is FullFat/`--mopac`-only (required=False); absent-on-standard-run
  → nullopt is the intended contract.

## Step 3 (2026-05-26) — scalar/vector family

**Source-of-truth rule for this step: the WRITER is definitive.** Each NPY's
column count / order / units / packing is whatever `src/*Result.cpp::WriteFeatures`
emits; `_catalog.py`, `_tensors.py`, `OBJECT_MODEL.md` are downstream and drift.
Verify every column against the writer line, not the catalog.

Files: `QtResultBlocks.h` additions (DsspScalars, DsspSs8, DsspHBondEnergy,
DsspChi, WaterShellCounts, HydrationShell, WaterPolarization, BondedEnergy,
GromacsEnergy) + the 9 group views QtDsspGroup, QtSasaGroup, QtWaterFieldGroup,
QtHydrationGroup, QtWaterPolarizationGroup, QtEeqGroup, QtBondedGroup,
QtGromacsGroup, QtApbsGroup.

Scrutinise (wide net — flag anything; we verify each lead):
- **GROMACS column count.** Block decodes 43; `GromacsEnergyResult.cpp:29-53`
  emits 43; fixture NPY is (1,43); catalog/gen.h say 42 (off-by-one, flagged for
  the team). Confirm the 43-index map (elec 0-2, bonded 3-8, vdw 9-11, thermo
  12-19, box 20-22, virial 23-31, pressure 32-40, T 41-42) and that virial/
  pressure `Mat3` fill is row-major XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ.
- **GROMACS is protein-axis** — `QtGromacsGroup::energy()` takes NO atomIdx,
  reads row(0). Verify it never indexes by atom.
- **dssp_chi** (DsspResult.cpp:283-321): interleaved per angle [cos,sin,exists],
  base=k*3, `forAngle(k)`→χ_{k+1}; missing angle → 0.0 (NOT NaN — the writer's
  own comment says NaN; the code writes 0.0).
- **dssp_backbone** (DsspResult.cpp:214-224): phi/psi are libdssp NEGATED-IUPAC
  (DsspResult.h:22-25); sasa is DSSP per-residue (≠ atom_sasa); ss_helix={H,G,I},
  ss_sheet={E,B}. Order phi,psi,sasa,ssHelix,ssSheet.
- **dssp_ss8** (DsspResult.cpp:234-255): one-hot col order H,G,I,E,B,T,S,C ==
  Types.h DsspCode ordinals; `dominant()`=argmax, `is(code)` bounds-checks.
- **dssp_hbond_energy**: acc0,acc1,don0,don1 (kcal/mol) — DsspResult.cpp:259-274.
- **EFG arrays** (water_efg, water_efg_first, apbs_efg): 5-col T2-only via
  `UnpackEfg` (writer PackT2 — WaterFieldResult.cpp:220/230, ApbsFieldResult.cpp
  :369). water_efield(_first) / apbs_E are 3-col xyz Vec3.
- **water_shell_counts** [nFirst,nSecond] (WaterFieldResult.cpp:255).
- **hydration_shell** [halfShellAsymmetry, meanWaterDipoleCos, nearestIonDist,
  nearestIonCharge] (HydrationShellResult.cpp:139-142); nearestIonDist is +inf
  when no ion in cutoff (:114-115) — `hasNearestIon()` guard.
- **water_polarization** 10-col [dipole(3), normal(3), asym, align, coher, count]
  (HydrationGeometryResult.cpp:141-150); normal duplicates SASA normal by design.
- **bonded_energy** 7-col [bond,angle,UB,proper,improper,cmap,total]
  (BondedEnergyResult.cpp:257-264); total=Σ6; FF-agnostic (UB/cmap 0 for ff14SB);
  reader deliberately does NOT propagate BondedEnergyResult.h's stale "CHARMM36m".
- **scalars**: atom_sasa (≠ DSSP sasa); sasa_normal (zero vector when buried);
  eeq_charges/eeq_cn — 1-col `row(atomIdx)[0]`.
- **absent→nullopt** at every accessor; `FromRow` decode matches column order.
- KNOWN cosmetic nit (logged, don't re-report): QtDsspGroup's `#include "Types.h"`
  is unused (transitive via QtResultBlocks.h); fix in post-review cleanup.

## Step 4 (2026-05-26) — AIMNet2 family

**Writer is definitive** (as Step 3). Files: `QtResultBlocks.h` `AIMNet2Embedding`
block + `QtAimnet2Group.h`. Writers: `AIMNet2Result.cpp` (WriteFeatures ~470-514)
+ `AIMNet2ChargeResponseGradientResult.cpp` (WriteFeatures ~214-248).

Scrutinise (wide net):
- **DTYPE: `aimnet2_aim` is float32** (`AIMNet2Result.cpp:492` WriteFloat32; fixture
  `<f4`) — the ONLY float32 array. The catalog has no dtype field, so the loader
  (#4) must read `<f4` from the NPY header and widen to double; the view then
  reads doubles. Confirm nothing decodes raw float32 bytes as double, and the
  comment flags the widen requirement.
- **`AIMNet2Embedding`** is a NON-OWNING view (const double* into the column, 256
  contiguous doubles), iterable (operator[]/size/begin/end), kDim=256 ==
  AIMNET2_AIM_DIMS. Lifetime = the snapshot's (transient use). Verify no copy and
  no retention assumption past the snapshot.
- **charges** (N,) `<f8` scalar, Hirshfeld (e); 1-col `row[0]`.
- **EFG total/aromatic/backbone** (N,5) T2 via `UnpackEfg`; mapping
  efg→AIMNet2EFG (`aimnet2_EFG_total_spherical`), efgAromatic→aromatic,
  efgBackbone→backbone (`AIMNet2Result.cpp:509-513`, PackT2). No sidechain EFG.
- **chargeResponseGradient** (N,3) Vec3 `dL/dr` (e²/Å, L=Σq²), parity-odd "1o",
  NOT a Buckingham polarisability; **chargeResponseGradientNorm** (N,) L2 norm.
  Verify the comment does not call it polarisability.
- **absent→nullopt** at every accessor; AIMNet2 is required + fails loud, so
  absence = the load path genuinely lacked it.
- FieldKind wiring: AIMNet2Charges/Aim/EFG/EFGAromatic/EFGBackbone/
  ChargeResponseGradient/ChargeResponseGradientScalar — cross-check vs gen.h.

## Step 5 (2026-05-26) — PlanarGeometry family

**Writer is definitive** (as Steps 3-4). File: `QtPlanarGeometryGroup.h` (all
scalars, no new blocks). Writer: `PlanarGeometryResult.cpp` (WriteFeatures
418-485) + `PlanarGeometryResult.h`.

Scrutinise (wide net):
- **MULTI-AXIS — first group with non-atom accessors.** Verify each accessor's
  index axis matches the writer's emission: pyramidalization=ATOM (`:427`, N),
  omegaActual/omegaDeviation/omegaIsXpro=RESIDUE (`:435/443/478`, ResidueCount),
  aromaticChi2=AROMATIC-RING (`:451`, AromaticRingCount), puckerQ/puckerTheta=
  SATURATED-RING (`:459/467`, SaturatedRingCount). Fixture sizes 846/54/54/15/1.
  Check the comment's "aromatic ordinal == global ring index while aromatic-first"
  and "saturated ordinal = global − AromaticRingCount" mapping.
- **DTYPE: omega_is_xpro is int8** (`PlanarGeometryResult.cpp:478` WriteInt8;
  fixture `|i1`). Loader widens → double; reader reads `row[0]!=0.0` → bool.
  Confirm the comment flags int8 + the loader-widen need.
- **omega NaN convention** (writer-definitive, header `:21-30`): NaN at
  chain-break / C-term / missing backbone cache ONLY. X→Pro gets its ACTUAL ω
  plus the omegaIsXpro mask — NOT NaN. (Catalog's "NaN at X-Pro" is stale.)
  Confirm the reader comment follows the writer, not the catalog.
- **Per-element NaN preserved**: ω/χ₂ accessors return optional{value} where the
  value may be NaN (a real "not applicable" sentinel) — NOT collapsed to nullopt
  (which means column-absent). Verify the distinction holds.
- **Units**: pyramidalization Å, omega/omega_deviation rad (deviation wrapped
  [−π,π]), aromatic_chi2 rad, pucker_Q Å, pucker_theta degrees [0,360).
- aromatic_chi2 honesty caveat (instantaneous χ₂, NOT flip kinetics; Akke &
  Weininger 2023). Cremer-Pople 1975 for pucker (Q, θ).
- absent→nullopt; FieldKind wiring vs gen.h (Pyramidalization/OmegaActual/
  OmegaDeviation/OmegaIsXpro/AromaticChi2/PuckerQ/PuckerTheta).

## Step 6 (2026-05-26) — Tripeptide + LarsenHBond reference-shielding family

**Writer is definitive** (as Steps 3-5). Files: `QtResultBlocks.h`
`TripeptideMethodTag` enum + `QtTripeptideGroup.h` + `QtLarsenHBondGroup.h`
(everything else reuses `UnpackSphericalTensor` + `Vec3` + scalars; no new
block). Writers: `TripeptideBackboneShieldingResult.cpp` (WriteFeatures
:306-378), `TripeptideNeighborShieldingResult.cpp` (:389-449),
`LarsenHBondShieldingResult.cpp` (:901-948); headers + `ConformationAtom.h`
:283-374 for the per-atom field defaults.

Scrutinise (wide net — flag anything; we verify each lead):
- **OPPOSITE SENTINELS — the headline.** Tripeptide writes **NaN** for "no
  match / no contribution" (writer inits the NPY buffer to NaN, writes the
  value only where `*_has_match`); Larsen writes **0.0** (the per-class
  tensors `PackFull9` UNCONDITIONALLY for every atom, default SphericalTensor
  = zero). Verify the reader comments + accessors keep both distinct from
  nullopt (= calculator absent this frame), and that the Larsen header does
  not claim NaN anywhere. Fixture cross-check: tripeptide_bb_* NaN on 38
  atoms; neighbor sum/prev/next NaN on 474/487/488; larsen_hbond_shielding
  657 all-zero rows == 657 count==0 rows, zero NaN.
- **Tripeptide method tag** (int8 → double widen): `TripeptideMethodTag`
  {NoMatch=0, Opbe=1, OrcaPbe=2} vs `MethodTagFromFrameType`
  (TripeptideBackboneShieldingResult.cpp:51) + ConformationAtom.h:283-291.
  Verify NoMatch(0) coincides with the NaN shielding rows (it does: 38==38).
  `hasBackboneMatch()` reads the int tag, not NaN — confirm that's sound.
- **Backbone fields** (all per-ATOM): bb_shielding (N,9) SphericalTensor;
  bb_residual_vec (N,3) Vec3 (ML feature — direction+magnitude load-bearing,
  feedback_residual_as_ml_feature); bb_match_distance (N,) = |residual|;
  bb_method_tag (N,) int8. Order + cols vs the writer blocks.
- **Neighbor fields** (Larsen 2015 Eq 3, per-ATOM): neighbor_shielding (N,9)
  is the SUM of i-1 + i+1 contributions (TripeptideNeighborShieldingResult.h
  :66-67); residual_vec_prev / _next (N,3) are per-direction, each NaN where
  THAT side did not contribute (so the three NaN counts differ — that's
  correct, not a bug). Verify prev↔VecPrev, next↔VecNext (no swap).
- **Larsen four-term decode** (Eq 5, per-ATOM, all N×9 SphericalTensor):
  shielding=total; 1pHB→shieldingPrimaryHB, 2pHB→shieldingSecondaryHB,
  1pHaB→shieldingPrimaryHaB, 2pHaB→shieldingSecondaryHaB. 1°=donor residue i,
  2°=acceptor residue i+1; HB=amide-H(HN) donor, HαB=Hα donor
  (LarsenHBondShieldingResult.h:85-93). Verify each accessor maps to the right
  FieldKind (no 1p/2p or HB/HaB transposition) vs the writer (:931-942).
- **Larsen water term**: `waterTerm` (N,) ppm = Δσ_w, 0.0 or up to 2.07 on a
  solvent-exposed amide H with zero H-bond candidates (header :12-14). Fixture
  range [0, 2.07], 833 zero. Never the full tensor.
- **Larsen Cβ diagnostic**: `diagnosticCBShielding` (N,9) SHOULD be ~0 (Larsen
  Table 2 Cβ all-false, LarsenHBondShieldingResult.h:75-77,99) — emitted to
  verify the parse→grid-load→rotate pipeline. Non-zero = methodological signal,
  NOT a bug. Fixture: 832/846 zero, the rest up to ±2.5 ppm. Confirm the
  comment frames it as a diagnostic, not a contribution.
- **Larsen count**: `count` (N,) int32 → double widen → int = larsen_hbond_n_pairs
  (LarsenHBondShieldingResult.cpp:945 WriteInt32). 0..6 in the fixture.
  `hasContribution()` = count>0; verify count==0 ⟺ zero shielding (657==657).
- **DTYPE for the loader (#4)**: tripeptide_bb_method_tag is `|i1`,
  larsen_hbond_count is `<i4` — both already in the verified non-`<f8`
  inventory; the views read the widened doubles. The 12 other arrays are `<f8`.
- **absent→nullopt** at every accessor; FieldKind wiring vs gen.h
  (TripeptideBBShielding/BBResidualVec/BBMatchDistance/BBMethodTag/
  NeighborShielding/NeighborResidualVecPrev/NeighborResidualVecNext;
  LarsenHBondShielding/1pHB/2pHB/1pHaB/2pHaB/DiagnosticCBShielding/WaterTerm/Count).

## Step 7 (2026-05-26) — MOPAC family (Core / Coulomb / McConnell)

**Writer is definitive** (as Steps 3-6). Files: `QtResultBlocks.h` additions
(`MopacScalars`, `MopacBondOrder`, `MopacGlobal`) + `QtMopacCoreGroup.h` +
`QtMopacCoulombGroup.h` + `QtMopacMcConnellGroup.h`. Writers: `MopacResult.cpp`
(WriteFeatures :482-535), `MopacCoulombResult.cpp` (:300-340),
`MopacMcConnellResult.cpp` (:300-339); header `MopacResult.h`. MOPAC is
FullFat `--mopac`-only (catalog `required=false`) → absent on standard runs →
nullopt is the intended contract.

Scrutinise (wide net — flag anything; we verify each lead):
- **NEW AXES — the headline.** This family introduces the Bond and Protein axes
  to the group views. (a) `mopac_bond_orders` is **(B,3) Bond-axis** — B is
  MOPAC's OWN unique-pair count (`bond_order_map_.size()`; 896 in the fixture)
  and is NOT the topology bond count, NOT parallel to QtProtein bonds; the
  bond-axis ordinal is unordered_map iteration order (arbitrary). Verify
  `bondOrder(bondIdx)` indexes by that ordinal and `bondOrderCount()` reads
  `NpyColumn.rows` (NOT AtomCount). (b) `mopac_global` is **(4,) Protein-axis**
  — verify `global()` takes no atom index and reads row(0); the on-disk shape is
  1-D `(4,)` NOT `(1,4)`, so flag that the #4 loader must map it to 1 row × 4
  cols (catalog cols=4), same as gromacs_energy.
- **MopacBondOrder decode**: `[atomA, atomB, wibergOrder]` with atomA/atomB cast
  double→size_t. Writer PairKey = min<<32|max, so **atomA < atomB** always
  (fixture: holds for all 896 rows, indices [0,845]). Verify the cast + the join
  semantics comment (atomA/atomB are the structure join, the ordinal is not).
- **MopacScalars decode**: `[charge, sPop, pPop, valency]` (MopacResult.cpp:494-503).
  col0 `charge` **duplicates** the `mopac_charges` array (fixture: exactly equal)
  — verify the comment says so. Units: charge e, populations electrons, valency
  = Wiberg Σ bond orders.
- **MopacGlobal decode**: `[heatOfFormation, dipole(3)]` (MopacResult.cpp:529-532).
  Units HoF **kcal/mol** (MopacResult.cpp:435), dipole **Debye**
  (MopacResult.cpp:212) — verify the comment's units against those writer lines.
- **MOPACCoulomb is an EXACT structural mirror of QtCoulombGroup**
  (MopacCoulombResult.cpp:312-338): shielding(9)/E(3)/efg_backbone(5)/
  efg_aromatic(5)/scalars(4). scalars = `[E_magnitude, E_bond_proj,
  E_backbone_frac, E_aromatic.norm()]` reusing `CoulombScalars` (the 3rd is a
  SIGNED projection, not a fraction). Verify NO field swap vs ff14SB Coulomb and
  the FieldKinds are the MOPAC* enumerators (not the ff14SB ones).
- **MOPACMcConnell is an EXACT structural mirror of QtMcConnellGroup**
  (MopacMcConnellResult.cpp:312-330): shielding(9)/category_T2(25)/scalars(6).
  Category order in the writer is `{backbone, sidechain, aromatic, CO_nearest,
  CN_nearest}` == `McConnellCategory` {0,1,2,3,4} — verify the 5×5 `c*5+m`
  decode + that scalars = `[co_sum, cn_sum, sidechain_sum, aromatic_sum,
  nearest_CO_dist, nearest_CN_dist]` (McConnellScalars order) with no swap.
- **Units in the catalog** (cross-check, writer wins on conflict): mopac_charges
  e; mopac_coulomb_shielding V/Å², _E V/Å, _efg_* V/Å²; mopac_mc_* Å⁻³.
- **absent→nullopt** at every accessor; FieldKind wiring vs gen.h (MOPACCharges/
  Scalars/BondOrders/Global; MOPACCoulombShielding/E/EFGBackbone/EFGAromatic/
  Scalars; MOPACMcShielding/McCategoryT2/McScalars). All arrays are `<f8`
  (no new non-`<f8` dtype here — bond indices are float64 on disk).

## Step 8 (2026-05-26) — Orca family (DFT quantum reference)

**Writer is definitive** (as Steps 3-7). File: `QtOrcaGroup.h` (3 accessors, no
new block — reuses `UnpackSphericalTensor`). Writer: `OrcaShieldingResult.cpp`
(parser :66-145, positional matchup :172-204, pack :235-254). **Fixture is the
single-pose baseline** — Orca is file-loaded, NOT in the trajectory fixture:
`baselines/old-system-baseline/orca-single_A0A062V9G2-WT_*/orca_{total,diamagnetic,paramagnetic}.npy`
(each (732,9) `<f8`).

Scrutinise (wide net):
- **3 accessors → right FieldKind, no swap**: total→OrcaTotal, diamagnetic→
  OrcaDiamagnetic, paramagnetic→OrcaParamagnetic; stems orca_total /
  orca_diamagnetic / orca_paramagnetic (gen.h:254-256). 9-col SphericalTensor
  via `UnpackSphericalTensor` ([T0, T1[3], T2[5]]) — matches the writer's
  `SphericalTensor::Decompose(Mat3)` → `PackFull9`.
- **Physics**: these are the DFT quantum-reference shielding (ppm), NOT a kernel
  (mechanism "quantum_reference") — the calibration target. total = diamagnetic +
  paramagnetic (Ramsey); fixture residual ~1.6e-3 ppm (print precision). Confirm
  the comment frames it as reference-not-kernel.
- **File-loaded / single-conformation only**: present only on `--orca`/`--mutant`
  runs, absent on `--trajectory` → nullopt on a trajectory snapshot. Confirm the
  comment says so and doesn't imply it's in the per-frame fixture.
- All three are `<f8` (no dtype subtlety). absent→nullopt at every accessor.

## How to report

`reader file:line | ground-truth file:line (or value) | expected | actual |
why it matters | confidence`, plus a corner-of-the-eye section. Precise
citations help most — we check each one.
