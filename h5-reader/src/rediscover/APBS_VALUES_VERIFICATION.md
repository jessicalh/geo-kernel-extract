# APBS Values Verification

Date: 2026-06-02. Scope: read-only verification of the producer APBS path and
existing emitted artifacts. No ORCA run, no extraction run, and no
`trajectory.h5` write was performed.

## Formal APBS Reference

- APBS molecular inputs formally carry atomic coordinates, partial charge, and
  atomic radius. The official PQR format lists `X Y Z Charge Radius`, with
  charge in electrons and radius in Å:
  <https://apbs.readthedocs.io/en/latest/formats/pqr.html>.
- APBS also states that pseudo-PDB inputs need a parameter file to provide
  charge and radius parameters:
  <https://apbs.readthedocs.io/en/latest/using/input/old/read.html>.
- For the molecular surface (`srfm mol`), APBS defines the dielectric and
  ion-accessibility coefficients from the biomolecular atom radii; the solvent
  free volume and inflated ion-exclusion volume are radius-dependent:
  <https://apbs.readthedocs.io/en/latest/using/input/old/elec/srfm.html>.
  Therefore, a uniform placeholder radius is not a harmless display parameter;
  it changes the dielectric boundary and the PB solution.
- Standard APBS parameter ranges line up with the producer defaults: `pdie` is
  usually 2-20 for biomolecules, `sdie` for water is usually 78-80, `ion`
  concentrations are in M, and common `dime` values include 65/97/129/161:
  <https://apbs.readthedocs.io/en/latest/using/input/old/elec/pdie.html>,
  <https://apbs.readthedocs.io/en/latest/using/input/old/elec/sdie.html>,
  <https://apbs.readthedocs.io/en/stable/using/input/old/elec/ion.html>,
  <https://apbs.readthedocs.io/en/nathan-docs/using/input/elec/dime.html>.
- APBS preferred potential units are `kBT/e`; at ~298 K this is about 25.6 mV:
  <https://apbs.readthedocs.io/en/latest/using/units.html>.
- Known protein local electric fields are commonly reported around
  `10^6-10^8 V/cm` (`0.01-1 V/Å`), with enzyme active-site examples around
  tens of MV/cm (`~0.1-1 V/Å`):
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC4238891/> and
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC5600505/>.

## Producer Path And Radii

The trajectory producer path is:

- `nmr_extract --trajectory`: `RunTrajectory` builds a `TrajectoryProtein` from
  the production directory, then opens `trajectory.h5` with truncate mode only
  when actually writing output (`../src/nmr_extract.cpp:231`,
  `../src/nmr_extract.cpp:243`, `../src/nmr_extract.cpp:267`,
  `../src/nmr_extract.cpp:279`, `../src/nmr_extract.cpp:285`,
  `../src/nmr_extract.cpp:288`).
- `TrajectoryProtein::BuildFromTrajectory` reads `production.tpr`, then builds
  protein and charges via `FullSystemReader::BuildProtein(..., Amber_ff14SB, ...)`
  (`../src/TrajectoryProtein.cpp:27`, `../src/TrajectoryProtein.cpp:30`,
  `../src/TrajectoryProtein.cpp:42`, `../src/TrajectoryProtein.cpp:92`,
  `../src/TrajectoryProtein.cpp:95`, `../src/TrajectoryProtein.cpp:102`).
- During per-frame execution, trajectory mode supplies `tp.Charges()` as the
  charge source (`../src/Trajectory.cpp:184`, `../src/Trajectory.cpp:185`);
  APBS is not skipped in the production per-frame extraction set
  (`../src/RunConfiguration.cpp:137`, `../src/RunConfiguration.cpp:141`,
  `../src/RunConfiguration.cpp:142`, `../src/RunConfiguration.cpp:155`,
  `../src/RunConfiguration.cpp:156`).

**Radii verdict: production trajectory APBS uses the uniform placeholder.**

- `FullSystemReader::BuildProtein` copies real TPR charges into
  `AtomChargeRadius.partial_charge`, but sets every `pb_radius` to
  `kCompatibilityPlaceholderPbRadiusAngstrom` and marks the row
  `PlaceholderPbRadius` (`../src/FullSystemReader.cpp:940`,
  `../src/FullSystemReader.cpp:941`, `../src/FullSystemReader.cpp:942`,
  `../src/FullSystemReader.cpp:943`).
- The placeholder constant is exactly `1.5 Å`; the source comment says this is
  a compatibility quarantine because the TPR path does not provide a defensible
  APBS/PB radius set yet (`../src/ChargeSource.h:94`,
  `../src/ChargeSource.h:98`).
- The charge table counts every non-`Matched` radius as non-authoritative
  (`../src/ForceFieldChargeTable.cpp:73`, `../src/ForceFieldChargeTable.cpp:76`).
- `ApbsFieldResult` explicitly warns when such radii are present and says the
  frame's `apbs_E / apbs_efg` are populated but not physically validated
  (`../src/ApbsFieldResult.cpp:128`, `../src/ApbsFieldResult.cpp:131`,
  `../src/ApbsFieldResult.cpp:134`, `../src/ApbsFieldResult.cpp:136`,
  `../src/ApbsFieldResult.cpp:137`, `../src/ApbsFieldResult.cpp:138`).
- Existing extraction logs show this happened for all 846 atoms in the available
  1P9J trajectory artifact (`../molprobity_runs/1P9J_5801/extract.log:12`).

There are producer paths that can carry real radii: the flat ff14SB param file
parser reads a `pb_radius` column (`../src/ChargeSource.cpp:69`,
`../src/ChargeSource.cpp:70`, `../src/ChargeSource.cpp:248`), and the AMBER
prmtop parser reads `%FLAG RADII` and assigns those PB radii
(`../src/ChargeSource.cpp:340`, `../src/ChargeSource.cpp:341`,
`../src/ChargeSource.cpp:356`, `../src/ChargeSource.cpp:367`,
`../src/ChargeSource.cpp:368`). Those are not what the trajectory TPR path uses.

## Charges, Grid, And APBS Parameters

Charges are real stored force-field charges, not a zero/stub fallback:

- The TPR path takes `tpr_atoms.atom[ai].q` as the partial charge
  (`../src/FullSystemReader.cpp:940`, `../src/FullSystemReader.cpp:941`).
- `ChargeAssignmentResult` projects the prepared table's `PartialChargeAt` and
  `PbRadiusAt` into the conformation atoms before APBS runs
  (`../src/ChargeAssignmentResult.cpp:35`, `../src/ChargeAssignmentResult.cpp:48`,
  `../src/ChargeAssignmentResult.cpp:51`, `../src/ChargeAssignmentResult.cpp:52`).
- `ApbsFieldResult` passes positions, charges, and radii to the C bridge
  (`../src/ApbsFieldResult.cpp:146`, `../src/ApbsFieldResult.cpp:147`,
  `../src/ApbsFieldResult.cpp:154`, `../src/ApbsFieldResult.cpp:156`,
  `../src/ApbsFieldResult.cpp:208`, `../src/ApbsFieldResult.cpp:211`).
- The bridge sets APBS `Vatom` position, charge, and radius directly
  (`../src/apbs_bridge.c:46`, `../src/apbs_bridge.c:49`,
  `../src/apbs_bridge.c:50`, `../src/apbs_bridge.c:51`).

The APBS solve parameters are otherwise sane for a protein PB calculation:

- Default grid/condition config: `grid_dim=161`, `fine_padding_A=40`,
  `fine_min_dim_A=40`, `coarse_padding_A=30`, `pdie=4`, `sdie=78.54`,
  `temperature=298.15 K`, `ionic_strength=0.15 M`
  (`../src/CalculatorConfig.cpp:87`, `../src/CalculatorConfig.cpp:88`,
  `../src/CalculatorConfig.cpp:89`, `../src/CalculatorConfig.cpp:90`,
  `../src/CalculatorConfig.cpp:91`, `../src/CalculatorConfig.cpp:92`,
  `../src/CalculatorConfig.cpp:93`, `../src/CalculatorConfig.cpp:94`,
  `../src/CalculatorConfig.cpp:95`).
- `ApbsFieldResult` computes bbox-derived fine/coarse dimensions from those
  settings (`../src/ApbsFieldResult.cpp:166`, `../src/ApbsFieldResult.cpp:171`,
  `../src/ApbsFieldResult.cpp:182`, `../src/ApbsFieldResult.cpp:187`,
  `../src/ApbsFieldResult.cpp:188`, `../src/ApbsFieldResult.cpp:191`,
  `../src/ApbsFieldResult.cpp:195`).
- The bridge uses a linearized PB solve with +1/-1 ions at the configured
  concentration and radii 0.95/1.81 Å (`../src/apbs_bridge.c:84`,
  `../src/apbs_bridge.c:85`, `../src/apbs_bridge.c:86`,
  `../src/apbs_bridge.c:88`, `../src/apbs_bridge.c:89`,
  `../src/apbs_bridge.c:135`, `../src/apbs_bridge.c:179`).
- It is a single manual/coarse-grid solve with SDH boundary conditions, not
  APBS two-level focusing (`../src/apbs_bridge.c:5`,
  `../src/apbs_bridge.c:6`, `../src/apbs_bridge.c:99`,
  `../src/apbs_bridge.c:107`, `../src/apbs_bridge.c:109`,
  `../src/apbs_bridge.c:124`, `../src/apbs_bridge.c:137`).

The grid and dielectric/ion choices are defensible. The radius model is not.

## Emitted-Value Sanity

Read-only stats from existing sidecars:

- `../tests/golden/blessed/fleet/frame_001`: `N=2480`, APBS `|E|` median
  `0.3958 V/Å`, p95 `0.8372 V/Å`, max `1.7828 V/Å`; APBS `|EFG|` median
  `0.4294 V/Å²`, p95 `0.9447 V/Å²`, max `1.6997 V/Å²`. Vacuum Coulomb `|E|`
  median/p95/max was `2.69/5.97/10.38 V/Å`; APBS/Coulomb `|E|` ratio median
  was `0.159`, consistent with solvent screening.
- `../baseline_features/P84477`: `N=351`, APBS `|E|` median `0.5701 V/Å`,
  p95 `1.3522 V/Å`, max `2.5693 V/Å`; APBS `|EFG|` median `0.9717 V/Å²`,
  p95 `2.6246 V/Å²`, max `4.3833 V/Å²`. Vacuum Coulomb `|E|` median/p95/max
  was `2.70/5.65/10.59 V/Å`; APBS/Coulomb `|E|` ratio median was `0.206`.
- Component signs are balanced around zero in both samples; no obvious global
  sign bias was seen.
- Surface-vs-buried separation from existing `atom_sasa.npy` was weak. For
  `frame_001`, `corr(SASA, |E|)=-0.065` and `corr(SASA, |EFG|)=-0.031`; for
  `P84477`, the correlations were approximately `-0.002` and `+0.014`. This is
  not enough to validate the dielectric boundary, and the placeholder radii
  make it especially non-authoritative.

Magnitude interpretation:

- These finite APBS values are in the broad physical range expected for local
  protein electrostatics. `0.4-1.4 V/Å` corresponds to `40-140 MV/cm`, which is
  high but not absurd for atom-local protein fields. EFG norms of order
  `0.4-3 V/Å²` are also plausible for Å-scale fixed-charge gradients.
- The values are therefore not "nonsense APBS output". They are plausible
  numerical outputs from APBS-like inputs.
- However, they are not physically expected APBS protein fields because APBS's
  dielectric/ion boundary is radius-defined, and trajectory mode gives every
  atom the same `1.5 Å` radius.

One additional handling caveat: `ApbsFieldResult` clamps the APBS-native
E-field before unit conversion (`../src/ApbsFieldResult.cpp:261`,
`../src/ApbsFieldResult.cpp:266`, `../src/ApbsFieldResult.cpp:276`,
`../src/ApbsFieldResult.cpp:287`). The `P84477` max `|E|=2.5693 V/Å` is exactly
`100 * 0.025693`, so at least one artifact appears clipped by that guard. This
is secondary to the radius problem but should be made explicit or revised before
quantitative coefficient claims.

## Units And Sign Convention

The units thread is mostly consistent:

- APBS grid potential is documented in the bridge header as `kT/e`; input
  positions/radii are Å and charges are elementary charges
  (`../src/apbs_bridge.h:13`, `../src/apbs_bridge.h:14`,
  `../src/apbs_bridge.h:15`, `../src/apbs_bridge.h:16`).
- Source computes `E = -grad(phi)` and EFG as `dE_i/dr_j`, then symmetrizes and
  trace-projects the EFG (`../src/ApbsFieldResult.cpp:72`,
  `../src/ApbsFieldResult.cpp:80`, `../src/ApbsFieldResult.cpp:86`,
  `../src/ApbsFieldResult.cpp:88`, `../src/ApbsFieldResult.cpp:99`,
  `../src/ApbsFieldResult.cpp:104`, `../src/ApbsFieldResult.cpp:106`,
  `../src/ApbsFieldResult.cpp:115`).
- The conversion `kT/(e*A) -> V/A` and `kT/(e*A^2) -> V/A^2` happens by
  multiplying both E and EFG by `KT_OVER_E_298K`
  (`../src/ApbsFieldResult.cpp:286`, `../src/ApbsFieldResult.cpp:287`,
  `../src/ApbsFieldResult.cpp:288`; constant at
  `../src/PhysicalConstants.h:67`, `../src/PhysicalConstants.h:69`).
- Per-frame NPY emission is `apbs_E` `(N,3)` and `apbs_efg` T2-only
  (`../src/ApbsFieldResult.cpp:358`, `../src/ApbsFieldResult.cpp:365`,
  `../src/ApbsFieldResult.cpp:368`, `../src/ApbsFieldResult.cpp:373`,
  `../src/ApbsFieldResult.cpp:375`).
- Dense H5 APBS E-field is written under
  `/trajectory/apbs_efield_time_series` with `units=V/Angstrom`
  (`../src/ApbsEfieldTimeSeriesTrajectoryResult.cpp:89`,
  `../src/ApbsEfieldTimeSeriesTrajectoryResult.cpp:96`,
  `../src/ApbsEfieldTimeSeriesTrajectoryResult.cpp:99`,
  `../src/ApbsEfieldTimeSeriesTrajectoryResult.cpp:112`,
  `../src/ApbsEfieldTimeSeriesTrajectoryResult.cpp:114`).
- Dense H5 APBS EFG is T2-only under `/trajectory/apbs_efg_time_series` with
  `units=V/Å^2` (`../src/ApbsEfgTimeSeriesTrajectoryResult.cpp:148`,
  `../src/ApbsEfgTimeSeriesTrajectoryResult.cpp:155`,
  `../src/ApbsEfgTimeSeriesTrajectoryResult.cpp:159`,
  `../src/ApbsEfgTimeSeriesTrajectoryResult.cpp:187`,
  `../src/ApbsEfgTimeSeriesTrajectoryResult.cpp:189`).
- Reader/catalog labels agree dimensionally: APBS E `V/A`, APBS EFG `V/A^2`
  (`src/model/QtApbsGroup.h:3`, `src/model/QtApbsGroup.h:8`,
  `src/model/QtApbsGroup.h:9`, `src/io/QtFieldCatalog.gen.h:252`,
  `src/io/QtFieldCatalog.gen.h:253`, `../python/nmr_extract/_catalog.py:361`,
  `../python/nmr_extract/_catalog.py:364`).

Minor cleanup: standardize `V/Angstrom`, `V/A`, and `V/Å` spelling across attrs
and catalogs. This is label hygiene, not the main artifact.

## Verdict

**Do not treat the APBS field/EFG rediscover failure as a mechanism failure
yet. In trajectory mode, it is a producer-side APBS input artifact.**

This is an input-verification finding; the final scientific attribution and
roadmap call remain with the lead.

APBS is being invoked with real stored ff14SB/AMBER-family charges and broadly
reasonable PB grid/dielectric/ionic settings. The emitted values are finite,
screened relative to vacuum Coulomb, sign-balanced, and roughly plausible in
magnitude. But APBS's documented PB behavior depends on atom radii for the
dielectric and ion-accessibility boundary, and the trajectory producer currently
feeds a uniform `1.5 Å` placeholder radius to every atom. That makes the PB
solution, E-field, and EFG physically non-authoritative.

Producer-side fix to schedule, not done here:

1. Wire authoritative per-atom PB radii for the trajectory source. Preferred
   options are AMBER/prmtop `%FLAG RADII` where available, or a generated
   ff14SB/PB radius table mapped through the same canonical atom/residue
   identity used for charges.
2. Treat `NonAuthoritativePbRadiusCount() > 0` as APBS QC failure for physical
   datasets: either do not emit APBS fields or emit with explicit invalid/QC
   provenance until radii are fixed.
3. Add radius provenance and APBS parameter attrs to the H5 groups.
4. Re-run extraction and only then revisit the rediscover field/EFG call.
5. Review the APBS-native E-field clamp and make clipping explicit in metadata
   or move/remove the guard.

IUPAC note: no IUPAC naming issue was needed for this verification; the relevant
path is the GROMACS/AMBER readback and APBS charge/radius payload.
