# Larsen H-bond Grounding

This is grounding only. It records what the producer and reference actually say so
the `larsen_hbond` rediscover relationship can be designed before it is wired.
No decisions are taken here.

## Working conclusion

The current rediscover stub should not be wired as written. It treats
`larsen_hbond` like a generic `T0+T2` source-sum over a dipolar geometric
kernel. Larsen/ProCS15 is different: the named literature law is an isotropic
ppm shielding lookup over H-bond geometry, and the library implementation is the
"ppm exception". The library's grid result can be useful as a ppm cross-check,
but it is already the calibrated Larsen signal. For law recovery, the spine must
emit the H-bond geometry and classes, not only consume the precomputed Larsen ppm
time series.

For the HN home stratum, the first rediscover target should be a scalar HN
shielding target, i.e. a scalar sigma de-circularisation on HN. Treat this as
scalar ppm, not as a T2 relationship, unless the lead explicitly opens a
separate tensor-grid study.

## Producer: what the library computes

Primary files read:

- `../src/LarsenHBondShieldingResult.cpp`
- `../src/LarsenHBondShieldingResult.h`
- `../src/LarsenHBondGrid.h`
- `../src/LarsenHBondGrid.cpp`
- `src/model/QtLarsenHBondGroup.h`

The calculator is `LarsenHBondShieldingResult`. It depends on
`SpatialIndexResult`, enumerates geometric H-bond candidates, queries
`LarsenHBondGrid`, rotates returned DFT-grid tensors into the protein lab frame,
and accumulates per-atom ppm shielding fields.

Detection in the library:

- Donors are backbone amide hydrogens (`HN`) and any alpha hydrogen (`H-alpha`,
  including both Gly HA atoms). This is narrower than the old rediscover
  "exchangeable H" stub. Hydroxyl hydrogens are not donor candidates in this
  calculator.
- Candidate acceptors are oxygen atoms only. The code skips same-residue
  contacts before classification.
- Candidate enumeration uses a padded spatial radius of `4.2 Angstrom`.
  This is not the physical law cutoff. It is a sweep cutoff chosen to cover the
  Larsen grid ranges.
- The geometric H-bond gate is `theta` in `[90.0, 180.0]` degrees. Here
  `theta` is the angle `donor_H - acceptor_O - acceptor_C`, with the vertex at
  the acceptor oxygen. It is not the generic `D-H-A` angle named in the old
  rediscover stub.
- The grid then gates by its own `r` axis. The code comments record Larsen scan
  ranges of `rOH` in `[1.5, 3.0] Angstrom` for the NMA/amide-H donor grids and
  `[1.8, 4.0] Angstrom` for the ALA/H-alpha donor grids.
- Geometry sent to the grid is `(r, theta, rho)`:
  `r = |donor_H - acceptor_O|`,
  `theta = angle(donor_H, acceptor_O, acceptor_C)`,
  `rho = dihedral(donor_H, acceptor_O, acceptor_C, acceptor_third)`.
- IUPAC note: `LarsenHBondGrid.h` says Larsen's filename rho sign is opposite
  to the IUPAC convention used by this loader. The parser keys the grid on
  IUPAC-sign rho, so a rediscover-side geometry attacher should compute the
  standard IUPAC-sign dihedral and move on.

Acceptor classes in the library:

- `BackboneCarbonyl`: backbone O, with C as `C'` and the third atom as `N(j+1)`.
- `SidechainCarbonyl`: Asn OD1 / Gln OE1, using the backbone/NMA acceptor grid
  as an explicit approximation. No separate sidechain amide acceptor grid exists.
- `HydroxylOxygen`: Ser OG, Thr OG1, Tyr OH. The bonded heavy atom is the
  acceptor-C analog and the bonded hydroxyl H is the third atom.
- `CarboxylateOxygen`: Asp/Glu sidechain carboxylate O or C-terminal carboxylate
  O. The other carboxylate O is the third atom. The code keeps only the closer
  O of a symmetric carboxylate pair for a given donor.

The library grid classes are:

- `AmideHydrogen x {BackboneCarbonyl, Hydroxyl, Carboxylate}` for `Delta sigma_HB`.
- `AlphaHydrogen x {BackboneCarbonyl, Hydroxyl, Carboxylate}` for
  `Delta sigma_HalphaB`.
- `SidechainCarbonyl` folds into `BackboneCarbonyl`.

What is emitted:

- The grid returns ppm shielding tensors, not an unscaled geometric kernel.
  The read-side mirror describes the four decomposition terms as per-atom,
  9-column `SphericalTensor`, units `ppm`.
- The producer writes four per-class Larsen H-bond time series:
  `larsen_hbond_1pHB_shielding_time_series`,
  `larsen_hbond_2pHB_shielding_time_series`,
  `larsen_hbond_1pHaB_shielding_time_series`, and
  `larsen_hbond_2pHaB_shielding_time_series`.
- Each per-class shielding time series is packed as `T0,T1,T2` in 9 components
  with units `ppm`.
- The producer also writes `larsen_hbond_water_term_time_series`, a scalar ppm
  time series, and `larsen_hbond_count_time_series`, an integer count.
- The water term is `2.07 ppm`, applied only to amide H atoms with zero
  geometric H-bond candidates under the library's geometric finder.
- The count is not a raw geometric candidate count. It is the number of
  Table-2 pair contributions accumulated on that atom; the diagnostic C-beta
  tensor does not increment it.

Where the shift lands:

- Contributions are accumulated on `ConformationAtom` fields, then decomposed
  to `SphericalTensor`.
- Primary terms land on atoms in the donor residue according to the library's
  `LarsenContribDispatch` table.
- Secondary terms use acceptor-side readouts when the acceptor is a backbone
  carbonyl with a valid `j+1` residue mapping. The acceptor readout `C` routes
  to the acceptor residue's own carbonyl carbon; `N`, `HN`, `CA`, and `HA` route
  to residue `j+1`.
- For HN scalar work, the important practical fact is that HN can receive
  Larsen H-bond ppm from the four per-class terms and can receive the scalar
  water term. The lead still needs to decide which of those terms are in the
  first HN de-circularised target/cross-check.

## Rediscover stub status

Current runtime status:

- `src/rediscover/main_extract.cpp` marks `larsen_hbond` as a fail-loud stub.
  Running it exits during validation with "data/decision is not ready, no zeros
  emitted".
- `src/rediscover/Catalog.h` and `src/rediscover/Catalog.cpp` have no Larsen
  `ArrayId` entries. The rediscover catalog currently exposes `KernelBs`,
  `KernelMc`, APBS fields, AIMNet2 fields, charges, positions, and DFT targets,
  but not the four Larsen ppm time series, the water term, or the count.

Current design-doc status:

- `src/rediscover/REDISCOVERY_MAP.md` still lists `larsen_hbond` as `T0+T2` and
  stubbed.
- `src/rediscover/STUB_LANGUAGE.md` proposes a placeholder relationship over
  `isExchangeable`, `hnFrame`, `near(Atoms, 3.5)`, `dipolarForm`, N/O acceptors,
  and `KernelBs` as a placeholder cross-check.
- That placeholder is stale relative to the producer. The producer uses
  HN/H-alpha donors, O acceptors, a 4.2 Angstrom spatial sweep, a 90-180 degree
  acceptor-angle gate, Larsen classes, and ppm grid lookup. It is not a
  Biot-Savart or McConnell-style kernel.

What wiring would need:

- A settled stratum. For the requested HN home stratum, this likely means a
  backbone amide HN-only first pass, not all `QtAtom::isExchangeable` hydrogens.
- A settled detection mode. If the lead chooses `GeometricLarsen`, rediscover
  should use the same geometry definitions as the library. DSSP/Kabsch-Sander is
  a different H-bond set and must remain a separate named choice if used at all.
- A source selector over candidate acceptor O atoms, with the producer's
  geometry and same-residue/carboxylate decisions either mirrored or explicitly
  changed by lead decision.
- A rediscover-side acceptor classification attacher. The classification lives in
  `../src/LarsenHBondShieldingResult.cpp`; the reader does not currently expose
  per-pair Larsen classes from H5.
- A source-row schema carrying at least donor class, acceptor class, donor atom,
  acceptor O, acceptor C, acceptor third, `r_HO`, `theta_HOC`, `rho_HOCX`,
  geometric disposition, and any class keys needed for pooling.
- An aggregate schema for class-pooled geometric features and, if desired, a
  fixed Larsen ppm cross-check.
- Catalog entries/readers for the existing Larsen H5 time series only if the
  cross-check consumes producer ppm directly.

Can the library ppm be consumed directly?

- For a producer cross-check: yes, after adding catalog support. For HN scalar
  cross-check, the plausible read is the T0/isotropic components of the four
  Larsen per-class ppm time series, plus the scalar water term, on HN atoms.
  The exact combination is a lead decision.
- For law recovery: no, not by itself. Consuming the library ppm directly gives
  the answer produced by the Larsen grid. It does not expose the geometric law
  input and would make the recovery circular.
- Therefore, the law-recovery substrate must emit geometric H-bond features
  spine-side. At minimum emit Larsen's `H...O` distance, acceptor-centered
  `H-O-C` angle, and `H-O-C-third` dihedral. If the lead wants generic H-bond
  descriptors like `D-A` distance or `D-H-A` angle, emit those as additional
  descriptors with explicit names rather than replacing the Larsen coordinates.

## Reference: Larsen/ProCS15 law

Primary reference read:

- `../references/larsen-2015-procs15-dft-chemical-shift-predictor.pdf`
- `../references-meta/larsen-2015-procs15-dft-chemical-shift-predictor-summary-qwen-test.txt`

The paper defines ProCS15 as an isotropic chemical shielding predictor. Its
overall scalar shielding decomposition is:

```text
sigma_i =
  sigma_BB_i
  + Delta sigma_BB_i-1
  + Delta sigma_BB_i+1
  + Delta sigma_HB_i
  + Delta sigma_HalphaB_i
  + Delta sigma_RC_i
  + Delta sigma_w_i
```

The H-bond terms are scalar shielding contributions:

```text
Delta sigma_HB_i =
  Delta sigma_1HB_i(r_HO, theta, rho)
  + Delta sigma_2HB_i(r_OH, theta_O, rho_O)

Delta sigma_HalphaB_i =
  Delta sigma_1HalphaB_i(r_HalphaO, theta, rho)
  + Delta sigma_2HalphaB_i(r_OHalpha, theta_O, rho_O)
```

The paper's H-bond scan details:

- `Delta sigma_HB`: `rOH` from `1.5` to `3.0 Angstrom` in `0.125 Angstrom`
  steps; `thetaH` from `180` to `90` degrees in `10` degree steps; `rhoH`
  over `-180` to `180` degrees.
- `Delta sigma_HalphaB`: `rOHalpha` from `1.8` to `4.0 Angstrom` in
  `0.2 Angstrom` steps; `thetaHalpha` from `180` to `90` degrees in
  `10` degree steps; `rhoHalpha` over `-180` to `180` degrees in `15` degree
  steps.
- The paper defines the angle around the acceptor oxygen using `H..O=C` or
  `H..O-C`, and the dihedral using `H..O=C-N`, `H..O=C-C`, or
  `H..O-C(..)HO`. This matches the library's `H-O-C-third` coordinate family.
- The water correction is `Delta sigma_w = 2.07 ppm` for amide protons that do
  not form H-bonds to other atoms in the protein structure.

Donor/acceptor classes from the reference, as operationalized by the library:

- Donors: amide H for `Delta sigma_HB`; H-alpha for `Delta sigma_HalphaB`.
- Acceptor chemical models: amide/carbonyl, hydroxyl/methanol-like, and
  carboxylate/acetate-like O environments. The library maps these to backbone
  carbonyl, sidechain carbonyl-as-backbone approximation, hydroxyl oxygen, and
  carboxylate oxygen.

The literature law is isotropic ppm shielding. The library stores and emits
full tensor readouts because the parsed DFT grids contain tensor data and the
project preserves tensor pathways, but the named ProCS15 law in the paper is a
scalar shielding contribution. Do not let the library's 9-column storage shape
turn the first Larsen rediscovery cell into a T2 claim by default.

## Separate sub-thread: ProCS15 as approximate DFT

Keep this separate from the r2SCAN target.

ProCS15 produces approximate shieldings by additive lookup/interpolation from
fragment DFT scans: OPBE/6-31G(d,p)//PM6 GIAO calculations on capped
tripeptides and small H-bond model systems, with CPCM solvation for the
backbone/sidechain grids. It then converts shieldings to shifts with per-atom
linear regression parameters. This is useful as a named literature approximation
and as a producer cross-check.

It is not the thesis DFT target. The rediscover target remains the actual
r2SCAN/ORCA shielding data used elsewhere in the project unless the lead opens a
separate "ProCS15 approximate-DFT" analysis. Do not mix ProCS15 scalar ppm
outputs into the r2SCAN target definition.

## Design questions for the lead

1. Detection criteria:
   Should the first cell use `GeometricLarsen` exactly as the library does
   (`4.2 Angstrom` spatial sweep, O acceptors, same-residue exclusion, 90-180
   degree acceptor-angle gate, grid-range miss handling), or should a different
   named H-bond detection be used?

2. Donor scope:
   Is the first rediscover cell HN-only, as the HN home-stratum prompt implies,
   or should it include H-alpha donors from the full library calculator? Should
   all other exchangeable hydrogens be explicitly out of scope for the first
   pass?

3. Acceptor classes:
   Should rediscover mirror the library's four acceptor classes exactly,
   including the sidechain-carbonyl-as-backbone approximation and carboxylate
   closer-O rule? Are N acceptors intentionally excluded for the Larsen cell?

4. Backbone-first scope:
   Should the first implementation restrict to backbone amide HN targets and
   backbone/sidechain O acceptors needed by the Larsen law, leaving non-HN
   donor frames and broader exchangeable-H chemistry for a later cell?

5. Scalar de-circularised target:
   What is the exact scalar target definition for HN? Candidate: DFT isotropic
   shielding sigma on HN after whatever non-Hbond de-circularisation the lead
   wants for this cell. This should be specified as scalar ppm, not T2.

6. Cross-check:
   Which producer ppm signal should be the fixed Larsen cross-check for HN?
   Candidate: sum the HN `T0` components of `1pHB`, `2pHB`, `1pHaB`, and
   `2pHaB`, then add `water_term`. The lead should decide whether H-alpha donor
   terms are included in the first HN cross-check and whether to also emit the
   four components separately.

7. Geometry features:
   Should the source rows emit only Larsen coordinates (`H...O`, `H-O-C`,
   `H-O-C-third`) or also generic H-bond descriptors (`D-A`, `D-H-A`) for
   interpretability? If both are emitted, their names and units must make the
   distinction explicit.

8. Law-recovery route:
   Is the objective to recover the scalar ProCS15 H-bond law from emitted
   geometry, or only to test whether the fixed producer Larsen ppm signal
   correlates with the HN DFT residual? These are different depths and should
   not share the same success criterion.

9. Approximate-DFT sub-thread:
   Should ProCS15 be documented only as a fixed literature cross-check, or should
   a separate analysis compare ProCS15 approximate shieldings against r2SCAN and
   experiment? If opened, keep it clearly labelled and separate from the
   r2SCAN-target rediscovery cell.
