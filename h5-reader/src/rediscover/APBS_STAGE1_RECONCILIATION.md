# APBS Stage 1 Reconciliation

Date: 2026-06-02. Branch: `h5-reader-pysr-spike`.

Scope: read-only reconciliation. No code edits, no ORCA run, no
`trajectory.h5` access. This report is the only file written. I did not use git
history. The active checkout has no `../learn` directory; hidden Claude
worktrees do contain `learn/` trees on other `worktree-agent-*` branches, but
those are not this branch and were not used.

## Short Verdict

The clean reconciliation is:

1. Stage 1 / mutation / ORCA mode used the prmtop path and therefore real
   AMBER `%FLAG RADII` PB radii.
2. The trajectory TPR path currently uses our flat `1.5 A` compatibility PB
   radius for every atom.
3. Therefore the earlier "trajectory APBS has placeholder radii" finding is a
   path-specific finding in our pipeline, not a vendor or APBS-calculator
   defect.
4. The placeholder radii make trajectory APBS non-final for quantitative PB
   boundary claims, but the available evidence does not prove that they make the
   fields/EFGs ruinous. The decisive test is one real-vs-placeholder APBS solve
   on the same structure/frame.
5. Stage 1 success proves the prmtop/APBS wiring can work without poisoning the
   720-protein ridge study. It does not, from the active artifacts available
   here, prove that APBS field/EFG was load-bearing in the R^2 = 0.818 model.

## 1. Stage 1 Radii Path

Stage 1 is described in the active Stage 2 plan as the settled 2026-04-10
mutation/static-snapshot calibration: per-element, per-atom-type ridge on 720
proteins, 446K atoms, 55 kernels, R^2 = 0.818
(`spec/stage2_pysr_campaign_2026-05-29.md:10-12`).

The current single-pose path confirms the prmtop route:

- `--orca` root expansion requires `{root}.xyz`, `{root}.prmtop`, and
  `{root}_nmr.out`; `--mutant` requires WT and ALA roots with the same
  expansion (`../src/Cli/Parse.cpp:68-72`, `../src/Cli/Parse.cpp:127-152`).
- `RunOrca` and `RunMutant` both call `BuildFromOrca`, then pass
  `build.charges.get()` into `OperationRunner` (`../src/nmr_extract.cpp:161-179`,
  `../src/nmr_extract.cpp:193-217`).
- `BuildFromOrca` rejects missing or unreadable prmtops as a hard load error
  (`../src/OrcaRunLoader.cpp:366-373`).
- `LoadWithPrmtop` stores `files.prmtop_path` in `ProteinBuildContext`
  (`../src/OrcaRunLoader.cpp:322-329`).
- The resolver branch for a non-empty `build_context.prmtop_path` is explicitly
  `PrmtopChargeSource`, never re-tleap (`../src/AmberChargeResolver.h:143-147`,
  `../src/AmberChargeResolver.cpp:343-352`).
- `PrmtopChargeSource` reads `%FLAG CHARGE` and `%FLAG RADII`, requires enough
  `RADII` rows for every protein atom, assigns `pb_radius = raw_radii[i]`, and
  marks the row `Matched` (`../src/ChargeSource.cpp:273-275`,
  `../src/ChargeSource.cpp:340-368`).
- `OperationRunner` runs `ApbsFieldResult` after `ChargeAssignmentResult` is
  present (`../src/OperationRunner.cpp:134-153`).

That confirms the working statement: **Stage 1 / mutation / ORCA path = real
prmtop PB radii** in the present code.

The trajectory path is different:

- `FullSystemReader` copies TPR partial charges but sets every `pb_radius` to
  `kCompatibilityPlaceholderPbRadiusAngstrom` and marks each row
  `PlaceholderPbRadius` (`../src/FullSystemReader.cpp:940-943`).
- `ApbsFieldResult` warns when non-authoritative radii are present: APBS still
  populates finite `apbs_E / apbs_efg`, but the PB boundary is
  placeholder-quality and not physically validated
  (`../src/ApbsFieldResult.cpp:128-138`).

So the corrected path split is: **Stage 1 = real radii; trajectory = placeholder
radii**.

IUPAC/naming note: the ORCA loader canonicalizes AMBER/prmtop atom names and
preserves cap-only `H1/H2/H3` semantics (`../src/OrcaRunLoader.cpp:203-219`).
That naming surface does not change the APBS radii conclusion.

## 2. Placeholder Effect

The prior report was right that a flat radius is not just display metadata:
APBS uses atom radii to define the dielectric and ion-accessibility boundary.
But "radius-dependent" is not the same as "ruinous."

What can be quantified without running APBS:

- A read-only spot check of the existing ORCA test prmtops gives real radii
  close to, but not identical to, 1.5 A. For
  `../tests/data/orca/A0A7C5FAR6_WT.prmtop`, `n=543`, min `0.8`, mean `1.470`,
  max `1.8`, with counts `0.8:5`, `1.3:258`, `1.5:54`, `1.55:45`,
  `1.7:176`, `1.8:5`. Mean absolute deviation from 1.5 A is `0.173 A`.
- For `../tests/data/orca/A0A7C5FAR6_ALA.prmtop`, `n=501`, min `0.8`, mean
  `1.465`, max `1.8`, with counts `0.8:4`, `1.3:244`, `1.5:53`,
  `1.55:42`, `1.7:153`, `1.8:5`. Mean absolute deviation from 1.5 A is
  `0.171 A`.
- The prior APBS values report found emitted trajectory-side values that are
  finite and broadly screened relative to vacuum Coulomb, not numerically
  nonsensical. Example values cited there: APBS `|E|` medians `0.3958` and
  `0.5701 V/A`, APBS/Coulomb `|E|` median ratios `0.159` and `0.206`, and APBS
  `|EFG|` medians `0.4294` and `0.9717 V/A^2`
  (`src/rediscover/APBS_VALUES_VERIFICATION.md:139-147`).

Physics reading:

- Charges and geometry remain the dominant source terms. The placeholder
  changes the dielectric/ion boundary around those sources.
- Many real radii differ from 1.5 A by about 0.2 A; a few differ more. That is
  a material boundary error, especially near solvent-exposed surfaces and
  high-curvature atom contacts.
- The emitted values being finite, screened, and in a plausible range argues
  against calling the field/EFG "garbage" solely from this audit.
- Therefore the honest statement is: **the placeholder effect is unknown in
  magnitude from current evidence. It could be modest to moderate; it has not
  been proven catastrophic.**

Decisive test, not run here:

Run one identical APBS case twice, once with prmtop `%FLAG RADII` and once with
flat `1.5 A`, then compare:

- `E` vector cosine and magnitude ratio per atom.
- `EFG` T2 component correlation, `|T2|` correlation, and magnitude ratio.
- Surface vs buried stratification.
- Whether the downstream Buckingham/EFG fits materially change.

Until that A/B exists, "trajectory APBS cannot make the rediscovery work because
of placeholder radii" is only a hypothesis. The stronger proven claim is that
our trajectory APBS inputs are not final PB inputs.

## 3. Was APBS Load-Bearing In Stage 1?

The active branch does not contain the requested `learn/stage1*` coefficient or
importance artifacts. The active Stage 2 plan names the intended rebuild files
as `transforms/stage1_rebuild_features.py` and `sdk/ridge_stage1.py`, but those
files are not present in this checkout (`spec/stage2_pysr_campaign_2026-05-29.md:510-514`).

What the active docs do establish:

- Stage 1 was a 55-kernel ridge study with R^2 = 0.818
  (`spec/stage2_pysr_campaign_2026-05-29.md:10-12`).
- The current reader glossary describes the ridge predictions as Stage 1
  outputs and names core additive groups such as ring-current, EFG,
  bond-anisotropy, quadrupole, dispersion, and H-bond. It separately says newer
  water/charges/SASA groups have zero Stage 1 weights and ride along for later
  calibration (`notes/H5_FIELD_GLOSSARY.md:1013-1027`).

What the active docs do not establish:

- No active artifact gives Stage 1 per-kernel coefficients or importances.
- No active artifact separates APBS E-field/EFG marginal contribution from other
  EFG/charge groups inside the 55-kernel Stage 1 model.
- Therefore Stage 1 success does **not** by itself prove "APBS was
  load-bearing." It proves the prmtop/APBS path did not break the settled
  calibration.

Later one-protein diagnostics are mixed and should not be retrofitted onto
Stage 1. For example, `src/rediscover/ENSEMBLE_LAYER3.md` shows field/Buckingham
and EFG carrying positive drop-one deltas in several 1P9J layer-3 rows, while
`src/rediscover/analysis/EFG_ARC_EVIDENCE.md` says corrected local-frame
APBS-EFG is not the current explanatory axis for DFT T2. That is useful Stage 2
context, not Stage 1 importance evidence.

## Reconciled Position

The prior framing should be revised from "producer-side APBS artifact,
scheduled fix" to:

> Our trajectory extraction path currently supplies placeholder PB radii to
> APBS, while our Stage 1 static prmtop path supplied real AMBER PB radii. APBS
> itself remains the same solver. The inconsistency is in our pipeline wiring
> and in our interpretation of which path produced which evidence.

Actionable consequences:

1. For trajectory APBS, wire real per-atom PB radii into our trajectory path, or
   explicitly mark the trajectory APBS field/EFG as placeholder-boundary data.
2. Do the one-frame real-vs-placeholder A/B before claiming the placeholder
   radii are the reason APBS field/EFG rediscovery fails.
3. Locate or regenerate the Stage 1 coefficient/importances artifact before
   claiming APBS was or was not load-bearing in R^2 = 0.818.
4. If Stage 1 APBS weights are meaningful, then the trajectory placeholder path
   becomes a prime suspect for trajectory-specific APBS weakness.
5. If Stage 1 APBS weights are near zero, then the field/EFG mechanism is
   probably weak or redundant in that ridge model, independent of the trajectory
   placeholder issue.

Bottom line: **own the split. APBS worked for us on the Stage 1 prmtop path; the
trajectory placeholder radii are our wiring/validation gap, not a vendor defect
and not proof that APBS itself is bad.**
