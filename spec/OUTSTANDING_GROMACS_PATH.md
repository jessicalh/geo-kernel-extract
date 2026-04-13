# Outstanding: GROMACS Trajectory Path

Tracking remaining work for the trajectory extraction path.
Fixed items are struck through with the commit that resolved them.

## Must fix

### Object model update
OBJECT_MODEL.md has "UPDATE PENDING (2026-04-12)" for GromacsProtein,
GromacsFrameHandler, GromacsFinalResult, SasaResult, AIMNet2Result,
extended DsspResult. The CONSTITUTION has a one-line exception for
free-standing conformations (done). Full OBJECT_MODEL documentation
of the streaming classes, their lifetimes, and the accumulation
pattern is needed.

### ~~WaterFieldResult: no GeometryChoice recording~~
~~WaterFieldResult computes a Coulomb kernel over water atoms but does
not record GeometryChoices.~~ — Done. Parameter summary record +
per-atom exceptional events (singularity guard, E-field clamp).
CalculatorId::WaterField added. 2026-04-13

### ~~HydrationShellResult: no GeometryChoice recording~~
~~Same issue.~~ — Done. Parameter summary record with cutoff values,
atom count, water count. CalculatorId::HydrationShell added. 2026-04-13

### ~~WaterFieldResult: no KernelFilterSet~~
~~Uses a raw distance cutoff (15A) without the KernelFilterSet pattern.~~
— Done. MinDistanceFilter replaces inline r²<0.01 guard. SelfSource
not applicable (water ≠ protein atoms). DipolarNearField not applicable
(point charges, no source extent). 2026-04-13

### ~~HydrationShellResult: no KernelFilterSet~~
~~Uses raw distance cutoffs without filter registration.~~ — Not
applicable. No geometric kernel evaluated (counting + cos averaging).
Distance cutoffs are physics boundaries, not singularity guards.
TOML-registered as named constants. 2026-04-13

### ~~TOML registration for solvent calculator parameters~~
~~The 15A water cutoff, 3.5A first-shell radius, 5.5A second-shell
radius, and 20A ion cutoff are physics constants that should be
registered.~~ — Done. 4 parameters: water_efield_cutoff (15.0 A),
water_first_shell_cutoff (3.5 A, shared), water_second_shell_cutoff
(5.5 A), hydration_ion_cutoff (20.0 A). In CalculatorConfig defaults,
calculator_params.toml, CALCULATOR_PARAMETER_API.md, and test
assertions. 2026-04-13

## Should fix

### SelectFrames is a placeholder
Currently selects frames where water_n_first or water_emag have high
variance, then truncates to max_frames by sorted index (keeps earliest,
not most informative). The session notes describe a more intentional
strategy: 40 by water diversity + 10 by Boltzmann weight. Refine in R
after inspecting the atom_catalog.csv from a real run.

### SDK test data for trajectory path
No SDK test data on disk for the new arrays (water_efield, water_efg,
hydration_shell, gromacs_energy, atom_sasa). The catalog entries and
loader wiring are in place and verified against /tmp/trajectory_test
output, but no persistent test fixture exists under tests/data/.
Generate by running a small trajectory extraction and committing the
output.

### validate_smoke.py does not check new arrays
The Python smoke validation script does not include the water/hydration/
energy filenames in its per-atom or tensor validation lists.

### spec/EXTRACTION_SDK.md stale
Says "53 registered arrays covering all 10 calculators". Actual count
is now 60 arrays covering 14 calculators (added SASA, WaterField,
HydrationShell, GromacsEnergy). Update the count and add descriptions
of the new calculators.

## Must fix: HydrationGeometry + EEQ trajectory integration

HydrationGeometryResult and EeqResult are implemented as single-frame
calculators but not yet wired into the trajectory streaming path.
Both need integration into:

- **GromacsFrameHandler** — call HydrationGeometryResult and EeqResult
  per frame, feed results into Welford accumulators
- **GromacsProteinAtom** — add accumulator fields for water
  polarisation columns (dipole vector, asymmetry, coherence, etc.)
  and EEQ charges + coordination number
- **H5 master file (WriteH5)** — add rollup columns for the new
  accumulators (mean + std per atom)
- **WriteCatalog** — add corresponding CSV columns via AllWelfords()
- **SDK trajectory loader** — expose the new H5 columns in
  `load_trajectory()` return type

## Done

- ~~Double accumulation in pass 2~~ — ProcessFrame takes accumulate flag
- ~~total_frames_ not reset on Reopen~~ — reset in Reopen()
- ~~nearest_ion_distance default 0.0~~ — changed to infinity
- ~~CONSTITUTION free-standing conformation exception~~ — one-liner added
- ~~7 NPY files missing from SDK catalog~~ — all registered
- ~~SDK wrapper classes missing~~ — WaterFieldGroup, HydrationGroup added
- ~~Loader wiring missing~~ — wired in load()
- ~~API.md missing new arrays~~ — documented
- ~~Coulomb marked required but skippable~~ — changed to required=False
- ~~SDK test_required_count threshold~~ — updated from 30 to 27
- ~~H5 master file writing~~ — HighFive integration, WriteH5 with positions (T,N,3), 45-column rollup, per-bond stats. 2026-04-13
- ~~WriteCatalog expansion~~ — 45 columns, data-driven via AllWelfords(). 2026-04-13
- ~~Hard-fail error handling~~ — AIMNet2, WaterField, HydrationShell, GromacsEnergy all kill run on failure. 2026-04-13
- ~~Open() takes RunOptions~~ — frame 0 uses same calculator set as all frames. No MOPAC on frame 0. 2026-04-13
- ~~Production CLI (nmr_extract --trajectory) scan opts~~ — APBS + DSSP + AIMNet2 enabled per frame. 2026-04-13
- ~~SDK load_trajectory~~ — python/nmr_extract/_trajectory.py, reads H5 master file. 2026-04-13
- ~~Expanded accumulators~~ — 45 Welfords (6 classical T0+T2, AIMNet2, APBS, water, DSSP, chi1-4, bond angles), 4 DeltaTrackers, TransitionCounters, per-bond length, per-frame positions. 2026-04-13
