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

### WaterFieldResult: no GeometryChoice recording
WaterFieldResult computes a Coulomb kernel over water atoms but does
not record GeometryChoices. Every calculator that evaluates a geometric
kernel must record what was included, excluded, and why — with entities,
roles, outcomes, and named numbers (distance, cutoff). Without this,
the UI cannot display the water field evaluation and the TOML cutoff
is not inspectable.

### HydrationShellResult: no GeometryChoice recording
Same issue. Half-shell asymmetry, dipole orientation, and ion distance
all involve geometric decisions (cutoff radii, shell boundaries) that
are not recorded.

### WaterFieldResult: no KernelFilterSet
Uses a raw distance cutoff (15A) without the KernelFilterSet pattern.
Should have MinDistanceFilter at minimum. The SelfSourceFilter is not
applicable (water atoms are not protein atoms), but DipolarNearFieldFilter
may be relevant for water molecules very close to protein atoms.

### HydrationShellResult: no KernelFilterSet
Uses raw distance cutoffs (3.5A first shell, 5.5A second shell, 20A
ion cutoff) without filter registration.

### TOML registration for solvent calculator parameters
The 15A water cutoff, 3.5A first-shell radius, 5.5A second-shell
radius, and 20A ion cutoff are physics constants that should be
registered in CALCULATOR_PARAMETER_API.md and accessible via TOML
override, matching the existing pattern for ring current calculators.

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
