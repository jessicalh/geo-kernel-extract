# Session 2 handoff — ApbsEfgTS + MOPAC bundle (6 TRs)

**Status as of master `a3cb51d` (pushed origin/master 2026-05-21):**

S1 closed. AIMNet2 trio (TR1/TR2/TR3) landed and codex-reviewed
through R0+R1. 9 commits + plan doc update. Tests green: 6/6 C++ on
the trio binary, 135/135 Python SDK, Integration1P9J on real
AIMNet2 backward shows 846/846 atoms populated.

S2 entry conditions met (templates + helpers in place). See
`spec/plan/13-tr-build-plan-2026-05-20.md` §Session 2 for the
authoritative 6-TR scope. This doc is the boot-up cheat sheet —
read alongside the plan, not instead of.

## S2 scope (verbatim from the plan)

| # | TR | Source | TR class | Cadence | Output | Clone-from |
|---|----|--------|----------|---------|--------|------------|
| 4 | `ApbsEfgTimeSeriesTrajectoryResult` | `ApbsFieldResult.apbs_efg_spherical` (T2) | FO | every frame | (N, T, 5) double | `BsShieldingTimeSeriesTrajectoryResult` |
| 5 | `MopacChargeWelfordTrajectoryResult` | `MopacResult` (charges, scalar/atom) | AV | sparse (every 20 ps) | 1 Welford channel/atom | `HBondCountWelfordTrajectoryResult` + sparse gate pattern from Tripeptide BB Residual |
| 6 | `MopacBondOrderWelfordTrajectoryResult` | `MopacResult` (bond orders) | AV | sparse | 1 Welford channel/bond (bonds.npy axis) | clone of #5 with bond axis |
| 7 | `MopacCoulombShieldingTimeSeriesTrajectoryResult` | `MopacCoulombResult` (T2) | FO | sparse | (N, T, 5) double | combines #4 (T2 TS) + #5 (sparse gate) |
| 8 | `MopacMcConnellShieldingTimeSeriesTrajectoryResult` | `MopacMcConnellResult` (T2) | FO | sparse | (N, T, 5) double | clone of #7 |
| 9 | `MopacVsFf14SbReconciliationTrajectoryResult` | `MopacCoulombResult` + `CoulombResult` (FF14SB) | FO | sparse (both attached) | (N, T) double `|cos(MOPAC_T2, FF14SB_T2)|` | new cross-source — its own canonical |

Codex annealing pass at session close.

## Patterns S1 left in place that S2 should reuse

**1. Full canonical Welford row.** `vector_{mean,std,m2,min,max,min_frame,max_frame}` +
`scalar_{mean,std,m2,min,max,min_frame,max_frame}` + `n_per_atom`.
WelfordFinalize called in `Finalize()` gated on
`source_attached_count_ > 0` (avoids segfault on never-Computed atoms).
See `src/AIMNet2ChargeResponseGradientWelfordTrajectoryResult.cpp` as the
reference; codex R1 F2 made this the universal pattern.

**2. Per-dataset units attrs (codex R1 F3).** mean/std/min/max in base
units; m2 in squared base units; *_frame in `frame_index`; n_per_atom
in `frame_count`. Group-level `units_vector` / `units_scalar` stay as
overall channel labels.

**3. "Absent, not faked" gate.** Source-attached gate `HasResult<T>()`
in Compute; on absence, TS NaN-fills float channels + sets
`source_attached_per_frame=0`. Welford skips the update entirely.
See R1 F1: source calculators that detect degenerate output (non-finite
gradient) should return nullptr, NOT attach a NaN result — the
HasResult gate is the single boundary.

**4. ProteinConformation::ForceAttachResultForTesting** test helper
landed in S1. Bypasses dependency check. Use only for synthetic-
accumulation tests where the deep dep chain is otherwise impossible
(see Integration1P9J vs synthetic coverage-gap pattern).

**5. Integration test on 1P9J.** The synthetic-positions
accumulation path FATALs at canonicalisation. Use `Trajectory::Run`
on the 1P9J fleet_amber fixture with `GTEST_SKIP() << "..."` when
model/CUDA unavailable. Mirror the pattern in
`tests/test_aimnet2_embedding_charge_response_gradient_time_series.cpp`
TEST(AIMNet2ChargeResponseGradientWelford, Integration1P9J).

**6. Anti-pattern grep self-check before commit.** See
`feedback_anti_pattern_grep_self_check`. Grep PATTERNS.md anti-
patterns on the diff before commit, then adversarial codex review.

## Source result audit (S2 pre-flight)

Before writing any S2 TR, verify these source ConformationResults
exist and have the expected output shape:

- `src/ApbsFieldResult.{h,cpp}` — does `apbs_efg_spherical` exist on
  the per-atom ConformationAtom? T2 SphericalTensor shape (5
  components)?
- `src/MopacResult.{h,cpp}` — per-atom charges (scalar) + per-bond
  bond orders (sparse). Bond axis from `bonds.npy` topology sidecar.
- `src/MopacCoulombResult.{h,cpp}` — T2 SphericalTensor per atom.
- `src/MopacMcConnellResult.{h,cpp}` — T2 SphericalTensor per atom.
- `src/CoulombResult.{h,cpp}` — FF14SB T2 SphericalTensor (the
  reconciliation reference).

If any source is missing or has a different shape than the plan
assumes, FLAG IMMEDIATELY — don't paper over. The plan was written
2026-05-20 without re-auditing source-side state; a stale assumption
about Mopac output shape would derail TR #5/#6 or #7/#8.

## Pre-flight script (≤30 min, then start coding)

```bash
# 1. Verify source results exist
grep -rn 'class ApbsFieldResult\|apbs_efg_spherical' src/ | head -5
grep -rn 'class MopacResult\|mopac_charge\|mopac_bond_order' src/ | head -5
grep -rn 'class MopacCoulombResult\|mopac_coulomb' src/ | head -5
grep -rn 'class MopacMcConnellResult\|mopac_mcconnell' src/ | head -5
grep -rn 'class CoulombResult' src/ | head -5

# 2. Look at the template TRs S1 leaves in place
ls src/BsShieldingTimeSeriesTrajectoryResult.* src/HBondCountWelfordTrajectoryResult.*
ls src/AIMNet2ChargeResponseGradientWelfordTrajectoryResult.*  # canonical full Welford row

# 3. Verify sparse Welford pattern (cadence handling) on Tripeptide BB
ls src/TripeptideBackboneResidualVecTimeSeriesTrajectoryResult.*
ls src/HBondCountWelfordTrajectoryResult.*

# 4. Verify topology sidecar bonds.npy shape (TR #6 bond axis)
ls data/topology/  # or wherever bonds.npy lives
python3 -c "import numpy as np; b = np.load('<bonds.npy path>'); print(b.dtype, b.shape)"
```

After pre-flight, proceed in order TR4 → TR5 → TR6 → TR7 → TR8 → TR9.
TR7/TR8 are clones with the same wiring; codex annealing at the end.

## Open questions to surface to user if found in pre-flight

1. **Mopac cadence default**: plan says "every 20 ps" — is that an
   existing JobSpec stride, or a new TR-side override? Check
   `src/RunOptions.cpp` for sparse-Mopac stride flags.
2. **Bond axis storage**: TR #6 Welford-per-bond — does the
   `WelfordMoments` array layout need a new sidecar emit, or does
   `bonds.npy` already make this trivial? Plan §Risk register flags
   this as low risk, "Topology sidecar already emits bonds.npy".
3. **MopacVsFf14SbReconciliation source-pairing**: gate strictly on
   BOTH attached; emit NaN on asymmetric attach. Already documented
   in plan §Risk register but verify against codex feedback when
   writing TR9.

## What's pushed and where

- Branch: `master` (origin/master at `a3cb51d`)
- Remote: `https://github.com/jessicalh/geo-kernel-extract.git`
- Tests green at HEAD; PR not needed (single-developer thesis tree).
- New session can `git pull` to be at the same state.
