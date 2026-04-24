# Cold-Read Review Follow-up — 2026-04-24 (evening)

A fresh-reader evaluation of whether the 2026-04-24 landing — compile-time
write-surface enforcement, shadow-orchestration removal, Layer 0 test
template — leaves the tree in a state where a new contributor can cleanly
add a `TrajectoryResult` subclass with three discipline tests. Companion
and style reference: `spec/cold_read_review_20260424.md` (earlier today).

---

## 1. Reading path actually taken

In INDEX order: `INDEX.md:56-107` → `TRAJECTORY_RESULT_PLAN_2026-04-24.md` →
`TRAJECTORY_WRITE_SURFACE_2026-04-24.md` → `PATTERNS.md §§13-18 (:200-319)`
+ addition `:1062-1293` (§14 enforcement at `:1143-1154`, §15 test-template
at `:1194-1202`) → `OBJECT_MODEL.md :2476-3003` →
`src/TrajectoryProtein.h:28-177` → `src/BondLengthStatsTrajectoryResult.{h,cpp}` →
`src/BsWelfordTrajectoryResult.{h,cpp}` →
`tests/test_gromacs_streaming.cpp:626-903`. All clear; two live PATTERNS /
OBJECT_MODEL sections remain a tax but the addition wins per stated policy.

---

## 2. TR addability: **yes, with one doc-drift caveat**

For `HmWelfordTrajectoryResult` (Plan §3 Layer 1 clone: swap
`bs_shielding_contribution` → `hm_shielding_contribution`, dep becomes
`HaighMallionResult`), the essentials are all present:

- ABC surface + Compute/Finalize/WriteH5Group contract:
  `OBJECT_MODEL.md:2616-2660` table + `TrajectoryResult.h` header.
- Factory + Seed-precedes-attach rationale:
  `BondLengthStatsTrajectoryResult.cpp:18-23` + `PATTERNS.md:1156-1167`.
- AV lifecycle, one-writer-per-field, accumulator-on-Result:
  `PATTERNS.md:1111-1141` + `BsWelfordTrajectoryResult.h:7-32`.
- Factory registration site: neither INDEX nor Plan names the specific
  function; a contributor reads `RunConfiguration.cpp:120-157` and
  figures it out. Gap flagged in prior review, not closed.
- `TrajectoryAtom` field extension: no doc tells the contributor "edit
  `TrajectoryAtom.h`"; inferable from `PATTERNS.md:1098-1103` but
  `pending_decisions_20260423.md` `AllWelfords` revival is the right
  fix (`TRAJECTORY_RESULT_PLAN_2026-04-24.md:144-146`).

**Critical read path:** INDEX block → PATTERNS §§14-15 addition →
OBJECT_MODEL addition (TrajectoryProtein + TrajectoryResult tables) →
one exemplar (BondLengthStats if no deps; BsWelford if CR dep) →
`tests/test_gromacs_streaming.cpp:626-903` as clone target. Plan and
WRITE_SURFACE are strategic context, not critical path for one clone.

**Doc-drift caveat.** `TRAJECTORY_RESULT_PLAN_2026-04-24.md:108-115`
Layer 0 calls for `tests/test_bond_length_stats_trajectory_result.cpp`
+ `tests/test_bs_welford_trajectory_result.cpp` as separate files.
Actual landed tests are inline in `tests/test_gromacs_streaming.cpp`
at `:626, 705, 764, 822`. INDEX.md:100 and PATTERNS.md:1202 correctly
point at `test_gromacs_streaming.cpp:BondLengthStats*`; the Plan is
stale. `test_bs_welford_trajectory_result.cpp` does not exist either.
Not blocking (INDEX lands the contributor in the right file) but
reconcile the Plan or extract the file.

---

## 3. Scope-discipline clarity: **strong, and surfaces early**

The instant-buffer-vs-time-buffer rule surfaces in three places,
each short, each re-grounded against its document's frame:

- `INDEX.md:92-96` — checklist bullet "Do not mutate
  `ProteinConformation`, `ConformationAtom`, or `Protein` through
  `tp`." Names mechanism in one sentence ("const-only public API;
  `Trajectory` is the sole friend with mutable access"). Links to
  WRITE_SURFACE. **First encounter on the nominal path.**
- `TRAJECTORY_RESULT_PLAN_2026-04-24.md:71-86` — 8-row buffer-role
  table with the user's verbatim 2026-04-24 rule quoted, plus the
  static-PDB test ("could a user get this quantity without a
  trajectory?").
- `PATTERNS.md:1143-1154` §14 addendum — one paragraph, names
  mechanism, points at spec doc.

Repetition is fine. A reader who stops at INDEX already knows not to
mutate; further reading adds motivation (Plan) and mechanism
(PATTERNS + OBJECT_MODEL). Minor terminology drift: the phrase
"instant-observation buffer" appears in WRITE_SURFACE and PATTERNS §14
addition, once in the Plan (`:73`), not in INDEX. Grep-searchable but
slightly inconsistent.

---

## 4. Enforcement visibility from docs alone: **mostly, with one gap**

Docs state *what* is enforced and *by what mechanism* clearly enough
to prevent the mistake: `INDEX.md:92-96`, `PATTERNS.md:1143-1154`
(§14 addition), `OBJECT_MODEL.md:2543`, and the WRITE_SURFACE diff
sketch at `:86-122` all name the private `MutableCanonicalConformation_()`
+ `friend class Trajectory` mechanism.

The gap: **no single table of the legitimate TR write surface.**
`MutableAtomAt`, `AdoptDenseBuffer<T>`, events bag, and
`traj.MutableSelections()` are scattered across `OBJECT_MODEL.md:2546,
2556, 2642-2644` and `PATTERNS.md:1150-1152`. By reading
`TrajectoryProtein.h:92-95, 142-156` a contributor sees them in one
place, and the header is on the critical path anyway — but a one-table
summary in the OBJECT_MODEL TrajectoryProtein section would close the
loop without requiring the header dive.

**Good on balance.** The rule is crisp enough at INDEX that the reader
won't try the forbidden path; the mechanism is described in enough
detail (private helper + friend) that opening the header confirms
rather than surprises. For trajectory scope the convention is that the
header is read anyway — the docs are not trying to substitute.

---

## 5. Test-template clarity: **clear, with one structural nit**

The four tests at `tests/test_gromacs_streaming.cpp:626-903` are clean
clone targets:

- `BondLengthStatsEndToEnd` (`:626-695`) — real fixture, physical
  plausibility (0.5 < mean < 3.0 Å), `any_moved` non-zero-std check.
- `BondLengthStatsFrame0Semantics` (`:705-761`) — `SetStride(99999)`
  dispatches only frame 0; AV one-sample invariants (min==mean==max,
  std==0, delta_n==0). Comment at `:717-719` names the technique.
- `BondLengthStatsFinalizeIdempotency` (`:764-819`) — call `Finalize`
  twice, compare field-by-field.
- `BondLengthStatsH5RoundTrip` (`:822-903`) — write temp H5, re-read,
  compare attributes + every dataset bit-for-bit.

Each is ~60-80 lines, three-section shape (Arrange = Config +
RunConfiguration; Act = `Trajectory::Run`; Assert = discipline
invariant). Clone recipe is visible: rename `BondLengthStats` →
`HmWelford`, swap `typeid(GeometryResult)` → `typeid(HaighMallionResult)`,
keep `skip_mopac = true`, shift assertions from `PerBond()` to
`tp.AtomAt(i).hm_t0_mean`.

**Structural nit.** Four tests in one file with `BondLengthStats`
prefix; for 40+ forthcoming TRs this file would balloon. Plan §3.1
called for per-TR files. Follow-up: either (a) extract to
`test_bond_length_stats_trajectory_result.cpp` (matches Plan);
(b) update Plan to match reality; or (c) keep current shape with a
clearer prefix convention. Most visible drift from the 2026-04-24
landing. Does not block a single clone.

---

## 6. Doc-vs-code drift after 2026-04-24 landing

1. **Header enforcement matches docs.** `TrajectoryProtein.h:28-34`
   has `friend class Trajectory`; `:61` is const-only
   `CanonicalConformation()`; `:73` is const-only `ProteinRef()`;
   `:173-177` is the private `MutableCanonicalConformation_()`. All
   match the WRITE_SURFACE sketch at `:97-122`. `Trajectory.cpp:178`
   is the single non-const caller.

2. **Layer 0 file location drift** — §5 above. INDEX + PATTERNS
   point at `test_gromacs_streaming.cpp`; Plan still says separate
   files.

3. **Plan sequencing step 2** (`:356-360`) calls for BondLengthStats
   then BsWelford test files. BondLengthStats landed inline;
   BsWelford has no test file at all. Either land the BsWelford test
   or revise the Plan.

4. **8 phases** — PATTERNS §15 addition, OBJECT_MODEL `:2712-2723`,
   and `Trajectory.cpp` comments all agree. Good.

5. **CROSS-RESULT READ markers.** Confirmed in prior review at
   `BsWelfordTrajectoryResult.h:7-17` (writer side) +
   `BsAnomalousAtomMarkerTrajectoryResult.{h,cpp}` (reader side).
   Not re-read; no reason to expect drift.

---

## 7. Recommended onboarding for `HmWelfordTrajectoryResult`

For "clone BsWelford against `hm_shielding_contribution`":

1. **`INDEX.md:56-107`** — 10-bullet checklist. 5 min.
2. **`src/TrajectoryResult.h`** — 102-line ABC header docstring. 5 min.
3. **`src/TrajectoryProtein.h:28-177`** — legitimate write surface
   (`MutableAtomAt`, `AdoptDenseBuffer<T>`, events bag) + `friend
   class Trajectory` comment at `:28-34`. 10 min.
4. **`src/BsWelfordTrajectoryResult.h` + `.cpp`** — direct clone
   target; header `:7-32` names every field written, read sources,
   AV lifecycle. 20 min.
5. **`PATTERNS.md §§14-18 addition (:1111-1293)`** — skim for the
   rules that bite: one writer per field, accumulator on TR,
   duplicate over chain, CROSS-RESULT READ marker discipline if
   needed. 15 min.
6. **`tests/test_gromacs_streaming.cpp:626-903`** — four tests;
   three clone nearly verbatim with name substitution. 15 min.
7. **`OBJECT_MODEL.md:2616-2686`** — TR ABC table + current TRs
   table. Referenced, not read end-to-end.

Skip for a single clone (saves ~1 hour): Plan end-to-end (read
only the Layer 1 table `:126-137`); WRITE_SURFACE end-to-end
(PATTERNS §14 addition has the facts); `pending_include_*`
(flagged non-authoritative); original PATTERNS §§13-18
(:200-319) and OBJECT_MODEL earlier trajectory sections — the
additions win per stated policy.

---

## Overall

A cold reader can confidently add `HmWelfordTrajectoryResult` with
its three discipline tests using only the INDEX block, BsWelford as
exemplar, the four BondLengthStats tests, and `TrajectoryProtein.h`.
The 2026-04-24 enforcement and test template closed the biggest gap
the earlier review flagged (no TR tests; no contract for what a TR
test should assert). Compile-time enforcement now prevents the
canonical mutate-conf0 mistake.

Two residuals: (a) Plan says "separate test files per TR" but tree
has four inline tests in `test_gromacs_streaming.cpp`. (b) BsWelford
Layer 0 test promised by Plan §3.1 item 2 is not landed. Neither is
blocking; both should be reconciled in a next-pass cleanup.
