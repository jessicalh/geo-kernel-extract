# TrajectoryProtein Write-Surface Enforcement — 2026-04-24

**Status:** pre-code design sketch. Locked this in before any of the
40+ future TRs lands, so the compiler enforces the buffer discipline
instead of humans + AIs having to remember it every clone.

Companion: `feedback_object_model_scope_discipline` (memory,
full behavioural rule) and
`spec/TRAJECTORY_RESULT_PLAN_2026-04-24.md` (the TR plan this
enforcement supports).

---

## The rule, restated

A `TrajectoryResult` runs per frame. It is given:

- `const ProteinConformation& conf` — this frame's snapshot, already
  populated by the `ConformationResult` pipeline. **Read-only from the
  TR's perspective.**
- `TrajectoryProtein& tp` — the time-spanning protein context. TR
  writes go here, but only on the **legitimate write surface** (TrajectoryAtom
  rollup, dense buffer transfer at Finalize, per-atom events).
- `Trajectory& traj` — process/run state. TR can push to
  `traj.MutableSelections()` for run-scope events.

A TR **must not** mutate `Protein`, `ProteinConformation`, or any
`ConformationAtom`. Those are the instant-observation buffer —
universal across static PDB and trajectory use cases. Mutation breaks:

- The one-way data flow invariant (PATTERNS "Nothing flows backward").
- The static-PDB use case (single-conformation consumers get different
  results depending on whether a trajectory TR is also attached).
- Testability (CR and TR can no longer test independently).
- The duplication-over-chaining discipline (PATTERNS §17).

---

## Current enforcement state

The compiler already does half:

```cpp
// src/TrajectoryResult.h:64-68
virtual void Compute(const ProteinConformation& conf,   // <-- const ✓
                     TrajectoryProtein& tp,              // <-- non-const
                     Trajectory& traj,
                     std::size_t frame_idx,
                     double time_ps) = 0;
```

Through `conf` a TR cannot mutate. But the non-const `tp` parameter
re-opens two paths:

1. **`TrajectoryProtein::ProteinRef()`** currently has two overloads:

   ```cpp
   // src/TrajectoryProtein.h:61-62
   Protein& ProteinRef() { return *protein_; }
   const Protein& ProteinRef() const { return *protein_; }
   ```

   A TR has non-const `tp`, so it calls the non-const overload and
   gets `Protein&`. It could then mutate the invariant Protein —
   adding a residue, changing a bond, editing an atom. Every Protein
   method that accepts this is a silent loaded gun in TR code.

2. **`TrajectoryProtein::CanonicalConformation()`** same shape:

   ```cpp
   // src/TrajectoryProtein.h:51-52
   ProteinConformation& CanonicalConformation();
   const ProteinConformation& CanonicalConformation() const;
   ```

   Non-const return reachable from TR. `conf0`'s `ConformationAtom`
   fields could be mutated through this path, bypassing the const
   protection on the per-frame `conf` parameter.

Everything else on `TrajectoryProtein` that a TR legitimately calls
(`MutableAtomAt`, `AdoptDenseBuffer<T>`, `Atoms`, `AtomCount`,
`ProteinId`, etc.) is fine by design.

---

## Proposed enforcement

Close both holes by making the two non-const overloads private +
granting `Trajectory` friend access. Trajectory::Run (the
orchestrator) is the only legitimate caller of the non-const paths —
it needs `CanonicalConformation()` non-const at Phase 6 so
`OperationRunner::Run` can attach CRs to conf0.

Diff sketch:

```cpp
// src/TrajectoryProtein.h (proposed)
class TrajectoryProtein {
    friend class Trajectory;   // <-- new: only Trajectory::Run calls mutable paths.

public:
    // ── Identity delegation (const-only to TR context) ───────────
    const Protein& ProteinRef() const { return *protein_; }
    // Removed: Protein& ProteinRef();

    // ── Canonical conformation (const-only to TR context) ────────
    const ProteinConformation& CanonicalConformation() const;
    // Removed: ProteinConformation& CanonicalConformation();

    // ── Legit TR write surface — unchanged ───────────────────────
    TrajectoryAtom& MutableAtomAt(size_t i) { return atoms_[i]; }
    template <typename T>
    void AdoptDenseBuffer(std::unique_ptr<DenseBuffer<T>> buffer,
                          std::type_index owner);
    // ... all other public methods unchanged.

private:
    // Non-const access for Trajectory::Run orchestration only.
    // Callable via friend; TRs cannot reach these.
    Protein& MutableProtein_() { return *protein_; }
    ProteinConformation& MutableCanonicalConformation_();
};
```

Callers affected:

- `src/Trajectory.cpp` — one Phase 6 call site currently uses
  `tp.CanonicalConformation()` expecting non-const:
  `OperationRunner::Run(tp.CanonicalConformation(), frame_opts)`.
  Updates to `OperationRunner::Run(tp.MutableCanonicalConformation_(), frame_opts)`.
- Any TR that currently calls `tp.ProteinRef()` non-const or
  `tp.CanonicalConformation()` non-const: must switch to the const
  path (usually already const-compatible) or is doing something it
  shouldn't.

A grep of the current tree (2026-04-24, post-catch-up) shows only
`Trajectory.cpp` uses `tp.CanonicalConformation()`; none of the six
landed TRs touch it. `tp.ProteinRef()` is called from several TRs
but only in read contexts (`.AtomAt(i)`, `.BondCount()`,
`.ResidueAt(r)`); the const overload covers those.

---

## Why friend, specifically

Alternatives considered:

- **Separate `TrajectoryWriteSurface` value type** passed to
  `Compute` instead of `TrajectoryProtein&`. Exposes only the four
  writable paths (MutableAtomAt, AdoptDenseBuffer, and the event/
  selections bag pointers). Cleaner at the call site but requires
  every TR to update Compute's signature. Bigger diff, bigger churn
  risk across 6 existing TRs.
- **Everything on `TrajectoryProtein` const-public, with mutation
  via free functions in an `internal::` namespace**. Works, but
  indirects the writing pattern away from the obvious `tp.Method()`
  form. Less readable.
- **Static assert on a concept**. Fragile, harder to maintain.

Friend is the smallest viable change. `Trajectory` is already the
orchestrator that legitimately mutates `TrajectoryProtein`; naming it
as the only class granted mutable-access privilege matches that role
exactly. No behaviour change, just visibility tightening.

---

## Rollout

Small surgical patch:

1. Edit `src/TrajectoryProtein.h` — privatise the two non-const
   accessors, rename to trailing-underscore private helpers, add
   `friend class Trajectory;`.
2. Edit `src/TrajectoryProtein.cpp` — rename any definition bodies
   that moved.
3. Edit `src/Trajectory.cpp` — one call site updates to the new
   private-helper name.
4. Build + smoke_tests — baseline unchanged (build clean + 28 NPY
   binary-diff drift per
   `project_smoke_test_bless_deferred_20260424`).

Reversible: if anything unexpected surfaces, the public const
overloads stay backward compatible with current const-callers; only
the non-const call sites would error, and those are the ones being
deliberately constrained.

---

## Negative-test demonstration (recommended, not shipped)

Sanity-check the enforcement once by dropping a throwaway TR stub
into a test build that deliberately tries each prohibited path and
verifying the compiler rejects it:

```cpp
// Throwaway — should NOT compile after enforcement lands.
void ProofTR::Compute(const ProteinConformation& conf,
                      TrajectoryProtein& tp,
                      Trajectory& traj,
                      std::size_t i, double t) {
    // (1) Mutate through conf — already blocked by const on conf.
    // conf.AtomAt(42).bs_t0_mean = 0.0;  // error: cannot assign

    // (2) Reach the invariant Protein mutably via tp.
    tp.ProteinRef().SomeMutatingMethod();  // error: ProteinRef is const-only

    // (3) Reach canonical conformation mutably via tp.
    tp.CanonicalConformation().AttachResult(nullptr);  // error: const
}
```

Each line should produce a compile error. Record the errors in the
doc if surprise; delete the stub after verifying.

---

## Relationship to PATTERNS.md + OBJECT_MODEL.md

This doc is pre-code. Once the patch lands, a one-paragraph note
appears in PATTERNS.md §14 addition ("`TrajectoryProtein` exposes only
const accessors to `Protein` / `ProteinConformation`; `Trajectory` is
the sole friend with mutable access, enforcing the instant-buffer
discipline at compile time").

OBJECT_MODEL.md trajectory-scope addition tables update the
`TrajectoryProtein` method rows: `ProteinRef()` and
`CanonicalConformation()` columns say *const-only from TR context;
Trajectory::Run uses private helpers via friend*.

Both doc updates ship in the same commit as the code change, per the
one-commit-per-topic convention.
