# Codex review brief — h5-reader buffering / conformation architecture (2026-05-26)

## What to review, and how to respond

You are an adversarial senior Qt/C++ architect. Review the design below —
**especially the one hard-to-reverse decision in §6.** Two independent Claude
reviews already ran (model-clarity lens + YAGNI lens) and converged; their
verdict is in §4. **Build on that — do not just repeat it.** Ground every claim
about existing code in the actual files (§8); never fabricate file:line
citations. Be specific, opinionated, and fair (not contrarian for its own sake).
This is design critique — do NOT modify any files.

The distinction the user cares about most: this is being reviewed *because it is
expensive to undo*. The user said, verbatim, "the conformation is our working
buffer and I want that codex review too on all of this — we can't come back."

## 0. Read the Qt house style FIRST (so your "simpler"/"better" uses the project's own patterns)

- `/home/jessica/.claude/skills/qt6-cpp/SKILL.md` — guiding principle: "crash-diagnosis infrastructure is priority #1 — when this goes wrong at 2am, can you read the dump and say exactly what happened?"
- `/home/jessica/.claude/skills/qt6-cpp/references/architecture.md` (error bus, startup/shutdown, signal/slot)
- `/home/jessica/.claude/skills/qt6-cpp/references/crash-diagnosis.md` (object census, QPointer, sanitizers — directly relevant to a QObject working-buffer + lifetime safety)
- `/home/jessica/.claude/skills/qt6-cpp/references/model-view.md` (large datasets feeding views)
- `/home/jessica/.claude/skills/qt6-cpp/templates/error-bus.h` (Result<T,E> + reporting), `connection-auditor.h` (ACONNECT), `thread-guard.h` (ASSERT_THREAD)
- NOTE: the skill's ENV specifics are Windows/MSVC; the PATTERNS are cross-platform and ARE this reader's house style.

## 1. The project

h5-reader is a standalone Qt6/VTK C++17 trajectory reader for an NMR-shielding
research pipeline. It reads a `trajectory.h5` + a 5-NPY topology sidecar +
per-frame calculator NPYs via its own typed boundary classes. It does NOT link
the physics library; it never writes H5. It is "the player, the interface, and
the exposed formal model" of the trajectory data.

Data shape: ~846 atoms x ~751 sampled frames in the reference fixture. The H5
holds dense per-calculator time-series (N, T, K) for ~56 groups. Separately, each
sampled frame also has ~60-80 per-atom calculator NPYs on disk
(`per_frame_npys/frame_NNNNNN/`) — the full-fidelity detail. A single-pose run
(no H5) has the same per-atom NPYs flat in a run root.

**Verified current state (check it yourself):** the H5 is EAGER-loaded in FULL
(~1.6 GB resident) in `QtProteinLoader::Load`; `QtConformation` owns it by
`unique_ptr`; `QtFrame` is a copyable value view that delegates every accessor to
the resident H5 buffers. The reader is **single-threaded** (no
std::thread/QThread workers; only the playback QTimer). `QtConformationSnapshot`
(the sparse per-frame store) + its ~25 typed group views EXIST and compile but
are **dead code** — no producer (loader), no owner (cache), no consumer (nothing
in `src/app` references a snapshot). The transient-use contract is already written
in `QtResultBlocks.h:429` ("do NOT retain a view across an LRU eviction").

## 2. Required capabilities / nice-to-haves

REQUIRED (must work):
1. Open + render a trajectory; scrub + play frames smoothly.
2. Park on a frame and inspect the FULL per-atom calculator detail.
3. Open a single-pose run (no H5; sidecar + flat per-atom NPYs); render + inspect; one conformation, no playback.

NICE-TO-HAVES (deferrable, must not cruft the basics):
- **Prefetch around the current frame** — the user has COMMITTED to this as the *immediate next increment* after the basics. Treat it as near-certain, not hypothetical.
- Window-scoped streaming so the H5 need not be fully resident (enables µs trajectories on laptops). The user explicitly wants "revisit H5-loaded-in-full" to stay reachable.
- Multi-scale navigation (µs -> ns window -> frame detail).

## 3. Decision history (so you build on it, not from scratch)

A "scoped-lifetime resource model" was proposed: every load DECLARES a scope
(App / Window=frame±n / Frame); consumers ACQUIRE by (field, frame) and get an
opaque HANDLE (observing by default, opt-in pin); scope MANAGERS are QObjects
owning edge signals + a worker; piles are (field x frame)-granular plain data; the
field catalog carries a default scope per field; "it just works" = the consumer
never knows whether data was resident, prefetched, or just read from disk.

Two adversarial Claude reviews shot most of that down (see §4). The team then
adopted a LEAN shape and the user refined it to: **the conformation itself is the
working buffer** (the cache + ready-signal live on the conformation, not a
separate store), with **prefetch committed as the next step**. §5 is that shape;
§6 is the irreversible question it raises.

## 4. The two Claude reviews' convergent verdict (build on this)

"Right idea, wrong spine, too early." Specifically (all verified against files):
- **"Scope" conflates two orthogonal axes** — data SOURCE (sidecar vs H5 vs per-frame NPY dir) and RETENTION policy (forever / sliding window / few-deep LRU). Don't fuse them into one enum; name the stores instead.
- **observing/pin is `weak_ptr`/`shared_ptr` verbatim** — pin = hold the shared_ptr (the LRU physically can't evict what you hold); observe = weak_ptr::lock() -> valid-or-null. No bespoke handle type needed; the standard library enforces it at zero API cost.
- **Opaque acquire-by-(field,frame) fights crash-diagnosis-#1** — selling opacity ("never knows resident vs disk vs evicted") as the payoff is exactly backwards: at 2am a wrongly-null frame is a question you WANT answered. A plain accessor returning `shared_ptr<const>` (null = not resident) + an ErrorBus line at the load seam is legible; the fog is the anti-feature.
- **(field x frame) granularity is fictional at the snapshot layer** — a `frame_NNNNNN/` dir loads as ONE snapshot; the embedding is a column *inside* it (`snap.column(FieldKind::AIMNet2Embedding)`), not independently acquired. Per-field granularity only matters for the deferred windowed-H5, and even there it's premature.
- **catalog default-scope = drift bait** — the catalog is generated (`DO NOT EDIT BY HAND`) from the Python contract, where scope has no meaning; it's a reader runtime policy, not an on-disk array property.
- **edge-signal triple -> one signal** — `snapshotReady(frame)`. "expired" = your weak_ptr went null; "about-to-leave-window" has no v1 consumer.
- **worker thread deferred** — but the user's prefetch commitment pulls async forward, so the API must be worker-ready now (same signature sync->async).
- **The cheapest forward seam** = ONE `shared_ptr<const Snapshot>`-returning accessor + ONE `snapshotReady` signal, and make it the *only* way to obtain a snapshot. Prefetch = call it for t±n. Async = miss-path moves to a worker, signature unchanged. Windowed-H5 later = flip `QtTrajectoryH5`'s buffers behind the existing `QtFrame` getter seam — a load-site change, not a rewrite. The indirection that makes the flip cheap ALREADY EXISTS (the getter seam); scopes are not needed to get it.
- **Minimal v1** = `FrameNpyLoader` (dir-agnostic producer) + a small frame-keyed `shared_ptr<const>` cache + the conformation base-class split (already designed in the single-pose note) + inspector wiring; leave eager-H5 / `QtFrame` / playback untouched.

One correction the verification turned up: the `Result<T,E>` the reviews lean on
is the *skill template's* pattern — NOT yet in the reader (today it uses
`ok`-bool structs like `QtLoadResult`). So adopting it is a choice, not a reuse.

## 5. The PROPOSED shape (review THIS — it incorporates §4 + the user's "conformation is the working buffer")

```cpp
// model/Conformation.h  — abstract base AND the working buffer.
//   QObject: it owns the snapshot cache + emits the ready signal (the manager role
//   lives here, per the user's "conformation is our working buffer"). The PILES
//   (QtConformationSnapshot) stay plain immutable-after-load data behind shared_ptr.
class Conformation : public QObject {
    Q_OBJECT
public:
    const QtProtein* protein() const;
    virtual std::size_t frameCount() const = 0;
    virtual double      timePicoseconds(std::size_t frame) const = 0;
    virtual Vec3        atomPosition(std::size_t frame, std::size_t atom) const = 0;  // rendering + overlays read this

    std::shared_ptr<const QtConformationSnapshot> snapshot(std::size_t frame) const;  // cached-or-null, never blocks
    void requestSnapshot(std::size_t frame);   // ensure loaded. v1: synchronous. next: post to worker. SAME signature.
signals:
    void snapshotReady(std::size_t frame);     // the one signal; also the prefetch/async seam
protected:
    virtual std::shared_ptr<const QtConformationSnapshot> loadSnapshot(std::size_t frame) = 0;  // subclass: where it comes from
    // cache: std::map<size_t, shared_ptr<const Snapshot>> + deque eviction order, small bound (LRU).
};

class TrajectoryConformation final : public Conformation {   // owns QtTrajectoryH5
    Vec3 atomPosition(...) const override;                   // from h5_->positions()
    std::shared_ptr<const QtConformationSnapshot> loadSnapshot(std::size_t f) override; // FrameNpyLoader on per_frame_npys/frame_{originalIndex(f)}/
    QtFrame frame(std::size_t t) const;                      // unchanged dense per-atom H5 view
};

class SingleConformation final : public Conformation {       // no H5; one preloaded snapshot
    std::size_t frameCount() const override { return 1; }
    Vec3 atomPosition(std::size_t /*0*/, std::size_t a) const override;   // from the pose snapshot's FieldKind::Pos column
    std::shared_ptr<const QtConformationSnapshot> loadSnapshot(std::size_t) override; // returns the one preloaded pose
};

// io/FrameNpyLoader  — directory-agnostic producer (trajectory frame dir OR single-pose run root).
//   Reads each NPY header (shape + dtype), widens the 7 non-<f8 arrays, cols -1->1, fails loud on truncation;
//   try/catch at the NPY seam -> ErrorBus::Report -> degraded state. Returns snapshot-or-error.
```

Wiring: `MoleculeScene` + overlays bind to `Conformation` (base) and read
`atomPosition(t,i)` (so single-pose renders + gets overlays via the same seam);
`QtProteinLoader` gains a `LoadPose(runRoot)` sibling to `Load(h5)`; an
open-directory action sniffs the dir (has `trajectory.h5` => trajectory; else
sidecar + flat NPYs => single pose).

## 6. THE hard-to-reverse decision — attack this hardest

**Is "the Conformation IS the working buffer" the right irreversible foundation —
or does it bite later?**

The conformation is a QObject that owns the snapshot LRU + emits `snapshotReady`,
with H5/single subclasses. The user is committed to this being where buffering
lives ("our working buffer"). Pressure-test it:

- When the H5 time-series ALSO becomes windowed/streamed (the deferred-but-wanted
  future), does the conformation become an overloaded do-everything buffer
  managing two different residency policies (snapshot LRU + H5 window)? Does
  putting BOTH on the conformation age well, or would a separate manager (the
  Claude model-clarity reviewer's "named store") have been the cleaner spine?
- The two Claude reviews SPLIT here: one wanted a separate QObject `SnapshotStore`;
  the other put the cache on the conformation base. The user chose
  conformation-as-buffer. Is that the right call, given it's expensive to undo?
- A QObject base with a per-subclass `loadSnapshot` virtual is the one place this
  design uses inheritance for behavior. Is that the right seam, or does the
  project's "direct named code, no pluggable indirection" rule frown on it? Is a
  non-virtual conformation that simply *holds* a small loader-configured cache
  cleaner?
- Thread-safety forward note (from a Claude reviewer, verify): the `shared_ptr`
  control block is atomic, but the cache CONTAINER (map+deque) is not — when the
  committed prefetch moves loading to a worker, the container needs a guard.
  Confirm the proposed shape makes that a one-mutex change, not a redesign.

Sub-questions: (a) prefetch/async — will moving the miss-path to a worker truly
need no API change AND no consumer change? (b) overlays reading positions via the
base `atomPosition()` seam so single-pose gets them — feasible, or genuine tech
debt to flag-and-defer for the overlays that re-evaluate kernels
(QtBiotSavartCalc/QtHaighMallionCalc)? (c) is leaving the eager-H5/QtFrame/playback
path UNTOUCHED in v1 actually safe given the base-class split changes
`QtConformation`'s type (consumers rebind to `Conformation`)?

## 7. Project ethos (honor it)

- VERBATIM: "No pluggable interfaces unless the user asks. Factories and
  abstract-base-class indirection are OFF by default. Direct named code is the norm."
- Qt discipline: every QObject in the object census; `ACONNECT` audited
  connections; `ASSERT_THREAD`; UDP log (port 9997) is the debug channel.
- The "no simplification bias" rule is PHYSICS/DATA only (never collapse the
  per-atom rank-2 tensor complexity — that's the thesis). It does NOT license
  code-architecture complexity. Push hard for CODE simplicity.
- DECIDED, do not relitigate: piles are plain immutable-after-load data behind
  shared_ptr/weak_ptr (no thread affinity -> lock-free worker->GUI handoff); the
  MANAGER is the QObject. The open question (§6) is only WHERE the manager lives
  (on the conformation vs separate).

## 8. Files to read (do not modify)

- `/shared/2026Thesis/nmr-shielding/h5-reader/notes/VISION_AND_PROGRESS.md`
- `/shared/2026Thesis/nmr-shielding/h5-reader/notes/SINGLE_POSE_AND_ORCA_DESIGN_2026-05-26.md`
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/model/QtConformationSnapshot.h`
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/model/QtConformation.h` (and `.cpp`)
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/model/QtFrame.h` (and `.cpp`)
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/io/QtTrajectoryH5.h`
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/io/QtProteinLoader.cpp`
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/app/MoleculeScene.cpp` (consumer that binds to the conformation)
- `/shared/2026Thesis/nmr-shielding/h5-reader/src/app/QtAtomInspectorDock.cpp` (the real "park + inspect detail" consumer)
- `/shared/2026Thesis/nmr-shielding/h5-reader/CLAUDE.md`

## Deliverable

1. Your verdict on §6 (the irreversible call) — keep / change, and exactly why, tied to the windowed-H5 future and the verbatim no-indirection rule.
2. Any place the §5 shape still over-builds or under-builds for the 3 required capabilities + the committed prefetch.
3. The single smallest change to §5 you'd make before code is written, if any.
4. Distinguish "the PLAN proposes X" from "the CODE currently does Y" throughout.
