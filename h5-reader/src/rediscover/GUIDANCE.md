# Rediscover — guidance

Top-level orientation for the rediscover work. The detailed class model is
`DESIGN.md`; current state lives in `NOW.md` and `STATE.md`.

This doc is meant to re-establish context for a session or agent coming in
cold. It states what is built, what is only designed, and what is open or
caveated — not a finished story.

It is a living document. It guides the sessions and agents working this
branch, and it should change as emergent issues appear — the "Open /
undecided" list below is the running log; add to it rather than smoothing
it over. Keep it honest and a little loose, not falsely crisp.

## Context anchor

Stage 2 of the thesis, on the single protein 1P9J. We consume `nmr_extract`
output (a trajectory H5 + per-frame data + ORCA DFT shielding); we do not
modify the extractor. This lives in the h5-reader tree on the experimental
branch `h5-reader-pysr-spike`.

## The emit boundary — a note for readers

**The emit work lives on the C++ model side.** The typed C++ model is the spine: it
computes the physics — geometry, the kernels, the spherical (T0/T1/T2) decomposition,
the fixed-literature kernel tensors, the DFT target — and EMITS it into the
per-(atom, frame) substrate. The Python side (the fitters, the symbolic distillation)
only CONSUMES that substrate and fits; it never recomputes the physics. So when a new
analysis or relationship needs a quantity the substrate does not yet carry, the fix is
to **extend the C++ emit** (the spine), not to compute it in Python. This is recurring,
expected, per-cell work — e.g. the per-source orientation vectors (ring normal / bond
axis) and the fixed-literature kernel T2 were emit-extensions on the C++ side, not
Python kernels. The discipline is in `analysis/PATTERNS.md`; the full per-cell workflow
is in `REDISCOVERY_MAP.md`.

## What this is, and is not

**Is:** a headless C++ extractor that turns the resident 1P9J data into a
per-(atom, frame) feature/target substrate, so we can study whether
classical shielding mechanisms — and, if they fall short, something else —
account for the DFT shielding.

**Is not:** production code, fleet tooling, a rewrite, or anything merged
back; and not a commitment to one fitting method. It is a one-shot
experiment to see whether the signal is there.

**Honest framing:** the goal is to *look*, not to prove. If the classical
form does not drop out, that is a finding about the kernel set, not a
failure of the exercise. Nothing here is claimed to work until it runs on
1P9J and we have looked at the output.

## Why this shape

- **The model is the spine.** The physics lives in the reader's typed C++
  model; the extractor queries it. We do not recompute geometry in Python —
  that path has failed repeatedly and produces unmaintainable code.
- **The substrate is method-agnostic because the fitter is open.** The
  problem is permutation-invariant over a variable source set,
  rotation-equivariant, tensor-valued (T2), and a local sum. Candidate
  fitters: ridge on the kernels, scalar symbolic regression (PySR, on T0),
  equivariant/tensor SR, and equivariant sum-pooling (e3nn-class). The
  extraction is the same whichever wins, so we build a general
  substrate and choose the fitter separately, later. The GNN/equivariant
  path is not foreclosed.

## The shape, in brief

1. Load one calcset into an immutable, all-frames-resident `RunData`
   (protein + topology + trajectory H5 + DFT frames + frame map).
2. Walk (target atom, DFT frame). Per record emit: identity + frame; the
   local source neighbourhood (sources + relative displacement vectors +
   identity, in a per-atom local frame); the full DFT shielding tensor
   target (raw 3×3 + library-basis decomposition + σ_iso T0); the
   producer's bare kernels as cross-checks. The source set is un-summed.
3. Write per-case output (CSV, padded source slots). A fitter — a separate,
   open choice — consumes it later.

## The two extractions

- **Ring current (Pople), aromatic ring-facing H.** Sources are nearby
  aromatic rings; their per-frame geometry is *read* from the H5
  ring-neighbourhood (distance, ρ, z, in-plane angle), so little is
  recomputed. Ring membership is the frame-0 snapshot the H5 stores — it is
  not "every ring that ever enters cutoff."
- **McConnell anisotropy, backbone amide HN.** Sources are nearby
  anisotropic bonds, discovered per frame with a nanoflann KD-tree over
  bond midpoints; geometry is computed fresh. The form is `(3cos²θ−1)/r³`
  about the **bond axis** (the `Δχ` lives in the parameter, not the kernel).

## The basis (the foundation, built first)

DFT and kernel tensors must share a component basis before any T2
comparison. We decompose the DFT **raw 3×3** in the **library's** T2 order
(`SphericalBasis::DecomposeLibrary`, matching `src/Types.cpp`), pinned by a
fixture (`3ẑẑᵀ−I` ⇒ `T2[2]=√6`). The reader's display-side
`SphericalDecomposition` uses a different order and is left for the viewer.

- `σ_iso` (T0 = trace/3) is rotation-invariant, so the scalar target is safe
  regardless of frame.
- T2 *component* comparison additionally requires the ORCA tensors to be in
  the same Cartesian frame as the H5 kernels. The loader does no rotation
  check, so **cross-frame T2 comparison is unverified until checked** —
  do not treat T2-component agreement as meaningful until that holds.
- A per-atom local frame expresses the neighbourhood vectors and the target
  tensor consistently, which is what an equivariant fitter needs.

## Historical Build Order

Historical dependency order for the first extractor; the all-atom substrate, partition filters, and Build 3 fit architecture have now landed.

1. `SphericalBasis::DecomposeLibrary` + fixture — **written** (not yet
   compiled/tested in this tree).
2. Per-atom local frame.
3. `RunData` + `RunLoader` + `FrameMap` (validating `frame_index_basis` and
   that frame counts agree).
4. The record + the CSV sink.
5. Ring-current extraction, then McConnell extraction.
6. Headless `main` + CMake target.
7. Run on 1P9J; inspect the output.

The fitter is not part of this build.

## Discipline

- Reuse the reader's typed model and loaders. Additive edits only: keep the
  DFT raw `Mat3` (the parser currently discards it); add the aromatic-H
  predicate (`IsBackboneAmideHydrogen` already exists). The GUI binary is
  untouched.
- Experimental, one-shot branch — no integration target. Review is for
  doability, physics, and coherence, not merge risk.
- Docs are detailed but truthful: no hyperbole, no partial truths; state
  built vs designed vs open vs caveated. Reviewed against overclaiming.
- Do not re-open the decided design; `DESIGN.md` is the detailed record.

## Status (truthful)

- **Built and run:** `h5reader_extract`, the relationship engine, all-atom substrate emit, partition filters, `CaseHunter`, Build 2 partition fit, and Build 3 per-type fit architecture are landed on `h5-reader-pysr-spike`.
- **Current substrate:** `/tmp/rediscover-runs/2026-06-03-per-atom-substrate-build1` is the live emitted substrate; current analysis result dirs include Build 2 and Build 3 under `/tmp/rediscover-runs`.
- **Current fit result:** Build 3 settled that total-T2 per-type-within is the trustworthy read; dia/para split is due diligence, not the primary target.
- **Open / undecided:** between-calculator network, equations table, statistical-position grading, and the future hierarchical/type-interaction model.

## Where the detail is

- Read first, for the vocabulary this doc assumes: `DESIGN.md`, the
  reader's `h5-reader/CLAUDE.md`, and the library `OBJECT_MODEL.md` for the
  kernel / H5 definitions.
- Class model and every type: `DESIGN.md`.
- Per-component/history docs sit in this directory; `NOW.md` is the live pointer.
