// Relationship — Layer-2 relationship combinators as CURRIED CLOSURES, bound
// into a named bundle, per SURFACE_DESIGN.md §"Layer 2 — relationship
// combinators" and §"Iterated closures (keep this — it's the heart)".
//
// The shape, faithfully: each combinator is a closure that captures its CONFIG
// once (currying); the Layer-3 engine ITERATES the resulting closures over the
// (atom, frame) index space. Building a relationship is partial application —
// `near(cloud, cutoff)` returns the configured `(atom,frame)->[src]` closure;
// `dispIn(frameF)` captures the frame; the relationship is just these curried
// closures bound together. NOTHING is threaded through the loop. This is what
// makes "a set of lambdas each with a little API" actually compose; the
// one-off IS this with every closure hand-inlined and un-curried.
//
// C++ reality (SURFACE_DESIGN): the cool levels are std::function capturing
// `const Body&` + config; only the innermost per-source attacher would be
// templated IF profiling bit (it has not; the run is a batch job). The body is
// captured by every closure, not passed at each call. The verbs (Layer 1) are
// the vocabulary these closures curry over.
//
// This is NOT a plugin framework (SURFACE_DESIGN §"On pluggable interfaces"):
// relationships are named bundles DEFINED IN CODE (ComposedRelationships.cpp),
// not discovered or registered. Virtual dispatch on the typed domain objects
// (QtRing virtuals, QtAtom predicates) is the allowed polymorphism and is used
// inside the attachers.

#pragma once

#include "AnalysisBody.h"
#include "LocalFrameBasis.h"
#include "RecordSink.h"
#include "RediscoverTypes.h"
#include "SpatialIndexSet.h"
#include "Verbs.h"

#include <QString>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <vector>

namespace h5reader::model {
struct QtAtom;
}

namespace h5reader::rediscover {

// ── The frame a LocalFrameFn produces ──────────────────────────────────────
// The basis + the chemistry-typed anchor atom whose centroid→atom direction
// fixed the x azimuth (aromatic-H: CG/CD2; -1 for HN). Recorded on the record.
struct FrameResult {
    LocalFrame frame;
    int32_t anchor_atom_index = -1;
};

// ── A raw source, before attachment ────────────────────────────────────────
// The uniform carrier a SourceSelector hands the engine: either a ring slot
// (the H5 frozen-membership backend) or a spatial SourceRef (the KD backend).
// The attachers branch on `kind` (the typed source kind). This is the "more
// lambdas, no special case" heterogeneity of SURFACE_DESIGN: a relationship's
// `selectors` is a LIST, flattened by the engine.
struct RawSource {
    SourceKind kind = SourceKind::Ring;
    // Ring backend (verbs::slots): the frozen-membership slot row.
    verbs::RingSlot ring;
    // KD backend (verbs::near): the spatial reference into a cloud.
    SourceRef ref;
};

// ── The aggregated-row result a Reducer folds the source set into ──────────
// Mirrors RecordSink::WriteAggregatedRow's parameters exactly, so the engine
// can hand the reducer's output straight to the untouched sink — the proof the
// composition reproduces the procedural oracle byte-for-byte.
struct AggregateResult {
    double sum_all = 0.0;            // Σ over all sources (incl. self/bonded)
    double sum_valid = 0.0;          // Σ over the producer-valid set
    int n_valid = 0;                 // count of the valid set
    std::vector<double> per_type;    // per-ring-type / per-bond-category sums (valid set)
    double cutoff_A = 0.0;           // recorded source-discovery cutoff (NaN for H5-backed)
};

// ── The curried closure types (the cool levels: std::function) ─────────────
// Stratum captures only the body; the rest capture body + their config.
using Stratum = std::function<std::vector<std::size_t>(const Body&)>;
using LocalFrameFn = std::function<FrameResult(const Body&, std::size_t atom, std::size_t frame)>;
using SourceSelector =
    std::function<std::vector<RawSource>(const Body&, std::size_t atom, std::size_t frame)>;
// An Attacher fills its part of an already-default-constructed SourceSlot from
// one raw source, in the recorded local frame. (Geometry / identity / kernel
// attachers compose by each writing their fields.)
using Attacher = std::function<void(const Body&, const AtomState&, const FrameResult&,
                                    const RawSource&, SourceSlot&)>;
// A Classifier writes the per-source validity / key onto the attached slot
// (SURFACE_DESIGN's SourceClassifier: `is_self_or_bonded ⊕ type key`). It runs
// AFTER the attachers and BEFORE the per-source rows are written, so the flag
// it sets lands on those rows — exactly as the procedural cell sets the flag
// inline before pushing the slot. May be empty (no classification needed). The
// `prep` closure runs once per (atom) to build any per-atom context the slot
// classifier needs (the own-ring/own-atom sets); it is the state the
// classifier closes over.
struct ClassifierContext {
    // Opaque per-atom scratch the classifier closure interprets. Kept as a
    // shared bag so the combinator stays a plain closure pair, no new type per
    // relationship. ring_current uses ownRings/ownAtoms; mcconnell needs none.
    std::vector<int32_t> ownRings;
    std::vector<int32_t> ownAtoms;
};
using ClassifierPrep = std::function<ClassifierContext(const Body&, std::size_t atom)>;
using Classifier =
    std::function<void(const Body&, const ClassifierContext&, SourceSlot&)>;
// A SourceFilter decides whether an attached slot is KEPT (true) or DROPPED
// (false) before it becomes a per-source row. This is where a selector's
// geometric reject lands — the procedural cell's `continue` on a degenerate
// bond axis / coincident midpoint becomes a filter, so no NaN source row is
// emitted (matching the cell, which never pushes those). Empty ⇒ keep all.
using SourceFilter = std::function<bool(const SourceSlot&)>;
using Reducer = std::function<AggregateResult(const Body&, std::size_t atom,
                                              const std::vector<SourceSlot>&)>;
using TargetFn =
    std::function<DftTarget(const Body&, std::size_t atom, std::size_t originalIndex,
                            const LocalFrame&)>;
// Reads the producer's bare per-atom kernel cross-check into the record (sets
// bare_kernel + bare_kernel_present). Empty closure ⇒ no cross-check.
using BareKernelFn = std::function<void(const Body&, std::size_t atom, std::size_t frame,
                                        NeighborhoodRecord&)>;

// ── A relationship: a named bundle of curried closures ─────────────────────
// `selectors` is a LIST (heterogeneous sources by design). `attachers` is a
// LIST (geometry ⊕ identity ⊕ kernel — each writes its fields). The schema is
// declared alongside (the columns the attachers + reducer emit), per
// SURFACE_DESIGN's per-relationship output carrier.
struct Relationship {
    QString name;
    FeatureSchema schema;
    Stratum stratum;
    LocalFrameFn frame_fn;
    std::vector<SourceSelector> selectors;
    std::vector<Attacher> attachers;
    ClassifierPrep classifier_prep;  // optional (empty ⇒ no per-source classify)
    Classifier classifier;           // optional (empty ⇒ no per-source classify)
    SourceFilter source_filter;      // optional (empty ⇒ keep every attached slot)
    Reducer reducer;
    TargetFn target_fn;
    BareKernelFn bare_kernel_fn;  // optional (empty ⇒ none)
};

// ── Curried verb builders (Layer-2 little-API factories) ───────────────────
// Each returns a configured closure: the config is baked in, the body + (atom,
// frame) arrive at iteration. This is the currying SURFACE_DESIGN insists on.

// Stratum = atomsWhere(pred): the typed stratum filter.
Stratum atomsWhere(std::function<bool(const model::QtAtom&)> pred);

// SourceSelector backends.
//   slotsBackend()            — the H5 ring-neighbourhood frozen-membership walk.
//   nearBackend(cloud, cutoff)— the KD radius query on a per-cloud tree.
SourceSelector slotsBackend();
SourceSelector nearBackend(CloudKind cloud, double cutoff);

}  // namespace h5reader::rediscover
