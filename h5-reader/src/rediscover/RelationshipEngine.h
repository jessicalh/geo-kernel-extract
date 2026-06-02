// RelationshipEngine — Layer-3: ONE pure, order-free loop that iterates the
// curried Layer-2 closures over the (atom, frame) index space and folds each
// record into the sink. SURFACE_DESIGN.md §"Layer 3 — the engine":
//
//   for (atom, frame) in index_space:        # order is FREE (pure) → cache choice
//       lf   = frame_fn(body, atom, frame)
//       src  = flatten( sel(body, atom, frame) for sel in selectors )
//       vals = [ attach(body, atom, frame, s, lf) for attach in attachers for s in src ]
//       rec  = reducer( SourceSet(src, vals) )
//       sink.fold(rec)
//
// Iteration = map over the index space + an inner map over sources + a fold
// (the sink). The schedule here is DFT-row-outer × stratum-atom-inner — the
// SAME traversal order the procedural cells used, chosen so per-frame KD trees
// and ring geometry are reused across atoms, AND so the emitted rows land in
// byte-identical order to the oracle (the gate compares files). Referential
// transparency means the order is a cache choice, not a correctness one.
//
// The engine is generic: it knows the closure protocol, not ring/McConnell. A
// relationship is data; running it is this loop. That is the whole point — the
// composed ring_current / mcconnell run through THIS, not a fourth hand walk.

#pragma once

#include "AnalysisBody.h"
#include "RecordSink.h"
#include "Relationship.h"

#include <cstddef>
#include <functional>

namespace h5reader::rediscover {

// ── The carrier seam (#29) ──────────────────────────────────────────────────
// The traversal (stratum → frame_fn → selectors flattened → attachers →
// classifier → source_filter → the attached SourceSlot set) is the ONE thing
// that does not vary between drivers; the CARRIER is the only thing that does.
// `RunTraversal` owns the canonical walk and, per (atom, frame), hands the
// carrier the fully-attached record + the FrameResult + the H5 row / original
// index — everything a reducer needs. THE REDUCER IS A PROPERTY OF THE CARRIER
// (the per-record closure folds the sources and streams them to its sink), not
// of the engine. This is the unification BroadBackbone documented as latent
// (#29): RunRelationship's old coupling to the scalar-sum RecordSink is gone;
// the scalar-sum sink, the broad reducer→BroadBackboneSink, and any future
// carrier are just different PerRecordSink closures over the SAME walk.
//
// It is NOT a plugin/ABC framework: the carrier is a plain std::function the
// driver supplies in code (feedback_no_abstractions). Referential transparency
// of the closures means the DFT-row-outer × stratum-atom-inner schedule is a
// cache choice; it is kept so emitted rows land in the oracle's byte order.
using PerRecordSink = std::function<void(std::size_t atom, std::size_t h5_row,
                                         std::size_t originalIndex, const FrameResult& frame,
                                         const NeighborhoodRecord& record)>;

// Walk one relationship's closures over the body's (atom, frame) index space,
// building the attached NeighborhoodRecord per case and handing it to `carrier`.
// Returns the number of cases emitted. The single definition of the traversal
// every SourceSlot carrier shares.
std::size_t RunTraversal(const Relationship& rel, const Body& body,
                         const PerRecordSink& carrier);

// Run one composed relationship over the body, streaming both row kinds to the
// sink. Returns the number of cases (per-(atom, frame) records) emitted —
// identical to the procedural cell's return so main_extract's logging matches.
// Thin carrier over RunTraversal: the scalar-sum reducer + the RecordSink.
std::size_t RunRelationship(const Relationship& rel, const Body& body, RecordSink& sink);

}  // namespace h5reader::rediscover
