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

namespace h5reader::rediscover {

// Run one composed relationship over the body, streaming both row kinds to the
// sink. Returns the number of cases (per-(atom, frame) records) emitted —
// identical to the procedural cell's return so main_extract's logging matches.
std::size_t RunRelationship(const Relationship& rel, const Body& body, RecordSink& sink);

}  // namespace h5reader::rediscover
