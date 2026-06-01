// Verbs — Layer-1 primitive access verbs for the rediscover functional API.
//
// These are the "little APIs" of SURFACE_DESIGN.md §"Producer side (C++) —
// primitive verbs": pure functions over the immutable Body that return views /
// small values and do NOT rebuild any state. Each verb is a THIN wrapper over
// the resident spine that is already built day-1 (Catalog, SpatialIndexSet,
// RingGeometryCache, TypedAtomIndex, TemporalIndex, the typed QtProtein /
// QtTopology). They EXPOSE the indexes as verbs; they never duplicate or
// regenerate them. (BROAD_BACKBONE_NEXT.md "USE the resident topology — do NOT
// regenerate it"; feedback_identity_from_chemistry_not_position.)
//
// The verbs are the vocabulary the Layer-2 combinators curry over. They take a
// `const Body&` explicitly here; the combinators capture it once and the engine
// iterates the resulting closures (SURFACE_DESIGN.md §"Iterated closures").
//
// Naming maps to the SURFACE_DESIGN verb table:
//   pos      -> conf.atomPosition (atom-major)
//   at       -> AtomState { pos } at (atom, frame)
//   window   -> TemporalIndex.range (the time axis; contiguous via atom-major)
//   value    -> Catalog.value / valueVec3 / valueTensor (uniform catalog read)
//   near     -> SpatialIndexSet.near (per-cloud KD radius query; cutoff REQUIRED)
//   ringSlots-> the H5 ring-neighbourhood slot walk (frozen frame-0 membership)
//   ringGeom -> RingGeometryCache.at (canonical-normal-flipped ring geometry)
//   atomOf   -> TypedAtomIndex.selectUnique (typed, collision-safe; no name scan)
//   selectAll-> TypedAtomIndex.select (the equivalence-class set)
//   ringsOf  -> QtTopology ring memberships for an atom's heavy parent
//
// No QObject; plain math + view returns. Reader owns H5: every H5 read goes
// through the Catalog / the resident Conformation, never a parallel path.

#pragma once

#include "AnalysisBody.h"
#include "Catalog.h"
#include "LocalFrameBasis.h"
#include "RingGeometryCache.h"
#include "SpatialIndexSet.h"
#include "TemporalIndex.h"
#include "TypedAtomIndex.h"

#include "../model/Types.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace h5reader::rediscover {

using model::Vec3;

// The per-(atom, frame) primitive state a verb hands back — currently just the
// lab-frame position, which is what every downstream attacher needs as the
// origin of the displacement vectors. Kept a struct so kernel handles can be
// added without touching call sites (SURFACE_DESIGN Layer-1 "AtomState").
struct AtomState {
    Vec3 pos = Vec3::Zero();
    std::size_t atom = 0;
    std::size_t frame = 0;
};

namespace verbs {

// pos(atom, frame) -> vec3 — atom-major lab position, straight off the resident
// Conformation via the Catalog (reader owns H5).
Vec3 pos(const Body& body, std::size_t atom, std::size_t frame);

// at(atom, frame) -> AtomState — the primitive bundle for one (atom, frame).
AtomState at(const Body& body, std::size_t atom, std::size_t frame);

// window(atom, frame, before[, after]) -> FrameWindow — the time axis as a
// no-copy contiguous span (atom-major storage IS the index; TemporalIndex does
// the arithmetic). The Δ / rate-of-change combinator folds over this.
FrameWindow window(const Body& body, std::size_t atom, std::size_t centerFrame,
                   std::size_t before, std::size_t after = 0);

// near(cloud, atom, frame, cutoff) -> [SourceRef] — KD radius query on the
// day-1 per-cloud trees. cutoff is REQUIRED at the call site (no hidden
// cutoffs; SURFACE_DESIGN Layer-1). The query point is the atom's own position.
std::vector<SourceRef> near(const Body& body, CloudKind cloud, std::size_t atom,
                            std::size_t frame, double cutoff);

// near(cloud, query, frame, cutoff) -> [SourceRef] — same, from an explicit
// query point (e.g. an anchor's position rather than the atom's).
std::vector<SourceRef> nearPoint(const Body& body, CloudKind cloud, const Vec3& query,
                                 std::size_t frame, double cutoff);

// value(arrayId, atom, frame[, slot][, comp]) -> scalar — the uniform catalog
// read. valueVec3 / valueTensor / valueT2 for the wider ranks.
double value(const Body& body, ArrayId id, std::size_t atom, std::size_t frame,
             int slot = -1, int comp = -1);
Vec3 valueVec3(const Body& body, ArrayId id, std::size_t atom, std::size_t frame);
model::SphericalTensor valueTensor(const Body& body, ArrayId id, std::size_t atom,
                                   std::size_t frame);
std::array<double, 5> valueT2(const Body& body, ArrayId id, std::size_t atom, std::size_t frame);
bool present(const Body& body, ArrayId id, std::size_t atom, std::size_t frame);

// ringSlots(atom, frame) -> the H5 ring-neighbourhood per-slot rows for this atom.
// The frozen frame-0 ring-membership backend (NOT a KD query): each used slot
// carries {distance, rho, z, in_plane_angle} at this frame + the absolute
// source ring id. -1 ring id ⇒ unused slot (skipped here). Returns the
// non-empty slots only.
struct RingSlot {
    int32_t ring_index = -1;
    double distance = 0.0;
    double rho = 0.0;
    double z = 0.0;
    double in_plane_angle = 0.0;
};
std::vector<RingSlot> ringSlots(const Body& body, std::size_t atom, std::size_t frame);

// ringGeom(ring, frame) -> canonical (normal-flipped) ring geometry, straight
// off the RingGeometryCache (built day-1; never recomputed here).
const model::RingGeometry& ringGeom(const Body& body, std::size_t ring, std::size_t frame);

// atomOf(scope, selector) -> the unique typed atom in `scope` matching the
// partial typed selector, or nullopt (with err) if 0/N match — the
// collision-safe, IUPAC-trap-safe lookup. selectAll returns the whole
// equivalence class.
std::optional<int32_t> atomOf(const Body& body, const std::vector<int32_t>& scope,
                              const TypedAtomSelector& selector, QString* err_out = nullptr);
std::vector<int32_t> selectAll(const Body& body, const std::vector<int32_t>& scope,
                               const TypedAtomSelector& selector);

// ringsOf(atom) -> the AROMATIC rings the atom's heavy parent belongs to
// (absolute ring ids). The "own ring" set the ring-current self/bonded
// exclusion uses. ownRingAtoms returns the union of those rings' atom indices
// (catches fused partners via a shared bridgehead). Both walk the resident
// topology reverse index, never re-derive connectivity.
std::vector<int32_t> ringsOf(const Body& body, std::size_t atom);
std::vector<int32_t> ownRingAtoms(const Body& body, std::size_t atom);

// ownAromaticRing(atom) -> the first aromatic ring the atom's heavy parent
// belongs to (its own ring, for the aromatic-H frame), or -1.
int ownAromaticRing(const Body& body, std::size_t atom);

}  // namespace verbs

}  // namespace h5reader::rediscover
