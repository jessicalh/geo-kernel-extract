// ComposedRelationships — ring_current and mcconnell re-expressed as composed
// Relationship bundles (Layer-2 curried closures) run through the Layer-3
// engine, per the brief deliverable #4 and SURFACE_DESIGN.md.
//
// These are the proof the functional API is real: each is BUILT from the verbs
// + the curried combinators, NOT a fourth hand-written walk, and each must
// reproduce its procedural-cell oracle byte-for-byte (the gate). The procedural
// cells (RingCurrentNeighborhood / McConnellNeighborhood) are KEPT as the
// reference oracle the composed path is validated against.
//
// The bundles are DEFINED IN CODE (named, not registered) — the allowed
// "named bundles" of SURFACE_DESIGN §"On pluggable interfaces", with virtual
// dispatch on the typed QtRing / QtAtom objects inside the attachers.

#pragma once

#include "Relationship.h"

namespace h5reader::rediscover {

// The Pople ring-current relationship: aromatic ring-facing H ← aromatic rings
// (slots backend), aromatic-H ring-normal frame (typed CG/CD2 anchor), the
// ring identity / dipolar / source-normal attachers, the self/bonded-aware
// sum_dipolar_all vs _producer_valid reducer, the DFT target.
Relationship MakeRingCurrentRelationship();

// The McConnell bond-anisotropy relationship: backbone amide HN ← anisotropic
// bonds (KD near backend, cutoff_A), HN amide-plane frame, the bond
// identity / axis / dipolar attachers, the plain per-category sum reducer, the
// DFT target. `cutoff_A` is the recorded source-discovery cutoff.
Relationship MakeMcConnellRelationship(double cutoff_A);

}  // namespace h5reader::rediscover
