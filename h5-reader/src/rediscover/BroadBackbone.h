// BroadBackbone — the FIRST composed relationship the functional API was built
// to enable: EVERY backbone atom (N / CA / C / O / amide-H / HA, every residue),
// ALL mechanisms at once (rings + anisotropic bonds + charge FIELD), against the
// per-backbone-atom DFT target. The thesis target (backbone chemical shifts).
//
// THE SHAPE (load-bearing, BROAD_BACKBONE_NEXT.md): broad_backbone is a composed
// `Relationship` — `atomsWhere(IsBackbone)` stratum + a backbone-class-dispatching
// LocalFrameFn + a LIST of `nearBackend` selectors (ring-centres / bond-midpoints
// / charge-sites — the GENERAL KD backends, NOT the aromatic-H `slots` backend) +
// per-mechanism attachers branching on the typed source kind. It is run through
// the SAME curried-closure protocol the engine uses for ring/mc — NOT a
// procedural BroadBackboneNeighborhood::extract() walk.
//
// The carrier differs from ring/mc (a heterogeneous two-kind bundle with the
// target-repeat FIX), so it streams to BroadBackboneSink via RunBroadBackbone —
// a sibling runner that drives the SAME `Relationship` closures + a BroadReducer.
// RelationshipEngine::RunRelationship is UNTOUCHED, so the ring/mc oracle parity
// holds. (This is the precise latent narrowness the broad case surfaced: the
// engine's RunRelationship is coupled to the scalar-sum RecordSink carrier;
// composing the heterogeneous broad carrier needs the same closures driven into
// a different sink — documented in HANDOFF / the report, not papered over.)

#pragma once

#include "AnalysisBody.h"
#include "BroadBackboneSink.h"
#include "Relationship.h"

#include <QString>

#include <cstddef>
#include <functional>

namespace h5reader::rediscover {

// The broad reducer folds the attached heterogeneous source set into the
// per-mechanism aggregated features + the local Coulomb field (BroadAggregate).
// State-carrying fold, same role as Relationship::reducer — but the broad
// output shape (per-mechanism + field, not scalar-sum) needs its own type.
using BroadReducer =
    std::function<BroadAggregate(const Body&, std::size_t atom, const FrameResult& frame,
                                 std::size_t h5_row, const std::vector<SourceSlot>&)>;

// A composed broad relationship: the shared `Relationship` bundle (stratum /
// frame_fn / selectors / attachers / source_filter / target_fn — the SAME
// curried closures the engine runs) PLUS the broad reducer + the recorded
// cutoffs. relationship_kind = source_sum.
struct BroadRelationship {
    Relationship rel;          // the composed closures (stratum/frame/selectors/attachers/target)
    BroadReducer broad_reducer;
    double ring_cutoff_A = 0.0;
    double bond_cutoff_A = 0.0;
    double charge_cutoff_A = 0.0;
    double mc_near_field_ratio = 0.0;
    QString charge_source;
};

// Compose broad_backbone. Cutoffs are required + recorded (no-hidden-cutoffs);
// charge_cutoff_A is the swept long-range field reach (6/10/12 Å). charge_source
// is ff14sb (the only loaded static source); exclude_residue drops the target's
// own residue from the charge sum (through-space, as charge_dipole).
BroadRelationship MakeBroadBackboneRelationship(double ring_cutoff_A, double bond_cutoff_A,
                                                double charge_cutoff_A,
                                                const QString& charge_source,
                                                bool exclude_residue,
                                                double mc_near_field_ratio);

// Run the composed broad relationship over the body, streaming the two-kind
// carrier to the sink. Drives the SAME `Relationship` closures as the engine's
// RunRelationship (stratum → frame_fn → selectors → attachers → source_filter),
// then the broad reducer → BroadBackboneSink::Write. Returns the case count.
std::size_t RunBroadBackbone(const BroadRelationship& brel, const Body& body,
                             BroadBackboneSink& sink);

}  // namespace h5reader::rediscover
