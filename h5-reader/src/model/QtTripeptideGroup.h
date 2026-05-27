// QtTripeptideGroup — read-side mirror of the SDK's Tripeptide group: the
// ProCS15 tripeptide-database DFT reference shielding (Larsen 2015, ProCS15)
// and its geometric residuals. One read-side group fronts the library's two
// per-atom ConformationResults that share the `Tripeptide` FieldGroup:
//
//   • Backbone  (TripeptideBackboneShieldingResult): the central residue's
//     own shielding, by nearest-pose lookup in the tripeptide DFT table and
//     Kabsch alignment onto the protein N/Cα/C.
//   • Neighbor  (TripeptideNeighborShieldingResult): Larsen 2015 Eq 3 — the
//     i±1 side-chain effect on residue i, via the "AXA-scan reuse" trick
//     (Δσ at the flanking ALA cap of the neighbour's AXA tripeptide vs AAA).
//     The stored tensor is the SUM of the i-1 and i+1 contributions.
//
// All fields are per-ATOM (NativeAxis::Atom). Shielding is a 9-col
// SphericalTensor (ppm), residuals are Vec3 (Å), match distance is a scalar
// (Å), method tag is the TripeptideMethodTag provenance enum.
//
// RESIDUAL-AS-ML-FEATURE: the residual vectors are the displacement between a
// matched tripeptide atom and the protein atom it was assigned to. Direction
// AND magnitude are load-bearing model inputs (the model learns the shielding
// correction from the geometric mismatch), not mere diagnostics — see the
// feedback_residual_as_ml_feature project rule.
//
// SENTINELS (verified against the fixture):
//   • nullopt              — the calculator did not run this frame
//                            ("absent, not faked").
//   • NaN within a present column — this atom had no match / no contribution.
//     Backbone fields are NaN on the 38 unmatched atoms; the neighbour sum is
//     NaN where NEITHER side contributed, and each per-direction residual is
//     NaN where THAT side did not contribute (so the prev / next / sum NaN
//     counts differ — 487 / 488 / 474 in the fixture). Callers test
//     std::isnan (scalars) or Vec3::hasNaN() (vectors); hasBackboneMatch()
//     gives a clean integer-tag gate that avoids poking at NaN.
//   • methodTag NoMatch(0) — the in-band integer twin of a NaN backbone row.

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtTripeptideGroup {
public:
    explicit QtTripeptideGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // ── Backbone (central-residue) shielding ──────────────────────────
    // tripeptide_bb_shielding (N,9), ppm. The matched tripeptide pose's σ,
    // rotated onto the protein frame. NaN row where the residue had no DB
    // match (use hasBackboneMatch() to gate).
    std::optional<SphericalTensor> backboneShielding(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::TripeptideBBShielding))
            return std::nullopt;
        return UnpackSphericalTensor(snap_->column(io::FieldKind::TripeptideBBShielding).row(atomIdx));
    }

    // tripeptide_bb_residual_vec (N,3), Å — matched-atom → protein-atom
    // displacement (ML feature: direction + magnitude). NaN where no match.
    std::optional<Vec3> backboneResidual(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::TripeptideBBResidualVec))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::TripeptideBBResidualVec).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

    // tripeptide_bb_match_distance (N,), Å — |residual| (post-alignment
    // matched-atom distance). Redundant with backboneResidual().norm() but
    // cheap. NaN where no match.
    std::optional<double> backboneMatchDistance(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::TripeptideBBMatchDistance))
            return std::nullopt;
        return snap_->column(io::FieldKind::TripeptideBBMatchDistance).row(atomIdx)[0];
    }

    // tripeptide_bb_method_tag (N,), int8→double — which DFT method produced
    // the matched pose. NoMatch(0) is the in-band no-match sentinel.
    std::optional<TripeptideMethodTag> backboneMethodTag(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::TripeptideBBMethodTag))
            return std::nullopt;
        return static_cast<TripeptideMethodTag>(
            static_cast<std::int8_t>(snap_->column(io::FieldKind::TripeptideBBMethodTag).row(atomIdx)[0]));
    }

    // Clean integer-tag gate: true iff the method-tag column is present and the
    // atom matched a tripeptide pose (tag != NoMatch). Equivalent to a non-NaN
    // backbone shielding row but reads the in-band int sentinel, not NaN.
    bool hasBackboneMatch(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::TripeptideBBMethodTag))
            return false;
        return snap_->column(io::FieldKind::TripeptideBBMethodTag).row(atomIdx)[0] != 0.0;
    }

    // ── Neighbor (i±1 side-chain) shielding — Larsen 2015 Eq 3 ─────────
    // tripeptide_neighbor_shielding (N,9), ppm. SUM of the Δσ_BB^{i-1} and
    // Δσ_BB^{i+1} contributions (Larsen treats them as independent additive
    // terms). NaN where neither neighbour contributed.
    std::optional<SphericalTensor> neighborShielding(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::TripeptideNeighborShielding))
            return std::nullopt;
        return UnpackSphericalTensor(snap_->column(io::FieldKind::TripeptideNeighborShielding).row(atomIdx));
    }

    // tripeptide_neighbor_residual_vec_prev (N,3), Å — residual from the i-1
    // direction. NaN where i-1 did not contribute.
    std::optional<Vec3> neighborResidualPrev(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::TripeptideNeighborResidualVecPrev))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::TripeptideNeighborResidualVecPrev).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

    // tripeptide_neighbor_residual_vec_next (N,3), Å — residual from the i+1
    // direction. NaN where i+1 did not contribute.
    std::optional<Vec3> neighborResidualNext(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::TripeptideNeighborResidualVecNext))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::TripeptideNeighborResidualVecNext).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
