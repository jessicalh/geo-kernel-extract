// QtWaterFieldGroup — read-side mirror of the SDK's explicit-solvent field
// group (WaterFieldResult). Same Coulomb kernel as CoulombResult but summed
// over explicit water charges (TIP3P: O = −0.834e, H = +0.417e). This is the
// explicit field that APBS approximates with a continuum dielectric — it
// carries water orientation, cavities, bridging/structural water.
//
//   water_efield (_first) (N,3)  V/Å   total / first-shell (<3.5 Å) E-field
//   water_efg    (_first) (N,5)  V/Å²  T2-only EFG (QtEfg; T0,T1 structural 0)
//   water_shell_counts    (N,2)  first / second-shell water-O counts
//
// Thin const view; nullopt = no explicit solvent this frame ("absent, not
// faked"). EFG is symmetric-traceless → T2-only, decoded via UnpackEfg.

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtWaterFieldGroup {
public:
    explicit QtWaterFieldGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // Total water E-field at this atom (V/Å), summed to the 15 Å cutoff.
    std::optional<Vec3> efield(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::WaterEfield))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::WaterEfield).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

    // First-shell-only (<3.5 Å) water E-field (V/Å).
    std::optional<Vec3> efieldFirst(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::WaterEfieldFirst))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::WaterEfieldFirst).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

    // Total water EFG (V/Å², T2-only).
    std::optional<QtEfg> efg(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::WaterEFG))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::WaterEFG).row(atomIdx));
    }

    // First-shell-only water EFG (V/Å², T2-only).
    std::optional<QtEfg> efgFirst(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::WaterEFGFirst))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::WaterEFGFirst).row(atomIdx));
    }

    // First / second-shell water-oxygen counts.
    std::optional<WaterShellCounts> shellCounts(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::WaterShellCounts))
            return std::nullopt;
        return WaterShellCounts::FromRow(snap_->column(io::FieldKind::WaterShellCounts).row(atomIdx));
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
