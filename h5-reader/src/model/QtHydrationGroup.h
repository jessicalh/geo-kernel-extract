// QtHydrationGroup — read-side mirror of the SDK's hydration-shell group
// (HydrationShellResult). Explicit-water packing geometry around each atom:
// half-shell asymmetry (UCSB method, COM reference), mean water-dipole
// orientation, and nearest-ion distance/charge. These encode the local
// dielectric environment that no geometry-only kernel can see.
//
//   hydration_shell (N,4) → HydrationShell { halfShellAsymmetry,
//                           meanWaterDipoleCos, nearestIonDist, nearestIonCharge }
//
// Thin const view; nullopt = no explicit solvent this frame ("absent, not
// faked"). nearestIonDist is +INFINITY when no ion is within cutoff — callers
// guard via HydrationShell::hasNearestIon().

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtHydrationGroup {
public:
    explicit QtHydrationGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<HydrationShell> shell(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::HydrationShell))
            return std::nullopt;
        return HydrationShell::FromRow(snap_->column(io::FieldKind::HydrationShell).row(atomIdx));
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
