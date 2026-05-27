// QtWaterPolarizationGroup — read-side mirror of the SDK's water-polarisation
// group (HydrationGeometryResult). First-shell water orientation expressed in
// the SASA-derived surface-normal frame (a proper normal, vs the cruder
// protein-COM direction used by QtHydrationGroup). Depends on SASA.
//
//   water_polarization (N,10) → WaterPolarization { dipole(Vec3), normal(Vec3),
//                               asymmetry, alignment, coherence, shellCount }
//
// The `normal` (cols 3-5) is the SASA outward normal, duplicated from
// QtSasaGroup::normal by design (it is the reference frame for interpreting
// the dipole). Thin const view; nullopt = absent this frame ("absent, not faked").

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtWaterPolarizationGroup {
public:
    explicit QtWaterPolarizationGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<WaterPolarization> polarization(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::WaterPolarization))
            return std::nullopt;
        return WaterPolarization::FromRow(snap_->column(io::FieldKind::WaterPolarization).row(atomIdx));
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
