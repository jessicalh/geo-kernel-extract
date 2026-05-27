// QtBondedGroup — read-side mirror of the SDK's bonded-energy group
// (BondedEnergyResult). Per-atom decomposition of the run's force-field bonded
// energy, evaluated from positions each frame, each interaction's energy split
// evenly among its participating atoms (kJ/mol).
//
//   bonded_energy (N,7) → BondedEnergy { bond, angle, ureyBradley, proper,
//                         improper, cmap, total(=Σ of the six) }
//
// Columns are force-field-agnostic; ureyBradley + cmap are zero for force
// fields lacking those terms (AMBER ff14SB). Thin const view; nullopt =
// bonded energy not computed this frame ("absent, not faked").

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtBondedGroup {
public:
    explicit QtBondedGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<BondedEnergy> energy(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::BondedEnergy))
            return std::nullopt;
        return BondedEnergy::FromRow(snap_->column(io::FieldKind::BondedEnergy).row(atomIdx));
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
