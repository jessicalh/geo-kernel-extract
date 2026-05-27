// QtGromacsGroup — read-side mirror of the SDK's GROMACS-energy group
// (GromacsEnergyResult). The .edr energy terms for this frame: PME
// electrostatics, bonded, LJ, thermodynamic state, box, virial + pressure
// tensors, per-group temperatures. Aggregate whole-system quantities the MD
// engine computed with explicit solvent — NOT per-atom.
//
//   gromacs_energy (1,43) → GromacsEnergy (PROTEIN-axis: one row per frame)
//
// PROTEIN-axis: energy() takes NO atom index (reads row 0). nullopt = no EDR
// energy for this frame ("absent, not faked"). The 43-vs-42-column catalog
// off-by-one is documented on the GromacsEnergy block; the loader takes the
// actual NPY shape as truth.

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"

#include "../io/QtFieldCatalog.gen.h"

#include <optional>

namespace h5reader::model {

class QtGromacsGroup {
public:
    explicit QtGromacsGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // Whole-frame GROMACS energy terms. Protein-axis: one row, no atom index.
    std::optional<GromacsEnergy> energy() const {
        if (!snap_->has(io::FieldKind::GromacsEnergy))
            return std::nullopt;
        return GromacsEnergy::FromRow(snap_->column(io::FieldKind::GromacsEnergy).row(0));
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
