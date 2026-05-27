// QtApbsGroup — read-side mirror of the SDK's APBS group (ApbsFieldResult).
// Solvated electrostatics from a linearised Poisson-Boltzmann solve: E-field
// and EFG obtained by central-difference of the APBS potential grid, converted
// from native kT/(e·Å) to V/Å (×0.025693 V at 298.15 K) so they are directly
// comparable to the vacuum Coulomb fields. APBS is the canonical solvated
// field (always on in production).
//
//   apbs_E   (N,3)  V/Å   solvated PB E-field
//   apbs_efg (N,5)  V/Å²  solvated PB EFG (T2-only, symmetric-traceless → QtEfg)
//
// Thin const view; nullopt = APBS did not run / failed this frame ("absent,
// not faked" — the library returns no result rather than substituting vacuum).

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtApbsGroup {
public:
    explicit QtApbsGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // Solvated Poisson-Boltzmann E-field (V/Å).
    std::optional<Vec3> E(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::APBSE))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::APBSE).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

    // Solvated Poisson-Boltzmann EFG (V/Å², T2-only).
    std::optional<QtEfg> efg(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::APBSEFG))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::APBSEFG).row(atomIdx));
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
