// QtSasaGroup — read-side mirror of the SDK's SASA group (SasaResult;
// Shrake-Rupley, ~92-point Fibonacci sphere, 1.4 Å water probe).
//
//   atom_sasa   (N,)   per-ATOM solvent-accessible surface area (Å²)
//   sasa_normal (N,3)  outward unit normal (mean direction of non-occluded
//                      test points); the ZERO vector for fully buried atoms.
//
// NOTE: atom_sasa is the per-ATOM SASA — a DIFFERENT quantity from the DSSP
// per-RESIDUE SASA carried in QtDsspGroup::backbone().sasa. Do not conflate.
// Thin const view; nullopt = SASA did not run this frame ("absent, not faked").

#pragma once

#include "QtConformationSnapshot.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtSasaGroup {
public:
    explicit QtSasaGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // Per-atom Shrake-Rupley SASA (Å²).
    std::optional<double> sasa(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::AtomSASA))
            return std::nullopt;
        return snap_->column(io::FieldKind::AtomSASA).row(atomIdx)[0];
    }

    // Outward surface normal (unit vector; zero vector when fully buried).
    std::optional<Vec3> normal(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::SASANormal))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::SASANormal).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
