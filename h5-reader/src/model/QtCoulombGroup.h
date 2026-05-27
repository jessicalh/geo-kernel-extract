// QtCoulombGroup — read-side mirror of the SDK's CoulombGroup. Vacuum-Coulomb
// E-field + EFG from ff14SB charges. Retired from production (APBS supersedes
// vacuum Coulomb), so present only in --mopac FullFat runs where it feeds the
// MOPAC-vs-FF14SB probe. On a standard run these columns are absent, so every
// accessor returns nullopt — "absent, not faked".
//
//   coulomb_shielding     (N, 9)  SphericalTensor, V/Å²
//   coulomb_E             (N, 3)  Vec3, V/Å
//   coulomb_efg_backbone  (N, 5)  QtEfg (T2-only), V/Å²
//   coulomb_efg_aromatic  (N, 5)  QtEfg (T2-only), V/Å²
//   coulomb_scalars       (N, 4)  E-field scalars, V/Å

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtCoulombGroup {
public:
    explicit QtCoulombGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<SphericalTensor> shielding(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::CoulombShielding))
            return std::nullopt;
        return UnpackSphericalTensor(snap_->column(io::FieldKind::CoulombShielding).row(atomIdx));
    }

    std::optional<Vec3> E(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::CoulombE))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::CoulombE).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

    std::optional<QtEfg> efgBackbone(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::CoulombEFGBackbone))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::CoulombEFGBackbone).row(atomIdx));
    }

    std::optional<QtEfg> efgAromatic(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::CoulombEFGAromatic))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::CoulombEFGAromatic).row(atomIdx));
    }

    std::optional<CoulombScalars> scalars(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::CoulombScalars))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::CoulombScalars).row(atomIdx);
        return CoulombScalars{r[0], r[1], r[2], r[3]};
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
