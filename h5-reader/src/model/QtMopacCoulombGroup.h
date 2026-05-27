// QtMopacCoulombGroup — read-side mirror of the SDK's MOPACCoulomb group
// (MopacCoulombResult). STRUCTURALLY IDENTICAL to QtCoulombGroup — same blocks,
// same column layout (MopacCoulombResult.cpp:312-338 mirrors CoulombResult) —
// but the E-field + EFG are built from MOPAC PM7 charges instead of ff14SB.
// The ff14SB-vs-MOPAC shielding/EFG residual is exactly the charge-source
// coordinate the FullFat `--mopac` probe exists to measure. Present only on
// that path; every accessor returns nullopt otherwise ("absent, not faked").
//
//   mopac_coulomb_shielding     (N, 9)  SphericalTensor, V/Å²
//   mopac_coulomb_E             (N, 3)  Vec3, V/Å
//   mopac_coulomb_efg_backbone  (N, 5)  QtEfg (T2-only), V/Å²
//   mopac_coulomb_efg_aromatic  (N, 5)  QtEfg (T2-only), V/Å²
//   mopac_coulomb_scalars       (N, 4)  CoulombScalars (the 3rd field is a
//                                       SIGNED projection V/Å, not a fraction)

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtMopacCoulombGroup {
public:
    explicit QtMopacCoulombGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<SphericalTensor> shielding(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACCoulombShielding))
            return std::nullopt;
        return UnpackSphericalTensor(snap_->column(io::FieldKind::MOPACCoulombShielding).row(atomIdx));
    }

    std::optional<Vec3> E(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACCoulombE))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::MOPACCoulombE).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

    std::optional<QtEfg> efgBackbone(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACCoulombEFGBackbone))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::MOPACCoulombEFGBackbone).row(atomIdx));
    }

    std::optional<QtEfg> efgAromatic(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACCoulombEFGAromatic))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::MOPACCoulombEFGAromatic).row(atomIdx));
    }

    std::optional<CoulombScalars> scalars(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACCoulombScalars))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::MOPACCoulombScalars).row(atomIdx);
        return CoulombScalars{r[0], r[1], r[2], r[3]};
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
