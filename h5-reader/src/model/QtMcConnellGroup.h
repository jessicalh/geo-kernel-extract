// QtMcConnellGroup — read-side mirror of the SDK's McConnellGroup. Bond
// magnetic-anisotropy shielding: the dominant T2 contributor (bonds are
// everywhere, steep 1/r³ decay, large Δχ for C=O). Full asymmetric tensor
// (T0 + T1 + T2, per the full McConnell M_ab, not just the dipolar kernel).
//
//   mc_shielding   (N, 9)   SphericalTensor, Å⁻³
//   mc_category_T2 (N, 25)  T2 per McConnell category (5 × 5)
//   mc_scalars     (N, 6)   CO/CN/sidechain/aromatic angular sums + nearest dists

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtMcConnellGroup {
public:
    explicit QtMcConnellGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<SphericalTensor> shielding(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::McShielding))
            return std::nullopt;
        return UnpackSphericalTensor(snap_->column(io::FieldKind::McShielding).row(atomIdx));
    }

    // T2 contribution decomposed by McConnell category (backbone / sidechain /
    // aromatic totals + nearest-CO / nearest-CN).
    std::optional<PerBondCategoryT2> categoryT2(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::McCategoryT2))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::McCategoryT2).row(atomIdx);
        PerBondCategoryT2 out;
        for (std::size_t c = 0; c < kMcConnellCategoryCount; ++c)
            for (std::size_t i = 0; i < 5; ++i)
                out.byCategory[c][i] = r[c * 5 + i];
        return out;
    }

    std::optional<McConnellScalars> scalars(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::McScalars))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::McScalars).row(atomIdx);
        return McConnellScalars{r[0], r[1], r[2], r[3], r[4], r[5]};
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
