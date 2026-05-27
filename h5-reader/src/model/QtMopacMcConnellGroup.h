// QtMopacMcConnellGroup — read-side mirror of the SDK's MOPACMcConnell group
// (MopacMcConnellResult). STRUCTURALLY IDENTICAL to QtMcConnellGroup — same
// blocks + column layout (MopacMcConnellResult.cpp:312-330 mirrors
// McConnellResult, including the category order backbone/sidechain/aromatic/
// CO-nearest/CN-nearest) — but the bond magnetic-anisotropy shielding is driven
// by MOPAC PM7 charges/bond-orders rather than ff14SB. Present only on the
// FullFat `--mopac` path; nullopt otherwise ("absent, not faked").
//
//   mopac_mc_shielding   (N, 9)   SphericalTensor, Å⁻³
//   mopac_mc_category_T2 (N, 25)  PerBondCategoryT2 (5 categories × 5 T2)
//   mopac_mc_scalars     (N, 6)   McConnellScalars (CO/CN/sidechain/aromatic
//                                 angular sums + nearest CO/CN distances)

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtMopacMcConnellGroup {
public:
    explicit QtMopacMcConnellGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<SphericalTensor> shielding(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACMcShielding))
            return std::nullopt;
        return UnpackSphericalTensor(snap_->column(io::FieldKind::MOPACMcShielding).row(atomIdx));
    }

    // T2 contribution decomposed by McConnell category (backbone / sidechain /
    // aromatic totals + nearest-CO / nearest-CN), 5 × 5 row-major.
    std::optional<PerBondCategoryT2> categoryT2(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACMcCategoryT2))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::MOPACMcCategoryT2).row(atomIdx);
        PerBondCategoryT2 out;
        for (std::size_t c = 0; c < kMcConnellCategoryCount; ++c)
            for (std::size_t i = 0; i < 5; ++i)
                out.byCategory[c][i] = r[c * 5 + i];
        return out;
    }

    std::optional<McConnellScalars> scalars(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACMcScalars))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::MOPACMcScalars).row(atomIdx);
        return McConnellScalars{r[0], r[1], r[2], r[3], r[4], r[5]};
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
