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
        static constexpr io::FieldKind kFields[kMcConnellCategoryCount] = {
            io::FieldKind::MOPACMcBackboneTotal,
            io::FieldKind::MOPACMcSidechainTotal,
            io::FieldKind::MOPACMcAromaticTotal,
            io::FieldKind::MOPACMcNearestCoT2,
            io::FieldKind::MOPACMcNearestCNT2,
        };
        PerBondCategoryT2 out;
        for (std::size_t c = 0; c < kMcConnellCategoryCount; ++c) {
            if (!snap_->has(kFields[c]))
                return std::nullopt;
            out.byCategory[c] =
                UnpackSphericalTensor(snap_->column(kFields[c]).row(atomIdx)).T2;
        }
        return out;
    }

    std::optional<McConnellScalars> scalars(std::size_t atomIdx) const {
        static constexpr io::FieldKind kFields[6] = {
            io::FieldKind::MOPACMcCoSum,
            io::FieldKind::MOPACMcCNSum,
            io::FieldKind::MOPACMcSidechainSum,
            io::FieldKind::MOPACMcAromaticSum,
            io::FieldKind::MOPACMcNearestCoDist,
            io::FieldKind::MOPACMcNearestCNDist,
        };
        double values[6] = {};
        for (std::size_t i = 0; i < 6; ++i) {
            if (!snap_->has(kFields[i]))
                return std::nullopt;
            values[i] = snap_->column(kFields[i]).row(atomIdx)[0];
        }
        return McConnellScalars{
            values[0], values[1], values[2], values[3], values[4], values[5]};
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
