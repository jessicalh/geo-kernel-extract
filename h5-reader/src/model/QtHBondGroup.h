// QtHBondGroup — read-side mirror of the SDK's HBondGroup geometric summary
// scalars.
//
//   hbond_scalars   (N, 4)  nearest_dist, 1/r³, count within 3.5 Å, McConnell scalar Σ

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtHBondGroup {
public:
    explicit QtHBondGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<HBondScalars> scalars(std::size_t atomIdx) const {
        const auto& col = snap_->column(io::FieldKind::HBondScalars);
        if (!col.present || col.cols < 4)
            return std::nullopt;
        const double* r = col.row(atomIdx);
        return HBondScalars{r[0], r[1], r[2], r[3]};
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
