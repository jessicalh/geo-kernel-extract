// QtMcConnellGroup — read-side mirror of the SDK's McConnellGroup. Bond
// magnetic-anisotropy shielding: the dominant T2 contributor (bonds are
// everywhere, steep 1/r³ decay, large Δχ for C=O). Full asymmetric tensor
// (T0 + T1 + T2, per the full McConnell M_ab, not just the dipolar kernel).
//
//   mc_shielding   (N, 9)   SphericalTensor, Å⁻³
//   mc_category_T2 (N, 25)  T2 per McConnell category (5 × 5)
//   mc_scalars     (N, 6)   CO/CN/sidechain/aromatic angular sums + nearest dists
//   mc_*_bo        (N, 9)   per-category bond-order producer tensors (forward sum)

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

    std::optional<SphericalTensor> tensor(io::FieldKind field,
                                          std::size_t atomIdx) const {
        if (!snap_->has(field))
            return std::nullopt;
        return UnpackSphericalTensor(snap_->column(field).row(atomIdx));
    }

    std::optional<SphericalTensor> producerBo(io::FieldKind producer,
                                              std::size_t atomIdx) const {
        return tensor(producer, atomIdx);
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
