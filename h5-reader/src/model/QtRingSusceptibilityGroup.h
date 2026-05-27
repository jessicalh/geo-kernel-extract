// QtRingSusceptibilityGroup — read-side mirror of the SDK's
// ring_susceptibility group. The aromatic ring as a point magnetic dipole
// at the ring centre (full McConnell form with b̂ → ring normal). Only a
// summed shielding tensor is emitted (ringchi_shielding); no per-ring-type
// decomposition, so this is a standalone group rather than a
// QtRingKernelGroup subclass.
//
//   ringchi_shielding (N, 9)  SphericalTensor, Å⁻³ (T0+T1+T2; asymmetric)

#pragma once

#include "QtConformationSnapshot.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtRingSusceptibilityGroup {
public:
    explicit QtRingSusceptibilityGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<SphericalTensor> shielding(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::RingChiShielding))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::RingChiShielding).row(atomIdx);
        SphericalTensor st;
        st.T0 = r[0];
        for (std::size_t i = 0; i < 3; ++i)
            st.T1[i] = r[1 + i];
        for (std::size_t i = 0; i < 5; ++i)
            st.T2[i] = r[4 + i];
        return st;
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
