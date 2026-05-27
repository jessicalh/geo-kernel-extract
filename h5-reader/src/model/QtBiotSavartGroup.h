// QtBiotSavartGroup — read-side mirror of the SDK's BiotSavartGroup
// (a RingKernelGroup + total_B + ring_counts). Johnson-Bovey double-loop
// ring-current model. Inherits shielding / perTypeT0 / perTypeT2 from
// QtRingKernelGroup and adds the two BS-only fields.
//
// Field provenance (python/nmr_extract/_catalog.py):
//   bs_shielding     (N, 9)   SphericalTensor, ppm·T/nA, sign σ=-dB/dB0
//   bs_per_type_T0   (N, 8)   isotropic per ring type
//   bs_per_type_T2   (N, 40)  T2 per ring type (8 x 5)
//   bs_total_B       (N, 3)   summed secondary B-field vector (Tesla)
//   bs_ring_counts   (N, 4)   aromatic-ring-centre counts in 3/5/8/12 Å shells

#pragma once

#include "QtResultBlocks.h"
#include "QtRingKernelGroup.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtBiotSavartGroup : public QtRingKernelGroup {
public:
    explicit QtBiotSavartGroup(const QtConformationSnapshot& snapshot)
        : QtRingKernelGroup(snapshot,
                            io::FieldKind::BSShielding,
                            io::FieldKind::BSPerTypeT0,
                            io::FieldKind::BSPerTypeT2) {}

    // Summed secondary magnetic field vector at this atom (Tesla).
    std::optional<Vec3> totalB(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::BSTotalB))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::BSTotalB).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }

    // Aromatic-ring-centre counts within 3/5/8/12 Å of this atom.
    std::optional<RingCounts> ringCounts(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::BSRingCounts))
            return std::nullopt;
        const double* r = snap_->column(io::FieldKind::BSRingCounts).row(atomIdx);
        return RingCounts{r[0], r[1], r[2], r[3]};
    }
};

}  // namespace h5reader::model
