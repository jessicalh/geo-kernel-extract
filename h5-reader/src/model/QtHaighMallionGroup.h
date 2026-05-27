// QtHaighMallionGroup — read-side mirror of the SDK's HaighMallion
// RingKernelGroup. Surface-integral ring-current model; same shape as
// BiotSavart minus total_B / ring_counts. (BS and HM T2 are ~parallel,
// cosine 0.999 across 279K atoms — a thesis finding the reader can show
// side-by-side via the shared QtRingKernelGroup accessors.)
//
//   hm_shielding   (N, 9)   SphericalTensor, Å⁻¹
//   hm_per_type_T0 (N, 8)
//   hm_per_type_T2 (N, 40)

#pragma once

#include "QtRingKernelGroup.h"

namespace h5reader::model {

class QtHaighMallionGroup : public QtRingKernelGroup {
public:
    explicit QtHaighMallionGroup(const QtConformationSnapshot& snapshot)
        : QtRingKernelGroup(snapshot,
                            io::FieldKind::HMShielding,
                            io::FieldKind::HMPerTypeT0,
                            io::FieldKind::HMPerTypeT2) {}
};

}  // namespace h5reader::model
