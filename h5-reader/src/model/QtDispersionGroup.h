// QtDispersionGroup — read-side mirror of the SDK's Dispersion
// RingKernelGroup. van der Waals C6 anisotropy over ring vertices,
// 1/r⁶ (isotropic) + 1/r⁸ (anisotropic); same shape as the other
// ring kernels.
//
//   disp_shielding   (N, 9)   SphericalTensor, Å⁻⁶
//   disp_per_type_T0 (N, 8)
//   disp_per_type_T2 (N, 40)

#pragma once

#include "QtRingKernelGroup.h"

namespace h5reader::model {

class QtDispersionGroup : public QtRingKernelGroup {
public:
    explicit QtDispersionGroup(const QtConformationSnapshot& snapshot)
        : QtRingKernelGroup(snapshot,
                            io::FieldKind::DispShielding,
                            io::FieldKind::DispPerTypeT0,
                            io::FieldKind::DispPerTypeT2) {}
};

}  // namespace h5reader::model
