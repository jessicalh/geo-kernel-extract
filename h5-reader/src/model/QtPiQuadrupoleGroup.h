// QtPiQuadrupoleGroup — read-side mirror of the SDK's PiQuadrupole
// RingKernelGroup. π-electron quadrupole EFG (Stone T-tensor), pure-T2
// physics, 1/r⁵ leading decay; same shape as the other ring kernels.
//
//   pq_shielding   (N, 9)   SphericalTensor, Å⁻⁵
//   pq_per_type_T0 (N, 8)   Buckingham A-term scalar per ring type
//   pq_per_type_T2 (N, 40)

#pragma once

#include "QtRingKernelGroup.h"

namespace h5reader::model {

class QtPiQuadrupoleGroup : public QtRingKernelGroup {
public:
    explicit QtPiQuadrupoleGroup(const QtConformationSnapshot& snapshot)
        : QtRingKernelGroup(snapshot,
                            io::FieldKind::PQShielding,
                            io::FieldKind::PQPerTypeT0,
                            io::FieldKind::PQPerTypeT2) {}
};

}  // namespace h5reader::model
