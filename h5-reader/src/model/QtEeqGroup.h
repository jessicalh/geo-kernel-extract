// QtEeqGroup — read-side mirror of the SDK's EEQ group (EeqResult).
// Geometry-dependent partial charges from extended electronegativity
// equilibration (Caldeweyher et al. 2019, J. Chem. Phys. 150, 154122,
// DOI 10.1063/1.5090222; D4 parameters, Ohno-Klopman Coulomb kernel). Unlike
// fixed ff14SB charges, these respond to conformation — one N×N solve/frame.
//
//   eeq_charges (N,)  partial charge (elementary charges, e)
//   eeq_cn      (N,)  coordination number (error-function counting; emitted
//                     for traceability of the CN-dependent EN shift)
//
// Thin const view; nullopt = EEQ did not run this frame ("absent, not faked").

#pragma once

#include "QtConformationSnapshot.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtEeqGroup {
public:
    explicit QtEeqGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // Geometry-dependent EEQ partial charge (e).
    std::optional<double> charge(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::EEQCharges))
            return std::nullopt;
        return snap_->column(io::FieldKind::EEQCharges).row(atomIdx)[0];
    }

    // Coordination number (intermediate, for traceability).
    std::optional<double> coordinationNumber(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::EEQCN))
            return std::nullopt;
        return snap_->column(io::FieldKind::EEQCN).row(atomIdx)[0];
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
