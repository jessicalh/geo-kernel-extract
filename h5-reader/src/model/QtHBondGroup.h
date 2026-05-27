// QtHBondGroup — read-side mirror of the SDK's HBondGroup. Dipolar H-bond
// shielding (full McConnell form with the coupling direction = D-H…A), plus
// the geometric summary scalars.
//
//   hbond_shielding (N, 9)  SphericalTensor, Å⁻³ (asymmetric; T0+T1+T2)
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

    std::optional<SphericalTensor> shielding(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::HBondShielding))
            return std::nullopt;
        return UnpackSphericalTensor(snap_->column(io::FieldKind::HBondShielding).row(atomIdx));
    }

    std::optional<HBondScalars> scalars(std::size_t atomIdx) const {
        const auto& col = snap_->column(io::FieldKind::HBondScalars);
        if (!col.present)
            return std::nullopt;
        // Writer-definitive shape guard: the 1P9J fixtures (2026-05-24) carry
        // hbond_scalars with 3 columns, NOT the catalog's 4. NPY columns are
        // positional in the writer's order (HBondResult.cpp WriteFeatures), so a
        // short emit drops TRAILING fields — read only what is present and leave
        // the rest at their 0.0 default. This is the difference between a benign
        // missing field and an out-of-bounds read of the next atom's row (the
        // last atom's r[3] would run off the buffer). Which 4th field the writer
        // dropped is a writer question, flagged not guessed.
        HBondScalars s;
        const double* r = col.row(atomIdx);
        if (col.cols > 0) s.nearest_dist = r[0];
        if (col.cols > 1) s.inv_d3 = r[1];
        if (col.cols > 2) s.count_3_5A = r[2];
        if (col.cols > 3) s.mcconnell_scalar = r[3];
        return s;
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
