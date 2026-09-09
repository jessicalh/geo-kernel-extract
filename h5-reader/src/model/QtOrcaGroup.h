// QtOrcaGroup — read-side mirror of the SDK's Orca group (OrcaShieldingResult):
// the DFT-computed NMR shielding tensors (ORCA r²SCAN/def2-SVP GIAO, ppm) parsed
// from an ORCA `_nmr.out`. This is the QUANTUM REFERENCE shielding — the
// calibration target the geometric kernels are fit against — NOT a geometric
// kernel itself (the rest of the reader's groups are kernels; this is the thing
// they predict). Catalog mechanism = "quantum_reference".
//
//   orca_total         (N, 9)  SphericalTensor, ppm — the total shielding σ
//   orca_diamagnetic   (N, 9)  SphericalTensor, ppm — Ramsey diamagnetic term
//   orca_paramagnetic  (N, 9)  SphericalTensor, ppm — Ramsey paramagnetic term
//
// σ = σ_dia + σ_para (Ramsey decomposition); verified total == dia + para to
// ~1.6e-3 ppm in the fixture (the residual is the `_nmr.out` tensor print
// precision). Each is a full asymmetric 3×3 σ decomposed to T0 (isotropic) +
// T1 (antisymmetric) + T2 (symmetric-traceless) via SphericalTensor::Decompose
// then PackFull9 (writer OrcaShieldingResult.cpp:235-254); read here with the
// shared UnpackSphericalTensor ([T0, T1[3], T2[5]]).
//
// These fields may be present in a single-pose extraction. Trajectory DFT
// tensors are loaded separately by DftShieldingStore, so a trajectory snapshot
// normally leaves these accessors empty.

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtOrcaGroup {
public:
    explicit QtOrcaGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // Total DFT shielding tensor σ (ppm).
    std::optional<SphericalTensor> total(std::size_t atomIdx) const {
        return tensor(io::FieldKind::OrcaTotal, atomIdx);
    }

    // Ramsey diamagnetic contribution (ppm).
    std::optional<SphericalTensor> diamagnetic(std::size_t atomIdx) const {
        return tensor(io::FieldKind::OrcaDiamagnetic, atomIdx);
    }

    // Ramsey paramagnetic contribution (ppm). total == diamagnetic + paramagnetic.
    std::optional<SphericalTensor> paramagnetic(std::size_t atomIdx) const {
        return tensor(io::FieldKind::OrcaParamagnetic, atomIdx);
    }

private:
    std::optional<SphericalTensor> tensor(io::FieldKind k, std::size_t atomIdx) const {
        if (!snap_->has(k))
            return std::nullopt;
        return UnpackSphericalTensor(snap_->column(k).row(atomIdx));
    }

    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
