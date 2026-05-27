// QtAimnet2Group — read-side mirror of the SDK's AIMNet2 group (AIMNet2Result +
// AIMNet2ChargeResponseGradientResult). AIMNet2 is a neural-network potential
// (libtorch, wB97M-trained, GPU) — per AIMNet2Result.h it yields per-atom
// Hirshfeld charges, a 256-d learned electronic-structure embedding, and a
// Coulomb EFG (same dipolar kernel as CoulombResult) decomposed into total /
// aromatic / backbone (no sidechain EFG, mirroring Coulomb). A companion Result
// emits the charge-response gradient d(Σ q²)/d(r) via autograd.
//
//   aimnet2_charges                          (N,)    Hirshfeld charge (e)
//   aimnet2_aim                              (N,256) embedding — float32 on disk
//   aimnet2_efg / _aromatic / _backbone      (N,5)   T2-only EFG (QtEfg)
//   aimnet2_charge_response_gradient         (N,3)   dL/dr (e²/Å), parity-odd
//   aimnet2_charge_response_gradient_scalar  (N,)    |dL/dr| (e²/Å)
//
// AIMNet2 is required in production and fails loud (CUDA-mandatory, no silent
// degradation — AIMNet2Result.h), so a nullopt here means the load path
// genuinely lacked it (e.g. a pre-2026-05-09 extraction for the gradient pair),
// still surfaced honestly per the "absent, not faked" contract.

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtAimnet2Group {
public:
    explicit QtAimnet2Group(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // Per-atom AIMNet2 Hirshfeld charge (e).
    std::optional<double> charge(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::AIMNet2Charges))
            return std::nullopt;
        return snap_->column(io::FieldKind::AIMNet2Charges).row(atomIdx)[0];
    }

    // 256-d learned electronic-structure embedding (non-owning view; see
    // AIMNet2Embedding — float32 on disk, widened to double by the loader).
    std::optional<AIMNet2Embedding> embedding(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::AIMNet2Aim))
            return std::nullopt;
        return AIMNet2Embedding{snap_->column(io::FieldKind::AIMNet2Aim).row(atomIdx)};
    }

    // Coulomb EFG from AIMNet2 charges (V/Å², T2-only): total.
    std::optional<QtEfg> efg(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::AIMNet2EFG))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::AIMNet2EFG).row(atomIdx));
    }
    // Aromatic-source EFG component.
    std::optional<QtEfg> efgAromatic(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::AIMNet2EFGAromatic))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::AIMNet2EFGAromatic).row(atomIdx));
    }
    // Backbone-source EFG component.
    std::optional<QtEfg> efgBackbone(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::AIMNet2EFGBackbone))
            return std::nullopt;
        return UnpackEfg(snap_->column(io::FieldKind::AIMNet2EFGBackbone).row(atomIdx));
    }

    // Charge-response gradient dL/dr with L = Σ_j q_j² (e²/Å). Parity-odd vector
    // (catalog irreps "1o"). NOT a Buckingham polarisability α = ∂μ/∂E
    // (AIMNet2ChargeResponseGradientResult.h).
    std::optional<Vec3> chargeResponseGradient(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::AIMNet2ChargeResponseGradient))
            return std::nullopt;
        const double* r =
            snap_->column(io::FieldKind::AIMNet2ChargeResponseGradient).row(atomIdx);
        return Vec3(r[0], r[1], r[2]);
    }
    // L2 norm of the charge-response gradient (e²/Å). The library writes this
    // companion scalar directly; it equals |chargeResponseGradient()| up to fp.
    std::optional<double> chargeResponseGradientNorm(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::AIMNet2ChargeResponseGradientScalar))
            return std::nullopt;
        return snap_->column(io::FieldKind::AIMNet2ChargeResponseGradientScalar).row(atomIdx)[0];
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
