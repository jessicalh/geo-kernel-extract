// QtLarsenHBondGroup — read-side mirror of the SDK's LarsenHBond group
// (LarsenHBondShieldingResult): per-atom hydrogen-bond shielding by Larsen
// 2015 ProCS15 Eq 5, computed via direct DFT grid lookup (LarsenHBondGrid;
// DOI 10.7717/peerj.1344). Runs side-by-side with the kernel-form HBond
// group — the grid-vs-kernel residual is itself a methodological coordinate
// (feedback_methods_accumulate).
//
// Larsen Eq 5 four-term decomposition (all per-ATOM, 9-col SphericalTensor, ppm):
//     Δσ_HB  = Δσ_1°HB  + Δσ_2°HB     (amide-H, i.e. HN, is the donor)
//     Δσ_HαB = Δσ_1°HαB + Δσ_2°HαB    (Hα is the donor)
//   1° = the effect on the DONOR residue i; 2° = the effect on the ACCEPTOR
//   residue i+1 (LarsenHBondShieldingResult.h:85-93 Term enum). `shielding`
//   is the summed total of all applicable terms.
//
// Plus a water term and a diagnostic:
//   • waterTerm: Δσ_w = 2.07 ppm isotropic, applied ONLY to amide H atoms
//     with ZERO geometric H-bond candidates (solvent-exposed amides, from
//     Larsen's NMA-water complex DFT). 0.0 elsewhere; never the full tensor.
//   • diagnosticCBShielding: Larsen Table 2 marks Cβ as unaffected by every
//     term, so this tensor SHOULD be ~0 — it is emitted to verify the
//     parse → grid-load → rotation pipeline. A non-zero Cβ is a METHODOLOGICAL
//     SIGNAL, not a bug (the fixture carries up to ±2.5 ppm on ~14 atoms).
//
// ML note (catalog is_feature): the four decomposition terms AND the water
// term are ML features; the summed total, the Cβ diagnostic, and the count are
// not (the model consumes the decomposition, not the redundant sum).
//
// SENTINELS: unlike QtTripeptideGroup (which writes NaN), Larsen writes 0.0
// for "no contribution" — the per-class tensors are packed UNCONDITIONALLY,
// so an atom outside Larsen's Table-2 dispatch carries a structural ZERO
// tensor, not NaN (fixture: 657 all-zero `shielding` rows == 657 count==0
// rows; no NaN anywhere). count()==0 is the in-band twin. nullopt still means
// the calculator did not run this frame ("absent, not faked").

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtLarsenHBondGroup {
public:
    explicit QtLarsenHBondGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // larsen_hbond_shielding (N,9), ppm — summed total of all applicable
    // Larsen H-bond terms. Structural-zero tensor where the atom has no
    // H-bond contribution (count()==0).
    std::optional<SphericalTensor> shielding(std::size_t atomIdx) const {
        return tensor(io::FieldKind::LarsenHBondShielding, atomIdx);
    }

    // ── The four Eq-5 decomposition terms (the ML features) ───────────
    // 1°HB: amide-H donor effect on the donor residue i (larsen_hbond_1pHB).
    std::optional<SphericalTensor> shieldingPrimaryHB(std::size_t atomIdx) const {
        return tensor(io::FieldKind::LarsenHBond1pHBShielding, atomIdx);
    }
    // 2°HB: amide-H donor effect on the acceptor residue i+1 (larsen_hbond_2pHB).
    std::optional<SphericalTensor> shieldingSecondaryHB(std::size_t atomIdx) const {
        return tensor(io::FieldKind::LarsenHBond2pHBShielding, atomIdx);
    }
    // 1°HαB: Hα donor effect on the donor residue i (larsen_hbond_1pHaB).
    std::optional<SphericalTensor> shieldingPrimaryHaB(std::size_t atomIdx) const {
        return tensor(io::FieldKind::LarsenHBond1pHaBShielding, atomIdx);
    }
    // 2°HαB: Hα donor effect on the acceptor residue i+1 (larsen_hbond_2pHaB).
    std::optional<SphericalTensor> shieldingSecondaryHaB(std::size_t atomIdx) const {
        return tensor(io::FieldKind::LarsenHBond2pHaBShielding, atomIdx);
    }

    // larsen_hbond_diagnostic_CB_shielding (N,9), ppm — the all-but-zero Cβ
    // probe (see header). Non-zero = methodological signal, not a contribution.
    std::optional<SphericalTensor> diagnosticCBShielding(std::size_t atomIdx) const {
        return tensor(io::FieldKind::LarsenHBondDiagnosticCBShielding, atomIdx);
    }

    // larsen_hbond_water_term (N,), ppm — Δσ_w isotropic (0.0, or up to
    // 2.07 ppm on a solvent-exposed amide H with no H-bond candidate).
    std::optional<double> waterTerm(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::LarsenHBondWaterTerm))
            return std::nullopt;
        return snap_->column(io::FieldKind::LarsenHBondWaterTerm).row(atomIdx)[0];
    }

    // larsen_hbond_count (N,), int32→double — the number of Table-2 H-bond pair
    // contributions accumulated on this atom (source field larsen_hbond_n_pairs;
    // the Cβ diagnostic does NOT inflate it — LarsenHBondShieldingResult.cpp:652).
    // 0 = none. NOT a count of all geometric H-bond candidates (those are a
    // separate result-level aggregate, not emitted per atom).
    std::optional<int> count(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::LarsenHBondCount))
            return std::nullopt;
        return static_cast<int>(snap_->column(io::FieldKind::LarsenHBondCount).row(atomIdx)[0]);
    }

    // Convenience: true iff the count column is present and the atom received ≥1
    // Table-2 pair contribution. The clean gate for "is the zero `shielding`
    // tensor real or absent" (count==0 ⟺ all-zero tensor, verified 657==657).
    // NOTE: this does NOT account for the water term — a solvent-exposed amide H
    // can have count==0 (→ false here) yet a non-zero waterTerm() (3 such atoms
    // in the fixture). Test waterTerm() separately for "any Larsen effect at all".
    bool hasContribution(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::LarsenHBondCount))
            return false;
        return snap_->column(io::FieldKind::LarsenHBondCount).row(atomIdx)[0] != 0.0;
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
