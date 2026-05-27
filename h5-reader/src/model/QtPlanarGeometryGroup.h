// QtPlanarGeometryGroup — read-side mirror of the SDK's PlanarGeometry group
// (PlanarGeometryResult). Position-only conformation geometry: the actual
// per-frame deviation-from-canonical of the substrate's planar/ring fields.
//
// UNLIKE every other result group, this one spans FOUR axes — its accessors
// take an atom / residue / aromatic-ring / saturated-ring index accordingly,
// NOT a uniform atomIdx. Index conventions (writer PlanarGeometryResult.cpp,
// header :95-127):
//   • pyramidalization   — per ATOM            (atomIdx ∈ [0, AtomCount))
//   • omega*             — per RESIDUE          (residueIdx ∈ [0, ResidueCount))
//   • aromaticChi2       — per AROMATIC ring    (aromatic-ring ordinal; == the
//                          QtProtein global ring index while rings are stored
//                          aromatic-first, ∈ [0, AromaticRingCount))
//   • pucker*            — per SATURATED ring   (saturated-ring ordinal;
//                          == globalRingIndex − AromaticRingCount,
//                          ∈ [0, SaturatedRingCount))
//
// Thin const view; nullopt = PlanarGeometry did not run this frame ("absent,
// not faked"). Per-element NaN in the ω / χ₂ vectors is a REAL "not applicable
// here" sentinel (distinct from column-absence) — callers check std::isnan.

#pragma once

#include "QtConformationSnapshot.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtPlanarGeometryGroup {
public:
    explicit QtPlanarGeometryGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // ── Per-atom ──────────────────────────────────────────────────────
    // Signed sp2 out-of-plane displacement (Å) of a planar-group atom from the
    // plane of its three bonded neighbours; CHARMM/AMBER improper-dihedral sign
    // (right-hand rule on the two lowest-atom-index neighbour vectors,
    // build-stable). 0.0 for non-planar atoms (PlanarGeometryResult.h:32-39).
    std::optional<double> pyramidalization(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::Pyramidalization))
            return std::nullopt;
        return snap_->column(io::FieldKind::Pyramidalization).row(atomIdx)[0];
    }

    // ── Per-residue (index = Protein residue index) ───────────────────
    // ω = Cα(i)–C(i)–N(i+1)–Cα(i+1) peptide-bond dihedral (radians). NaN at a
    // chain-break boundary, the C-terminus, or where the i+1 backbone cache is
    // missing. X→Pro bonds get their ACTUAL ω (cis/trans is real signal) and are
    // flagged by omegaIsXpro — NOT NaN-filled (PlanarGeometryResult.h:21-30; the
    // catalog's "NaN at X-Pro" is stale). Check std::isnan for the sentinel.
    std::optional<double> omegaActual(std::size_t residueIdx) const {
        if (!snap_->has(io::FieldKind::OmegaActual))
            return std::nullopt;
        return snap_->column(io::FieldKind::OmegaActual).row(residueIdx)[0];
    }

    // ω − π wrapped to [−π, π] (radians) — deviation from a planar trans amide.
    // NaN wherever omegaActual is NaN.
    std::optional<double> omegaDeviation(std::size_t residueIdx) const {
        if (!snap_->has(io::FieldKind::OmegaDeviation))
            return std::nullopt;
        return snap_->column(io::FieldKind::OmegaDeviation).row(residueIdx)[0];
    }

    // True where the bond into residue i+1 is X→Pro (1 else 0; 0 at C-term).
    // int8 on disk (PlanarGeometryResult.cpp:478 WriteInt8) → widened to double
    // by the loader; read as a mask. Use it to interpret an omega row as
    // cis/trans Pro isomerism rather than amide non-planarity.
    std::optional<bool> omegaIsXpro(std::size_t residueIdx) const {
        if (!snap_->has(io::FieldKind::OmegaIsXpro))
            return std::nullopt;
        return snap_->column(io::FieldKind::OmegaIsXpro).row(residueIdx)[0] != 0.0;
    }

    // ── Per-aromatic-ring (index = aromatic-ring ordinal) ─────────────
    // Parent residue's χ₂ dihedral (radians): Cα–Cβ–Cγ–Cδ (PHE/TYR),
    // Cα–Cβ–Cγ–Nδ1 (HIS), Cα–Cβ–Cγ–Cδ1 (TRP). The ring-flip observable per Akke
    // & Weininger 2023 (J. Phys. Chem. B 127, 591). HONESTY CAVEAT: a single
    // frame gives the INSTANTANEOUS χ₂, not flip kinetics — hence the name
    // aromatic_chi2, not ring_flip_state. NaN if the parent has no χ₂.
    std::optional<double> aromaticChi2(std::size_t aromaticRingIdx) const {
        if (!snap_->has(io::FieldKind::AromaticChi2))
            return std::nullopt;
        return snap_->column(io::FieldKind::AromaticChi2).row(aromaticRingIdx)[0];
    }

    // ── Per-saturated-ring (index = saturated-ring ordinal) ───────────
    // Cremer-Pople (1975, J. Am. Chem. Soc. 97, 1354) puckering amplitude Q (Å)
    // and phase angle θ (degrees, [0, 360)). 5-rings (Pro pyrrolidine) only;
    // θ mod 72° gives the envelope/twist endo/exo classification.
    std::optional<double> puckerQ(std::size_t saturatedRingIdx) const {
        if (!snap_->has(io::FieldKind::PuckerQ))
            return std::nullopt;
        return snap_->column(io::FieldKind::PuckerQ).row(saturatedRingIdx)[0];
    }
    std::optional<double> puckerTheta(std::size_t saturatedRingIdx) const {
        if (!snap_->has(io::FieldKind::PuckerTheta))
            return std::nullopt;
        return snap_->column(io::FieldKind::PuckerTheta).row(saturatedRingIdx)[0];
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
