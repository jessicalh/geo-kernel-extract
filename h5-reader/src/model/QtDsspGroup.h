// QtDsspGroup — read-side mirror of the SDK's DSSP group (DsspResult; libdssp,
// Joosten 2011 / Kabsch & Sander 1983). Secondary structure + backbone
// geometry + sidechain dihedrals, all per-atom BROADCAST from the residue.
//
//   dssp_backbone     (N,5)  phi, psi (NEGATED-IUPAC rad), DSSP-residue sasa, ssHelix, ssSheet
//   dssp_ss8          (N,8)  8-class one-hot → DsspCode
//   dssp_hbond_energy (N,4)  acceptor0/1, donor0/1 (kcal/mol)
//   dssp_chi          (N,12) χ1-4 × {cos, sin, exists}, interleaved
//
// Thin const view over one QtConformationSnapshot; every accessor is
// std::optional<> — nullopt means DSSP did not run this frame ("absent, not
// faked"). phi/psi are libdssp's negated-IUPAC convention; negate for IUPAC.

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtDsspGroup {
public:
    explicit QtDsspGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    std::optional<DsspScalars> backbone(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::DSSPBackbone))
            return std::nullopt;
        return DsspScalars::FromRow(snap_->column(io::FieldKind::DSSPBackbone).row(atomIdx));
    }

    std::optional<DsspSs8> ss8(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::DSSPSs8))
            return std::nullopt;
        return DsspSs8::FromRow(snap_->column(io::FieldKind::DSSPSs8).row(atomIdx));
    }

    std::optional<DsspHBondEnergy> hbondEnergy(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::DSSPHBondEnergy))
            return std::nullopt;
        return DsspHBondEnergy::FromRow(snap_->column(io::FieldKind::DSSPHBondEnergy).row(atomIdx));
    }

    std::optional<DsspChi> chi(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::DSSPChi))
            return std::nullopt;
        return DsspChi::FromRow(snap_->column(io::FieldKind::DSSPChi).row(atomIdx));
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
