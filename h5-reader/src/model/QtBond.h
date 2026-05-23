// QtBond — one covalent bond, mirrored from bonds.npy (9 typed fields).
//
// Substrate gives us order, category, plus four convenience flags
// (is_rotatable, is_aromatic, is_peptide, is_backbone) typed at the
// upstream — overlays that need these don't recompute from the
// category.
//
// Per-frame geometry (length, midpoint, direction) is conformation-
// dependent and lives on QtFrame, not here.

#pragma once

#include "Types.h"

#include <cstdint>

namespace h5reader::model {

struct QtBond {
    int32_t bondIndex = -1;   // row index in bonds.npy
    int32_t atomIndexA = -1;  // into QtProtein.atoms()
    int32_t atomIndexB = -1;
    BondOrder order = BondOrder::Unknown;
    BondCategory category = BondCategory::Unknown;

    // Typed flags from the substrate (sidecar bonds.npy). Mirror
    // src/Bond.h's IsPeptideBond / IsBackbone / IsAromatic predicates
    // but populated from the bond's own typed metadata, not computed.
    bool isRotatable = false;
    bool isAromatic = false;
    bool isPeptide = false;
    bool isBackbone = false;

    // Discriminator queries (still typed, for ergonomic call sites)
    bool IsDisulfide() const { return category == BondCategory::Disulfide; }
    bool IsPeptideCO() const { return category == BondCategory::PeptideCO; }
    bool IsPeptideCN() const { return category == BondCategory::PeptideCN; }
};

}  // namespace h5reader::model
