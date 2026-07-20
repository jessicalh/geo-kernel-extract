#pragma once

#include "QtFieldCatalog.gen.h"

#include <cstddef>

namespace h5reader::io {

// MOPACDirect contains both compact per-atom scientific values and large
// restart/debug storage. These six compact fields are inputs to the F001
// MOPAC-rich shielding model and therefore belong in a Reader snapshot.
inline constexpr bool IsMopacModelInput(FieldKind kind) noexcept {
    switch (kind) {
    case FieldKind::MOPACChargesFullPrecision:
    case FieldKind::MOPACBondValenciesFullPrecision:
    case FieldKind::MOPACAtomSPopulation:
    case FieldKind::MOPACAtomPPopulation:
    case FieldKind::MOPACAtomDPopulation:
    case FieldKind::MOPACLewisBondCount:
        return true;
    default:
        return false;
    }
}

// A frame snapshot carries calculator results used by the Reader and its
// bundled models. Topology has its own typed sidecar, SpatialIndex is producer
// machinery, and raw MOPAC restart/coefficient storage is not a Reader field.
inline constexpr bool ShouldLoadFrameField(FieldKind kind) noexcept {
    if (kind == FieldKind::Count || kind == FieldKind::AtomsCategoryInfo)
        return false;

    const FieldGroup group = kFieldCatalog[static_cast<std::size_t>(kind)].group;
    switch (group) {
    case FieldGroup::Topology:
    case FieldGroup::SpatialIndex:
        return false;
    case FieldGroup::MOPACDirect:
        return IsMopacModelInput(kind);
    default:
        return true;
    }
}

}  // namespace h5reader::io
