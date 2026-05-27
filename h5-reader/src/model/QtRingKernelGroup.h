// QtRingKernelGroup — shared base for the through-space ring-current /
// ring-kernel groups (BiotSavart, HaighMallion, PiQuadrupole, Dispersion),
// mirroring the SDK's RingKernelGroup. Each carries the same three per-atom
// fields — a full shielding SphericalTensor (T0+T1+T2) plus per-aromatic-
// ring-type T0 and T2 decompositions — differing only in which catalog
// FieldKind backs them. Subclasses pass their FieldKinds in the ctor;
// QtBiotSavartGroup extends this with total_B + ring_counts.
//
// Thin, const, copyable view over one QtConformationSnapshot. Every
// accessor is std::optional<> — nullopt means the source calculator did
// not run this frame ("absent, not faked"); a present value of 0.0 is a
// real measurement.

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"
#include "Types.h"

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtRingKernelGroup {
public:
    QtRingKernelGroup(const QtConformationSnapshot& snapshot,
                      io::FieldKind shieldingKind,
                      io::FieldKind perTypeT0Kind,
                      io::FieldKind perTypeT2Kind)
        : snap_(&snapshot),
          shieldingKind_(shieldingKind),
          perTypeT0Kind_(perTypeT0Kind),
          perTypeT2Kind_(perTypeT2Kind) {}

    // Full ring-current shielding tensor at this atom (T0 + T1 + T2).
    std::optional<SphericalTensor> shielding(std::size_t atomIdx) const {
        return tensorAt(shieldingKind_, atomIdx);
    }

    // Isotropic (T0) contribution decomposed by aromatic ring type (8).
    std::optional<PerRingTypeT0> perTypeT0(std::size_t atomIdx) const {
        if (!snap_->has(perTypeT0Kind_))
            return std::nullopt;
        const double* r = snap_->column(perTypeT0Kind_).row(atomIdx);
        PerRingTypeT0 out;
        for (std::size_t t = 0; t < kAromaticRingTypeCount; ++t)
            out.byType[t] = r[t];
        return out;
    }

    // T2 (5-component) contribution decomposed by aromatic ring type (8x5).
    std::optional<PerRingTypeT2> perTypeT2(std::size_t atomIdx) const {
        if (!snap_->has(perTypeT2Kind_))
            return std::nullopt;
        const double* r = snap_->column(perTypeT2Kind_).row(atomIdx);
        PerRingTypeT2 out;
        for (std::size_t t = 0; t < kAromaticRingTypeCount; ++t)
            for (std::size_t i = 0; i < 5; ++i)
                out.byType[t][i] = r[t * 5 + i];
        return out;
    }

protected:
    // Decode a 9-column row into a SphericalTensor (layout [T0, T1[3], T2[5]]).
    std::optional<SphericalTensor> tensorAt(io::FieldKind kind, std::size_t atomIdx) const {
        if (!snap_->has(kind))
            return std::nullopt;
        const double* r = snap_->column(kind).row(atomIdx);
        SphericalTensor st;
        st.T0 = r[0];
        for (std::size_t i = 0; i < 3; ++i)
            st.T1[i] = r[1 + i];
        for (std::size_t i = 0; i < 5; ++i)
            st.T2[i] = r[4 + i];
        return st;
    }

    const QtConformationSnapshot* snap_ = nullptr;
    io::FieldKind shieldingKind_;
    io::FieldKind perTypeT0Kind_;
    io::FieldKind perTypeT2Kind_;
};

}  // namespace h5reader::model
