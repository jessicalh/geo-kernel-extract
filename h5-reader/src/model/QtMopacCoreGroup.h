// QtMopacCoreGroup — read-side mirror of the SDK's MOPACCore group
// (MopacResult, PM7+MOZYME semiempirical via libmopac). The QM charge source
// for the FullFat `--mopac` trajectory path, reconciled against ff14SB; absent
// on every other run (→ nullopt, "absent, not faked").
//
// THE FAMILY THAT INTRODUCES THE BOND AND PROTEIN AXES to the reader's group
// views — its three accessors span three different axes:
//   • charge / scalars  — per ATOM       (atomIdx ∈ [0, AtomCount))
//   • bondOrder         — per MOPAC BOND  (bondIdx ∈ [0, bondOrderCount()); a
//                         SPARSE axis that is NOT the topology bond axis — its
//                         ordinal is arbitrary hash order, so join back to
//                         structure via MopacBondOrder.atomA/atomB)
//   • global            — per PROTEIN     (one row for the whole frame; no index)
//
//   mopac_charges     (N,)    per-atom Mulliken charge (e)
//   mopac_scalars     (N, 4)  MopacScalars [charge, sPop, pPop, valency]
//   mopac_bond_orders (B, 3)  MopacBondOrder [atomA, atomB, wibergOrder]
//   mopac_global      (4,)    MopacGlobal [heatOfFormation, dipole]

#pragma once

#include "QtConformationSnapshot.h"
#include "QtResultBlocks.h"  // MopacScalars / MopacBondOrder / MopacGlobal (+ Vec3 via Types.h)

#include "../io/QtFieldCatalog.gen.h"

#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtMopacCoreGroup {
public:
    explicit QtMopacCoreGroup(const QtConformationSnapshot& snapshot) : snap_(&snapshot) {}

    // ── Per-atom ──────────────────────────────────────────────────────
    // Mulliken partial charge (e). Identical value to scalars().charge.
    std::optional<double> charge(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACCharges))
            return std::nullopt;
        return snap_->column(io::FieldKind::MOPACCharges).row(atomIdx)[0];
    }

    // Electronic summary [charge, sPop, pPop, valency].
    std::optional<MopacScalars> scalars(std::size_t atomIdx) const {
        if (!snap_->has(io::FieldKind::MOPACScalars))
            return std::nullopt;
        return MopacScalars::FromRow(snap_->column(io::FieldKind::MOPACScalars).row(atomIdx));
    }

    // ── Per-MOPAC-bond (sparse axis) ──────────────────────────────────
    // Number of MOPAC-reported bond-order rows (the Bond-axis length B). Iterate
    // bondOrder(0 .. bondOrderCount()-1). FIRST reader accessor to read a column's
    // row count — the Bond axis differs from AtomCount, so it relies on the loader
    // (#4) populating NpyColumn.rows from the NPY shape (846 atoms but 896 MOPAC
    // bonds in the fixture). Returns 0 both when MOPAC was absent AND when it ran
    // but reported no bonds (the writer emits an explicit empty (0,3) array —
    // MopacResult.cpp:522-525); test snapshot.has(FieldKind::MOPACBondOrders)
    // directly to distinguish the two.
    std::size_t bondOrderCount() const {
        if (!snap_->has(io::FieldKind::MOPACBondOrders))
            return 0;
        return static_cast<std::size_t>(snap_->column(io::FieldKind::MOPACBondOrders).rows);
    }

    // One MOPAC bond-order row by Bond-axis ordinal ∈ [0, bondOrderCount()).
    std::optional<MopacBondOrder> bondOrder(std::size_t bondIdx) const {
        if (!snap_->has(io::FieldKind::MOPACBondOrders))
            return std::nullopt;
        return MopacBondOrder::FromRow(snap_->column(io::FieldKind::MOPACBondOrders).row(bondIdx));
    }

    // ── Per-protein (frame-global) ────────────────────────────────────
    // PM7 heat of formation + molecular dipole. Protein-axis: no atom index,
    // reads the single frame row.
    std::optional<MopacGlobal> global() const {
        if (!snap_->has(io::FieldKind::MOPACGlobal))
            return std::nullopt;
        return MopacGlobal::FromRow(snap_->column(io::FieldKind::MOPACGlobal).row(0));
    }

private:
    const QtConformationSnapshot* snap_ = nullptr;
};

}  // namespace h5reader::model
