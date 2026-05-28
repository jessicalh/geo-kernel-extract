// QtConformationSnapshot — one sampled frame's full per-atom calculator
// state: the per-frame analog of the SDK's `Protein`.
//
// Two linked layers (see notes/H5_READER_REWRITE_DESIGN + the
// project_h5reader_formal_model_design memory): the dense H5 drives
// animation/time-series; THIS sparse, full-fidelity snapshot is loaded
// on demand when the user parks on a sampled frame, joined to the H5 by
// frame index. Identity and topology stay on QtProtein (shared, loaded
// once) — the snapshot holds only the per-frame *calculator* groups,
// exactly the ~60-80 NPYs in per_frame_npys/frame_NNNNNN/ (the same
// ConformationResult::WriteAllFeatures payload as extraction modes 1-4).
//
// Storage is the no-strings load boundary's downstream side: a dense
// per-FieldKind array of raw columns, filled by FrameNpyLoader (task #4)
// after FindFieldByStem() resolves each filename to its typed FieldKind.
// The typed group VIEWS (QtBiotSavartGroup, ...) interpret those columns
// into SphericalTensor / Vec3 / the typed blocks. The views are where
// "objects answer questions about themselves"; this class is the store.
//
// Held only as current-frame source data by the conformation. Long-lived
// dashboard history belongs in strip buffers after observers sample the frame.

#pragma once

#include "../io/QtFieldCatalog.gen.h"

#include <array>
#include <cstddef>
#include <vector>

namespace h5reader::model {

class QtProtein;

// One loaded NPY array: row-major doubles, rows x cols. cols comes from
// the catalog (a catalog cols of -1 — variable / 1-D — resolves to 1 at
// load). `present == false` means the source calculator did not run for
// this frame (e.g. MOPAC absent, or a conditional source) — the reader
// contract is "absent, not faked", so a group view returns nullopt.
struct NpyColumn {
    bool present = false;
    int rows = 0;
    int cols = 0;
    std::vector<double> data;  // size == rows * cols, row-major

    const double* row(std::size_t r) const {
        return data.data() + r * static_cast<std::size_t>(cols);
    }
};

class QtConformationSnapshot {
public:
    QtConformationSnapshot(const QtProtein* protein, std::size_t frameIndex, double timePs)
        : protein_(protein), frameIndex_(frameIndex), timePs_(timePs) {}

    const QtProtein* protein() const { return protein_; }
    std::size_t frameIndex() const { return frameIndex_; }  // original XTC index (frame_NNNNNN)
    double timePs() const { return timePs_; }

    // Boundary store — indexed by FieldKind ordinal. Filled by the loader;
    // read by the typed group views.
    bool has(io::FieldKind k) const { return column(k).present; }
    const NpyColumn& column(io::FieldKind k) const {
        return columns_[static_cast<std::size_t>(k)];
    }
    NpyColumn& mutableColumn(io::FieldKind k) {
        return columns_[static_cast<std::size_t>(k)];
    }

    // Typed group views are constructed directly on a snapshot, e.g.
    //   QtBiotSavartGroup bs(snapshot);  bs.shielding(atomIdx);
    // (Constructed-on rather than returned-by-method to keep this header
    // free of the group headers; the snapshot is pure storage.)

private:
    const QtProtein* protein_ = nullptr;
    std::size_t frameIndex_ = 0;
    double timePs_ = 0.0;
    std::array<NpyColumn, static_cast<std::size_t>(io::FieldKind::Count)> columns_;
};

}  // namespace h5reader::model
