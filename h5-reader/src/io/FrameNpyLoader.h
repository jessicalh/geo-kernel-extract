// FrameNpyLoader — directory-agnostic producer of one QtConformationSnapshot
// from a directory of per-atom calculator NPYs.
//
// One mechanism serves both a trajectory frame directory and a single-pose
// extraction root.
//
// Per NPY it resolves the filename stem -> typed FieldKind exactly once
// (FindFieldByStem, the no-strings boundary), reads BOTH shape and dtype from
// the NPY header, and widens fields selected by FrameFieldPolicy into the
// snapshot's per-FieldKind NpyColumn store. Fixed catalog widths and Atom or
// Residue topology-axis lengths are enforced before exposing a column. A 1-D array is shaped by
// its NativeAxis
// (Atom/Residue/Ring -> N x 1; Protein -> 1 x K).
//
// "Absent, not faked": a calculator NPY that is missing simply leaves its
// FieldKind column absent (a group view returns nullopt). A malformed NPY logs
// loud via ErrorBus and is skipped (that column absent), but the snapshot is
// still returned — a partial snapshot is success. The loader returns null only
// when the directory itself is unreadable / has no recognised NPYs.

#pragma once

#include <QString>
#include <cstddef>
#include <memory>

namespace h5reader::model {
class QtProtein;
class QtConformationSnapshot;
}  // namespace h5reader::model

namespace h5reader::io {

class FrameNpyLoader {
public:
    // dir: the directory of per-atom NPYs. protein: the topology spine (for
    // axis-length cross-checks). frameIndex / timePs: stamped onto the snapshot
    // (original XTC index + time; 0 for a single pose). Returns null on a
    // directory-level failure; a partial snapshot otherwise.
    static std::shared_ptr<h5reader::model::QtConformationSnapshot>
    LoadSnapshotDir(const QString& dir, const h5reader::model::QtProtein* protein, std::size_t frameIndex, double timePs);
};

}  // namespace h5reader::io
