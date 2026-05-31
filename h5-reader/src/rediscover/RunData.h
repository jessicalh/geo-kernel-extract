// RunData / RunLoader / FrameMap / DftFrameSet — the immutable, all-frames-
// resident carrier for one 1P9J calcset, and the loader that builds it.
//
// RunData owns the typed protein spine, the H5-backed trajectory conformation
// (positions + kernel TS + ring-neighbourhood TS, all frames resident), the
// DFT shielding frames keyed by original (trr) frame index, and the frame map
// (H5 row → original index → DFT target). It is built once by RunLoader::Load
// and read-only thereafter.
//
// Reuse, not reinvention (DESIGN.md "Reuse"): the load path is
// QtProteinLoader::LoadRunPath (which resolves the `.LGS`, builds QtProtein +
// TrajectoryConformation from the sidecar + trajectory.h5) plus a walk of
// CalcsetManifest.dft.frames[] through DftShieldingLoader::LoadAndValidate.
// No file discovery — every path comes from the manifest.

#pragma once

#include "../io/CalcsetManifest.h"
#include "../model/Conformation.h"
#include "../model/DftShielding.h"
#include "../model/QtProtein.h"
#include "../model/TrajectoryConformation.h"

#include <QString>

#include <cstddef>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

namespace h5reader::rediscover {

// ── DftFrameSet ───────────────────────────────────────────────────────────
// The DFT shielding frames, keyed by ORIGINAL (trr) frame index — the key the
// `.LGS` dft.frames[] and the H5 frame_indices share. A missing original index
// is an honest gap (the campaign is partial), never a faked value.
class DftFrameSet {
public:
    void Insert(std::size_t originalIndex,
                std::shared_ptr<const model::DftShieldingFrame> frame) {
        byOriginal_[originalIndex] = std::move(frame);
    }

    bool Has(std::size_t originalIndex) const { return byOriginal_.count(originalIndex) != 0; }

    const model::DftShieldingFrame* Frame(std::size_t originalIndex) const {
        auto it = byOriginal_.find(originalIndex);
        return it == byOriginal_.end() ? nullptr : it->second.get();
    }

    // The single atom's shielding for (atom, originalIndex), or nullptr.
    const model::DftAtomShielding* AtomShielding(std::size_t atom, std::size_t originalIndex) const {
        const model::DftShieldingFrame* fr = Frame(originalIndex);
        if (!fr || atom >= fr->atoms.size()) return nullptr;
        return &fr->atoms[atom];
    }

    std::size_t frameCount() const { return byOriginal_.size(); }

private:
    std::unordered_map<std::size_t, std::shared_ptr<const model::DftShieldingFrame>> byOriginal_;
};

// ── FrameMap ──────────────────────────────────────────────────────────────
// H5 row → original (trr) frame index, and the subset of rows that have a DFT
// target. Validates frame_index_basis and frame-count agreement at build time
// (fail-loud per DESIGN.md). The H5 row→original mapping comes from the
// TrajectoryConformation (which reuses TrajectoryFrameMap::OriginalIndex).
class FrameMap {
public:
    // expected_basis: the manifest's trajectory.frame_index_basis (the
    // substrate expects "trr_frame_index"). ok==false + error set on mismatch.
    static std::optional<FrameMap> Build(const model::TrajectoryConformation& traj,
                                         const DftFrameSet& dft,
                                         const QString& frame_index_basis,
                                         QString* err_out);

    std::size_t frameCount() const { return originalByRow_.size(); }
    std::size_t originalIndex(std::size_t row) const { return originalByRow_[row]; }

    // Rows (sorted) that have a DFT target — the rows the case loop walks.
    const std::vector<std::size_t>& dftRows() const { return dftRows_; }

private:
    std::vector<std::size_t> originalByRow_;  // (n_frames,)
    std::vector<std::size_t> dftRows_;        // rows with a DFT target
};

// ── RunData ───────────────────────────────────────────────────────────────
struct RunData {
    std::unique_ptr<model::QtProtein> protein;
    std::unique_ptr<model::Conformation> conformation;  // owns the H5 trajectory
    io::CalcsetManifest manifest;
    DftFrameSet dft;
    FrameMap frameMap;

    // Typed downcast convenience: non-null because RunLoader rejects
    // non-trajectory calcsets.
    const model::TrajectoryConformation* trajectory() const {
        return conformation ? conformation->asTrajectory() : nullptr;
    }
    const io::QtTrajectoryH5* h5() const {
        const auto* t = trajectory();
        return t ? t->h5() : nullptr;
    }
};

// ── RunLoader ─────────────────────────────────────────────────────────────
class RunLoader {
public:
    // Load one trajectory calcset (directory holding the single `.LGS`, or the
    // `.LGS` path). Returns std::nullopt + writes a message to err_out on any
    // failure (not a trajectory kind, H5/sidecar load failure, basis mismatch,
    // frame-count disagreement). A partial DFT campaign is NOT a failure —
    // missing frames are gaps.
    static std::optional<RunData> Load(const QString& calcset_path, QString* err_out);
};

}  // namespace h5reader::rediscover
