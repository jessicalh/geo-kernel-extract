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
#include "../io/QtFieldCatalog.gen.h"
#include "../model/Conformation.h"
#include "../model/DftShielding.h"
#include "../model/QtProtein.h"
#include "../model/TrajectoryConformation.h"

#include <QString>

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <optional>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace h5reader::rediscover {

// ── DftFrameSet ───────────────────────────────────────────────────────────
// The DFT shielding frames, keyed by ORIGINAL (trr) frame index — the key the
// `.LGS` dft.frames[] and the H5 frame_indices share. A missing original index
// is an honest gap (the campaign is partial), never a faked value.
class DftFrameSet {
public:
    class SelectionReadGuard {
    public:
        explicit SelectionReadGuard(const DftFrameSet& dft) : dft_(dft) {
            dft_.EnterSelectionReadGuard();
        }
        ~SelectionReadGuard() { dft_.LeaveSelectionReadGuard(); }

        SelectionReadGuard(const SelectionReadGuard&) = delete;
        SelectionReadGuard& operator=(const SelectionReadGuard&) = delete;

    private:
        const DftFrameSet& dft_;
    };

    void Insert(std::size_t originalIndex,
                std::shared_ptr<const model::DftShieldingFrame> frame) {
        byOriginal_[originalIndex] = std::move(frame);
    }

    bool Has(std::size_t originalIndex) const {
        AssertReadable();
        return byOriginal_.count(originalIndex) != 0;
    }

    const model::DftShieldingFrame* Frame(std::size_t originalIndex) const {
        AssertReadable();
        auto it = byOriginal_.find(originalIndex);
        return it == byOriginal_.end() ? nullptr : it->second.get();
    }

    // The single atom's shielding for (atom, originalIndex), or nullptr.
    const model::DftAtomShielding* AtomShielding(std::size_t atom, std::size_t originalIndex) const {
        AssertReadable();
        auto it = byOriginal_.find(originalIndex);
        const model::DftShieldingFrame* fr = it == byOriginal_.end() ? nullptr : it->second.get();
        if (!fr || atom >= fr->atoms.size()) return nullptr;
        return &fr->atoms[atom];
    }

    std::size_t frameCount() const {
        AssertReadable();
        return byOriginal_.size();
    }

private:
    void EnterSelectionReadGuard() const { ++selectionReadGuardDepth_; }
    void LeaveSelectionReadGuard() const {
        if (selectionReadGuardDepth_ > 0) --selectionReadGuardDepth_;
    }
    void AssertReadable() const {
        if (selectionReadGuardDepth_ > 0)
            throw std::runtime_error("DFT target read during CaseHunter selection");
    }

    std::unordered_map<std::size_t, std::shared_ptr<const model::DftShieldingFrame>> byOriginal_;
    mutable int selectionReadGuardDepth_ = 0;
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
    static FrameMap Static(std::size_t originalIndex, bool hasDft);

    std::size_t frameCount() const { return originalByRow_.size(); }
    std::size_t originalIndex(std::size_t row) const { return originalByRow_[row]; }

    // Rows (sorted) that have a DFT target — the rows the case loop walks.
    const std::vector<std::size_t>& dftRows() const { return dftRows_; }

private:
    std::vector<std::size_t> originalByRow_;  // (n_frames,)
    std::vector<std::size_t> dftRows_;        // rows with a DFT target
};

// ── RunData ───────────────────────────────────────────────────────────────
enum class PoseKind : int { Trajectory = 0, Static = 1 };

struct StaticNpyArray {
    QString stem;
    QString path;
    std::size_t rows = 0;
    std::size_t cols = 0;
    bool frameVarying = false;
    std::size_t atomsPerFrame = 0;
    std::size_t frameCount = 1;
    QString dtype_descr;
    std::vector<double> values;
    std::vector<float> floatValues;

    bool empty() const { return rows == 0 || cols == 0; }
    bool hasRow(std::size_t row) const { return row < rows; }
    std::size_t rowFor(std::size_t atom, std::size_t frame = 0) const {
        return frameVarying ? frame * atomsPerFrame + atom : atom;
    }
    double value(std::size_t row, std::size_t col = 0) const {
        if (row >= rows || col >= cols || values.empty()) return 0.0;
        return values[row * cols + col];
    }
    const double* rowData(std::size_t row) const {
        if (row >= rows || values.empty()) return nullptr;
        return values.data() + row * cols;
    }
    const float* floatRowData(std::size_t row) const {
        if (row >= rows || floatValues.empty()) return nullptr;
        return floatValues.data() + row * cols;
    }
};

struct RunData {
    RunData() = default;
    RunData(RunData&&) noexcept = default;
    RunData& operator=(RunData&&) noexcept = default;
    RunData(const RunData&) = delete;
    RunData& operator=(const RunData&) = delete;

    std::unique_ptr<model::QtProtein> protein;
    std::unique_ptr<model::Conformation> conformation;  // owns the H5 trajectory
    io::CalcsetManifest manifest;
    DftFrameSet dft;
    FrameMap frameMap;
    std::unordered_map<int, StaticNpyArray> producerArrays;

    // Typed downcast convenience: non-null because RunLoader rejects
    // non-trajectory calcsets.
    const model::TrajectoryConformation* trajectory() const {
        return conformation ? conformation->asTrajectory() : nullptr;
    }
    const io::QtTrajectoryH5* h5() const {
        const auto* t = trajectory();
        return t ? t->h5() : nullptr;
    }
    PoseKind poseKind() const { return trajectory() ? PoseKind::Trajectory : PoseKind::Static; }
    double timePs(std::size_t row) const {
        const auto* t = trajectory();
        if (t) return t->timePicoseconds(row);
        return std::numeric_limits<double>::quiet_NaN();
    }
    const StaticNpyArray* producerArray(io::FieldKind kind) const {
        auto it = producerArrays.find(static_cast<int>(kind));
        return it == producerArrays.end() ? nullptr : &it->second;
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
