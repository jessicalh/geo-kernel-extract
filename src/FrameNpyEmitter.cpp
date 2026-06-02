#include "FrameNpyEmitter.h"

#include "CategoryInfoProjection.h"
#include "ConformationResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "TopologySidecar.h"

#include <filesystem>
#include <string>

namespace nmr {

namespace fs = std::filesystem;

// ============================================================================
// Singleton state -- file-local. Configure / OnFrame / Reset are the
// only access surface; nothing else needs to see this struct.
// ============================================================================

namespace {

struct EmitterState {
    bool active = false;
    bool sidecars_written = false;
    const Protein* protein = nullptr;
    FrameNpyEmitter::Config config{};
};

EmitterState& State() {
    static EmitterState s;
    return s;
}

std::string FrameDirName(std::size_t frame_idx) {
    std::string n = std::to_string(frame_idx);
    if (n.size() < 6)
        n.insert(0, 6 - n.size(), '0');
    return "frame_" + n;
}

}  // namespace


// ============================================================================
// FrameNpyEmitter -- public surface
// ============================================================================

void FrameNpyEmitter::Configure(const Protein& protein, Config config) {
    auto& s = State();
    s.protein = &protein;
    s.config = std::move(config);
    s.sidecars_written = false;
    s.active = !s.config.output_dir.empty();

    if (!s.active) {
        OperationLog::Warn("FrameNpyEmitter", "Configure called with empty output_dir; emitter remains inactive");
        return;
    }

    OperationLog::Info(LogCalcOther,
                       "FrameNpyEmitter",
                       "configured: dir=" + s.config.output_dir.string());
}


void FrameNpyEmitter::OnFrame(const ProteinConformation& conf, std::size_t frame_idx, double time_ps) {
    auto& s = State();
    if (!s.active)
        return;
    if (!s.protein)
        return;

    // Per-protein sidecars (atoms_category_info, topology) are written on
    // the first OnFrame. Configure runs before tp.Seed in trajectory
    // mode, so the Protein is not yet finalized at Configure time — no
    // bonds, no rings, no LegacyAmberTopology. By the first OnFrame the
    // pipeline has finalized.
    if (!s.sidecars_written) {
        fs::create_directories(s.config.output_dir);
        const std::string parent = s.config.output_dir.string();
        const int cat = CategoryInfoProjection::WriteFeatures(*s.protein, parent);
        const int topo = TopologySidecar::WriteFeatures(*s.protein, parent);
        OperationLog::Info(
            LogCalcOther,
            "FrameNpyEmitter",
            "sidecars: category_info=" + std::to_string(cat) + " topology=" + std::to_string(topo) + " at " + parent);
        s.sidecars_written = true;
    }

    // No stride/window gate: Trajectory::Run only calls OnFrame for the
    // frames the single --stride dispatched.
    const fs::path frame_dir = s.config.output_dir / FrameDirName(frame_idx);
    fs::create_directories(frame_dir);
    const int arrays = ConformationResult::WriteAllFeatures(conf, frame_dir.string());
    OperationLog::Info(LogCalcOther,
                       "FrameNpyEmitter",
                       "frame " + std::to_string(frame_idx) + " t=" + std::to_string(time_ps) + "ps wrote "
                           + std::to_string(arrays) + " arrays to " + frame_dir.string());
}


void FrameNpyEmitter::Reset() {
    auto& s = State();
    s.active = false;
    s.sidecars_written = false;
    s.protein = nullptr;
    s.config = Config{};
}


bool FrameNpyEmitter::IsActive() {
    return State().active;
}

}  // namespace nmr
