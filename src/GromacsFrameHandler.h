#pragma once
//
// GromacsFrameHandler: XTC/TPR stream reader with PBC fix.
//
// Pure reader. Open mounts the XTC stream and builds the PBC fixer
// from the TPR. ReadNextFrame advances one frame — reads raw XTC,
// applies PBC fix to the protein slice, and splits the full-system
// coordinates into protein positions + SolventEnvironment. Accessors
// expose the read state for Trajectory::Run to orchestrate the
// per-frame pipeline.
//
// Does NOT create ProteinConformations (that's tp.Seed for frame 0
// and tp.TickConformation for frames 1..N, called by Trajectory::Run).
// Does NOT run per-frame ConformationResults (that's OperationRunner::
// Run, called by Trajectory::Run). Does NOT write to Trajectory's env
// (that's Trajectory::Run populating env from the handler's
// accessors + the preloaded EDR lookup each frame). Does NOT know
// about TrajectoryResults.
//
// Name stays GROMACS-specific: this is the format-specific reader.
// Trajectory-scope orchestration lives in Trajectory::Run.
//

#include "FullSystemReader.h"
#include "SolventEnvironment.h"
#include "Types.h"
#include "xtc_reader.h"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class TrajectoryProtein;

class GromacsFrameHandler {
public:
    explicit GromacsFrameHandler(TrajectoryProtein& tp);

    // Mount XTC stream + build PBC fixer from TPR. Sanity-checks that
    // the TPR atom count matches the wrapped Protein's topology. Does
    // NOT read a frame.
    bool Open(const std::string& xtc_path, const std::string& tpr_path);

    // Advance one frame. Reads XTC, applies PBC fix, splits into
    // protein positions + SolventEnvironment. Populates internal
    // state readable via accessors. Returns false at EOF.
    bool ReadNextFrame();

    // Advance one frame without extracting. Returns false at EOF.
    // After Skip, accessors (ProteinPositions/Solvent) return stale
    // data from the last successful ReadNextFrame — only Index()
    // updates. Intended for stride-based frame dropping.
    bool Skip();

    // Reopen XTC from the start. The next ReadNextFrame reads frame 0.
    bool Reopen();

    // ── Per-frame accessors (valid after a successful ReadNextFrame) ─

    const std::vector<Vec3>& ProteinPositions() const {
        return protein_positions_;
    }
    const SolventEnvironment& Solvent() const { return solvent_; }

    // Position in the XTC of the most recently advanced frame
    // (counting both ReadNextFrame and Skip). Undefined before any
    // advance.
    std::size_t Index() const { return current_index_; }

    // Simulation time of the most recently ReadNextFrame'd frame.
    double Time() const { return last_frame_time_; }

    bool HasRead() const { return has_read_; }

    const std::string& error() const { return error_; }

private:
    TrajectoryProtein& tp_;

    XtcStreamReader reader_;

    std::vector<Vec3> protein_positions_;
    SolventEnvironment solvent_;

    double last_frame_time_ = 0.0;
    std::size_t current_index_ = 0;   // valid iff has_read_
    bool has_read_ = false;

    std::string error_;
};

}  // namespace nmr
