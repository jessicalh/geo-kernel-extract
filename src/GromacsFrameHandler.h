#pragma once
//
// GromacsFrameHandler: reads GROMACS trajectory frames one at a time,
// runs calculators, feeds results into GromacsProtein accumulators.
//
// Owns: XTC stream, MoleculeWholer (PBC fix).
// Borrows: GromacsProtein (for accumulation + sys_reader + protein).
//
// Frame lifecycle per frame:
//   1. Read full-system XTC frame (streaming, one at a time)
//   2. Extract protein float coords, PBC fix via MoleculeWholer
//   3. Put fixed coords back, split via FullSystemReader
//   4. Create free-standing ProteinConformation
//   5. Run calculators via OperationRunner
//   6. gp.AccumulateFrame(conf, frame_idx)
//   7. Conformation dies — memory freed
//
// First frame is special: finalizes protein (bond detection) and
// creates conformation 0 (permanent, lives in Protein's vector).
//

#include "GromacsProtein.h"
#include "OperationRunner.h"
#include "xtc_reader.h"
#include "pbc_whole.h"
#include <memory>
#include <string>

namespace nmr {

class GromacsFrameHandler {
public:
    // Construct with the GromacsProtein that accumulates results.
    explicit GromacsFrameHandler(GromacsProtein& gp);

    // Open trajectory for streaming. Creates MoleculeWholer from TPR.
    // Reads first frame to finalize protein (bond detection, conformation 0).
    // Frame 0 runs with the provided opts (same calculator set as all
    // subsequent frames). Returns false on error (check error()).
    bool Open(const std::string& xtc_path, const std::string& tpr_path,
              const RunOptions& opts = RunOptions());

    // Process the next frame: read from XTC, PBC fix, split, create
    // conformation, run calculators, accumulate. Returns false at EOF.
    // If output_dir is non-empty, writes NPY to output_dir/frame_NNNN/.
    // Set accumulate=false in pass 2 to avoid double-counting frames
    // that were already accumulated during pass 1.
    bool Next(const RunOptions& opts, const std::string& output_dir = "",
              bool accumulate = true);

    // Skip one frame without processing. Returns false at EOF.
    bool Skip();

    // Reopen XTC from the start for pass 2.
    bool Reopen();

    // Frame index of the frame most recently processed by Next().
    size_t frame_index() const { return frame_idx_; }

    // Total frames seen after a complete pass.
    size_t total_frames() const { return total_frames_; }

    const std::string& error() const { return error_; }

private:
    GromacsProtein& gp_;
    XtcStreamReader reader_;
    std::unique_ptr<MoleculeWholer> wholer_;
    size_t frame_idx_ = 0;
    size_t total_frames_ = 0;
    bool first_frame_done_ = false;
    std::string error_;

    // Process one full-system XTC frame. Common path for first and
    // subsequent frames.
    bool ProcessFrame(const XtcFrame& frame, const RunOptions& opts,
                      const std::string& output_dir, bool accumulate);
};

}  // namespace nmr
