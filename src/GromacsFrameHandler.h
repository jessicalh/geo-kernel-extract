#pragma once
//
// GromacsFrameHandler: reads GROMACS trajectory, creates conformations,
// runs calculators, writes results.
//
// Owns the XTC reader and PBC fixer. Creates free-standing
// ProteinConformations that point at the GromacsProtein's Protein
// for topology but are never added to its conformations_ vector.
//
// Each frame: read XTC → PBC fix → nm→A → create conformation →
// run OperationRunner → write NPY → log → conformation dies.
//

#include "GromacsProtein.h"
#include "OperationRunner.h"
#include "xtc_reader.h"
#include "pbc_whole.h"
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class GromacsFrameHandler {
public:
    // Construct with the GromacsProtein (topology source) and run options.
    GromacsFrameHandler(GromacsProtein& gp, const RunOptions& opts);

    // Open a trajectory. Loads all XTC frames from a walker directory,
    // applies PBC fix. Returns false on error (check error()).
    bool OpenTrajectory(const std::string& xtc_path,
                        const std::string& tpr_path);

    // Process up to max_frames frames. Creates free-standing conformations,
    // runs geometry calculators, writes NPY to temp_dir/frame_{i}/.
    // Returns number of frames processed.
    // Frame 0 is skipped here — it's already in the Protein via AddMDFrame.
    size_t ProcessFrames(size_t max_frames,
                         const std::string& temp_dir);

    // Number of frames available in the trajectory.
    size_t FrameCount() const { return frames_.size(); }

    const std::string& error() const { return error_; }

private:
    GromacsProtein& gp_;
    RunOptions opts_;
    std::vector<XtcFrame> frames_;
    std::string error_;
};

}  // namespace nmr
