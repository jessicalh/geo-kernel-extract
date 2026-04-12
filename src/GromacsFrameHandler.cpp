#include "GromacsFrameHandler.h"
#include "ProteinConformation.h"
#include "ConformationResult.h"
#include "OperationLog.h"

#include <filesystem>

namespace fs = std::filesystem;

namespace nmr {

GromacsFrameHandler::GromacsFrameHandler(
        GromacsProtein& gp, const RunOptions& opts)
    : gp_(gp), opts_(opts)
{}


bool GromacsFrameHandler::OpenTrajectory(
        const std::string& xtc_path,
        const std::string& tpr_path) {

    OperationLog::Scope scope("GromacsFrameHandler::OpenTrajectory", xtc_path);

    // Read all XTC frames
    try {
        frames_ = read_all_xtc_frames(xtc_path);
    } catch (const std::exception& e) {
        error_ = e.what();
        return false;
    }

    if (frames_.empty()) {
        error_ = "no frames in " + xtc_path;
        return false;
    }

    // Verify atom count matches protein
    size_t protein_atoms = gp_.protein().AtomCount();
    if (static_cast<size_t>(frames_[0].natoms) != protein_atoms) {
        error_ = "atom count mismatch: protein=" +
                 std::to_string(protein_atoms) +
                 " xtc=" + std::to_string(frames_[0].natoms);
        return false;
    }

    // PBC fix — make protein whole in each frame
    try {
        MoleculeWholer wholer(tpr_path);
        for (auto& frame : frames_)
            wholer.make_whole(frame);
    } catch (const std::exception& e) {
        error_ = std::string("PBC fix failed: ") + e.what();
        return false;
    }

    OperationLog::Info(LogCalcOther, "GromacsFrameHandler::OpenTrajectory",
        std::to_string(frames_.size()) + " frames, " +
        std::to_string(protein_atoms) + " atoms");

    return true;
}


size_t GromacsFrameHandler::ProcessFrames(
        size_t max_frames,
        const std::string& temp_dir) {

    const size_t N = gp_.protein().AtomCount();
    size_t n_processed = 0;
    size_t n_to_process = std::min(max_frames, frames_.size());

    for (size_t fi = 0; fi < n_to_process; ++fi) {
        const XtcFrame& frame = frames_[fi];

        // Convert nm → Angstrom and build position vector
        std::vector<Vec3> positions(N);
        for (size_t ai = 0; ai < N; ++ai) {
            positions[ai] = Vec3(
                static_cast<double>(frame.x[ai * 3 + 0]) * 10.0,
                static_cast<double>(frame.x[ai * 3 + 1]) * 10.0,
                static_cast<double>(frame.x[ai * 3 + 2]) * 10.0);
        }

        // Create free-standing conformation — NOT in the Protein's vector
        auto conf = std::make_unique<MDFrameConformation>(
            &gp_.protein(), std::move(positions),
            0,                                    // walker
            static_cast<double>(frame.time),      // time_ps
            1.0,                                  // weight (default, COLVAR later)
            0.0,                                  // rmsd_nm
            0.0);                                 // rg_nm

        // Set frame time for GromacsEnergyResult matching
        RunOptions frame_opts = opts_;
        frame_opts.frame_time_ps = static_cast<double>(frame.time);

        // Run calculators
        RunResult rr = OperationRunner::Run(*conf, frame_opts);

        // Write features to temp directory
        std::string frame_dir = temp_dir + "/frame_" + std::to_string(fi);
        // Use POSIX mkdir — fs::create_directories resolves to libtorch's
        // broken std::filesystem stub when linked against libtorch.
        std::string mkdir_cmd = "mkdir -p " + frame_dir;
        (void)system(mkdir_cmd.c_str());
        int n_files = ConformationResult::WriteAllFeatures(*conf, frame_dir);

        // Record the path
        gp_.AddFramePath(frame_dir);

        OperationLog::Info(LogCalcOther, "GromacsFrameHandler::ProcessFrames",
            "frame " + std::to_string(fi) +
            " t=" + std::to_string(frame.time) + "ps" +
            " atoms=" + std::to_string(N) +
            " files=" + std::to_string(n_files) +
            " attached=[" + [&]{
                std::string s;
                for (size_t ri = 0; ri < rr.attached.size(); ++ri) {
                    if (ri) s += ", ";
                    s += rr.attached[ri];
                }
                return s;
            }() + "]");

        // Conformation dies here — memory freed
        ++n_processed;
    }

    return n_processed;
}

}  // namespace nmr
