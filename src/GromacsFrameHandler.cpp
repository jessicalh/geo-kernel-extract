#include "GromacsFrameHandler.h"
#include "ProteinConformation.h"
#include "ConformationResult.h"
#include "OperationLog.h"

#include <filesystem>

namespace fs = std::filesystem;

namespace nmr {

GromacsFrameHandler::GromacsFrameHandler(GromacsProtein& gp)
    : gp_(gp)
{}


bool GromacsFrameHandler::Open(
        const std::string& xtc_path,
        const std::string& tpr_path,
        const RunOptions& caller_opts) {

    OperationLog::Scope scope("GromacsFrameHandler::Open", xtc_path);

    // Open streaming reader
    if (!reader_.Open(xtc_path)) {
        error_ = "cannot open XTC: " + xtc_path;
        return false;
    }

    // PBC fixer from TPR (trimmed to protein topology, verbatim
    // from fes-sampler — do not modify)
    try {
        wholer_ = std::make_unique<MoleculeWholer>(tpr_path);
    } catch (const std::exception& e) {
        error_ = std::string("MoleculeWholer: ") + e.what();
        return false;
    }

    // Verify protein atom count agreement between MoleculeWholer
    // (which trims TPR to protein blocks) and FullSystemReader
    // (which parsed atom ranges from the same TPR).
    const auto& topo = gp_.sys_reader().Topology();
    if (static_cast<size_t>(wholer_->natoms()) != topo.protein_count) {
        error_ = "protein atom count mismatch: MoleculeWholer=" +
                 std::to_string(wholer_->natoms()) +
                 " FullSystemReader=" +
                 std::to_string(topo.protein_count);
        return false;
    }

    // Read and process first frame — finalizes protein
    XtcFrame first_frame;
    if (!reader_.ReadNext(first_frame)) {
        error_ = "cannot read first frame from " + xtc_path;
        return false;
    }

    // PBC fix on protein portion of full-system frame
    const size_t pstart = topo.protein_start;
    const size_t pcount = topo.protein_count;
    std::vector<float> protein_coords(
        first_frame.x.begin() + pstart * 3,
        first_frame.x.begin() + (pstart + pcount) * 3);
    wholer_->make_whole(protein_coords, first_frame.box);

    // Write fixed coords back into full frame
    std::copy(protein_coords.begin(), protein_coords.end(),
              first_frame.x.begin() + pstart * 3);

    // Split into protein positions + solvent
    std::vector<Vec3> protein_pos;
    SolventEnvironment solvent;
    if (!gp_.sys_reader().ExtractFrame(first_frame.x, protein_pos, solvent)) {
        error_ = "failed to extract first frame";
        return false;
    }

    // Finalize protein: bond detection + conformation 0 + accumulators
    gp_.FinalizeProtein(protein_pos, static_cast<double>(first_frame.time));

    // Run calculators on conformation 0 with the SAME opts as all frames.
    // Charges come from GromacsProtein (TPR) regardless of caller opts.
    auto& conf0 = gp_.protein().ConformationAt(0);
    RunOptions init_opts = caller_opts;
    init_opts.charge_source = gp_.charges();
    init_opts.solvent = &solvent;
    init_opts.frame_time_ps = static_cast<double>(first_frame.time);
    RunResult rr = OperationRunner::Run(conf0, init_opts);
    if (!rr.Ok()) {
        error_ = "frame 0 calculator failed: " + rr.error;
        OperationLog::Error("GromacsFrameHandler", error_);
        return false;
    }

    // Accumulate frame 0
    gp_.AccumulateFrame(conf0, 0, static_cast<double>(first_frame.time));

    first_frame_done_ = true;
    frame_idx_ = 0;
    total_frames_ = 1;

    OperationLog::Info(LogCalcOther, "GromacsFrameHandler::Open",
        std::to_string(reader_.natoms()) + " atoms/frame, " +
        std::to_string(pcount) + " protein, " +
        std::to_string(topo.water_count) + " water, " +
        std::to_string(topo.ion_count) + " ions");

    return true;
}


bool GromacsFrameHandler::Next(
        const RunOptions& opts, const std::string& output_dir,
        bool accumulate) {

    XtcFrame frame;
    if (!reader_.ReadNext(frame)) {
        total_frames_ = frame_idx_ + 1;
        return false;
    }

    ++frame_idx_;
    ++total_frames_;

    return ProcessFrame(frame, opts, output_dir, accumulate);
}


bool GromacsFrameHandler::Skip() {
    if (!reader_.Skip()) {
        total_frames_ = frame_idx_ + 1;
        return false;
    }
    ++frame_idx_;
    ++total_frames_;
    return true;
}


bool GromacsFrameHandler::Reopen() {
    if (!reader_.Reopen()) {
        error_ = "failed to reopen XTC";
        return false;
    }
    frame_idx_ = 0;
    total_frames_ = 1;  // frame 0 already counted

    // Skip first frame — already processed during Open
    XtcFrame discard;
    if (!reader_.ReadNext(discard)) {
        error_ = "cannot re-read first frame";
        return false;
    }

    return true;
}


bool GromacsFrameHandler::ProcessFrame(
        const XtcFrame& frame, const RunOptions& opts,
        const std::string& output_dir, bool accumulate) {

    const auto& topo = gp_.sys_reader().Topology();
    const size_t pstart = topo.protein_start;
    const size_t pcount = topo.protein_count;

    // PBC fix on protein coordinates (in nm, in-place on a copy)
    std::vector<float> protein_coords(
        frame.x.begin() + pstart * 3,
        frame.x.begin() + (pstart + pcount) * 3);
    wholer_->make_whole(protein_coords, frame.box);

    // Build a modified frame with fixed protein coords for ExtractFrame
    std::vector<float> fixed_xyz = frame.x;
    std::copy(protein_coords.begin(), protein_coords.end(),
              fixed_xyz.begin() + pstart * 3);

    // Split into protein positions (Angstroms) + solvent
    std::vector<Vec3> protein_pos;
    SolventEnvironment solvent;
    if (!gp_.sys_reader().ExtractFrame(fixed_xyz, protein_pos, solvent)) {
        error_ = "frame " + std::to_string(frame_idx_) + " extract failed";
        return false;
    }

    // Create free-standing conformation (NOT in Protein's vector)
    ProteinConformation conf(&gp_.protein(), std::move(protein_pos));

    // Run calculators
    RunOptions frame_opts = opts;
    frame_opts.solvent = &solvent;
    frame_opts.frame_time_ps = static_cast<double>(frame.time);
    RunResult rr = OperationRunner::Run(conf, frame_opts);

    if (!rr.Ok()) {
        error_ = "frame " + std::to_string(frame_idx_) +
                 " calculator failed: " + rr.error;
        OperationLog::Error("GromacsFrameHandler", error_);
        return false;
    }

    // Accumulate into GromacsProtein (skip in pass 2 — already counted)
    if (accumulate)
        gp_.AccumulateFrame(conf, static_cast<int>(frame_idx_),
                            static_cast<double>(frame.time));

    // Write NPY if requested
    if (!output_dir.empty()) {
        char frame_dir[512];
        std::snprintf(frame_dir, sizeof(frame_dir), "%s/frame_%04zu",
                      output_dir.c_str(), frame_idx_);
        fs::create_directories(frame_dir);
        int n_files = ConformationResult::WriteAllFeatures(conf, frame_dir);
        gp_.AddFramePath(frame_dir);

        OperationLog::Info(LogCalcOther, "GromacsFrameHandler",
            "frame " + std::to_string(frame_idx_) +
            " t=" + std::to_string(frame.time) + "ps" +
            " files=" + std::to_string(n_files));
    }

    // Conformation dies here — memory freed
    return true;
}

}  // namespace nmr
