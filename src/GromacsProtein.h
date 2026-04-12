#pragma once
//
// GromacsProtein: trajectory manager. NOT a feature extractor.
//
// This is an adapter around a real Protein + its trajectory-level state.
// Feature extraction happens in ConformationResult subclasses
// (WaterFieldResult, HydrationShellResult, etc.) which take a
// ProteinConformation + SolventEnvironment and know nothing about
// this class.
//
// Owns the immutable Protein (topology + conformation 0 via AddMDFrame).
// Conformation 0 is permanent — it feeds the .h5 master file.
// All other frames are free-standing conformations created by
// GromacsFrameHandler, processed, accumulated, then freed.
//
// Holds trajectory-specific state:
//   - per-atom Welford accumulators (GromacsProteinAtom) updated by
//     GromacsFrameHandler after each frame via AccumulateFrame()
//   - charges from TPR
//   - FullSystemReader (topology + frame splitting)
//
// Does NOT read XTC. Does NOT create conformations. Does NOT run
// calculators. Those are GromacsFrameHandler's job.
//

#include "Protein.h"
#include "GromacsProteinAtom.h"
#include "GromacsEnsembleLoader.h"
#include "FullSystemReader.h"

#include <memory>
#include <string>
#include <vector>

namespace nmr {

class GromacsProtein {
public:
    GromacsProtein();
    ~GromacsProtein();

    // ── Build ────────────────────────────────────────────────────

    // Build from fleet paths (pre-extracted PDB poses).
    bool Build(const FleetPaths& paths);

    // Build from full-system trajectory (TPR only).
    // Reads topology for atom ranges (protein/water/ion), builds
    // Protein from TPR. Does NOT finalize — FinalizeProtein() must
    // be called after the first XTC frame provides positions.
    bool BuildFromTrajectory(const std::string& tpr_path);

    // Finalize protein construction using first frame positions.
    // Calls FinalizeConstruction (bond detection) + AddMDFrame
    // (conformation 0) + InitAtoms (accumulators).
    // Called by GromacsFrameHandler after reading the first frame.
    void FinalizeProtein(std::vector<Vec3> first_frame_positions,
                         double time_ps);

    // ── Identity ─────────────────────────────────────────────────

    Protein& protein() { return *protein_; }
    const Protein& protein() const { return *protein_; }
    ChargeSource* charges() { return charges_.get(); }
    const std::string& protein_id() const { return protein_id_; }
    int net_charge() const { return net_charge_; }
    const std::string& error() const { return error_; }
    const std::vector<std::string>& pose_names() const { return pose_names_; }

    // ── Trajectory infrastructure ────────────────────────────────

    // FullSystemReader: topology + frame splitting. Borrowed by
    // GromacsFrameHandler for ExtractFrame on each frame.
    const FullSystemReader& sys_reader() const { return sys_reader_; }

    // ── Per-atom trajectory accumulators ─────────────────────────

    GromacsProteinAtom& AtomAt(size_t i) { return atoms_[i]; }
    const GromacsProteinAtom& AtomAt(size_t i) const { return atoms_[i]; }
    size_t AtomCount() const { return atoms_.size(); }

    // Accumulate per-atom stats from one frame's computed conformation.
    // Called by GromacsFrameHandler after running calculators.
    void AccumulateFrame(const ProteinConformation& conf, int frame_idx);

    // Write the per-atom trajectory catalog CSV.
    void WriteCatalog(const std::string& path) const;

    // Select frames for full extraction from accumulated stats.
    // Returns sorted, deduplicated frame indices.
    std::vector<size_t> SelectFrames(size_t max_frames) const;

    // ── Frame manifest (fleet path) ──────────────────────────────

    const std::vector<std::string>& frame_paths() const { return frame_paths_; }
    void AddFramePath(const std::string& path) { frame_paths_.push_back(path); }

private:
    void InitAtoms();

    std::unique_ptr<Protein> protein_;
    std::unique_ptr<ChargeSource> charges_;
    FullSystemReader sys_reader_;
    std::vector<GromacsProteinAtom> atoms_;
    std::vector<std::string> frame_paths_;
    std::vector<std::string> pose_names_;
    std::string protein_id_;
    std::string error_;
    int net_charge_ = 0;
};

}  // namespace nmr
