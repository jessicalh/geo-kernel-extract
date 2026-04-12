#pragma once
//
// GromacsProtein: trajectory manager. NOT a feature extractor.
//
// This is trajectory lifecycle management + accumulation, not physics.
// Feature extraction happens in ConformationResult subclasses
// (WaterFieldResult, HydrationShellResult, etc.) which take a
// ProteinConformation + SolventEnvironment and know nothing about
// this class.
//
// Owns the immutable Protein (topology + conformation 0 via AddMDFrame).
// Free-standing conformations are created from the Protein's topology
// but never added to its conformations_ vector — they live for one
// frame and die.
//
// Holds trajectory-specific state:
//   - frame manifest (where frame outputs were written)
//   - per-atom Welford accumulators (GromacsProteinAtom) updated by
//     GromacsFrameHandler after each frame
//   - charges from TPR
//
// The accumulators survive across all frames and are written to
// atom_catalog.csv by GromacsFinalResult at the end.
//
// Future: rename to TrajectoryProtein / TrajectoryProteinAtom when
// the GROMACS-specific loader is factored out.
//

#include "Protein.h"
#include "GromacsProteinAtom.h"
#include "GromacsEnsembleLoader.h"
#include "FullSystemReader.h"
#include "OperationRunner.h"

// XtcFrame is defined in xtc_reader.h which includes xdrfile.h (C).
// xdrfile.h conflicts with GROMACS/torch C++ headers, so we can't
// include it here. The frame storage is opaque — managed in .cpp.
struct XtcFrame;
#include <memory>
#include <set>
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

    // Build from full-system trajectory (TPR + reference PDB).
    // Does NOT read the XTC — frames come later via ScanFrame/ExtractFrame.
    bool BuildFromTrajectory(const std::string& tpr_path,
                             const std::string& ref_pdb_path);

    // ── Identity ─────────────────────────────────────────────────

    Protein& protein() { return *protein_; }
    const Protein& protein() const { return *protein_; }
    ChargeSource* charges() { return charges_.get(); }
    const std::string& protein_id() const { return protein_id_; }
    int net_charge() const { return net_charge_; }
    const std::string& error() const { return error_; }
    const std::vector<std::string>& pose_names() const { return pose_names_; }

    // ── Per-atom trajectory accumulators ─────────────────────────

    GromacsProteinAtom& AtomAt(size_t i) { return atoms_[i]; }
    const GromacsProteinAtom& AtomAt(size_t i) const { return atoms_[i]; }
    size_t AtomCount() const { return atoms_.size(); }

    // ── Trajectory processing ────────────────────────────────────

    // Load all frames from a full-system XTC. Call after BuildFromTrajectory.
    bool LoadTrajectory(const std::string& xtc_path);
    size_t FrameCount() const;

    // Scan one frame: run lightweight calculators, accumulate stats.
    // No NPY output. Returns false on error.
    bool ScanFrame(size_t frame_idx, const RunOptions& scan_opts);

    // Extract one frame: run full calculators, write NPY to output_dir.
    // Returns number of NPY files written, or -1 on error.
    int ExtractFrame(size_t frame_idx, const RunOptions& extract_opts,
                     const std::string& output_dir);

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
    void AccumulateFrame(const ProteinConformation& conf, int frame_idx);

    std::unique_ptr<Protein> protein_;
    std::unique_ptr<ChargeSource> charges_;
    FullSystemReader sys_reader_;
    struct FrameStorage;                  // pimpl for XtcFrame vector
    std::unique_ptr<FrameStorage> frames_;
    std::vector<GromacsProteinAtom> atoms_;
    std::vector<std::string> frame_paths_;
    std::vector<std::string> pose_names_;
    std::string protein_id_;
    std::string error_;
    int net_charge_ = 0;
};

}  // namespace nmr
