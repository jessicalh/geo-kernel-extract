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
#include "GromacsRunContext.h"

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

    // Build from full-system trajectory directory.
    // dir_path must contain md.tpr, md.xtc, md.edr.
    // Reads TPR once for topology (atom ranges, charges, bonded parameters).
    // Reads EDR (all frames, preloaded for O(1) per-frame lookup).
    // Does NOT finalize — FinalizeProtein() must be called after the
    // first XTC frame provides positions for bond detection.
    bool BuildFromTrajectory(const std::string& dir_path);

    // The run context: bonded params, EDR, cursor state.
    // Lifetime = GromacsProtein lifetime.
    // Mutable accessor for GromacsFrameHandler to advance the cursor.
    const GromacsRunContext& run_context() const { return run_ctx_; }
    GromacsRunContext& run_context_mut() { return run_ctx_; }

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

    // Convenience delegates to run_context().
    const BondedParameters& bonded_params() const { return run_ctx_.bonded_params; }

    // ── Per-atom trajectory accumulators ─────────────────────────

    GromacsProteinAtom& AtomAt(size_t i) { return atoms_[i]; }
    const GromacsProteinAtom& AtomAt(size_t i) const { return atoms_[i]; }
    size_t AtomCount() const { return atoms_.size(); }

    // ── Per-bond trajectory accumulators ─────────────────────────

    GromacsProteinBond& BondAt(size_t i) { return bonds_accum_[i]; }
    const GromacsProteinBond& BondAt(size_t i) const { return bonds_accum_[i]; }
    size_t BondAccumCount() const { return bonds_accum_.size(); }

    // Accumulate per-atom stats from one frame's computed conformation.
    // Called by GromacsFrameHandler after running calculators.
    void AccumulateFrame(const ProteinConformation& conf, int frame_idx,
                         double time_ps = 0.0);

    // Write the per-atom trajectory catalog CSV.
    void WriteCatalog(const std::string& path) const;

    // Write master HDF5 file: positions, rollup stats, bonds.
    void WriteH5(const std::string& path) const;

    // Select frames for full extraction from accumulated stats.
    // Returns sorted, deduplicated frame indices.
    std::vector<size_t> SelectFrames(size_t max_frames) const;

    // ── Stored frame positions (for movies + training) ────────────

    size_t StoredFrameCount() const { return n_stored_frames_; }
    const std::vector<double>& FramePositions() const { return frame_positions_; }
    const std::vector<double>& FrameTimes() const { return frame_times_; }

    // ── Frame manifest (fleet path) ──────────────────────────────

    const std::vector<std::string>& frame_paths() const { return frame_paths_; }
    void AddFramePath(const std::string& path) { frame_paths_.push_back(path); }

private:
    void InitAtoms();

    std::unique_ptr<Protein> protein_;
    std::unique_ptr<ChargeSource> charges_;
    int net_charge_ = 0;
    GromacsRunContext run_ctx_;
    FullSystemReader sys_reader_;
    std::vector<GromacsProteinAtom> atoms_;
    std::vector<GromacsProteinBond> bonds_accum_;

    // Per-frame positions: (frame_idx, atom_idx) → Vec3.
    // Stored flat as [frame0_atom0_x, frame0_atom0_y, frame0_atom0_z,
    //                 frame0_atom1_x, ...].  N_atoms * 3 doubles per frame.
    // These enable VTK movie playback and waveform model training.
    std::vector<double> frame_positions_;   // flat (T * N * 3)
    std::vector<double> frame_times_;       // (T,) in ps
    size_t n_stored_frames_ = 0;

    std::vector<std::string> frame_paths_;
    std::vector<std::string> pose_names_;
    std::string protein_id_;
    std::string error_;
};

}  // namespace nmr
