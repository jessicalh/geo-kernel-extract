#pragma once
//
// GromacsFrameHandler: TRR/TPR stream reader with PBC fix.
//
// Pure reader. Open mounts the TRR stream and builds the PBC fixer
// from the TPR. ReadNextFrame advances one frame — reads raw TRR
// (positions + velocities + box), applies PBC fix to the protein
// slice, and splits the full-system coordinates into protein positions,
// protein velocities, and SolventEnvironment. Accessors expose the
// read state for Trajectory::Run to orchestrate the per-frame pipeline.
//
// TRR (vs the prior XTC reader) carries velocities and box matrix
// directly. Forces are not currently captured (our production .mdp
// uses nstfout=0). The reader is built on libgromacs's gmx_trr_*
// API — same library we already link for TPR parsing — rather than
// vendoring xdrfile_trr.
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

#include <Eigen/Dense>
#include <cstddef>
#include <memory>
#include <string>
#include <vector>

// Forward declare libgromacs file handle so the gromacs headers stay
// in GromacsFrameHandler.cpp (where they already live for trrio/pbc).
struct t_fileio;

namespace nmr {

class TrajectoryProtein;

class GromacsFrameHandler {
public:
    explicit GromacsFrameHandler(TrajectoryProtein& tp);
    ~GromacsFrameHandler();

    GromacsFrameHandler(const GromacsFrameHandler&) = delete;
    GromacsFrameHandler& operator=(const GromacsFrameHandler&) = delete;

    // Mount TRR stream + use the protein-only mtop owned by the
    // FullSystemReader (already parsed at TrajectoryProtein build
    // time) for PBC fixing. Sanity-checks that the TRR atom count
    // matches the wrapped Protein + solvent + ions.  trr_path here
    // is the production.trr path; tpr_path is unused (kept for
    // backward call-site compatibility) since the TPR was already
    // parsed.
    bool Open(const std::string& trr_path, const std::string& tpr_path);

    // Advance one frame. Reads positions + velocities + box from TRR,
    // converts to Å / Å/ps, PBC-fixes protein, splits into protein
    // positions + protein velocities + SolventEnvironment. Populates
    // internal state readable via accessors. Returns false at EOF.
    bool ReadNextFrame();

    // Advance one frame without extracting. Returns false at EOF.
    // After Skip, the per-frame accessors return stale data from the
    // last successful ReadNextFrame — only Index() updates. Intended
    // for stride-based frame dropping.
    bool Skip();

    // Reopen the TRR from the start. The next ReadNextFrame reads
    // frame 0.
    bool Reopen();

    // ── Per-frame accessors (valid after a successful ReadNextFrame) ─

    const std::vector<Vec3>& ProteinPositions() const {
        return protein_positions_;
    }
    const std::vector<Vec3>& ProteinVelocities() const {
        return protein_velocities_;
    }
    const Eigen::Matrix3d& BoxMatrix() const { return box_matrix_; }
    const SolventEnvironment& Solvent() const { return solvent_; }

    // Position in the TRR of the most recently advanced frame
    // (counting both ReadNextFrame and Skip). Undefined before any
    // advance.
    std::size_t Index() const { return current_index_; }

    // Simulation time of the most recently ReadNextFrame'd frame.
    double Time() const { return last_frame_time_; }

    bool HasRead() const { return has_read_; }

    const std::string& error() const { return error_; }

private:
    TrajectoryProtein& tp_;

    // Libgromacs TRR handle. Opaque pointer; lifetime matches Open/Close.
    t_fileio* trr_fio_ = nullptr;
    std::string trr_path_;
    int natoms_ = 0;

    // Per-frame buffers (full-system, in nm/ps as TRR stores them).
    std::vector<float> raw_x_;       // natoms × 3
    std::vector<float> raw_v_;       // natoms × 3 (empty if TRR has no velocities)
    bool                trr_has_v_ = false;

    // Per-frame extracted state (in Å / Å/ps).
    std::vector<Vec3>  protein_positions_;
    std::vector<Vec3>  protein_velocities_;  // empty when TRR carries no v
    Eigen::Matrix3d    box_matrix_ = Eigen::Matrix3d::Zero();
    SolventEnvironment solvent_;

    double last_frame_time_ = 0.0;
    std::size_t current_index_ = 0;   // valid iff has_read_
    bool has_read_ = false;

    std::string error_;
};

}  // namespace nmr
