#pragma once
//
// Trajectory: process entity representing one traversal of a
// trajectory source (XTC + TPR + EDR). Peer of TrajectoryProtein in
// the trajectory-scope object model but distinct in role:
//
//   TrajectoryProtein = the protein as observed (physical-world).
//   Trajectory        = the run that observed it (process).
//
// Holds both process-role state (cursor, handler) during the active
// phase and record-role state (frame times, indices, selections)
// after Run completes.
//
// Run(tp, config, session, extras, output_dir) drives the traversal
// in eight named phases. The per-frame body is uniform for frame 0
// and frames 1..N; the only asymmetry is that frame 0 goes through
// tp.Seed (canonical conformation, seated permanently on Protein) and
// frames 1..N go through tp.TickConformation (ephemeral conformation
// owned by the iteration).
//
//   1. Open the handler (mount XTC + build PBC fixer from TPR).
//   2. Read frame 0 and tp.Seed — finalize Protein topology, create
//      the canonical ProteinConformation (conf0), allocate
//      TrajectoryAtoms.
//   3. Attach TrajectoryResults from config.TrajectoryResultFactories()
//      + extras. Now that Protein is finalized, factories see bonds,
//      rings, and the initialized TrajectoryAtoms vector.
//   4. Validate dependencies (TrajectoryResult-type vs attached set,
//      ConformationResult-type vs config's required set) and
//      caller-supplied resources (AIMNet2 model if required).
//   5. Build the per-frame RunOptions template (config's per-frame
//      opts + tp's charges, bonded params, and session's AIMNet2).
//   6. Frame 0: update env from handler + EDR; run OperationRunner on
//      the canonical conformation; tp.DispatchCompute; record.
//   7. Per-frame loop: skip stride-1; read next frame; tickify;
//      update env; run OperationRunner; tp.DispatchCompute; record.
//   8. tp.FinalizeAllResults. Selections are NOT collected here —
//      any TrajectoryResult that emits them pushes directly to
//      traj.MutableSelections() during its own Compute or Finalize.
//
// EDR frames are preloaded in the constructor so Phase 6 already has
// them available. Bonded parameters live on TrajectoryProtein
// (topology-scope, set at tp.BuildFromTrajectory). Cursor state lives
// inside the handler during Run only.
//

#include "GromacsEnergyResult.h"   // GromacsEnergy
#include "RecordBag.h"
#include "SelectionRecord.h"
#include "SolventEnvironment.h"
#include "errors.h"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace HighFive { class File; }

namespace nmr {

class TrajectoryProtein;
class RunConfiguration;
class Session;
class TrajectoryResult;
class GromacsFrameHandler;


// Per-frame environment stash owned by Trajectory. Single-slot,
// overwritten each frame. Carries what a frame brings into the
// per-frame pipeline beyond atomic coordinates: solvent (water + ions),
// the current EDR energy pointer, frame index + time. Named so
// TrajectoryResults can read the current frame's environment during
// dispatch without going through RunOptions; handlers write into it
// at each read.
//
// Multi-frame environment history (e.g. water dipoles per frame over
// the whole trajectory) is NOT stored here — that's per-handler
// buffers-from-ctor, deposited from this env's current slot each frame.
struct TrajectoryEnv {
    SolventEnvironment solvent;
    const GromacsEnergy* current_energy = nullptr;
    std::size_t current_frame_idx = 0;
    double current_frame_time = 0.0;
};

class Trajectory {
public:
    Trajectory(std::filesystem::path xtc_path,
               std::filesystem::path tpr_path,
               std::filesystem::path edr_path);
    ~Trajectory();

    Trajectory(const Trajectory&) = delete;
    Trajectory& operator=(const Trajectory&) = delete;

    // Drive the traversal. Updates tp in place. Post-Run, tp holds
    // finalized trajectory-scope data and this Trajectory holds the
    // run record (frame times, selections, env last-frame state).
    //
    // Inputs:
    //   tp        — the trajectory-scope protein model, already built
    //                via tp.BuildFromTrajectory.
    //   config    — the named run shape (PerFrameExtractionSet etc.).
    //                Owns per-frame RunOptions base, handler factories,
    //                required ConformationResult types, stride,
    //                whether AIMNet2 is mandatory.
    //   session   — process-wide resources. Supplies AIMNet2 model
    //                (if the config requires it), acts as the single
    //                place resources live.
    //   extras    — optional ad-hoc TrajectoryResults (tests, custom
    //                analyses). Attached after config.Factories().
    //   output_dir — where trajectory-scope output will be written
    //                at end. Held as a record field; the writer
    //                (traj.WriteH5 / tp.WriteH5) consumes it.
    //
    // Returns kOk on success, or a non-zero ErrorCode with a matching
    // OperationLog::Error diagnostic.
    Status Run(TrajectoryProtein& tp,
               const RunConfiguration& config,
               const Session& session,
               std::vector<std::unique_ptr<TrajectoryResult>> extras = {},
               std::filesystem::path output_dir = {});

    // ── Record fields (valid post-Run) ──────────────────────────

    bool IsComplete() const { return state_ == State::Complete; }
    size_t FrameCount() const { return frame_count_; }
    double TotalTimePs() const;
    const std::vector<double>& FrameTimes() const { return frame_times_; }
    const std::vector<size_t>& FrameIndices() const { return frame_indices_; }

    // ── Selection bag (run-scope event stream) ──────────────────
    //
    // TrajectoryResults push directly during Compute or Finalize.
    // Reducer TrajectoryResults read the bag via ByKind<...>() at
    // Finalize and push their reduced sets back under their own
    // kind. Top-level H5 emission walks selections_.Kinds() and emits
    // one group per kind.

    RecordBag<SelectionRecord>& MutableSelections() { return selections_; }
    const RecordBag<SelectionRecord>& Selections() const {
        return selections_;
    }

    // ── Source paths ────────────────────────────────────────────

    const std::filesystem::path& XtcPath() const { return xtc_path_; }
    const std::filesystem::path& TprPath() const { return tpr_path_; }
    const std::filesystem::path& EdrPath() const { return edr_path_; }

    // ── EDR energy lookup (O(log T) by time) ────────────────────

    bool HasEdr() const { return !edr_frames_.empty(); }
    const GromacsEnergy* EnergyAtTime(double time_ps) const;

    // ── Per-frame environment stash (single-slot) ──────────────

    // The current frame's environment: solvent, energy pointer, frame
    // index + time. Written by the handler each frame before
    // OperationRunner::Run; read by calculators through RunOptions,
    // and by TrajectoryResults directly during tp.DispatchCompute.
    // Overwritten each frame; handlers that want per-frame environment
    // across the trajectory copy from this slot into their own
    // buffers-from-ctor during Compute.
    TrajectoryEnv& MutableEnv() { return env_; }
    const TrajectoryEnv& Env() const { return env_; }

    // Output directory recorded on the run. Used by the writer after
    // Run completes. Empty path means "no writer output requested."
    const std::filesystem::path& OutputDir() const { return output_dir_; }

    // ── H5 emission: /trajectory/ group ─────────────────────────
    //
    // /trajectory/source/        — attributes: xtc_path, tpr_path, edr_path
    // /trajectory/frames/time_ps — (T,) float64
    // /trajectory/frames/original_index — (T,) uint64
    // /trajectory/selections/... — selection records
    void WriteH5(HighFive::File& file) const;

private:
    // EDR preload (called from constructor — happens before any
    // traversal starts so Phase 3 has data).
    bool LoadEdr(const std::filesystem::path& edr_path);

    enum class State { Constructed, Running, Complete };
    State state_ = State::Constructed;

    // Source paths.
    std::filesystem::path xtc_path_;
    std::filesystem::path tpr_path_;
    std::filesystem::path edr_path_;

    // Preloaded EDR energy frames (time-sorted).
    std::vector<GromacsEnergy> edr_frames_;

    // Post-Run record fields.
    size_t frame_count_ = 0;
    std::vector<double> frame_times_;
    std::vector<size_t> frame_indices_;
    RecordBag<SelectionRecord> selections_;
    std::filesystem::path output_dir_;

    // Per-frame env stash (single-slot; overwritten each frame).
    TrajectoryEnv env_;

    // During-Run state.
    std::unique_ptr<GromacsFrameHandler> handler_;
};

}  // namespace nmr
