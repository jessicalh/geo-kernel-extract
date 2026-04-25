#include "GromacsFrameHandler.h"
#include "TrajectoryProtein.h"
#include "OperationLog.h"

#include <algorithm>
#include <exception>

namespace nmr {


GromacsFrameHandler::GromacsFrameHandler(TrajectoryProtein& tp)
    : tp_(tp) {}


// ── Open ─────────────────────────────────────────────────────────
//
// Mount the XTC stream. The TPR was already parsed by
// TrajectoryProtein::BuildFromTrajectory (which called
// FullSystemReader::ReadTopology), so the protein slice and the
// trimmed protein-only mtop used for PBC are already on hand. The
// pre-2026-04-25 cross-check between MoleculeWholer's independent
// reparse and FullSystemReader's slice is gone — both are now the
// same source of truth, with composition predicates rather than
// moltype-name string matching.

bool GromacsFrameHandler::Open(const std::string& xtc_path,
                               const std::string& /*tpr_path*/) {
    OperationLog::Scope scope("GromacsFrameHandler::Open", xtc_path);

    if (!reader_.Open(xtc_path)) {
        error_ = "cannot open XTC: " + xtc_path;
        return false;
    }

    const auto& topo = tp_.SysReader().Topology();

    has_read_ = false;
    current_index_ = 0;

    OperationLog::Info(LogCalcOther, "GromacsFrameHandler::Open",
        std::to_string(reader_.natoms()) + " atoms/frame, " +
        std::to_string(topo.protein_count) + " protein, " +
        std::to_string(topo.water_count) + " water, " +
        std::to_string(topo.ion_count) + " ions");

    return true;
}


// ── ReadNextFrame ────────────────────────────────────────────────
//
// Read one XTC frame, PBC-fix protein coords, split full-system into
// protein positions + SolventEnvironment. Populates internal state.

bool GromacsFrameHandler::ReadNextFrame() {
    XtcFrame frame;
    if (!reader_.ReadNext(frame)) return false;

    const auto& topo = tp_.SysReader().Topology();
    const std::size_t pstart = topo.protein_start;
    const std::size_t pcount = topo.protein_count;

    // PBC fix on protein slice. FullSystemReader owns the protein-only
    // mtop and pbcType; passing the frame's box is correct under NPT.
    std::vector<float> protein_coords(
        frame.x.begin() + pstart * 3,
        frame.x.begin() + (pstart + pcount) * 3);
    if (!tp_.SysReader().MakeProteinWhole(protein_coords, frame.box)) {
        error_ = std::string("MakeProteinWhole failed at index ") +
                 std::to_string(has_read_ ? current_index_ + 1 : 0);
        return false;
    }

    // Write fixed coords back.
    std::vector<float> fixed_xyz = frame.x;
    std::copy(protein_coords.begin(), protein_coords.end(),
              fixed_xyz.begin() + pstart * 3);

    // Split into protein positions (Å) + solvent.
    protein_positions_.clear();
    solvent_ = {};
    if (!tp_.SysReader().ExtractFrame(fixed_xyz, protein_positions_, solvent_)) {
        error_ = std::string("frame extract failed at index ") +
                 std::to_string(has_read_ ? current_index_ + 1 : 0);
        return false;
    }

    last_frame_time_ = static_cast<double>(frame.time);
    current_index_ = has_read_ ? current_index_ + 1 : 0;
    has_read_ = true;
    return true;
}


// ── Skip ─────────────────────────────────────────────────────────
//
// Advance the XTC cursor without extracting. Used for stride-based
// frame dropping: skip stride-1 frames between each ReadNextFrame.

bool GromacsFrameHandler::Skip() {
    if (!reader_.Skip()) return false;
    current_index_ = has_read_ ? current_index_ + 1 : 0;
    has_read_ = true;  // cursor has moved past a real frame
    return true;
}


// ── Reopen ───────────────────────────────────────────────────────
//
// Reset the XTC cursor to the start. Per the static path's legacy
// expectation, Protein is not re-finalized — the caller has already
// done that on the first pass and should treat conf0 as valid.

bool GromacsFrameHandler::Reopen() {
    if (!reader_.Reopen()) {
        error_ = "failed to reopen XTC";
        return false;
    }
    has_read_ = false;
    current_index_ = 0;
    return true;
}

}  // namespace nmr
