#include "GromacsFrameHandler.h"
#include "TrajectoryProtein.h"
#include "OperationLog.h"

#include "gromacs/fileio/trrio.h"
#include "gromacs/utility/vectypes.h"

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>

namespace nmr {

// nm → Å conversion (TRR stores positions / velocities in nm; our
// per-frame buffers are in Å / Å/ps to match every other ConformationAtom
// position field in the system).
static constexpr float NM_TO_ANGSTROM = 10.0f;


GromacsFrameHandler::GromacsFrameHandler(TrajectoryProtein& tp)
    : tp_(tp) {}


GromacsFrameHandler::~GromacsFrameHandler() {
    if (trr_fio_) {
        gmx_trr_close(trr_fio_);
        trr_fio_ = nullptr;
    }
}


// ── Open ─────────────────────────────────────────────────────────
//
// Mount the TRR stream via libgromacs gmx_trr_open. The TPR was
// already parsed by TrajectoryProtein::BuildFromTrajectory (which
// called FullSystemReader::ReadTopology), so the protein slice and
// the trimmed protein-only mtop used for PBC fixing are already on
// hand. Atom-count sanity check: TRR header reports total system
// atoms, which must match topology.total_atoms.

bool GromacsFrameHandler::Open(const std::string& trr_path,
                               const std::string& /*tpr_path*/) {
    OperationLog::Scope scope("GromacsFrameHandler::Open", trr_path);

    if (trr_fio_) {
        gmx_trr_close(trr_fio_);
        trr_fio_ = nullptr;
    }

    trr_path_ = trr_path;
    trr_fio_ = gmx_trr_open(trr_path, "r");
    if (!trr_fio_) {
        error_ = "cannot open TRR: " + trr_path;
        return false;
    }

    // Probe the first frame's header to learn natoms and whether
    // velocities are present. We do this with a peek then rewind via
    // Reopen — gmx_trr's stream model doesn't natively peek-and-rewind,
    // so instead we rely on the system topology's total atom count
    // (already parsed from the TPR) and treat the first ReadNextFrame
    // as the source of truth for natoms/v presence.
    const auto& topo = tp_.SysReader().Topology();
    natoms_ = static_cast<int>(topo.total_atoms);

    raw_x_.assign(static_cast<std::size_t>(natoms_) * 3, 0.0f);
    raw_v_.clear();
    trr_has_v_ = false;

    has_read_ = false;
    current_index_ = 0;

    OperationLog::Info(LogCalcOther, "GromacsFrameHandler::Open",
        std::to_string(natoms_) + " atoms/frame, " +
        std::to_string(topo.protein_count) + " protein, " +
        std::to_string(topo.water_count) + " water, " +
        std::to_string(topo.ion_count) + " ions");

    return true;
}


// ── ReadNextFrame ────────────────────────────────────────────────
//
// Read one TRR frame: positions + velocities + box. PBC-fix protein,
// split full-system into protein positions, protein velocities, and
// SolventEnvironment. Populates internal state.
//
// First-frame note: TRR per-frame headers can have x_size / v_size /
// f_size set independently. We allocate the velocity buffer on
// first-read if v_size is non-zero, then keep it allocated for the
// duration. Forces are skipped (NULL pointer; our nstfout=0 anyway).

bool GromacsFrameHandler::ReadNextFrame() {
    if (!trr_fio_) return false;

    int64_t step = 0;
    real    t = 0.0f;
    real    lambda = 0.0f;
    rvec    box_rv[3];
    int     na_frame = natoms_;

    // Pre-fill the velocity buffer with NaN sentinels before every read.
    // libgromacs's gmx_trr_read_frame returns success and leaves the v
    // buffer untouched if the TRR frame carries no velocities (nstvout=0
    // / XTC-like TRR). The unified read API gives us no other channel
    // for "this frame had v_size > 0", so we sentinel-fill, then check
    // post-read whether the sentinel was overwritten. Production MD
    // (nstvout=5000) overwrites; non-velocity TRRs leave NaN, telling
    // us trr_has_v_ is false for that frame.
    const float nan_sentinel = std::numeric_limits<float>::quiet_NaN();
    if (raw_v_.size() != static_cast<std::size_t>(natoms_) * 3) {
        raw_v_.assign(static_cast<std::size_t>(natoms_) * 3, nan_sentinel);
    } else {
        std::fill(raw_v_.begin(), raw_v_.end(), nan_sentinel);
    }

    auto* x_ptr = reinterpret_cast<rvec*>(raw_x_.data());
    auto* v_ptr = reinterpret_cast<rvec*>(raw_v_.data());

    if (!gmx_trr_read_frame(trr_fio_, &step, &t, &lambda,
                            box_rv, &na_frame,
                            x_ptr, v_ptr, /*f=*/nullptr)) {
        return false;
    }

    if (na_frame != natoms_) {
        error_ = "TRR frame natoms mismatch: header=" +
                 std::to_string(na_frame) + " expected=" +
                 std::to_string(natoms_);
        return false;
    }

    // Sentinel-survival check: if the first velocity slot is still NaN
    // after the read, libgromacs didn't write velocities (this TRR has
    // no velocity stream). Otherwise gmx_trr_read_frame overwrote our
    // sentinels with real values.
    trr_has_v_ = !raw_v_.empty() && !std::isnan(raw_v_[0]);

    // Box matrix in Å (TRR stores in nm).
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            box_matrix_(i, j) = static_cast<double>(box_rv[i][j]) * 10.0;
        }
    }

    // PBC fix on protein slice. FullSystemReader owns the protein-only
    // mtop and pbcType; passing the frame's box is correct under NPT.
    const auto& topo = tp_.SysReader().Topology();
    const std::size_t pstart = topo.protein_start;
    const std::size_t pcount = topo.protein_count;

    std::vector<float> protein_coords(
        raw_x_.begin() + pstart * 3,
        raw_x_.begin() + (pstart + pcount) * 3);

    // FullSystemReader::MakeProteinWhole signature takes a 3×3 float box.
    float box_for_pbc[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            box_for_pbc[i][j] = box_rv[i][j];
        }
    }
    if (!tp_.SysReader().MakeProteinWhole(protein_coords, box_for_pbc)) {
        error_ = std::string("MakeProteinWhole failed at index ") +
                 std::to_string(has_read_ ? current_index_ + 1 : 0);
        return false;
    }

    // Write fixed coords back over the protein slice for ExtractFrame.
    std::vector<float> fixed_xyz(raw_x_.begin(), raw_x_.end());
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

    // Carve protein velocities (nm/ps → Å/ps). Empty if TRR had no v.
    protein_velocities_.clear();
    if (trr_has_v_) {
        protein_velocities_.reserve(pcount);
        for (std::size_t a = 0; a < pcount; ++a) {
            const std::size_t base = (pstart + a) * 3;
            protein_velocities_.emplace_back(
                static_cast<double>(raw_v_[base + 0]) * 10.0,
                static_cast<double>(raw_v_[base + 1]) * 10.0,
                static_cast<double>(raw_v_[base + 2]) * 10.0);
        }
    }

    last_frame_time_ = static_cast<double>(t);
    current_index_ = has_read_ ? current_index_ + 1 : 0;
    has_read_ = true;
    return true;
}


// ── Skip ─────────────────────────────────────────────────────────
//
// Advance the TRR cursor without extracting. Used for stride-based
// frame dropping. Calls gmx_trr_read_frame with NULL output pointers
// for x, v, f — libgromacs reads the header and skips the data.

bool GromacsFrameHandler::Skip() {
    if (!trr_fio_) return false;

    int64_t step = 0;
    real    t = 0.0f;
    real    lambda = 0.0f;
    rvec    box_rv[3];
    int     na_frame = natoms_;

    if (!gmx_trr_read_frame(trr_fio_, &step, &t, &lambda,
                            box_rv, &na_frame,
                            /*x=*/nullptr, /*v=*/nullptr, /*f=*/nullptr)) {
        return false;
    }

    current_index_ = has_read_ ? current_index_ + 1 : 0;
    has_read_ = true;  // cursor has moved past a real frame
    return true;
}


// ── Reopen ───────────────────────────────────────────────────────
//
// Reset the TRR cursor to the start. Per the static path's legacy
// expectation, Protein is not re-finalized — the caller has already
// done that on the first pass and should treat conf0 as valid.

bool GromacsFrameHandler::Reopen() {
    if (trr_fio_) {
        gmx_trr_close(trr_fio_);
        trr_fio_ = nullptr;
    }
    trr_fio_ = gmx_trr_open(trr_path_, "r");
    if (!trr_fio_) {
        error_ = "failed to reopen TRR: " + trr_path_;
        return false;
    }
    has_read_ = false;
    current_index_ = 0;
    return true;
}

}  // namespace nmr
