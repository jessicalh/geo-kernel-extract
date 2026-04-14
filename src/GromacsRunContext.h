#pragma once
//
// GromacsRunContext: the complete state of one GROMACS trajectory run.
//
// Preloaded data (read once from TPR + EDR at build time):
//   - Bonded interaction parameters for per-atom energy decomposition
//   - All EDR energy frames for O(1) per-frame lookup
//
// Cursor state (advances as the frame handler iterates):
//   - Current frame index and time
//   - Current frame's energy (looked up from preloaded EDR)
//
// Owned by GromacsProtein. Borrowed by GromacsFrameHandler (which
// advances the cursor) and by RunAnalysis (which reads from it to
// populate RunOptions per frame).
//

#include "BondedEnergyResult.h"
#include "GromacsEnergyResult.h"

#include <cstddef>
#include <string>
#include <vector>

namespace nmr {

class FullSystemReader;  // forward — only needed in Build impl

struct GromacsRunContext {

    // ── Build ────────────────────────────────────────────────────
    // Extracts bonded params via the caller's FullSystemReader (already
    // parsed from the TPR — no re-read). Loads EDR if path non-empty.
    // Returns false on error (check error).
    bool Build(FullSystemReader& reader,
               const std::string& tpr_path,
               const std::string& edr_path = "");

    // ── Preloaded data ──────────────────────────────────────────

    BondedParameters bonded_params;
    std::vector<GromacsEnergy> edr_frames;  // sorted by time

    bool HasBondedParams() const { return !bonded_params.interactions.empty(); }
    bool HasEdr() const { return !edr_frames.empty(); }

    // O(log T) lookup by simulation time. Returns nullptr if no EDR.
    const GromacsEnergy* EnergyAtTime(double time_ps) const;

    // ── Cursor state (advanced by frame handler) ────────────────

    size_t frame_index = 0;
    double frame_time_ps = 0.0;

    // Current frame's energy from preloaded EDR. Updated by AdvanceFrame.
    // nullptr if no EDR loaded.
    const GromacsEnergy* current_energy = nullptr;

    // Advance cursor to a new frame. Updates current_energy from EDR.
    void AdvanceFrame(size_t idx, double time_ps);

    // ── Error ───────────────────────────────────────────────────
    std::string error;

private:
    bool LoadEdr(const std::string& edr_path);
};

}  // namespace nmr
