#include "GromacsEnergyResult.h"
#include "NpyWriter.h"

#include <cassert>

namespace nmr {

// ── Compute (preloaded) ─────────────────────────────────────────
// Data already looked up by Trajectory::EnergyAtTime.

std::unique_ptr<GromacsEnergyResult> GromacsEnergyResult::Compute(
        ProteinConformation& conf,
        const GromacsEnergy& energy) {

    auto result = std::make_unique<GromacsEnergyResult>();
    result->energy_ = energy;
    return result;
}


// ── WriteFeatures ──────────────────────────────────────────────

int GromacsEnergyResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {

    // Write as a (1, COLS) array: one row of per-frame scalars.
    // Column order matches the GromacsEnergy struct layout.
    std::vector<double> row = {
        // Electrostatic (3)
        energy_.coulomb_sr, energy_.coulomb_recip, energy_.coulomb_14,
        // Bonded (6)
        energy_.bond, energy_.angle, energy_.urey_bradley,
        energy_.proper_dih, energy_.improper_dih, energy_.cmap_dih,
        // VdW (3)
        energy_.lj_sr, energy_.lj_14, energy_.disper_corr,
        // Thermodynamic state (8)
        energy_.potential, energy_.kinetic, energy_.total_energy,
        energy_.enthalpy, energy_.temperature, energy_.pressure,
        energy_.volume, energy_.density,
        // Box (3)
        energy_.box_x, energy_.box_y, energy_.box_z,
        // Virial tensor (9)
        energy_.vir[0], energy_.vir[1], energy_.vir[2],
        energy_.vir[3], energy_.vir[4], energy_.vir[5],
        energy_.vir[6], energy_.vir[7], energy_.vir[8],
        // Pressure tensor (9)
        energy_.pres[0], energy_.pres[1], energy_.pres[2],
        energy_.pres[3], energy_.pres[4], energy_.pres[5],
        energy_.pres[6], energy_.pres[7], energy_.pres[8],
        // Per-group temperature (2)
        energy_.T_protein, energy_.T_non_protein,
    };

    NpyWriter::WriteFloat64(output_dir + "/gromacs_energy.npy",
                            row.data(), 1, static_cast<int>(row.size()));
    return 1;
}

}  // namespace nmr
