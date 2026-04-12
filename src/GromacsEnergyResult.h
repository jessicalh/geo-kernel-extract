#pragma once
//
// GromacsEnergyResult: per-frame energy terms from the GROMACS simulation.
//
// Reads the .edr energy file and extracts the frame matching the
// conformation's time (ps). Stores the GROMACS-computed electrostatic
// energy (PME: short-range + reciprocal), LJ, potential, temperature,
// and pressure components.
//
// These are aggregate (whole-system) quantities, not per-atom. They
// characterise the electrostatic environment that GROMACS computed
// during the actual MD simulation with explicit solvent and PME.
//
// Use case: compare GROMACS electrostatic energy per frame against
// APBS per-atom fields to quantify what information the trajectory
// carries that our calculators do or don't capture.
//
// Dependencies: none (reads from external file, not from conformation).
//

#include "ConformationResult.h"
#include "ProteinConformation.h"

#include <string>
#include <vector>

namespace nmr {

// Per-frame energy terms extracted from GROMACS .edr file.
struct GromacsEnergy {
    double time_ps       = 0.0;   // frame time (ps)
    double coulomb_sr    = 0.0;   // PME real-space Coulomb (kJ/mol)
    double coulomb_recip = 0.0;   // PME reciprocal-space Coulomb (kJ/mol)
    double coulomb_14    = 0.0;   // 1-4 intramolecular Coulomb (kJ/mol)
    double lj_sr         = 0.0;   // Lennard-Jones short-range (kJ/mol)
    double potential     = 0.0;   // total potential energy (kJ/mol)
    double temperature   = 0.0;   // system temperature (K)
    double pressure      = 0.0;   // system pressure (bar)
    double volume        = 0.0;   // box volume (nm^3)

    // Total Coulomb = SR + reciprocal (excludes 1-4 which is bookkeeping)
    double CoulombTotal() const { return coulomb_sr + coulomb_recip; }
};


class GromacsEnergyResult : public ConformationResult {
public:
    std::string Name() const override { return "GromacsEnergyResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return {};  // reads from file, no conformation dependencies
    }

    // Factory: read .edr file, find frame nearest to target_time_ps.
    // Returns nullptr if .edr cannot be read or no matching frame found.
    static std::unique_ptr<GromacsEnergyResult> Compute(
        ProteinConformation& conf,
        const std::string& edr_path,
        double target_time_ps);

    const GromacsEnergy& Energy() const { return energy_; }

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    GromacsEnergy energy_;
};

}  // namespace nmr
