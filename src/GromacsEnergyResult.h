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
// All energies in kJ/mol (GROMACS native), all else in stated units.
struct GromacsEnergy {
    double time_ps       = 0.0;   // frame time (ps)

    // ── Electrostatic ──────────────────────────────────────────
    double coulomb_sr    = 0.0;   // PME real-space Coulomb (kJ/mol)
    double coulomb_recip = 0.0;   // PME reciprocal-space Coulomb (kJ/mol)
    double coulomb_14    = 0.0;   // 1-4 intramolecular Coulomb (kJ/mol)

    // ── Bonded (internal strain) ───────────────────────────────
    double bond          = 0.0;   // bond stretching (kJ/mol)
    double angle         = 0.0;   // harmonic angle bending (kJ/mol)
    double urey_bradley  = 0.0;   // Urey-Bradley 1-3 distance (kJ/mol)
    double proper_dih    = 0.0;   // proper dihedral (kJ/mol)
    double improper_dih  = 0.0;   // improper dihedral — planarity (kJ/mol)
    double cmap_dih      = 0.0;   // CMAP backbone correction (kJ/mol)

    // ── Van der Waals ──────────────────────────────────────────
    double lj_sr         = 0.0;   // Lennard-Jones short-range (kJ/mol)
    double lj_14         = 0.0;   // 1-4 LJ (kJ/mol)
    double disper_corr   = 0.0;   // long-range dispersion correction (kJ/mol)

    // ── Thermodynamic state ────────────────────────────────────
    double potential     = 0.0;   // total potential energy (kJ/mol)
    double kinetic       = 0.0;   // total kinetic energy (kJ/mol)
    double total_energy  = 0.0;   // potential + kinetic (kJ/mol)
    double enthalpy      = 0.0;   // H = E + pV (kJ/mol)
    double temperature   = 0.0;   // system temperature (K)
    double pressure      = 0.0;   // scalar pressure (bar)
    double volume        = 0.0;   // box volume (nm^3)
    double density       = 0.0;   // system density (kg/m^3)

    // ── Box dimensions (NPT cell breathing) ────────────────────
    double box_x         = 0.0;   // (nm)
    double box_y         = 0.0;   // (nm)
    double box_z         = 0.0;   // (nm)

    // ── Virial tensor (internal stress, 3x3 symmetric) ─────────
    double vir[9]        = {};    // XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ (kJ/mol)

    // ── Pressure tensor (3x3) ──────────────────────────────────
    double pres[9]       = {};    // XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ (bar)

    // ── Per-group temperature ──────────────────────────────────
    double T_protein     = 0.0;   // protein group temperature (K)
    double T_non_protein = 0.0;   // solvent+ions temperature (K)

    // Total Coulomb = SR + reciprocal (excludes 1-4 which is bookkeeping)
    double CoulombTotal() const { return coulomb_sr + coulomb_recip; }
};


class GromacsEnergyResult : public ConformationResult {
public:
    std::string Name() const override { return "GromacsEnergyResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return {};  // reads from file, no conformation dependencies
    }

    // Factory: create from preloaded energy data (from GromacsRunContext).
    static std::unique_ptr<GromacsEnergyResult> Compute(
        ProteinConformation& conf,
        const GromacsEnergy& energy);

    const GromacsEnergy& Energy() const { return energy_; }

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    GromacsEnergy energy_;
};

}  // namespace nmr
