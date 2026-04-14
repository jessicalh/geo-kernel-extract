#pragma once
//
// BondedEnergyResult: per-atom bonded energy decomposition from
// GROMACS CHARMM36m force field parameters.
//
// Reads bonded interaction lists and parameters (extracted from TPR
// at build time). For each frame, evaluates the force field energy
// functions from conformation positions and accumulates per-atom.
//
// Energy types: bond stretch, angle bend, Urey-Bradley, proper
// dihedral, improper dihedral, CMAP correction. Each stored as a
// per-atom double (kJ/mol, split evenly among participating atoms).
//
// Dependencies: GeometryResult (needs positions, but that's always
// present). The bonded parameters come through RunOptions, not from
// another ConformationResult.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace nmr {

// One bonded interaction extracted from the GROMACS TPR.
// Atom indices are protein-local (0-based, matching our Protein).
struct BondedInteraction {
    enum Type : uint8_t {
        Bond,           // harmonic bond stretch
        Angle,          // harmonic angle bend
        UreyBradley,    // Urey-Bradley 1-3 distance term
        ProperDih,      // periodic proper dihedral
        ImproperDih,    // harmonic improper dihedral
        CMAP            // dihedral energy correction map
    };

    Type type;
    // Atom indices (2 for bond, 3 for angle/UB, 4 for dihedral, 5 for CMAP)
    size_t atoms[5] = {};
    int n_atoms = 0;

    // Parameters (meaning depends on type):
    // Bond:       p[0]=r0 (nm), p[1]=k (kJ/mol/nm^2)
    // Angle:      p[0]=theta0 (rad), p[1]=k (kJ/mol/rad^2)
    // UreyBradley: p[0]=r13_0 (nm), p[1]=k_ub (kJ/mol/nm^2)
    // ProperDih:  p[0]=phi0 (rad), p[1]=k (kJ/mol), p[2]=multiplicity
    // ImproperDih: p[0]=phi0 (rad), p[1]=k (kJ/mol/rad^2)
    // CMAP:       p[0]=cmap_type_index (into cmap grid array)
    double p[3] = {};
};

// All bonded interactions for the protein, extracted from TPR once.
// Stored on GromacsProtein, passed to BondedEnergyResult via RunOptions.
struct BondedParameters {
    std::vector<BondedInteraction> interactions;

    // CMAP grids: each grid is (grid_spacing x grid_spacing) doubles.
    // Indexed by BondedInteraction::p[0] for CMAP type.
    int cmap_grid_spacing = 0;
    std::vector<std::vector<double>> cmap_grids;  // [type][phi_idx * spacing + psi_idx]
};


class BondedEnergyResult : public ConformationResult {
public:
    std::string Name() const override { return "BondedEnergyResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return {};  // only needs positions (always present)
    }

    // Factory: compute per-atom bonded energies from positions + parameters.
    // params must outlive this call but is not stored.
    static std::unique_ptr<BondedEnergyResult> Compute(
        ProteinConformation& conf,
        const BondedParameters& params);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

    // Per-atom accumulated energies (kJ/mol, split among participating atoms)
    const std::vector<double>& BondEnergy()     const { return bond_energy_; }
    const std::vector<double>& AngleEnergy()    const { return angle_energy_; }
    const std::vector<double>& UBEnergy()       const { return ub_energy_; }
    const std::vector<double>& ProperDihEnergy() const { return proper_energy_; }
    const std::vector<double>& ImproperDihEnergy() const { return improper_energy_; }
    const std::vector<double>& CmapEnergy()     const { return cmap_energy_; }
    const std::vector<double>& TotalBonded()    const { return total_bonded_; }

private:
    std::vector<double> bond_energy_;      // (N,) per atom
    std::vector<double> angle_energy_;     // (N,)
    std::vector<double> ub_energy_;        // (N,)
    std::vector<double> proper_energy_;    // (N,)
    std::vector<double> improper_energy_;  // (N,)
    std::vector<double> cmap_energy_;      // (N,)
    std::vector<double> total_bonded_;     // (N,) sum of all above
};

}  // namespace nmr
