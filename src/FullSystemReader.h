#pragma once
//
// FullSystemReader: reads full-system GROMACS topology and trajectory
// frames, splitting atoms into protein + solvent.
//
// The TPR contains the topology for ALL atoms (protein, water, ions).
// The full-system .xtc contains positions for ALL atoms at each frame.
// This reader:
//   1. Reads the TPR to identify atom ranges (protein, water, ions)
//   2. Reads a full-system .xtc frame
//   3. Splits positions: protein → Vec3 vector, solvent → SolventEnvironment
//
// The protein positions feed into MDFrameConformation as usual.
// The SolventEnvironment feeds into WaterFieldResult and friends.
//
// Coordinates are in Angstroms (converted from nm).
//

#include "SolventEnvironment.h"
#include "BondedEnergyResult.h"

#include <Eigen/Dense>
#include <string>
#include <vector>

namespace nmr {

// Topology layout extracted from TPR.
struct SystemTopology {
    // Atom ranges in the full-system frame.
    size_t protein_start = 0;
    size_t protein_count = 0;
    size_t water_O_start = 0;   // first water oxygen index
    size_t water_count = 0;     // number of water MOLECULES (3 atoms each)
    size_t ion_start = 0;
    size_t ion_count = 0;
    size_t total_atoms = 0;

    // Per-water-molecule charge (from first water in topology).
    double water_O_charge = 0.0;
    double water_H_charge = 0.0;

    // Per-ion charges and elements (one per ion).
    std::vector<double> ion_charges;
    std::vector<int>    ion_atomic_numbers;
};


class FullSystemReader {
public:
    // Read the TPR to build the SystemTopology.
    // Returns false on error (check error()).
    bool ReadTopology(const std::string& tpr_path);

    // Given a full-system XTC frame (all atoms, in nm),
    // extract protein positions (in Angstroms) and SolventEnvironment.
    // protein_positions will have protein_count entries.
    bool ExtractFrame(const std::vector<float>& full_frame_xyz,
                      std::vector<Vec3>& protein_positions,
                      SolventEnvironment& solvent) const;

    // Extract bonded interaction parameters for the protein from the TPR.
    // Must be called after ReadTopology(). Populates params with all
    // bond, angle, UB, proper dihedral, improper dihedral, and CMAP
    // interactions for the protein atoms. Atom indices are protein-local
    // (0-based, same as our Protein).
    bool ExtractBondedParameters(const std::string& tpr_path,
                                 BondedParameters& params) const;

    const SystemTopology& Topology() const { return topo_; }
    const std::string& error() const { return error_; }

private:
    SystemTopology topo_;
    std::string error_;
};

}  // namespace nmr
