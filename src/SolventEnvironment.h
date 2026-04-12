#pragma once
//
// SolventEnvironment: water and ion positions + charges from a
// full-system GROMACS trajectory frame.
//
// Built by reading the full-system .xtc (all atoms) and the TPR
// topology (which identifies protein, water, ion atoms).  The
// protein positions go to the existing calculators.  This struct
// holds everything else — the explicit solvent that APBS
// approximates with a continuum.
//
// Coordinates are in Angstroms (converted from GROMACS nm).
// Charges are in elementary charges (e).
//

#include <Eigen/Dense>
#include <vector>

namespace nmr {

using Vec3 = Eigen::Vector3d;

// One water molecule: O position + two H positions + charges.
struct WaterMolecule {
    Vec3   O_pos;
    Vec3   H1_pos;
    Vec3   H2_pos;
    double O_charge;   // typically -0.834e (TIP3P) or similar
    double H_charge;   // typically +0.417e (TIP3P)

    // Dipole moment vector (from O toward midpoint of H's, scaled by charges)
    Vec3 Dipole() const {
        Vec3 mid_H = 0.5 * (H1_pos + H2_pos);
        // μ = Σ q_i * r_i (relative to O)
        return O_charge * Vec3::Zero()
             + H_charge * (H1_pos - O_pos)
             + H_charge * (H2_pos - O_pos);
    }
};

// One ion: position + charge + element.
struct Ion {
    Vec3   pos;
    double charge;     // +1 for Na+, -1 for Cl-, etc.
    int    atomic_number;  // 11=Na, 17=Cl, 19=K, etc.
};

// Full solvent environment for one trajectory frame.
struct SolventEnvironment {
    std::vector<WaterMolecule> waters;
    std::vector<Ion> ions;

    // All solvent atom positions for spatial indexing.
    // water_O_positions[i] corresponds to waters[i].
    std::vector<Vec3> water_O_positions;

    bool Empty() const { return waters.empty() && ions.empty(); }
    size_t WaterCount() const { return waters.size(); }
    size_t IonCount() const { return ions.size(); }
};

}  // namespace nmr
