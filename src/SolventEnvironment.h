#pragma once
// Coordinates are in Angstroms (converted from GROMACS nm).
// Charges are in elementary charges (e).

#include <Eigen/Dense>
#include <vector>

namespace nmr {

using Vec3 = Eigen::Vector3d;

struct WaterMolecule {
    Vec3   O_pos;
    Vec3   H1_pos;
    Vec3   H2_pos;
    double O_charge;
    double H_charge;

    // Dipole moment vector: μ = Σ q_i (r_i − r_O). O's own term is zero;
    // each H contributes q_H (r_H − r_O).
    Vec3 Dipole() const {
        return H_charge * (H1_pos - O_pos)
             + H_charge * (H2_pos - O_pos);
    }
};

struct Ion {
    Vec3   pos;
    double charge;
    int    atomic_number;
};

struct SolventEnvironment {
    std::vector<WaterMolecule> waters;
    std::vector<Ion> ions;

    // water_O_positions[i] corresponds to waters[i].
    std::vector<Vec3> water_O_positions;

    bool Empty() const { return waters.empty() && ions.empty(); }
    size_t WaterCount() const { return waters.size(); }
    size_t IonCount() const { return ions.size(); }
};

}  // namespace nmr
