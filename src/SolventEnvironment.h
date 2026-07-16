#pragma once
//
// SolventEnvironment: water and ion positions + charges from a
// full-system GROMACS trajectory frame.
//
// Built by splitting a full-system GROMACS trajectory frame with the
// TPR topology that identifies protein, water, and ion atoms. The
// protein positions go to the existing calculators. This struct holds
// everything else for explicit-solvent calculators.
//
// Coordinates are in Angstroms (converted from GROMACS nm).
// Charges are in elementary charges (e).
//

#include <Eigen/Dense>
#include <memory>
#include <vector>

enum class PbcType : int;

namespace nmr {

namespace detail {
struct SolventPbcState;
}

using Vec3 = Eigen::Vector3d;

// One water molecule: O position + two H positions + charges.
struct WaterMolecule {
    Vec3   O_pos;
    Vec3   H1_pos;
    Vec3   H2_pos;
    double O_charge;   // typically -0.834e (TIP3P) or similar
    double H_charge;   // typically +0.417e (TIP3P)

    // Dipole moment vector: μ = Σ q_i (r_i − r_O). O's own term is zero;
    // each H contributes q_H (r_H − r_O).
    Vec3 Dipole() const {
        return H_charge * (H1_pos - O_pos)
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

    // Attach the authoritative TPR PBC type and this frame's full TRR cell.
    // The matrix uses the producer convention: lattice vectors are columns,
    // coordinates and cell are both in Angstroms.
    void SetPeriodicCell(const Eigen::Matrix3d& box_matrix,
                         PbcType pbc_type);
    bool HasPeriodicCell() const { return static_cast<bool>(pbc_state_); }
    const Eigen::Matrix3d& BoxMatrix() const;

    // GROMACS minimum-image x1-x2 displacement. Without a periodic cell
    // (the legal static/synthetic path), this is the ordinary difference.
    Vec3 MinimumImageDisplacement(const Vec3& x1, const Vec3& x2) const;

    // Put an extraction-whole water in the image whose oxygen is nearest
    // target while preserving its already-canonical H-O vectors.
    WaterMolecule WaterAtNearestImage(const WaterMolecule& water,
                                      const Vec3& target) const;

private:
    std::shared_ptr<const detail::SolventPbcState> pbc_state_;
};

}  // namespace nmr
