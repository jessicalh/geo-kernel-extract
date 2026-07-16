#include "SolventEnvironment.h"

#include "gromacs/mdtypes/md_enums.h"
#include "gromacs/pbcutil/pbc.h"

namespace nmr {

namespace detail {

struct SolventPbcState {
    SolventPbcState(const Eigen::Matrix3d& box_matrix, PbcType pbc_type)
        : box_matrix(box_matrix) {
        matrix gromacs_box;
        for (int i = 0; i < DIM; ++i) {
            for (int j = 0; j < DIM; ++j) {
                // Eigen stores lattice vectors as columns; GROMACS stores
                // the same vectors as rows.
                gromacs_box[i][j] = static_cast<real>(box_matrix(j, i));
            }
        }
        set_pbc(&pbc, pbc_type, gromacs_box);
    }

    Eigen::Matrix3d box_matrix;
    t_pbc pbc{};
};

}  // namespace detail


void SolventEnvironment::SetPeriodicCell(
        const Eigen::Matrix3d& box_matrix,
        PbcType pbc_type) {
    pbc_state_ = std::make_shared<const detail::SolventPbcState>(
        box_matrix, pbc_type);
}


const Eigen::Matrix3d& SolventEnvironment::BoxMatrix() const {
    static const Eigen::Matrix3d no_box = Eigen::Matrix3d::Zero();
    return pbc_state_ ? pbc_state_->box_matrix : no_box;
}


Vec3 SolventEnvironment::MinimumImageDisplacement(
        const Vec3& x1,
        const Vec3& x2) const {
    if (!pbc_state_) return x1 - x2;

    dvec gmx_x1 = {x1.x(), x1.y(), x1.z()};
    dvec gmx_x2 = {x2.x(), x2.y(), x2.z()};
    dvec gmx_dx;
    pbc_dx_d(&pbc_state_->pbc, gmx_x1, gmx_x2, gmx_dx);
    return Vec3(gmx_dx[0], gmx_dx[1], gmx_dx[2]);
}


WaterMolecule SolventEnvironment::WaterAtNearestImage(
        const WaterMolecule& water,
        const Vec3& target) const {
    if (!pbc_state_) return water;

    WaterMolecule imaged = water;
    imaged.O_pos = target + MinimumImageDisplacement(water.O_pos, target);
    imaged.H1_pos = imaged.O_pos + (water.H1_pos - water.O_pos);
    imaged.H2_pos = imaged.O_pos + (water.H2_pos - water.O_pos);
    return imaged;
}

}  // namespace nmr
