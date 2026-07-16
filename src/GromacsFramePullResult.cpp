#include "GromacsFramePullResult.h"
#include "NpyWriter.h"

namespace nmr {

std::unique_ptr<GromacsFramePullResult> GromacsFramePullResult::Compute(
        ProteinConformation& /*conf*/,
        const std::vector<Vec3>* velocities,
        const Eigen::Matrix3d* box_matrix) {

    // Catch-all rule: if the frame yielded nothing, attach nothing.
    // The opts gate in OperationRunner skips the call when both
    // pointers are null; defend the contract here too.
    if (!velocities && !box_matrix) {
        return nullptr;
    }

    auto result = std::make_unique<GromacsFramePullResult>();
    if (velocities) result->velocities_ = *velocities;
    if (box_matrix) result->box_matrix_ = *box_matrix;
    return result;
}


int GromacsFramePullResult::WriteFeatures(
        const ProteinConformation& /*conf*/,
        const std::string& output_dir) const {
    if (!HasBoxMatrix()) return 0;

    // One frame-level row. Components are the exact applied TRR cell in
    // Eigen storage convention: a 3x3 matrix whose columns are lattice
    // vectors, flattened in row-major order for the NPY boundary.
    double row[9];
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            row[i * 3 + j] = box_matrix_(i, j);
    NpyWriter::WriteFloat64(output_dir + "/gromacs_box.npy", row, 1, 9);
    return 1;
}

}  // namespace nmr
