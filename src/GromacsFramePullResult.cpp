#include "GromacsFramePullResult.h"

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

}  // namespace nmr
