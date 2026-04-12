#include "GromacsProtein.h"
#include "OperationLog.h"

#include <filesystem>

namespace nmr {

bool GromacsProtein::Build(const FleetPaths& paths) {
    auto result = BuildFromGromacs(paths);
    if (!result.ok) {
        error_ = result.error;
        return false;
    }

    protein_ = std::move(result.protein);
    charges_ = std::move(result.charges);
    net_charge_ = result.net_charge;
    pose_names_ = std::move(result.pose_names);

    protein_id_ = std::filesystem::path(
        paths.sampled_poses_dir).filename().string();

    OperationLog::Info(LogCalcOther, "GromacsProtein::Build",
        protein_id_ + ": " +
        std::to_string(protein_->AtomCount()) + " atoms, " +
        std::to_string(protein_->ResidueCount()) + " residues, " +
        std::to_string(protein_->ConformationCount()) + " poses loaded");

    return true;
}

}  // namespace nmr
