#pragma once
//
// GromacsProtein: wraps a real Protein for trajectory processing.
//
// Owns the immutable Protein (topology + conformation 0 via AddMDFrame).
// Streaming conformations are free-standing — created from the Protein's
// topology but never added to its conformations_ vector.
//
// Holds trajectory-specific state: frame manifest, accumulation buffers.
// GromacsFrameHandler processes each frame.
// GromacsFinalResult writes the output.
//

#include "Protein.h"
#include "GromacsEnsembleLoader.h"
#include "OperationRunner.h"
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class GromacsProtein {
public:
    // Build from fleet paths. Returns false on error (check error()).
    bool Build(const FleetPaths& paths);

    // The real protein. Conformation 0 is in its vector.
    Protein& protein() { return *protein_; }
    const Protein& protein() const { return *protein_; }

    // Charges from the TPR — caller attaches via ChargeAssignmentResult.
    ChargeSource* charges() { return charges_.get(); }

    // Frame manifest — paths where frames were written.
    const std::vector<std::string>& frame_paths() const { return frame_paths_; }
    void AddFramePath(const std::string& path) { frame_paths_.push_back(path); }

    // Protein identity.
    const std::string& protein_id() const { return protein_id_; }

    // Error from Build.
    const std::string& error() const { return error_; }

    // Net charge from TPR.
    int net_charge() const { return net_charge_; }

    // Pose names from PDB filenames (fes-sampler naming).
    const std::vector<std::string>& pose_names() const { return pose_names_; }

private:
    std::unique_ptr<Protein> protein_;
    std::unique_ptr<ChargeSource> charges_;
    std::vector<std::string> frame_paths_;
    std::vector<std::string> pose_names_;
    std::string protein_id_;
    std::string error_;
    int net_charge_ = 0;
};

}  // namespace nmr
