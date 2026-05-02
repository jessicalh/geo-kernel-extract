#pragma once
//
// GromacsFramePullResult: catch-all ConformationResult for per-frame
// data the GROMACS frame pull yields.
//
// Source-available capture at the trajectory-frame boundary. One read,
// all of it. The CR exists because trying to "split this into pieces"
// (one CR per data-source-substream) turns into nonsense over time;
// one CR holds whatever the frame gave us, and any consumer that
// needs a piece of it pulls from here.
//
// Currently captured:
//   - per-atom velocities       (Å/ps; empty when the load path has
//                                no velocity stream — XTC-only legacy,
//                                TRR with nstvout=0, non-trajectory load)
//   - simulation box matrix     (Å; zero when the load path doesn't
//                                carry a periodic cell)
//
// Future scope (intentional, named so the catch-all stays one thing):
// solvent slice and EDR-row-at-frame-time will likely fold here. They
// currently live in `Trajectory::env_` with their own consumers
// (WaterFieldResult, HydrationShellResult, HydrationGeometryResult,
// GromacsEnergyResult). Folding is a separate migration.
//
// Calculator-side rule: per-frame data that varies by load path goes
// here, never as direct fields on ProteinConformation. PDB and GROMACS
// are both legal load paths; ConformationAtom / ProteinConformation
// should not carry fields that one path leaves empty. See memory
// `feedback_capture_at_the_boundary`.
//
// Dependencies: none. Reads from RunOptions pointers populated from
// TrajectoryEnv each frame; not a calculator.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "Types.h"

#include <Eigen/Core>
#include <memory>
#include <string>
#include <vector>

namespace nmr {

class GromacsFramePullResult : public ConformationResult {
public:
    std::string Name() const override { return "GromacsFramePullResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Factory: stash whatever the frame yielded. Either pointer may be
    // null (XTC-only legacy, non-trajectory load); the result records
    // empty / zero in that case. Returns nullptr if both are null —
    // there is nothing to capture and the OperationRunner gate
    // shouldn't have called us.
    static std::unique_ptr<GromacsFramePullResult> Compute(
        ProteinConformation& conf,
        const std::vector<Vec3>* velocities,
        const Eigen::Matrix3d* box_matrix);

    // Empty when the load path doesn't carry velocities (XTC-only,
    // TRR with nstvout=0).
    const std::vector<Vec3>& Velocities() const { return velocities_; }
    bool HasVelocities() const { return !velocities_.empty(); }

    // Zero matrix when the load path doesn't carry a periodic cell.
    const Eigen::Matrix3d& BoxMatrix() const { return box_matrix_; }
    bool HasBoxMatrix() const { return !box_matrix_.isZero(); }

private:
    std::vector<Vec3> velocities_;
    Eigen::Matrix3d box_matrix_ = Eigen::Matrix3d::Zero();
};

}  // namespace nmr
