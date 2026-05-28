// SingleConformation — one pose, no H5. The single-conformation sibling of
// TrajectoryConformation under the shared Conformation base (NOT a faked
// one-frame trajectory — "subclass, don't hack").
//
// Backs a --orca / --mutant / --pdb single-pose run: sidecar + manifest →
// QtProtein (loaded today, unchanged) + the run-root per-atom NPYs → one
// QtConformationSnapshot. There is nothing to animate (frameCount() == 1) and
// no trajectory.h5; positions come from the snapshot's Pos column rather than a
// dense per-TR buffer. The reader never writes H5 — the single conformation is
// constructed in memory.

#pragma once

#include "Conformation.h"

#include <cstddef>
#include <memory>

namespace h5reader::model {

class QtProtein;
class QtConformationSnapshot;

class SingleConformation final : public Conformation {
    Q_OBJECT

public:
    // `pose` is the run-root snapshot (loaded by FrameNpyLoader on the run
    // directory). Held for the conformation's lifetime because a single-pose
    // run has exactly one source frame.
    SingleConformation(const QtProtein* protein, std::shared_ptr<const QtConformationSnapshot> pose);
    ~SingleConformation() override;

    std::size_t frameCount() const override { return 1; }
    double timePicoseconds(std::size_t) const override { return 0.0; }
    Vec3 atomPosition(std::size_t frame, std::size_t atomIdx) const override;

protected:
    std::shared_ptr<const QtConformationSnapshot> loadSnapshot(std::size_t frame) override;

private:
    std::shared_ptr<const QtConformationSnapshot> pose_;
};

}  // namespace h5reader::model
