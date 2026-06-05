// TransformedConformation — decorator over a Conformation that applies a
// per-frame rigid-body transform (R, T) at the atomPosition() seam.
//
// Implements the upstream "data transform" layer described in
// feedback_viewer_two_layers_transform_and_camera. The raw GROMACS-output
// trajectory has 6 rigid-body degrees of freedom relative to the
// simulation box (3 translation + 3 rotation) which both drift over MD
// time. None of this is removed before display. The decorator lets a
// downstream consumer (renderer, picker, overlays) see a stabilised
// frame: RMSD-fit to a reference using either all atoms or a typed subset
// (the default UI subset is the backbone).
//
// Architectural shape (per the memory entry's prescription): the wrapper
// is a Conformation itself, holding a Conformation* inner (non-owning).
// All Conformation virtuals delegate to inner EXCEPT atomPosition(frame,
// atom) which applies the per-frame transform. Consumers that hold a
// Conformation* (MoleculeScene, MeasurementOverlay, picker, REST
// /positions) see the wrapped one; polymorphism does the rest.
//
// Per-frame transforms are cached as a whole trajectory sequence. The
// reference is the converged iterative mean of the selected fit atoms; each
// frame's raw Kabsch rotation is then optionally smoothed. Translation is
// always derived from the same-frame fit centroid so the fit-set centroid maps
// exactly to a constant reference anchor before atomPosition() applies
// R * raw + T. The cache invalidates when the transform mode, reference anchor,
// or smoothing window changes. Cache access is not thread-safe but the reader
// is single-threaded on the GUI thread; ASSERT_THREAD guards entry points.
//
// PBC unwrap: deliberately NOT implemented in this decorator. The
// canonical PBC unwrap (fes-sampler's pbc_whole.h via do_pbc_mtop)
// requires libgromacs which h5-reader does not link by policy
// (CLAUDE.md: standalone, never links the library). Per
// feedback_pbc_verbatim the rule is "port verbatim or skip cleanly".
// We skip cleanly here. The RMSD fit modes still deliver most of the
// stabilisation value on a trajectory whose PBC unwrap was already
// done at extraction time (the typical case for 1P9J and friends).
//
// ReaderMainWindow sets the startup mode before the scene builds:
// FitSubset over the typed backbone subset, seeded/anchored at frame 0. That
// keeps the displayed molecule stationary on open while preserving all internal
// motion because the transform is rigid.

#pragma once

#include "Conformation.h"
#include "Types.h"

#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

namespace h5reader::model {

class QtConformationSnapshot;
class TrajectoryConformation;

class TransformedConformation final : public Conformation {
    Q_OBJECT

public:
    static constexpr int kDefaultStabilisationWindow = 0;

    enum class Mode {
        FitReference,     // Kabsch fit to iterative mean using all atoms
        FitSubset,        // Kabsch fit to iterative mean using provided subset atom indices
    };

    // `inner` is the underlying Conformation (TrajectoryConformation or
    // SingleConformation). NOT owned — the wrapper does not extend its
    // lifetime. Caller (ReaderMainWindow) must outlive both. `parent` is
    // the Qt parent (ReaderMainWindow) so destruction order is determined.
    explicit TransformedConformation(Conformation* inner, QObject* parent = nullptr);
    ~TransformedConformation() override;

    // ----- Conformation seam (delegated to inner unless noted) -----
    std::size_t frameCount() const override;
    double      timePicoseconds(std::size_t frame) const override;
    std::size_t originalFrameIndex(std::size_t frame) const override;
    const TrajectoryConformation* asTrajectory() const override;

    // The ONE virtual this decorator actually decorates: applies the
    // per-frame rigid-body transform to the inner conformation's raw
    // position. R * raw + T, from the precomputed transform sequence.
    Vec3 atomPosition(std::size_t frame, std::size_t atomIdx) const override;

    // ----- Transform control -----
    Mode mode() const { return mode_; }
    std::size_t referenceFrame() const { return referenceFrame_; }
    const std::vector<std::size_t>& subsetAtoms() const { return subsetAtoms_; }
    int stabilisationWindow() const { return stabilisationWindow_; }

    // Switch transform mode. `referenceFrame` seeds the iterative mean and
    // anchors the displayed fit-set centroid for both fit modes.
    // `subsetAtoms` is used only by FitSubset. Bumps the generation counter and clears the
    // per-frame cache atomically; emits transformChanged() so consumers
    // can request a re-render. ASSERT_THREAD(this).
    void setMode(Mode mode,
                 std::size_t referenceFrame = 0,
                 std::vector<std::size_t> subsetAtoms = {});

    // Set the symmetric smoothing half-width in frames. 0 means no temporal
    // smoothing: atomPosition() uses the raw per-frame Kabsch rotation with
    // centroid-pinned translation. Non-zero windows smooth only rotation;
    // translation is still re-derived from the current fit-set centroid.
    void setStabilisationWindow(int halfWidth);

    // Convenience for the harness / REST: build a backbone-only subset
    // for FitSubset by walking QtProtein's atoms and selecting those
    // with QtAtom::IsBackbone() == true.
    static std::vector<std::size_t> BackboneSubset(const QtProtein& protein);

signals:
    // Emitted after setMode() finishes. ReaderMainWindow connects this
    // to refreshCurrentFrame so the scene re-evaluates atom positions
    // for the new transform without waiting for the next playback tick.
    void transformChanged();

protected:
    // Delegates to inner. The base Conformation API contract requires the
    // snapshot facade on this object (so REST /selection/instrument etc.
    // continue to work via the wrapper), but the actual loader lives on
    // the inner subclass — we forward unchanged.
    std::shared_ptr<const QtConformationSnapshot> loadSnapshot(std::size_t frame) override;

private:
    // (rotation, translation). atomPosition returns R * raw + T.
    struct Transform3D {
        Mat3 R = Mat3::Identity();
        Vec3 T = Vec3::Zero();
    };

    struct FrameFit {
        Transform3D transform;
        Vec3 currentCentroid = Vec3::Zero();
    };

    // Rebuild the whole-trajectory cache after mode/reference/window changes.
    void rebuildTransformCache();
    void rebuildReferenceMean();

    std::vector<Vec3> fitPositions(std::size_t frame) const;

    // Compute the unsmoothed transform for `frame` from the inner
    // conformation's RAW positions using fitAtomIndices_.
    FrameFit computeRawFrameFit(std::size_t frame) const;
    Transform3D computeRawTransform(std::size_t frame) const;

    std::vector<FrameFit> computeRawTransformSequence() const;
    std::vector<Transform3D> smoothTransformSequence(const std::vector<FrameFit>& raw) const;

    // Kabsch algorithm — compute the rotation+translation that minimises
    // sum of squared distances between `current` atoms and the reference's
    // same atoms. Both vectors hold the SAME N atoms in the SAME order.
    // Returns identity transform if N < 3 (degenerate / underdetermined).
    static Transform3D KabschFit(const std::vector<Vec3>& current,
                                 const std::vector<Vec3>& reference);

    Conformation* inner_ = nullptr;

    Mode mode_ = Mode::FitReference;
    std::size_t referenceFrame_ = 0;
    std::vector<std::size_t> subsetAtoms_;
    std::vector<std::size_t> fitAtomIndices_;
    int stabilisationWindow_ = kDefaultStabilisationWindow;

    // One display transform per frame, rebuilt on setMode() /
    // setStabilisationWindow().
    std::vector<Transform3D> transformCache_;
    mutable std::uint64_t generation_ = 0;

    // Cached converged-mean reference positions for FitReference / FitSubset,
    // in fitAtomIndices_ order and centred around zero. referenceCentroid_
    // is the constant world anchor for the fit-set centroid.
    mutable std::vector<Vec3> referencePositions_;
    mutable Vec3 referenceCentroid_ = Vec3::Zero();
};

}  // namespace h5reader::model
