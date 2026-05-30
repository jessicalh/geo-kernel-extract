// TransformedConformation — decorator over a Conformation that applies a
// per-frame rigid-body transform (R, T) at the atomPosition() seam.
//
// Implements the upstream "data transform" layer described in
// feedback_viewer_two_layers_transform_and_camera. The raw GROMACS-output
// trajectory has 6 rigid-body degrees of freedom relative to the
// simulation box (3 translation + 3 rotation) which both drift over MD
// time. None of this is removed before display. The decorator lets a
// downstream consumer (renderer, picker, overlays) see a stabilised
// frame: COM-centred, RMSD-fit to a reference, or RMSD-fit on a subset
// (e.g. backbone).
//
// Architectural shape (per the memory entry's prescription): the wrapper
// is a Conformation itself, holding a Conformation* inner (non-owning).
// All Conformation virtuals delegate to inner EXCEPT atomPosition(frame,
// atom) which applies the per-frame transform. Consumers that hold a
// Conformation* (MoleculeScene, MeasurementOverlay, picker, REST
// /positions) see the wrapped one; polymorphism does the rest.
//
// Per-frame transforms are cached. The cache invalidates when the
// transform mode changes (setMode bumps a generation counter and clears
// the table). Cache lookups are not thread-safe but the reader is single-
// threaded on the GUI thread; ASSERT_THREAD guards entry points.
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
// Default mode is `Identity` — same atomPosition(frame, atom) as the
// inner Conformation, no allocation, no cache use. Behaviour at
// startup is therefore unchanged.

#pragma once

#include "Conformation.h"
#include "Types.h"

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <vector>

namespace h5reader::model {

class QtConformationSnapshot;
class TrajectoryConformation;

class TransformedConformation final : public Conformation {
    Q_OBJECT

public:
    enum class Mode {
        Identity,         // R=I, T=0 — equivalent to "no transform"
        CenterCom,        // R=I, T=-COM(frame, all_atoms)
        FitReference,     // Kabsch fit of frame against reference frame, all atoms
        FitSubset,        // Kabsch fit using only the provided subset atom indices
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
    // position. R * raw + T, computed on demand and cached per frame.
    Vec3 atomPosition(std::size_t frame, std::size_t atomIdx) const override;

    // ----- Transform control -----
    Mode mode() const { return mode_; }
    std::size_t referenceFrame() const { return referenceFrame_; }
    const std::vector<std::size_t>& subsetAtoms() const { return subsetAtoms_; }

    // Switch transform mode. `referenceFrame` is used by FitReference /
    // FitSubset (ignored by Identity / CenterCom). `subsetAtoms` is used
    // only by FitSubset. Bumps the generation counter and clears the
    // per-frame cache atomically; emits transformChanged() so consumers
    // can request a re-render. ASSERT_THREAD(this).
    void setMode(Mode mode,
                 std::size_t referenceFrame = 0,
                 std::vector<std::size_t> subsetAtoms = {});

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

    // Compute the transform for `frame` from the inner conformation's
    // RAW positions. The implementation dispatches on mode_.
    Transform3D computeTransform(std::size_t frame) const;

    // Kabsch algorithm — compute the rotation+translation that minimises
    // sum of squared distances between `current` atoms and the reference's
    // same atoms. Both vectors hold the SAME N atoms in the SAME order.
    // Returns identity transform if N < 3 (degenerate / underdetermined).
    static Transform3D KabschFit(const std::vector<Vec3>& current,
                                 const std::vector<Vec3>& reference);

    Conformation* inner_ = nullptr;

    Mode mode_ = Mode::Identity;
    std::size_t referenceFrame_ = 0;
    std::vector<std::size_t> subsetAtoms_;

    // Per-frame cache keyed by frame index. Mutable because atomPosition
    // is const. Cleared on setMode() via generation_ bump.
    mutable std::unordered_map<std::size_t, Transform3D> transformCache_;
    mutable std::uint64_t generation_ = 0;

    // Cached reference positions for FitReference / FitSubset — set once
    // per setMode() so we don't re-read the inner conformation for every
    // frame. The size matches subsetAtoms_ for FitSubset, or the protein
    // atom count for FitReference. Empty for Identity / CenterCom.
    mutable std::vector<Vec3> referencePositions_;
};

}  // namespace h5reader::model
