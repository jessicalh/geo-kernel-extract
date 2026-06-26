// TensorGhostTrail -- a transient "ghost trail" of a single atom's shielding
// tensor across its last N measured frames, for heroshot figures. Each ghost is
// the SAME shared TensorGlyphActor (ovaloid + principal-axis arrows) the live
// CSA glyph uses, drawn at that frame's REAL atom position with its REAL
// PAS-frame eigendecomposition -- no interpolation, every angle true -- and
// faded by the caller (newest opaque, oldest faint) so the re-orientation reads
// "from the side" as a stack of translucent tensors telling the story over time.
//
// Heroshot-only: owned by the REST layer (RestServer), built and cleared on
// demand by POST /heroshot/ghost_trail and /heroshot/clear. It does NOT live in
// the reader UI and does NOT move playback; the reader's own single live glyph
// is untouched. Pure renderer: the caller supplies each ghost's lab-frame
// eigendecomposition (via ReaderMainWindow::probeAtomCsa) and its fade.

#pragma once

#include "TensorGlyphActor.h"

#include "../model/Types.h"

#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <array>
#include <cstddef>
#include <memory>
#include <vector>

namespace h5reader::app {

class TensorGhostTrail {
public:
    // One ghost: a tensor's lab-frame eigendecomposition at one frame plus the
    // fade applied to it. Fields mirror TensorGlyphActor::show's arguments.
    struct Sample {
        model::Vec3 center = model::Vec3::Zero();      // atom display position that frame
        std::array<double, 3> principalValues{{0.0, 0.0, 0.0}};  // ascending sigma11..33
        model::Mat3 pasAxes = model::Mat3::Identity();  // lab-frame eigenvector columns
        double iso = 0.0;                               // isotropic reference (= trace/3)
        double opacity = 1.0;                           // 1 = newest/opaque, ->0 = oldest/faint
        TensorGlyphActor::Style style;                  // per-ghost figure styling
    };

    explicit TensorGhostTrail(vtkSmartPointer<vtkRenderer> sceneRenderer);
    ~TensorGhostTrail();

    TensorGhostTrail(const TensorGhostTrail&) = delete;
    TensorGhostTrail& operator=(const TensorGhostTrail&) = delete;

    // Replace the trail with these ghosts (index 0 drawn first). Grows the pool
    // of glyph actors to match; any actors beyond the new length are cleared.
    // An empty vector clears the whole trail.
    void show(const std::vector<Sample>& samples);
    void clear();
    std::size_t size() const { return visible_; }

private:
    vtkSmartPointer<vtkRenderer> renderer_;  // scene renderer (depth-peeled), one ref held
    std::vector<std::unique_ptr<TensorGlyphActor>> glyphs_;  // one full glyph per ghost
    std::size_t visible_ = 0;                // currently-shown ghost count
};

}  // namespace h5reader::app
