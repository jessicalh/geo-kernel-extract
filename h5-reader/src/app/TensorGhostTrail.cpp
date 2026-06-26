#include "TensorGhostTrail.h"

#include "TensorGlyphActor.h"

#include <utility>

namespace h5reader::app {

TensorGhostTrail::TensorGhostTrail(vtkSmartPointer<vtkRenderer> sceneRenderer)
    : renderer_(std::move(sceneRenderer)) {}

// Out-of-line so the unique_ptr<TensorGlyphActor> elements are destroyed here,
// where TensorGlyphActor is a complete type (it is forward-declared in the .h).
TensorGhostTrail::~TensorGhostTrail() = default;

void TensorGhostTrail::show(const std::vector<Sample>& samples) {
    // Grow the actor pool to cover the request -- each TensorGlyphActor owns its
    // own VTK pipeline, so there is one per ghost. The pool only ever grows;
    // surplus actors are hidden (cleared) below and reused on the next show().
    while (glyphs_.size() < samples.size())
        glyphs_.push_back(std::make_unique<TensorGlyphActor>(renderer_));

    for (std::size_t i = 0; i < samples.size(); ++i) {
        const Sample& s = samples[i];
        glyphs_[i]->show(s.center, s.principalValues, s.pasAxes, s.iso, s.opacity, s.style);
    }
    for (std::size_t i = samples.size(); i < glyphs_.size(); ++i)
        glyphs_[i]->clear();

    visible_ = samples.size();
}

void TensorGhostTrail::clear() {
    for (auto& g : glyphs_)
        if (g) g->clear();
    visible_ = 0;
}

}  // namespace h5reader::app
