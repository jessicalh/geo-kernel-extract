#include "CsaTensorOverlay.h"

#include "TensorGlyphActor.h"

#include <array>
#include <utility>

namespace h5reader::app {

CsaTensorOverlay::CsaTensorOverlay(vtkSmartPointer<vtkRenderer> sceneRenderer,
                                   vtkSmartPointer<vtkRenderer> /*hudRenderer*/,
                                   QObject* parent)
    : QObject(parent),
      glyph_(std::make_unique<TensorGlyphActor>(std::move(sceneRenderer))) {}

CsaTensorOverlay::~CsaTensorOverlay() = default;

void CsaTensorOverlay::show(const model::Vec3& atomPos,
                            const model::CsaShape& shape,
                            const std::optional<model::Mat3>& /*molecularAxes*/) {
    if (!shape.valid) {
        glyph_->clear();
        return;
    }
    // CSA convention: principal_values ascending sigma11<=sigma22<=sigma33, with
    // pas_axes columns following the same order; iso = sigma_iso. The shared
    // actor draws the index-coloured principal-axis arrows.
    glyph_->show(atomPos,
                 {shape.principal_values[0], shape.principal_values[1], shape.principal_values[2]},
                 shape.pas_axes,
                 shape.sigma_iso);
}

void CsaTensorOverlay::clear() {
    glyph_->clear();
}

void CsaTensorOverlay::setVisible(bool on) {
    glyph_->setVisible(on);
}

bool CsaTensorOverlay::isActive() const {
    return glyph_->isActive();
}

}  // namespace h5reader::app
