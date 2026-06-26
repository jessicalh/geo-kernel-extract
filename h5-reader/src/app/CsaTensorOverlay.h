// CsaTensorOverlay -- atom-centered NMR shielding-tensor glyph for the viewer's
// CSA view. A thin ADAPTER: it computes nothing about how a tensor is drawn --
// that lives in the shared TensorGlyphActor, the ONE representation every
// tensor in the scene uses. This class only feeds the shielding tensor's
// eigendecomposition (CsaShape: principal values, PAS axes, sigma_iso) to that
// actor, so CSA and the bond-orientation tensor render identically --
// consistent, not ad hoc.
//
// NO in-scene text: the numbers (iso/span/skew/eta, per-axis principal values)
// and the colour key live in the Atom Info panel (QtAtomInspectorDock); the
// scene stays pure geometry. PURE RENDERER: the controller (ReaderMainWindow)
// owns the DFT store + selection + frame, computes the CsaShape, feeds it via
// show(), and mirrors the numbers into the Atom Info panel.

#pragma once

#include "../model/CsaShape.h"
#include "../model/Types.h"

#include <QObject>

#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <memory>
#include <optional>

namespace h5reader::app {

class TensorGlyphActor;

class CsaTensorOverlay final : public QObject {
    Q_OBJECT
public:
    explicit CsaTensorOverlay(vtkSmartPointer<vtkRenderer> sceneRenderer,
                              vtkSmartPointer<vtkRenderer> hudRenderer,
                              QObject* parent = nullptr);
    ~CsaTensorOverlay() override;

    // Draw the shielding tensor's sigma_11/22/33 principal-axis arrows at
    // atomPos via the shared TensorGlyphActor. molecularAxes is accepted for
    // signature stability but unused (the PAS axes carry the orientation).
    void show(const model::Vec3& atomPos,
              const model::CsaShape& shape,
              const std::optional<model::Mat3>& molecularAxes);
    void clear();
    void setVisible(bool on);
    bool isActive() const;

private:
    std::unique_ptr<TensorGlyphActor> glyph_;
};

}  // namespace h5reader::app
