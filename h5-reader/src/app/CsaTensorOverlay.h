// CsaTensorOverlay -- atom-centered 3-D glyph for the viewer's standalone
// tensor view. Draws sigma's PAS ellipsoid (oriented by its OWN eigenvectors,
// sized by the principal shielding deviations) together with the chemical
// molecular-frame axes, so the directional story the stats report ("sigma_33
// along C=O", "ring current loads the ring-normal axis") is visible on a bare
// extractor run -- no rediscover pass required.
//
// A scene overlay on the layer-1 overlay renderer (always legible, like
// SceneRevealOverlay). PURE RENDERER: the controller (ReaderMainWindow) owns
// the DFT store + selection + frame, computes the CsaShape and molecular axes,
// and feeds them via show(); the overlay only draws. No DFT / model-load
// coupling here -- mirrors how SceneRevealOverlay is fed by its controller.

#pragma once

#include "../model/CsaShape.h"
#include "../model/Types.h"

#include <QObject>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkUnsignedCharArray.h>

#include <optional>

namespace h5reader::app {

class CsaTensorOverlay final : public QObject {
    Q_OBJECT
public:
    explicit CsaTensorOverlay(vtkSmartPointer<vtkRenderer> renderer, QObject* parent = nullptr);
    ~CsaTensorOverlay() override;

    // Draw sigma's PAS ellipsoid at atomPos, plus (when present) the molecular
    // frame axes. Replaces any prior glyph. molecularAxes columns are the x, y,
    // z directors in lab coordinates; std::nullopt for an unframed atom (the
    // ellipsoid still shows the intrinsic PAS orientation).
    void show(const model::Vec3& atomPos,
              const model::CsaShape& shape,
              const std::optional<model::Mat3>& molecularAxes);
    void clear();
    void setVisible(bool on);
    bool isActive() const { return active_; }

private:
    void ensureActors();

    vtkSmartPointer<vtkRenderer> renderer_;

    // PAS ellipsoid (unit sphere -> scaled+rotated+translated by a vtkTransform).
    vtkSmartPointer<vtkSphereSource>            ellipsoidSource_;
    vtkSmartPointer<vtkTransform>               ellipsoidTransform_;
    vtkSmartPointer<vtkTransformPolyDataFilter> ellipsoidFilter_;
    vtkSmartPointer<vtkActor>                   ellipsoidActor_;

    // Molecular-frame axes: three colored line segments (x=R, y=G, z=B).
    vtkSmartPointer<vtkPolyData>          axesData_;
    vtkSmartPointer<vtkPoints>            axesPoints_;
    vtkSmartPointer<vtkCellArray>         axesLines_;
    vtkSmartPointer<vtkUnsignedCharArray> axesColors_;
    vtkSmartPointer<vtkActor>             axesActor_;

    bool actorsBuilt_ = false;
    bool active_ = false;
    bool hasAxes_ = false;
};

}  // namespace h5reader::app
