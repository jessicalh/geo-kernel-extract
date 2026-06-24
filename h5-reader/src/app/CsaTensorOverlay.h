// CsaTensorOverlay -- atom-centered NMR shielding-tensor glyph for the viewer's
// standalone CSA view. ONE mode (focused for submission): the surface is the
// TensorView OVALOID (a unit sphere deformed radially by r ~ |dev(n)|, the
// anisotropy shielding surface; it pinches toward the dev(n)=0 cone so the
// shielded/deshielded lobes separate -- sign by shape AND colour, the
// NMR-community "more correct" form). The descriptive element is the three
// PRINCIPAL AXES (sigma_11/22/33) drawn as double-headed arrows along the PAS
// eigenvectors, index-coloured (amber/teal/violet).
//
// NO in-scene text: billboard labels + a corner readout on the molecule were
// fiddly and cluttered the scene (text overlays on molecules are notoriously
// hard). The numbers (iso/span/skew/eta, the per-axis principal values) and the
// colour key now live in the Atom Info panel (QtAtomInspectorDock); the scene
// stays pure geometry. The translucent surface + opaque arrows draw on the
// depth-peeled MAIN renderer so they compose seamlessly with the molecule (like
// the field-grid isosurfaces). PURE RENDERER: the controller (ReaderMainWindow)
// owns the DFT store + selection + frame, computes the CsaShape, feeds it via
// show(), and mirrors the numbers into the Atom Info panel.

#pragma once

#include "../model/CsaShape.h"
#include "../model/Types.h"

#include <QObject>

#include <vtkActor.h>
#include <vtkArrowSource.h>
#include <vtkColorTransferFunction.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <array>
#include <optional>

namespace h5reader::app {

class CsaTensorOverlay final : public QObject {
    Q_OBJECT
public:
    explicit CsaTensorOverlay(vtkSmartPointer<vtkRenderer> sceneRenderer,
                              vtkSmartPointer<vtkRenderer> hudRenderer,
                              QObject* parent = nullptr);
    ~CsaTensorOverlay() override;

    // Draw the ovaloid surface + sigma_11/22/33 principal-axis arrows at atomPos
    // (PURE GEOMETRY -- no in-scene text). The numeric shape (iso/span/skew/eta)
    // and the colour key live in the Atom Info panel, not on the molecule.
    // Replaces any prior glyph. molecularAxes is accepted for signature stability
    // but no longer drawn (the surface + the PAS arrows carry the orientation).
    void show(const model::Vec3& atomPos,
              const model::CsaShape& shape,
              const std::optional<model::Mat3>& molecularAxes);
    void clear();
    void setVisible(bool on);
    bool isActive() const { return active_; }

private:
    void ensureActors();
    void hideAll();

    vtkSmartPointer<vtkRenderer> renderer_;     // scene renderer (depth-peeled): surface + arrows
    vtkSmartPointer<vtkRenderer> hudRenderer_;  // overlay renderer: sigma labels + readout

    // Ovaloid surface pipeline: sphere -> (deep copy + radial deform + per-vertex
    // sign scalar) -> transform (orient onto PAS, translate) -> normals -> actor.
    vtkSmartPointer<vtkSphereSource>            ovaloidSource_;
    vtkSmartPointer<vtkPolyData>                glyphLocal_;
    vtkSmartPointer<vtkTransform>               glyphTransform_;
    vtkSmartPointer<vtkTransformPolyDataFilter> glyphFilter_;
    vtkSmartPointer<vtkPolyDataNormals>         glyphNormals_;
    vtkSmartPointer<vtkColorTransferFunction>   glyphLut_;   // diverging sign map
    vtkSmartPointer<vtkActor>                   glyphActor_;

    // Principal-axis arrows (the descriptive element): one shared arrow source +
    // mapper, six actors (each PAS axis double-headed, +/- director), index-coloured.
    vtkSmartPointer<vtkArrowSource>          arrowSource_;
    vtkSmartPointer<vtkPolyDataMapper>       arrowMapper_;
    std::array<vtkSmartPointer<vtkActor>, 6> arrowActors_;

    bool actorsBuilt_ = false;
    bool active_ = false;
};

}  // namespace h5reader::app
