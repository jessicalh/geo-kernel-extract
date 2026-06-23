// CsaTensorOverlay -- atom-centered superquadric tensor glyph for the viewer's
// standalone CSA view. Follows the Schultz & Kindlmann superquadric tensor-glyph
// best practice (IEEE TVCG 2010): the glyph's SHAPE is sharpened by the
// deviatoric anisotropy so axial (cigar), planar/rhombic (disc) and isotropic
// (sphere) tensors are distinguishable at a glance -- the disambiguation a plain
// ellipsoid cannot give -- sized by |sigma_i - sigma_iso|, oriented onto the
// principal-axis system, and SIGN-COLOURED per surface point (shielded vs
// deshielded directions). sigma_11/22/33 tip labels and a corner readout
// (iso / span / eta / skew) carry the exact numbers. No arrow thicket: the
// surface carries orientation + magnitude, the way TensorView / Kindlmann
// glyphs do. Visible on a bare extractor run (DFT, no rediscover pass).
//
// The translucent glyph draws on the depth-peeled MAIN renderer, so it composes
// seamlessly with the molecule (order-independent transparency, the same path
// the field-grid isosurfaces and occupancy shells use); the sigma_11/22/33
// labels + corner readout draw on the always-on-top overlay renderer so they
// stay readable regardless of depth. PURE RENDERER: the controller
// (ReaderMainWindow) owns the DFT store + selection + frame, computes the
// CsaShape, and feeds it via show(); the overlay only draws.

#pragma once

#include "../model/CsaShape.h"
#include "../model/Types.h"

#include <QObject>

#include <vtkActor.h>
#include <vtkBillboardTextActor3D.h>
#include <vtkColorTransferFunction.h>
#include <vtkCornerAnnotation.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataNormals.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkSuperquadricSource.h>
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

    // Draw the superquadric tensor glyph + sigma_11/22/33 labels + corner readout
    // at atomPos. Replaces any prior glyph. molecularAxes is accepted for
    // signature stability with the controller but no longer drawn (the surface
    // carries the orientation); pass std::nullopt or anything.
    void show(const model::Vec3& atomPos,
              const model::CsaShape& shape,
              const std::optional<model::Mat3>& molecularAxes);
    void clear();
    void setVisible(bool on);
    bool isActive() const { return active_; }

    // Three glyph styles the controller cycles (re-feeding show() after a change):
    //  * Superquadric (Kindlmann, default) -- shape sharpens by anisotropy type.
    //  * Ovaloid (TensorView, NMR-community "more correct") -- a sphere deformed
    //    radially by r ~ |dev(n)|, sign-coloured, pinching at the zero-crossing.
    //  * Ellipsoid (the popular smooth form) -- a roundness-1 superquadric.
    // Enum ORDER is the cycle order + the REST/int index (0/1/2); do not reorder.
    enum class GlyphStyle { Superquadric, Ovaloid, Ellipsoid };
    void setStyle(GlyphStyle s) { style_ = s; }
    GlyphStyle style() const { return style_; }

private:
    void ensureActors();
    void hideAll();

    vtkSmartPointer<vtkRenderer> renderer_;     // scene renderer (depth-peeled): the glyph
    vtkSmartPointer<vtkRenderer> hudRenderer_;  // overlay renderer: sigma labels + readout

    // Superquadric glyph pipeline: source -> (deep copy + per-vertex sign scalar)
    // -> transform (scale by |dev|, orient onto PAS, translate to atom) -> actor.
    vtkSmartPointer<vtkSuperquadricSource>      glyphSource_;
    vtkSmartPointer<vtkSphereSource>            ovaloidSource_;  // base for the ovaloid
    vtkSmartPointer<vtkPolyData>                glyphLocal_;
    vtkSmartPointer<vtkTransform>               glyphTransform_;
    vtkSmartPointer<vtkTransformPolyDataFilter> glyphFilter_;
    vtkSmartPointer<vtkPolyDataNormals>         glyphNormals_;  // re-derive after deform
    vtkSmartPointer<vtkColorTransferFunction>   glyphLut_;   // diverging sign map
    vtkSmartPointer<vtkActor>                   glyphActor_;

    // sigma_11/22/33 tip labels (billboard = always camera-facing).
    std::array<vtkSmartPointer<vtkBillboardTextActor3D>, 3> axisLabels_;

    // Corner readout: iso / span / eta / skew + the colour key.
    vtkSmartPointer<vtkCornerAnnotation> readout_;

    bool actorsBuilt_ = false;
    bool active_ = false;
    GlyphStyle style_ = GlyphStyle::Superquadric;
};

}  // namespace h5reader::app
