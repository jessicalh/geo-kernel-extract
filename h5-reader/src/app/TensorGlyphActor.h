// TensorGlyphActor -- the single, shared scene representation for a symmetric
// 3x3 tensor. ONE form, drawn for EVERY tensor the viewer shows (CSA shielding,
// bond-orientation order, ...) so the representation stays consistent, not ad
// hoc. Usable-mode default is principal-axis arrows only; the older deformed
// ovaloid surface remains an explicit style option for figure experiments.
//
// Pure renderer: the caller supplies the eigendecomposition ALREADY resolved
// into the lab frame -- centre, the three principal values, their eigenvector
// columns, and the isotropic reference (= trace/3). No tensor-specific semantics
// live here (a shielding tensor and an order tensor differ only in what the
// caller computes and feeds). No in-scene text -- the numbers belong in a Qt
// panel. Extracted verbatim from CsaTensorOverlay so the CSA glyph is unchanged.

#pragma once

#include "../model/Types.h"

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

namespace h5reader::app {

class TensorGlyphActor {
public:
    struct Style {
        double ovaloidScale = 1.0;
        double arrowLengthScale = 1.0;
        double arrowWidthScale = 1.0;
        double surfaceOpacity = 0.50;
        double arrowOpacity = 1.0;
        bool showSurface = false;
        bool showArrows = true;
        std::array<bool, 3> showAxes{{true, true, true}};
    };

    explicit TensorGlyphActor(vtkSmartPointer<vtkRenderer> sceneRenderer);
    ~TensorGlyphActor();

    TensorGlyphActor(const TensorGlyphActor&) = delete;
    TensorGlyphActor& operator=(const TensorGlyphActor&) = delete;

    // Draw the principal-axis arrows for the symmetric tensor whose
    // lab-frame eigendecomposition is (principalValues[i] along pasAxes.col(i)),
    // centred at `center`, isotropic reference `iso` (= trace/3). The ovaloid
    // radius along direction n is ~ |dev(n)| where dev = value - iso, so an
    // isotropic tensor collapses toward a point and an anisotropic one bulges
    // along its dominant deviation. Arrows are coloured by principal-value index
    // (0/1/2 -> amber/teal/violet). Replaces any prior glyph.
    // `opacity` (default 1.0) scales the arrows and any explicitly enabled
    // surface so a caller can fade it (e.g. a ghost trail's older frames).
    void show(const model::Vec3& center,
              const std::array<double, 3>& principalValues,
              const model::Mat3& pasAxes,
              double iso,
              double opacity = 1.0);
    void show(const model::Vec3& center,
              const std::array<double, 3>& principalValues,
              const model::Mat3& pasAxes,
              double iso,
              double opacity,
              const Style& style);
    void clear();
    void setVisible(bool on);
    bool isActive() const { return active_; }

private:
    void ensureActors();
    void hideAll();

    vtkSmartPointer<vtkRenderer> renderer_;  // scene renderer (depth-peeled)

    // Ovaloid surface pipeline: sphere -> (deep copy + radial deform + per-vertex
    // sign scalar) -> transform (orient onto PAS, translate) -> normals -> actor.
    vtkSmartPointer<vtkSphereSource>            ovaloidSource_;
    vtkSmartPointer<vtkPolyData>                glyphLocal_;
    vtkSmartPointer<vtkTransform>               glyphTransform_;
    vtkSmartPointer<vtkTransformPolyDataFilter> glyphFilter_;
    vtkSmartPointer<vtkPolyDataNormals>         glyphNormals_;
    vtkSmartPointer<vtkColorTransferFunction>   glyphLut_;   // diverging sign map
    vtkSmartPointer<vtkActor>                   glyphActor_;

    // Principal-axis arrows: one shared arrow source + mapper, six actors (each
    // PAS axis double-headed, +/- director), index-coloured.
    vtkSmartPointer<vtkArrowSource>          arrowSource_;
    vtkSmartPointer<vtkPolyDataMapper>       arrowMapper_;
    std::array<vtkSmartPointer<vtkActor>, 6> arrowActors_;
    std::array<bool, 6> arrowSlotVisible_{{false, false, false, false, false, false}};
    bool surfaceVisible_ = false;

    bool actorsBuilt_ = false;
    bool active_ = false;
};

}  // namespace h5reader::app
