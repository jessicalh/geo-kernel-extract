// QtFieldGridOverlay — volumetric BS/HM butterfly isosurfaces.
//
// For each aromatic ring we sample the closed-form Biot-Savart (and/or
// Haigh-Mallion) kernel on a structured 3D grid centred on the ring,
// feed the scalar field into vtkImageData, and run vtkContourFilter to
// extract isosurfaces at ±threshold. The "butterfly" pattern above/
// below the ring is the iconic chemist picture of ring-current shielding.
//
// Per-frame update: on setFrame(t) the ring geometry (center, normal,
// radius) and vertex positions come from the H5 slab and the atom
// position slab, kernel re-eval rebuilds scalar arrays, contour
// filters rerun on Modified(). Expensive but bounded — at 3 rings ×
// 20³ grid it's ~10-30 ms per frame on our target hardware.
//
// Kernel evaluators are free functions in h5reader::calculators —
// thread-safe today even though we call them on the GUI thread. See
// the threading discussion in memory project_viewer_hardwon_lessons.

#pragma once

#include "../model/QtConformation.h"
#include "../model/QtProtein.h"

#include <QObject>

#include <vtkActor.h>
#include <vtkContourFilter.h>
#include <vtkFloatArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkImageData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>

#include <vector>

namespace h5reader::app {

enum class FieldGridMode {
    BiotSavart    = 0,
    HaighMallion  = 1,
    Sum           = 2,   // BS + HM — redundant by design (cos~0.999)
};

class QtFieldGridOverlay final : public QObject {
    Q_OBJECT

public:
    explicit QtFieldGridOverlay(
        vtkSmartPointer<vtkRenderer>                  renderer,
        vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
        QObject* parent = nullptr);
    ~QtFieldGridOverlay() override;

    void Build(const model::QtProtein&      protein,
               const model::QtConformation& conformation);

public slots:
    // Recompute per-ring scalar grids from kernel re-eval at frame t,
    // then rerun the contour filters. One Render() is issued by the
    // MoleculeScene owner after all overlays have updated.
    void setFrame(int t);

    void setMode(FieldGridMode mode);
    void setThresholdPpm(double threshold);
    void setOpacity(double opacity);

    // Master visibility; also split shielded / deshielded toggles.
    void setVisible(bool visible);
    void setShieldedVisible(bool visible);
    void setDeshieldedVisible(bool visible);

private:
    struct RingGrid {
        vtkSmartPointer<vtkImageData>      imageData;
        vtkSmartPointer<vtkFloatArray>     scalars;
        vtkSmartPointer<vtkTrivialProducer> producer;
        vtkSmartPointer<vtkContourFilter>  contourShielded;
        vtkSmartPointer<vtkContourFilter>  contourDeshielded;
        vtkSmartPointer<vtkActor>          actorShielded;    // sky blue, T0 < -threshold
        vtkSmartPointer<vtkActor>          actorDeshielded;  // coral,    T0 > +threshold
    };

    // Rebuild the scalar field for one ring at the given frame. Uses
    // the active FieldGridMode to decide BS / HM / Sum.
    void RecomputeRingScalars(size_t ringIdx, int t);

    // Apply current threshold to each contour filter (Modified()).
    void UpdateThresholds();

    void ApplyActorStyling(RingGrid& rg);

    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    const model::QtProtein*                       protein_      = nullptr;
    const model::QtConformation*                  conformation_ = nullptr;
    std::vector<RingGrid>                         rings_;

    FieldGridMode mode_              = FieldGridMode::BiotSavart;
    double        thresholdPpm_       = 0.10;    // ppm
    double        opacity_            = 0.40;
    bool          visible_            = false;   // off by default — user enables
    bool          shieldedVisible_    = true;
    bool          deshieldedVisible_  = true;
};

}  // namespace h5reader::app
