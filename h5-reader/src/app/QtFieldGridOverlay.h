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

#include "../model/Conformation.h"
#include "../model/QtProtein.h"

#include <QObject>
#include <QPointer>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkContourFilter.h>
#include <vtkFloatArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkImageData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>

#include <optional>
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

    void Build(const model::QtProtein& protein,
               model::Conformation&    conformation);

    size_t ringCount() const { return rings_.size(); }
    std::optional<size_t> visibleRing() const { return visibleRing_; }

public slots:
    // Recompute per-ring scalar grids from kernel re-eval at frame t,
    // then rerun the contour filters. One Render() is issued by the
    // MoleculeScene owner after all overlays have updated.
    void setFrame(int t);

    void setMode(FieldGridMode mode);
    void setThresholdPpm(double threshold);
    void setOpacity(double opacity);

    // Two orthogonal butterfly knobs (re-evaluate the scalar grid).
    void setGaussianExtent(double sigmaA);   // radial taper σ in Å — lobe reach
    void setGaussianPeak(double amplitude);  // amplitude gain (1.0 = true physics)

    // Master visibility; also split shielded / deshielded toggles.
    void setVisible(bool visible);
    void setShieldedVisible(bool visible);
    void setDeshieldedVisible(bool visible);
    void setNullConeVisible(bool visible);
    void setNullConeOpacity(double opacity);
    void setNullConeLength(double lengthA);
    void setVisibleRing(std::optional<size_t> ringIdx);

private:
    struct RingGrid {
        vtkSmartPointer<vtkImageData>      imageData;
        vtkSmartPointer<vtkFloatArray>     scalars;
        vtkSmartPointer<vtkTrivialProducer> producer;
        vtkSmartPointer<vtkContourFilter>  contourShielded;
        vtkSmartPointer<vtkContourFilter>  contourDeshielded;
        vtkSmartPointer<vtkActor>          actorShielded;    // sky blue, T0 < -threshold
        vtkSmartPointer<vtkActor>          actorDeshielded;  // coral,    T0 > +threshold
        vtkSmartPointer<vtkPoints>         nullConePoints;
        vtkSmartPointer<vtkCellArray>      nullConePolys;
        vtkSmartPointer<vtkPolyData>       nullConePoly;
        vtkSmartPointer<vtkActor>          actorNullCone;
    };

    // Rebuild the scalar field for one ring at the given frame. Uses
    // the active FieldGridMode to decide BS / HM / Sum.
    void RecomputeRingScalars(size_t ringIdx, int t);

    // Apply current threshold to each contour filter (Modified()).
    void UpdateThresholds();

    bool RingIsAllowed(size_t ringIdx) const;
    void RecomputeVisibleRings();
    void UpdateNullConeGeometry(size_t ringIdx, RingGrid& rg, const model::RingGeometry& geo);
    void UpdateVisibleNullCones();
    void ApplyActorStyling(size_t ringIdx, RingGrid& rg);

    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;
    const model::QtProtein*                       protein_      = nullptr;
    QPointer<model::Conformation>                 conformation_;
    std::vector<RingGrid>                         rings_;

    FieldGridMode mode_              = FieldGridMode::BiotSavart;
    int           currentFrame_      = 0;
    double        thresholdPpm_       = 0.30;    // ppm — dominant-zone default
                                                 // (~6–7.5 Å reach; fits the
                                                 // 9 Å grid). Runtime-tunable
                                                 // via setThresholdPpm().
    // Physics-preserving Gaussian radial taper on the kernel field:
    //   T0_shown(p) = gaussianPeak_ · T0_kernel(p) · exp(-r²/(2·σ²))
    // r = |p - ring centre|. σ (gaussianExtentA_) bounds the lobe REACH by
    // smoothly killing the unbounded 1/r³ tail, so the isosurface always
    // closes; gaussianPeak_ is an amplitude gain (1.0 keeps the true near-ring
    // ppm). Both runtime-tunable (POST /field/gaussian) — the two orthogonal
    // "extent" and "peak" controls.
    double        gaussianExtentA_    = 9.0;    // σ, Å — reach (taper width)
    double        gaussianPeak_       = 1.0;    // amplitude gain (1 = physics)
    double        opacity_            = 0.40;
    double        nullConeOpacity_    = 0.12;
    double        nullConeLengthA_    = 4.5;
    std::optional<size_t> visibleRing_;         // nullopt = all rings
    bool          visible_            = false;   // off by default — user enables
    bool          nullConeVisible_    = false;
    bool          shieldedVisible_    = true;
    bool          deshieldedVisible_  = true;
};

}  // namespace h5reader::app
