// QtAtomTrajectoryOverlay -- selected-atom local trajectory envelope.
//
// Replacement for the old whole-trajectory Occupancy UI. For the focused atom,
// sample the selected atom across the displayed trajectory, estimate an
// occupation density with the tested OccupancyShellsMath KDE/HDR core, and draw
// two translucent nested shells. Changing frames moves the atom through the
// same envelope; it does not rebuild a different local shape. This keeps the
// useful visual form -- an implied volume/arc of motion -- without claiming
// frame-to-frame interpolation or physics that we do not yet have.
//
// ORCA shielding samples are .LGS-backed and may be logged for the current
// frame, but they do not define the envelope geometry and do not live in
// trajectory.h5.

#pragma once

#include "../model/Conformation.h"
#include "../model/DftShieldingStore.h"
#include "../model/QtProtein.h"
#include "OccupancyShellsMath.h"

#include <QObject>
#include <QPointer>

#include <vtkActor.h>
#include <vtkContourFilter.h>
#include <vtkDoubleArray.h>
#include <vtkImageData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkTrivialProducer.h>

#include <array>
#include <cstddef>
#include <optional>
#include <unordered_map>

namespace h5reader::model {
class AtomSelection;
}

namespace h5reader::app {

class QtAtomTrajectoryOverlay final : public QObject {
    Q_OBJECT

public:
    explicit QtAtomTrajectoryOverlay(vtkSmartPointer<vtkRenderer> renderer,
                                     QObject* parent = nullptr);
    ~QtAtomTrajectoryOverlay() override;

    void Build(const model::QtProtein& protein,
               model::Conformation& conformation);

    void setSelection(model::AtomSelection* selection);
    void setDftStore(model::DftShieldingStore* store);

public slots:
    void setFrame(int frame);
    void onFocusChanged(std::size_t atomIdx);
    void onSelectionCleared();
    void onTransformChanged();
    void setVisible(bool on);

signals:
    void rebuildStarted(int frameCount);
    void rebuildFinished(int frameCount, int dftSamples, int loadMs);

private:
    struct Shell {
        vtkSmartPointer<vtkContourFilter> contour;
        vtkSmartPointer<vtkActor> actor;
        double fraction = 0.0;
        double opacity = 0.3;
    };

    void rebuild();
    void hideEnvelope();
    void applyShellStyling(std::optional<double> trendDelta, double trendScale);
    std::optional<double> sampleOrcaT0(std::size_t frame, std::size_t atom);
    void clearScalarCacheForAtom(std::size_t atom);

    vtkSmartPointer<vtkRenderer> renderer_;
    vtkSmartPointer<vtkImageData> imageData_;
    vtkSmartPointer<vtkDoubleArray> scalars_;
    vtkSmartPointer<vtkTrivialProducer> producer_;
    std::array<Shell, 2> shells_;

    const model::QtProtein* protein_ = nullptr;
    QPointer<model::Conformation> conformation_;
    QPointer<model::AtomSelection> selection_;
    QPointer<model::DftShieldingStore> dftStore_;
    std::optional<std::size_t> cachedAtom_;
    std::unordered_map<std::size_t, std::optional<double>> orcaT0ByOriginal_;

    math::OccupancyConfig cfg_;
    int currentFrame_ = 0;
    int windowStart_ = -1;
    int windowEnd_ = -1;
    bool visible_ = false;
    bool dirty_ = false;
    bool hasEnvelope_ = false;
};

}  // namespace h5reader::app
