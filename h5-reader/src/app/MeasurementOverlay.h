// MeasurementOverlay marks an ordered AtomSelection with up to four coloured
// spheres and a connecting polyline. Positions are read again on every frame,
// so the geometry remains attached to the atoms during playback. Numeric values
// are shown in the measurement dock rather than over the molecule.
//
// MoleculeScene owns the overlay. Build() runs once, setFrame() updates it, and
// the scene performs rendering after all overlays have advanced. VTK state is
// changed only on the GUI thread.

#pragma once

#include "../model/Conformation.h"
#include "../model/QtProtein.h"

#include <QObject>
#include <QPointer>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>

#include <array>
#include <cstddef>

namespace h5reader::model {
class AtomSelection;
}

namespace h5reader::app {

class MeasurementOverlay final : public QObject {
    Q_OBJECT
public:
    explicit MeasurementOverlay(
        vtkSmartPointer<vtkRenderer>                  renderer,
        vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
        QObject*                                      parent = nullptr);
    ~MeasurementOverlay() override;

    void Build(const model::QtProtein& protein,
               model::Conformation&    conformation);

    // Bind the selection whose atoms this overlay marks. Wired by
    // ReaderMainWindow once the scene + selection both exist.
    void setSelection(model::AtomSelection* selection);

    bool isVisible() const { return visible_; }

public:
    // Selection membership changed — re-evaluate which spheres show, at the
    // current frame. No Render (the scene refreshes the frame right after).
    void onSelectionChanged();

    // Reposition the visible spheres to frame t's atom positions.
    void setFrame(int t);

    // Show or hide the whole overlay.
    void setVisible(bool on);

    // Instrument mode — a marker preset for the harness work. When on,
    // sphere colours switch to a CPK-distinct fixed table (magenta, spring
    // green, deep pink, vivid violet), opacity goes to 1.0, and radius to
    // 1.5 Å. Designed so a Python harness can locate the marker via
    // connected-component blob analysis on a snapshot PNG: the colours are
    // outside every CPK element colour, so a hue threshold finds the marker
    // without confusing it with the molecule. When off, restore Okabe-Ito
    // colours + the default opacity/radius. No Render — the caller (REST
    // handler today) issues scene_->requestRender() after, matching the
    // setVisible() flow at ReaderMainWindow.cpp:591.
    //
    // `focusOnly` (default false): when true AND `on` is true, all four
    // sphere actors get the slot-0 magenta colour and ONLY the focus-slot
    // sphere is rendered (the others are SetVisibility(0)), preventing another
    // marker from occluding the focus marker. Disabling instrument mode
    // (on=false) restores the multi-slot view regardless of focusOnly.
    void setInstrumentMode(bool on, bool focusOnly = false);

private:
    void applyFrame(int t);
    void applyInstrumentModeToActors();

    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow_;

    static constexpr std::size_t kMaxSpheres = 4;  // == AtomSelection::kMaxAtoms
    std::array<vtkSmartPointer<vtkSphereSource>, kMaxSpheres> spheres_;
    std::array<vtkSmartPointer<vtkActor>,        kMaxSpheres> actors_;

    // Connecting polyline through the ordered selected atoms (slot order).
    // Lines only — the measured value lives in the strip chart readout.
    vtkSmartPointer<vtkPolyData>  lineData_;
    vtkSmartPointer<vtkPoints>    linePoints_;
    vtkSmartPointer<vtkCellArray> lineCells_;
    vtkSmartPointer<vtkActor>     lineActor_;

    const model::QtProtein*         protein_ = nullptr;
    QPointer<model::Conformation>   conformation_;
    QPointer<model::AtomSelection>  selection_;
    bool visible_        = true;
    int  lastFrame_      = 0;
    bool instrumentMode_ = false;
    bool focusOnly_      = false;
};

}  // namespace h5reader::app
