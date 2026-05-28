// SceneRevealOverlay -- transient molecule highlight requested by dashboard strips.
//
// This overlay is intentionally separate from MeasurementOverlay. Measurement
// follows AtomSelection; reveal follows pinned strip bindings and must not
// mutate the user's current atom selection. MoleculeScene owns it and fans
// setFrame() through it like the other overlays.

#pragma once

#include "../model/Conformation.h"
#include "../model/QtProtein.h"
#include "../model/SignalDictionary.h"

#include <QObject>
#include <QPointer>

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>

#include <cstddef>
#include <vector>

namespace h5reader::app {

class SceneRevealOverlay final : public QObject {
    Q_OBJECT
public:
    explicit SceneRevealOverlay(vtkSmartPointer<vtkRenderer> renderer,
                                QObject* parent = nullptr);
    ~SceneRevealOverlay() override;

    void Build(const model::QtProtein& protein,
               model::Conformation& conformation);

    const std::vector<std::size_t>& activeAtoms() const { return activeAtoms_; }
    bool isActive() const { return active_; }

public slots:
    void reveal(const model::SignalBinding& binding, int frame);
    void clear();
    void setFrame(int t);

private:
    void ensureSphereCount(std::size_t count);
    std::vector<std::size_t> atomsForBinding(const model::SignalBinding& binding) const;
    void applyFrame(int t);

    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkSphereSource>> spheres_;
    std::vector<vtkSmartPointer<vtkActor>> sphereActors_;

    vtkSmartPointer<vtkPolyData> lineData_;
    vtkSmartPointer<vtkPoints> linePoints_;
    vtkSmartPointer<vtkCellArray> lineCells_;
    vtkSmartPointer<vtkActor> lineActor_;

    const model::QtProtein* protein_ = nullptr;
    QPointer<model::Conformation> conformation_;
    std::vector<std::size_t> activeAtoms_;
    std::vector<std::size_t> lineAtoms_;
    bool active_ = false;
    int lastFrame_ = 0;
};

}  // namespace h5reader::app
