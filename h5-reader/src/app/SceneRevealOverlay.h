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
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <array>
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

    // L-3a (2026-05-29): per-bond-vector Mat3 ellipsoid glyph.
    // `tensor` is row-major 9 doubles (the bond_orientation_tensor row
    // from QtReorientationalDynamics). The ellipsoid is positioned at
    // the bond midpoint each frame via applyFrame; principal axes +
    // radii are computed by TensorGlyphMath::decomposeSymmetric3x3.
    // Pre-condition: Build() has been called for the active protein +
    // conformation. The caller (controller / main window) is
    // responsible for resolving the tail/head atom pair from the
    // BondVectorAnchor; this method takes the resolved atoms directly
    // so the overlay stays decoupled from the iRED/Reorient H5 lookup
    // path that lookupBondVector already centralises.
    void revealTensor(std::size_t tailAtom,
                      std::size_t headAtom,
                      const std::array<double, 9>& tensor,
                      int frame);
    void clearTensor();

private:
    void ensureSphereCount(std::size_t count);
    std::vector<std::size_t> atomsForBinding(const model::SignalBinding& binding) const;
    void applyFrame(int t);
    void applyTensorFrame(int t);

    vtkSmartPointer<vtkRenderer> renderer_;
    std::vector<vtkSmartPointer<vtkSphereSource>> spheres_;
    std::vector<vtkSmartPointer<vtkActor>> sphereActors_;

    vtkSmartPointer<vtkPolyData> lineData_;
    vtkSmartPointer<vtkPoints> linePoints_;
    vtkSmartPointer<vtkCellArray> lineCells_;
    vtkSmartPointer<vtkActor> lineActor_;

    // L-3a tensor glyph state.
    vtkSmartPointer<vtkSphereSource>           tensorSphere_;
    vtkSmartPointer<vtkTransform>              tensorTransform_;
    vtkSmartPointer<vtkTransformPolyDataFilter> tensorFilter_;
    vtkSmartPointer<vtkActor>                  tensorActor_;
    bool tensorActive_ = false;
    std::size_t tensorTail_ = 0;
    std::size_t tensorHead_ = 0;
    std::array<double, 9> tensorData_{};

    const model::QtProtein* protein_ = nullptr;
    QPointer<model::Conformation> conformation_;
    std::vector<std::size_t> activeAtoms_;
    std::vector<std::size_t> lineAtoms_;
    bool active_ = false;
    int lastFrame_ = 0;
};

}  // namespace h5reader::app
