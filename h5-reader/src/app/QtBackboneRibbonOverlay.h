// QtBackboneRibbonOverlay — backbone ribbon via vtkProteinRibbonFilter.
//
// Adapts ui/src/BackboneRibbonOverlay.cpp for trajectory animation.
// Build() constructs the input polydata once; setFrame(t) updates CA
// positions (and DSSP codes, since the H5 carries per-frame DSSP) and
// reruns the ribbon pipeline.
//
// Per-residue colouring is typed (AminoAcid enum), never string-based.
// Secondary structure comes from dssp/ss8 per frame — the ribbon
// shape and coloring can therefore evolve as helices form or break
// during MD.

#pragma once

#include "../model/QtConformation.h"
#include "../model/QtProtein.h"

#include <QObject>

#include <vtkActor.h>
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkProteinRibbonFilter.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkStringArray.h>
#include <vtkUnsignedCharArray.h>

#include <vector>

namespace h5reader::app {

class QtBackboneRibbonOverlay final : public QObject {
    Q_OBJECT

public:
    explicit QtBackboneRibbonOverlay(
        vtkSmartPointer<vtkRenderer> renderer,
        QObject* parent = nullptr);
    ~QtBackboneRibbonOverlay() override;

    // Build the ribbon pipeline from the protein identity + frame-0
    // positions + DSSP codes. Must be called once before setFrame().
    void Build(const model::QtProtein&      protein,
               const model::QtConformation& conformation);

public slots:
    // Update CA positions (and other backbone atom positions) and DSSP
    // codes, then rerun vtkProteinRibbonFilter. Relatively expensive
    // (spline subdivision) but well within the frame budget at typical
    // protein sizes (~2-5 ms on 300-atom proteins).
    void setFrame(int t);

    // Show/hide the ribbon actor.
    void setVisible(bool visible);

private:
    void UpdateInputArrays(int t);
    void ApplyResidueColors();

    vtkSmartPointer<vtkRenderer>                  renderer_;
    vtkSmartPointer<vtkPoints>                    points_;
    vtkSmartPointer<vtkStringArray>               atomTypes_;
    vtkSmartPointer<vtkIdTypeArray>               atomType_;
    vtkSmartPointer<vtkIdTypeArray>               residue_;
    vtkSmartPointer<vtkUnsignedCharArray>         chain_;
    vtkSmartPointer<vtkUnsignedCharArray>         ss_;
    vtkSmartPointer<vtkUnsignedCharArray>         ssBegin_;
    vtkSmartPointer<vtkUnsignedCharArray>         ssEnd_;
    vtkSmartPointer<vtkUnsignedCharArray>         ishetatm_;
    vtkSmartPointer<vtkPolyData>                  inputPd_;
    vtkSmartPointer<vtkProteinRibbonFilter>       ribbon_;
    vtkSmartPointer<vtkActor>                     actor_;

    const model::QtProtein*      protein_      = nullptr;
    const model::QtConformation* conformation_ = nullptr;
    int                          subdivideFactor_ = 20;
    bool                         visible_      = true;

    // Cached per-segment residue lists for coloring — rebuilt each
    // setFrame because SS (and therefore segment boundaries) may change.
    struct Segment { std::vector<model::AminoAcid> residues; };
    std::vector<Segment> segments_;
};

}  // namespace h5reader::app
