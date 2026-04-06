#include "PeptideBondOverlay.h"

#include "Protein.h"
#include "ProteinConformation.h"
#include "Bond.h"

#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkTubeFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkNew.h>
#include <vtkLine.h>

using namespace nmr;

namespace {

// Build one actor from many line segments, all the same color
vtkSmartPointer<vtkActor> makeBatchedTubes(
    const std::vector<std::pair<Vec3,Vec3>>& segments,
    double r, double g, double b,
    double radius, double opacity)
{
    if (segments.empty()) return nullptr;

    vtkNew<vtkPoints> pts;
    vtkNew<vtkCellArray> lines;

    for (const auto& [start, end] : segments) {
        vtkIdType id0 = pts->InsertNextPoint(start.data());
        vtkIdType id1 = pts->InsertNextPoint(end.data());
        vtkNew<vtkLine> line;
        line->GetPointIds()->SetId(0, id0);
        line->GetPointIds()->SetId(1, id1);
        lines->InsertNextCell(line);
    }

    vtkNew<vtkPolyData> pd;
    pd->SetPoints(pts);
    pd->SetLines(lines);

    vtkNew<vtkTubeFilter> tube;
    tube->SetInputData(pd);
    tube->SetRadius(radius);
    tube->SetNumberOfSides(8);

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(tube->GetOutputPort());

    auto actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(r, g, b);
    actor->GetProperty()->SetOpacity(opacity);
    actor->SetForceTranslucent(true);
    actor->GetProperty()->SetBackfaceCulling(true);
    return actor;
}

} // anonymous namespace

PeptideBondOverlay::PeptideBondOverlay(
    vtkSmartPointer<vtkRenderer> renderer,
    const Protein& protein,
    const ProteinConformation& conf)
    : renderer_(renderer)
{
    buildActors(protein, conf);
}

PeptideBondOverlay::~PeptideBondOverlay() {
    for (auto& a : bbActors_) renderer_->RemoveActor(a);
    for (auto& a : scActors_) renderer_->RemoveActor(a);
}

void PeptideBondOverlay::buildActors(
    const Protein& protein,
    const ProteinConformation& conf)
{
    std::vector<std::pair<Vec3,Vec3>> peptideSegments, scSegments, otherSegments;

    for (size_t i = 0; i < protein.BondCount(); ++i) {
        const Bond& bond = protein.BondAt(i);
        Vec3 posA = conf.AtomAt(bond.atom_index_a).Position();
        Vec3 posB = conf.AtomAt(bond.atom_index_b).Position();
        auto seg = std::make_pair(posA, posB);

        if (bond.IsPeptideBond()) {
            peptideSegments.push_back(seg);
        } else if (bond.category == BondCategory::SidechainOther ||
                   bond.category == BondCategory::SidechainCO) {
            scSegments.push_back(seg);
        } else {
            otherSegments.push_back(seg);
        }
    }

    // Peptide bonds — red
    if (auto a = makeBatchedTubes(peptideSegments, 0.9, 0.2, 0.2, tubeRadius_, 0.7)) {
        renderer_->AddActor(a);
        bbActors_.push_back(a);
    }
    // Other backbone bonds — blue
    if (auto a = makeBatchedTubes(otherSegments, 0.2, 0.3, 0.9, tubeRadius_, 0.7)) {
        renderer_->AddActor(a);
        bbActors_.push_back(a);
    }
    // Sidechain bonds — orange
    if (auto a = makeBatchedTubes(scSegments, 0.9, 0.6, 0.1, tubeRadius_ * 0.8, 0.6)) {
        renderer_->AddActor(a);
        scActors_.push_back(a);
    }
}

void PeptideBondOverlay::setVisible(bool visible) {
    for (auto& a : bbActors_) a->SetVisibility(visible ? 1 : 0);
    for (auto& a : scActors_) a->SetVisibility(visible ? 1 : 0);
}

void PeptideBondOverlay::setShowSidechain(bool show) {
    for (auto& a : scActors_) a->SetVisibility(show ? 1 : 0);
}
