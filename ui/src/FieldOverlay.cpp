#include "FieldOverlay.h"
#include "Colormap.h"
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkArrowSource.h>
#include <vtkGlyph3D.h>
#include <vtkTrivialProducer.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkNew.h>
#include <cmath>
#include <algorithm>

FieldOverlay::FieldOverlay(vtkSmartPointer<vtkRenderer> renderer)
    : renderer_(renderer)
{
}

FieldOverlay::~FieldOverlay() {
    clear();
}

void FieldOverlay::setScalarField(
    const std::vector<nmr::Vec3>& positions,
    const std::vector<double>& values,
    double sphereRadius,
    double opacity)
{
    clearScalar();
    if (positions.empty()) return;

    vtkNew<vtkPoints> points;
    vtkNew<vtkFloatArray> scalars;
    scalars->SetName("FieldValue");
    scalars->SetNumberOfComponents(1);

    double vmin = 1e30, vmax = -1e30;
    for (size_t i = 0; i < positions.size(); ++i) {
        points->InsertNextPoint(positions[i].data());
        float val = (float)values[i];
        scalars->InsertNextValue(val);
        if (val < vmin) vmin = val;
        if (val > vmax) vmax = val;
    }

    vtkNew<vtkPolyData> pd;
    pd->SetPoints(points);
    pd->GetPointData()->SetScalars(scalars);

    // Wrap in TrivialProducer for proper pipeline connection
    vtkNew<vtkTrivialProducer> producer;
    producer->SetOutput(pd);

    vtkNew<vtkSphereSource> sphere;
    sphere->SetRadius(sphereRadius);
    sphere->SetThetaResolution(8);
    sphere->SetPhiResolution(8);

    vtkNew<vtkGlyph3D> glyph;
    glyph->SetInputConnection(producer->GetOutputPort());
    glyph->SetSourceConnection(sphere->GetOutputPort());
    glyph->SetScaleModeToDataScalingOff();
    glyph->SetColorModeToColorByScalar();
    glyph->Update();

    auto ctf = makeDivergingCTF(vmin, vmax);
    double absMax = std::max(std::abs(vmin), std::abs(vmax));
    if (absMax < 1e-10) absMax = 1.0;

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(glyph->GetOutputPort());
    mapper->SetLookupTable(ctf);
    mapper->SetScalarRange(-absMax, absMax);
    mapper->SetScalarModeToUsePointData();

    scalarActor_ = vtkSmartPointer<vtkActor>::New();
    scalarActor_->SetMapper(mapper);
    scalarActor_->GetProperty()->SetOpacity(opacity);
    scalarActor_->SetForceTranslucent(true);
    renderer_->AddActor(scalarActor_);

    // Scalar bar (color legend)
    auto bar = vtkSmartPointer<vtkScalarBarActor>::New();
    bar->SetLookupTable(ctf);
    bar->SetNumberOfLabels(5);
    bar->SetWidth(0.08);
    bar->SetHeight(0.4);
    bar->SetPosition(0.02, 0.3);
    bar->GetLabelTextProperty()->SetFontSize(10);
    bar->GetLabelTextProperty()->SetColor(0.9, 0.9, 0.9);
    bar->GetTitleTextProperty()->SetColor(0.9, 0.9, 0.9);
    scalarBarActor_ = bar;
    renderer_->AddActor2D(bar);
}

void FieldOverlay::setArrows(
    const std::vector<nmr::Vec3>& positions,
    const std::vector<nmr::Vec3>& vectors,
    double scale,
    double opacity)
{
    clearArrows();
    if (positions.empty()) return;

    vtkNew<vtkPoints> points;
    vtkNew<vtkFloatArray> vecs;
    vecs->SetName("EField");
    vecs->SetNumberOfComponents(3);

    for (size_t i = 0; i < positions.size(); ++i) {
        double mag = vectors[i].norm();
        if (mag < 1e-10) continue;
        points->InsertNextPoint(positions[i].data());
        float v[3] = {(float)vectors[i](0), (float)vectors[i](1), (float)vectors[i](2)};
        vecs->InsertNextTuple(v);
    }

    if (points->GetNumberOfPoints() == 0) return;

    vtkNew<vtkPolyData> pd;
    pd->SetPoints(points);
    pd->GetPointData()->SetVectors(vecs);

    vtkNew<vtkTrivialProducer> producer;
    producer->SetOutput(pd);

    vtkNew<vtkArrowSource> arrow;
    arrow->SetTipResolution(8);
    arrow->SetShaftResolution(8);

    vtkNew<vtkGlyph3D> glyph;
    glyph->SetInputConnection(producer->GetOutputPort());
    glyph->SetSourceConnection(arrow->GetOutputPort());
    glyph->SetVectorModeToUseVector();
    glyph->SetScaleModeToDataScalingOff();  // fixed size — vectors are pre-normalized
    glyph->SetScaleFactor(scale);
    glyph->OrientOn();
    glyph->Update();

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(glyph->GetOutputPort());
    mapper->ScalarVisibilityOff();

    arrowActor_ = vtkSmartPointer<vtkActor>::New();
    arrowActor_->SetMapper(mapper);
    arrowActor_->GetProperty()->SetColor(0.95, 0.6, 0.1);
    arrowActor_->GetProperty()->SetOpacity(opacity);
    arrowActor_->SetForceTranslucent(true);
    renderer_->AddActor(arrowActor_);
}

void FieldOverlay::clearScalar() {
    if (scalarActor_) {
        renderer_->RemoveActor(scalarActor_);
        scalarActor_ = nullptr;
    }
    if (scalarBarActor_) {
        renderer_->RemoveActor2D(scalarBarActor_);
        scalarBarActor_ = nullptr;
    }
}

void FieldOverlay::clearArrows() {
    if (arrowActor_) {
        renderer_->RemoveActor(arrowActor_);
        arrowActor_ = nullptr;
    }
}

void FieldOverlay::clear() {
    clearScalar();
    clearArrows();
}

void FieldOverlay::setScalarVisible(bool visible) {
    if (scalarActor_) scalarActor_->SetVisibility(visible ? 1 : 0);
    if (scalarBarActor_) scalarBarActor_->SetVisibility(visible ? 1 : 0);
}

void FieldOverlay::setArrowsVisible(bool visible) {
    if (arrowActor_) arrowActor_->SetVisibility(visible ? 1 : 0);
}
