#include "EllipsoidGlyph.h"
#include "Colormap.h"
#include <vtkPoints.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkTensorGlyph.h>
#include <vtkTrivialProducer.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkScalarBarActor.h>
#include <vtkTextProperty.h>
#include <vtkNew.h>
#include <cmath>
#include <algorithm>

EllipsoidGlyph::EllipsoidGlyph(vtkSmartPointer<vtkRenderer> renderer)
    : renderer_(renderer)
{
}

EllipsoidGlyph::~EllipsoidGlyph() {
    clear();
}

void EllipsoidGlyph::clear() {
    if (actor_) {
        renderer_->RemoveActor(actor_);
        actor_ = nullptr;
    }
    if (barActor_) {
        renderer_->RemoveActor2D(barActor_);
        barActor_ = nullptr;
    }
}

void EllipsoidGlyph::setVisible(bool visible) {
    if (actor_) actor_->SetVisibility(visible ? 1 : 0);
    if (barActor_) barActor_->SetVisibility(visible ? 1 : 0);
}

void EllipsoidGlyph::setData(
    const std::vector<nmr::Vec3>& positions,
    const std::vector<nmr::Mat3>& tensors,
    const std::vector<double>& isoValues,
    double scale,
    double opacity)
{
    clear();
    if (positions.empty()) return;

    int N = static_cast<int>(positions.size());

    vtkNew<vtkPoints> points;
    points->SetNumberOfPoints(N);

    // 9-component tensor array (row-major for VTK)
    // Mat3 is Eigen column-major, so we transpose when packing
    vtkNew<vtkDoubleArray> tensorArray;
    tensorArray->SetName("Tensors");
    tensorArray->SetNumberOfComponents(9);
    tensorArray->SetNumberOfTuples(N);

    // Scalar array for coloring by isotropic value
    vtkNew<vtkFloatArray> scalars;
    scalars->SetName("IsoValue");
    scalars->SetNumberOfComponents(1);
    scalars->SetNumberOfTuples(N);

    double vmin = 1e30, vmax = -1e30;

    for (int i = 0; i < N; ++i) {
        points->SetPoint(i, positions[i].data());

        // Transpose Mat3 (column-major Eigen) to row-major for VTK
        const nmr::Mat3& t = tensors[i];
        double row_major[9] = {
            t(0,0), t(0,1), t(0,2),
            t(1,0), t(1,1), t(1,2),
            t(2,0), t(2,1), t(2,2)
        };
        tensorArray->SetTuple(i, row_major);

        float val = static_cast<float>(isoValues[i]);
        scalars->SetValue(i, val);
        if (val < vmin) vmin = val;
        if (val > vmax) vmax = val;
    }

    vtkNew<vtkPolyData> pd;
    pd->SetPoints(points);
    pd->GetPointData()->SetTensors(tensorArray);
    pd->GetPointData()->SetScalars(scalars);

    vtkNew<vtkTrivialProducer> producer;
    producer->SetOutput(pd);

    // Sphere source for glyph shape
    vtkNew<vtkSphereSource> sphere;
    sphere->SetThetaResolution(16);
    sphere->SetPhiResolution(16);

    // Tensor glyph filter — VTK handles eigendecomposition internally
    vtkNew<vtkTensorGlyph> tensorGlyph;
    tensorGlyph->SetInputConnection(producer->GetOutputPort());
    tensorGlyph->SetSourceConnection(sphere->GetOutputPort());
    tensorGlyph->SetExtractEigenvalues(true);
    tensorGlyph->SetScaleFactor(scale);
    tensorGlyph->ColorGlyphsOn();
    tensorGlyph->SetColorModeToScalars();
    tensorGlyph->Update();

    // Diverging colormap
    auto ctf = makeDivergingCTF(vmin, vmax);
    double absMax = std::max(std::abs(vmin), std::abs(vmax));
    if (absMax < 1e-10) absMax = 1.0;

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputConnection(tensorGlyph->GetOutputPort());
    mapper->SetLookupTable(ctf);
    mapper->SetScalarRange(-absMax, absMax);
    mapper->SetScalarModeToUsePointData();

    actor_ = vtkSmartPointer<vtkActor>::New();
    actor_->SetMapper(mapper);
    actor_->GetProperty()->SetOpacity(opacity);
    actor_->GetProperty()->SetInterpolationToPhong();
    actor_->SetForceTranslucent(true);
    renderer_->AddActor(actor_);

    // Color bar
    auto bar = vtkSmartPointer<vtkScalarBarActor>::New();
    bar->SetLookupTable(ctf);
    bar->SetNumberOfLabels(5);
    bar->SetWidth(0.08);
    bar->SetHeight(0.4);
    bar->SetPosition(0.02, 0.3);
    bar->GetLabelTextProperty()->SetFontSize(10);
    bar->GetLabelTextProperty()->SetColor(0.9, 0.9, 0.9);
    bar->GetTitleTextProperty()->SetColor(0.9, 0.9, 0.9);
    bar->SetTitle("ppm");
    barActor_ = bar;
    renderer_->AddActor2D(bar);
}
