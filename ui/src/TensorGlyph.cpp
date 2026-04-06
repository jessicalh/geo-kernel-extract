#include "TensorGlyph.h"
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPointData.h>
#include <vtkTriangle.h>
#include <vtkNew.h>
#include <cmath>

TensorGlyph::TensorGlyph(vtkSmartPointer<vtkRenderer> renderer)
    : renderer_(renderer)
{
}

TensorGlyph::~TensorGlyph() {
    clear();
}

void TensorGlyph::clear() {
    for (auto& actor : actors_)
        renderer_->RemoveActor(actor);
    actors_.clear();
}

void TensorGlyph::setVisible(bool visible) {
    for (auto& actor : actors_)
        actor->SetVisibility(visible ? 1 : 0);
}

vtkSmartPointer<vtkActor> TensorGlyph::createGlyph(
    const nmr::Vec3& center,
    const nmr::SphericalTensor& st,
    double scale,
    double opacity)
{
    if (st.T2Magnitude() < 1e-6)
        return nullptr;

    vtkNew<vtkPoints> points;
    vtkNew<vtkCellArray> triangles;
    vtkNew<vtkUnsignedCharArray> colors;
    colors->SetNumberOfComponents(3);
    colors->SetName("Colors");

    // Generate sphere vertices
    int nTheta = THETA_RES;
    int nPhi = PHI_RES;

    // Store point indices in a grid
    std::vector<std::vector<int>> grid(nTheta + 1, std::vector<int>(nPhi));

    for (int i = 0; i <= nTheta; ++i) {
        double theta = M_PI * i / nTheta;
        for (int j = 0; j < nPhi; ++j) {
            double phi = 2.0 * M_PI * j / nPhi;

            // f(theta,phi) = sum_m T2_m * Y_2^m(theta,phi)
            // Real spherical harmonics (Condon-Shortley convention):
            double ct = std::cos(theta), st2 = std::sin(theta);
            double cp = std::cos(phi), sp = std::sin(phi);
            double c2p = std::cos(2*phi), s2p = std::sin(2*phi);
            // Y_2^{-2} = 0.5*sqrt(15/pi)*sin^2(theta)*sin(2phi)
            // Y_2^{-1} = 0.5*sqrt(15/pi)*sin(2theta)*sin(phi)
            // Y_2^{0}  = 0.25*sqrt(5/pi)*(3cos^2(theta)-1)
            // Y_2^{1}  = -0.5*sqrt(15/pi)*sin(2theta)*cos(phi)
            // Y_2^{2}  = 0.5*sqrt(15/pi)*sin^2(theta)*cos(2phi)
            double s15p = 0.5 * std::sqrt(15.0 / M_PI);
            double s5p  = 0.25 * std::sqrt(5.0 / M_PI);
            double f = st.T2[0] * s15p * st2 * st2 * s2p     // m=-2
                     + st.T2[1] * s15p * std::sin(2*theta) * sp  // m=-1
                     + st.T2[2] * s5p * (3*ct*ct - 1)         // m=0
                     + st.T2[3] * (-s15p) * std::sin(2*theta) * cp // m=1
                     + st.T2[4] * s15p * st2 * st2 * c2p;     // m=2
            double r = std::abs(f) * scale;

            double x = center(0) + r * std::sin(theta) * std::cos(phi);
            double y = center(1) + r * std::sin(theta) * std::sin(phi);
            double z = center(2) + r * std::cos(theta);

            int id = (int)points->InsertNextPoint(x, y, z);
            grid[i][j] = id;

            // Color: red for positive (deshielded), blue for negative (shielded)
            if (f > 0) {
                unsigned char c[3] = {200, 50, 50};
                colors->InsertNextTypedTuple(c);
            } else {
                unsigned char c[3] = {50, 50, 200};
                colors->InsertNextTypedTuple(c);
            }
        }
    }

    // Create triangles
    for (int i = 0; i < nTheta; ++i) {
        for (int j = 0; j < nPhi; ++j) {
            int j1 = (j + 1) % nPhi;

            vtkNew<vtkTriangle> t1;
            t1->GetPointIds()->SetId(0, grid[i][j]);
            t1->GetPointIds()->SetId(1, grid[i + 1][j]);
            t1->GetPointIds()->SetId(2, grid[i + 1][j1]);
            triangles->InsertNextCell(t1);

            vtkNew<vtkTriangle> t2;
            t2->GetPointIds()->SetId(0, grid[i][j]);
            t2->GetPointIds()->SetId(1, grid[i + 1][j1]);
            t2->GetPointIds()->SetId(2, grid[i][j1]);
            triangles->InsertNextCell(t2);
        }
    }

    vtkNew<vtkPolyData> polyData;
    polyData->SetPoints(points);
    polyData->SetPolys(triangles);
    polyData->GetPointData()->SetScalars(colors);

    vtkNew<vtkPolyDataMapper> mapper;
    mapper->SetInputData(polyData);

    vtkNew<vtkActor> actor;
    actor->SetMapper(mapper);
    actor->GetProperty()->SetOpacity(opacity);
    actor->GetProperty()->SetInterpolationToPhong();

    return actor;
}

void TensorGlyph::setData(
    const std::vector<nmr::Vec3>& positions,
    const std::vector<nmr::SphericalTensor>& tensors,
    double scale,
    double opacity)
{
    clear();

    for (size_t i = 0; i < positions.size(); ++i) {
        auto actor = createGlyph(positions[i], tensors[i], scale, opacity);
        if (actor) {
            renderer_->AddActor(actor);
            actors_.push_back(actor);
        }
    }
}
