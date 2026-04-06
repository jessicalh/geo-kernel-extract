#include "RingCurrentOverlay.h"

#include "Protein.h"
#include "ProteinConformation.h"
#include "Ring.h"
#include "Residue.h"

#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkNew.h>
#include <vtkPolygon.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkArrowSource.h>
#include <cmath>

using namespace nmr;

RingCurrentOverlay::RingCurrentOverlay(
    vtkSmartPointer<vtkRenderer> renderer,
    const Protein& protein,
    const ProteinConformation& conf)
    : renderer_(renderer)
{
    buildActors(protein, conf);
}

RingCurrentOverlay::~RingCurrentOverlay() {
    for (auto& actor : actors_)
        renderer_->RemoveActor(actor);
}

// Color by parent residue type (the ring knows its residue)
static void ringColor(AminoAcid type, double& r, double& g, double& b) {
    switch (type) {
        case AminoAcid::PHE: r = 0.3;  g = 0.9;  b = 0.3;  break;
        case AminoAcid::TYR: r = 0.2;  g = 0.85; b = 0.85; break;
        case AminoAcid::TRP: r = 0.85; g = 0.3;  b = 0.85; break;
        default:             r = 0.95; g = 0.9;  b = 0.25; break; // HIS variants
    }
}

void RingCurrentOverlay::buildActors(
    const Protein& protein,
    const ProteinConformation& conf)
{
    for (size_t i = 0; i < protein.RingCount(); ++i) {
        const Ring& ring = protein.RingAt(i);
        const RingGeometry& geo = conf.ring_geometries[i];

        if (geo.radius < 0.01) continue;

        AminoAcid resType = protein.ResidueAt(ring.parent_residue_index).type;

        // Generate polygon vertices from center + normal + radius
        constexpr int nVerts = 12;
        Vec3 n = geo.normal.normalized();

        Vec3 arbitrary = (std::abs(n.x()) < 0.9) ? Vec3(1, 0, 0) : Vec3(0, 1, 0);
        Vec3 u = n.cross(arbitrary).normalized();
        Vec3 v = n.cross(u);

        // Translucent filled polygon — two copies offset above/below ring plane
        for (int side = 0; side < 2; ++side) {
            double sign = (side == 0) ? 1.0 : -1.0;
            Vec3 offset = n * (0.15 * sign);

            vtkNew<vtkPoints> pts;
            vtkNew<vtkCellArray> polys;
            vtkNew<vtkPolygon> polygon;
            polygon->GetPointIds()->SetNumberOfIds(nVerts);
            for (int vi = 0; vi < nVerts; ++vi) {
                double angle = 2.0 * M_PI * vi / nVerts;
                Vec3 p = geo.center + geo.radius * (u * std::cos(angle) + v * std::sin(angle)) + offset;
                pts->InsertNextPoint(p.data());
                polygon->GetPointIds()->SetId(vi, vi);
            }
            polys->InsertNextCell(polygon);

            vtkNew<vtkPolyData> polyData;
            polyData->SetPoints(pts);
            polyData->SetPolys(polys);

            vtkNew<vtkPolyDataMapper> mapper;
            mapper->SetInputData(polyData);

            vtkNew<vtkActor> actor;
            actor->SetMapper(mapper);

            double r, g, b;
            ringColor(resType, r, g, b);
            actor->GetProperty()->SetColor(r, g, b);
            actor->GetProperty()->SetOpacity(0.35);
            actor->GetProperty()->LightingOff();
            actor->GetProperty()->SetEdgeVisibility(1);
            actor->GetProperty()->SetEdgeColor(r * 0.7, g * 0.7, b * 0.7);
            actor->GetProperty()->SetLineWidth(2.0);

            renderer_->AddActor(actor);
            actors_.push_back(actor);
        }

        // Small arrow along ring normal indicating current direction
        vtkNew<vtkArrowSource> arrow;
        arrow->SetTipResolution(16);
        arrow->SetShaftResolution(16);
        arrow->SetTipLength(0.3);
        arrow->SetTipRadius(0.12);
        arrow->SetShaftRadius(0.04);

        vtkNew<vtkTransform> xform;
        xform->PostMultiply();

        Vec3 up = geo.normal;
        Vec3 xAxis(1, 0, 0);
        Vec3 rotAxis = xAxis.cross(up);
        double rotAngle = std::acos(std::clamp(xAxis.dot(up), -1.0, 1.0))
                          * 180.0 / M_PI;

        xform->Scale(1.2, 1.2, 1.2);
        if (rotAxis.norm() > 1e-6) {
            rotAxis.normalize();
            xform->RotateWXYZ(rotAngle, rotAxis.x(), rotAxis.y(), rotAxis.z());
        } else if (xAxis.dot(up) < 0) {
            xform->RotateZ(180.0);
        }
        Vec3 arrowBase = geo.center + geo.normal * 0.3;
        xform->Translate(arrowBase.x(), arrowBase.y(), arrowBase.z());

        vtkNew<vtkTransformPolyDataFilter> xformFilter;
        xformFilter->SetInputConnection(arrow->GetOutputPort());
        xformFilter->SetTransform(xform);

        vtkNew<vtkPolyDataMapper> arrowMapper;
        arrowMapper->SetInputConnection(xformFilter->GetOutputPort());

        vtkNew<vtkActor> arrowActor;
        arrowActor->SetMapper(arrowMapper);

        double r, g, b;
        ringColor(resType, r, g, b);
        arrowActor->GetProperty()->SetColor(r, g, b);
        arrowActor->GetProperty()->SetOpacity(0.9);

        renderer_->AddActor(arrowActor);
        actors_.push_back(arrowActor);
    }
}

void RingCurrentOverlay::setVisible(bool visible) {
    for (auto& actor : actors_)
        actor->SetVisibility(visible ? 1 : 0);
}

void RingCurrentOverlay::setTubeRadius(double radius) {
    tubeRadius_ = radius;
}
