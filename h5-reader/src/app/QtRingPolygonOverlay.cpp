#include "QtRingPolygonOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/QtFrame.h"
#include "../model/QtResidue.h"

#include <QLoggingCategory>

#include <vtkArrowSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolygon.h>
#include <vtkProperty.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>

#include <cmath>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cRings, "h5reader.overlay.rings")

// Number of vertices sampled around each ring polygon. 12 is the
// library-viewer value; gives a visibly round polygon at normal zoom.
constexpr int kPolygonVertices = 12;

// Offset of the translucent "above / below" polygons from the ring
// plane, in Ångströms. 0.15 Å is the library-viewer value.
constexpr double kPolygonOffsetA = 0.15;

// Polygon colour by amino-acid type. Typed dispatch, not string-match.
// Palette matches the existing ui/src/RingCurrentOverlay.cpp.
void RingRgbForAminoAcid(model::AminoAcid aa, unsigned char rgb[3]) {
    switch (aa) {
        case model::AminoAcid::PHE: rgb[0] = 77;  rgb[1] = 230; rgb[2] = 77;  break;
        case model::AminoAcid::TYR: rgb[0] = 51;  rgb[1] = 217; rgb[2] = 217; break;
        case model::AminoAcid::TRP: rgb[0] = 217; rgb[1] = 77;  rgb[2] = 217; break;
        case model::AminoAcid::HIS: rgb[0] = 242; rgb[1] = 230; rgb[2] = 64;  break;
        default:                    rgb[0] = 200; rgb[1] = 200; rgb[2] = 100; break;
    }
}

}  // namespace

QtRingPolygonOverlay::QtRingPolygonOverlay(
    vtkSmartPointer<vtkRenderer>                  renderer,
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
    QObject* parent)
    : QObject(parent),
      renderer_(std::move(renderer)),
      renderWindow_(std::move(renderWindow))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtRingPolygonOverlay"));
}

QtRingPolygonOverlay::~QtRingPolygonOverlay() {
    for (const auto& ra : rings_) {
        if (ra.actorAbove) renderer_->RemoveActor(ra.actorAbove);
        if (ra.actorBelow) renderer_->RemoveActor(ra.actorBelow);
        if (ra.arrowActor) renderer_->RemoveActor(ra.arrowActor);
    }
}

void QtRingPolygonOverlay::Build(const model::QtProtein&      protein,
                                  const model::QtConformation& conformation) {
    ASSERT_THREAD(this);

    if (protein_ == &protein && conformation_ == &conformation && !rings_.empty())
        return;

    for (const auto& ra : rings_) {
        if (ra.actorAbove) renderer_->RemoveActor(ra.actorAbove);
        if (ra.actorBelow) renderer_->RemoveActor(ra.actorBelow);
        if (ra.arrowActor) renderer_->RemoveActor(ra.arrowActor);
    }
    rings_.clear();

    protein_      = &protein;
    conformation_ = &conformation;

    const size_t n_rings = protein.ringCount();
    rings_.reserve(n_rings);

    for (size_t ri = 0; ri < n_rings; ++ri) {
        const auto& ring = protein.ring(ri);
        const int residueIndex = ring.parentResidueIndex;
        const model::AminoAcid aa =
            residueIndex >= 0 && static_cast<size_t>(residueIndex) < protein.residueCount()
                ? protein.residue(residueIndex).aminoAcid
                : model::AminoAcid::Unknown;

        RingActor ra;
        RingRgbForAminoAcid(aa, ra.rgb);

        for (int side = 0; side < 2; ++side) {
            auto points = vtkSmartPointer<vtkPoints>::New();
            points->SetNumberOfPoints(kPolygonVertices);
            for (int v = 0; v < kPolygonVertices; ++v) {
                points->SetPoint(v, 0.0, 0.0, 0.0);   // filled in setFrame
            }

            auto polygon = vtkSmartPointer<vtkPolygon>::New();
            polygon->GetPointIds()->SetNumberOfIds(kPolygonVertices);
            for (int v = 0; v < kPolygonVertices; ++v)
                polygon->GetPointIds()->SetId(v, v);

            auto polys = vtkSmartPointer<vtkCellArray>::New();
            polys->InsertNextCell(polygon);

            auto pd = vtkSmartPointer<vtkPolyData>::New();
            pd->SetPoints(points);
            pd->SetPolys(polys);

            auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mapper->SetInputData(pd);

            auto actor = vtkSmartPointer<vtkActor>::New();
            actor->SetMapper(mapper);
            const double r = ra.rgb[0] / 255.0;
            const double g = ra.rgb[1] / 255.0;
            const double b = ra.rgb[2] / 255.0;
            actor->GetProperty()->SetColor(r, g, b);
            actor->GetProperty()->SetOpacity(0.35);
            actor->GetProperty()->LightingOff();
            actor->GetProperty()->SetEdgeVisibility(1);
            actor->GetProperty()->SetEdgeColor(r * 0.7, g * 0.7, b * 0.7);
            actor->GetProperty()->SetLineWidth(2.0);
            actor->SetForceTranslucent(true);

            renderer_->AddActor(actor);
            if (side == 0) { ra.pointsAbove = points; ra.polyAbove = pd; ra.actorAbove = actor; }
            else           { ra.pointsBelow = points; ra.polyBelow = pd; ra.actorBelow = actor; }
        }

        // Arrow along the ring normal. Geometry is static except for
        // the transform, which we update per frame via the actor's
        // user transform rather than re-running the pipeline.
        auto arrow = vtkSmartPointer<vtkArrowSource>::New();
        arrow->SetTipResolution(16);
        arrow->SetShaftResolution(16);
        arrow->SetTipLength(0.3);
        arrow->SetTipRadius(0.12);
        arrow->SetShaftRadius(0.04);

        auto arrowMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        arrowMapper->SetInputConnection(arrow->GetOutputPort());

        auto arrowActor = vtkSmartPointer<vtkActor>::New();
        arrowActor->SetMapper(arrowMapper);
        arrowActor->GetProperty()->SetColor(ra.rgb[0] / 255.0,
                                             ra.rgb[1] / 255.0,
                                             ra.rgb[2] / 255.0);
        arrowActor->GetProperty()->SetOpacity(0.9);
        renderer_->AddActor(arrowActor);
        ra.arrowActor = arrowActor;

        rings_.push_back(ra);
    }

    qCInfo(cRings).noquote()
        << "Built" << rings_.size() << "ring polygon actors";

    setFrame(0);
}

void QtRingPolygonOverlay::UpdateRingActor(
    RingActor& ra, const model::RingGeometry& geo) {
    if (geo.radius < 1e-6) {
        ra.actorAbove->SetVisibility(0);
        ra.actorBelow->SetVisibility(0);
        ra.arrowActor->SetVisibility(0);
        return;
    }
    if (visible_) {
        ra.actorAbove->SetVisibility(1);
        ra.actorBelow->SetVisibility(1);
        ra.arrowActor->SetVisibility(1);
    }

    const auto basis = model::OrthoBasisFromNormal(geo.normal);
    const model::Vec3& n = basis.n;
    const model::Vec3& u = basis.u;
    const model::Vec3& v = basis.v;

    for (int side = 0; side < 2; ++side) {
        const double sign = side == 0 ? +1.0 : -1.0;
        const model::Vec3 offset = n * (kPolygonOffsetA * sign);

        auto* points = side == 0 ? ra.pointsAbove.Get() : ra.pointsBelow.Get();
        auto* poly   = side == 0 ? ra.polyAbove.Get()   : ra.polyBelow.Get();

        for (int i = 0; i < kPolygonVertices; ++i) {
            const double a = 2.0 * M_PI * i / kPolygonVertices;
            const model::Vec3 p =
                geo.center + geo.radius * (u * std::cos(a) + v * std::sin(a)) + offset;
            points->SetPoint(i, p.x(), p.y(), p.z());
        }
        points->Modified();
        poly->Modified();
    }

    // Arrow transform: align +x (the arrow source's default axis) with
    // the ring normal, translate to 0.3 Å above the centre, scale 1.2x.
    auto xform = vtkSmartPointer<vtkTransform>::New();
    xform->PostMultiply();
    const model::Vec3 up = geo.normal;
    const model::Vec3 xAxis(1, 0, 0);
    model::Vec3 rotAxis = xAxis.cross(up);
    const double rotAngle =
        std::acos(std::clamp(xAxis.dot(up), -1.0, 1.0)) * 180.0 / M_PI;
    xform->Scale(1.2, 1.2, 1.2);
    if (rotAxis.norm() > 1e-6) {
        rotAxis.normalize();
        xform->RotateWXYZ(rotAngle, rotAxis.x(), rotAxis.y(), rotAxis.z());
    } else if (xAxis.dot(up) < 0) {
        xform->RotateZ(180.0);
    }
    const model::Vec3 arrowBase = geo.center + geo.normal * 0.3;
    xform->Translate(arrowBase.x(), arrowBase.y(), arrowBase.z());
    ra.arrowActor->SetUserTransform(xform);
}

void QtRingPolygonOverlay::setFrame(int t) {
    ASSERT_THREAD(this);
    if (!protein_ || !conformation_) return;
    if (t < 0 || static_cast<size_t>(t) >= conformation_->frameCount()) return;

    const auto& frame = conformation_->frame(static_cast<size_t>(t));
    for (size_t ri = 0; ri < rings_.size(); ++ri) {
        const auto geo = frame.ringGeometry(ri);
        UpdateRingActor(rings_[ri], geo);
    }
    // MoleculeScene will invoke Render() after all overlays finish — we
    // don't Render here to avoid double renders per frame.
}

void QtRingPolygonOverlay::setVisible(bool visible) {
    ASSERT_THREAD(this);
    visible_ = visible;
    for (auto& ra : rings_) {
        if (ra.actorAbove) ra.actorAbove->SetVisibility(visible ? 1 : 0);
        if (ra.actorBelow) ra.actorBelow->SetVisibility(visible ? 1 : 0);
        if (ra.arrowActor) ra.arrowActor->SetVisibility(visible ? 1 : 0);
    }
}

}  // namespace h5reader::app
