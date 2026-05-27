#include "MeasurementOverlay.h"

#include "../model/AtomSelection.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <QLoggingCategory>

#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cMeas, "h5reader.overlay.measurement")

// A touch smaller than QtSelectionOverlay's 1.0 Å so the slot colour wraps the
// atom without hiding the ball-and-stick sphere underneath; translucent so
// overlapping selections stay readable.
constexpr double kSphereRadiusA = 0.85;
constexpr double kOpacity       = 0.50;

// Connecting polyline: a neutral, near-white line that reads against any slot
// colour without competing with the spheres. Width is in screen pixels.
constexpr double kLineRgb[3]  = {0.92, 0.92, 0.92};
constexpr double kLineWidth   = 2.5;
constexpr double kLineOpacity = 0.90;
}  // namespace

MeasurementOverlay::MeasurementOverlay(
    vtkSmartPointer<vtkRenderer>                  renderer,
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
    QObject*                                      parent)
    : QObject(parent),
      renderer_(std::move(renderer)),
      renderWindow_(std::move(renderWindow))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("MeasurementOverlay"));
}

MeasurementOverlay::~MeasurementOverlay() {
    for (auto& actor : actors_) {
        if (actor) renderer_->RemoveActor(actor);
    }
    if (lineActor_) renderer_->RemoveActor(lineActor_);
}

void MeasurementOverlay::Build(const model::QtProtein& protein,
                                model::Conformation&    conformation) {
    ASSERT_THREAD(this);
    if (protein_ == &protein && conformation_ == &conformation && actors_[0])
        return;

    for (auto& actor : actors_) {
        if (actor) renderer_->RemoveActor(actor);
        actor = nullptr;
    }
    if (lineActor_) {
        renderer_->RemoveActor(lineActor_);
        lineActor_ = nullptr;
    }

    protein_      = &protein;
    conformation_ = &conformation;

    // One sphere+actor per slot, each fixed to its Okabe-Ito slot colour and
    // hidden until the selection populates the slot. Reusing a fixed set of
    // four actors (rather than creating/destroying per pick) keeps the VTK
    // pipeline stable across selection changes.
    for (std::size_t s = 0; s < kMaxSpheres; ++s) {
        spheres_[s] = vtkSmartPointer<vtkSphereSource>::New();
        spheres_[s]->SetRadius(kSphereRadiusA);
        spheres_[s]->SetPhiResolution(16);
        spheres_[s]->SetThetaResolution(16);

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(spheres_[s]->GetOutputPort());

        const std::array<double, 3> c = model::AtomSelection::SlotColorRgb(s);
        actors_[s] = vtkSmartPointer<vtkActor>::New();
        actors_[s]->SetMapper(mapper);
        actors_[s]->GetProperty()->SetColor(c[0], c[1], c[2]);
        actors_[s]->GetProperty()->SetOpacity(kOpacity);
        actors_[s]->SetVisibility(0);
        renderer_->AddActor(actors_[s]);
    }

    // Connecting polyline (created once, points + cell refilled per frame). A
    // single vtkPolyData with one polyline cell through the selected atoms in
    // slot order; hidden until the selection defines at least a distance.
    linePoints_ = vtkSmartPointer<vtkPoints>::New();
    lineCells_  = vtkSmartPointer<vtkCellArray>::New();
    lineData_   = vtkSmartPointer<vtkPolyData>::New();
    lineData_->SetPoints(linePoints_);
    lineData_->SetLines(lineCells_);

    auto lineMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    lineMapper->SetInputData(lineData_);

    lineActor_ = vtkSmartPointer<vtkActor>::New();
    lineActor_->SetMapper(lineMapper);
    lineActor_->GetProperty()->SetColor(kLineRgb[0], kLineRgb[1], kLineRgb[2]);
    lineActor_->GetProperty()->SetLineWidth(kLineWidth);
    lineActor_->GetProperty()->SetOpacity(kLineOpacity);
    lineActor_->SetVisibility(0);
    renderer_->AddActor(lineActor_);
}

void MeasurementOverlay::setSelection(model::AtomSelection* selection) {
    ASSERT_THREAD(this);
    selection_ = selection;
    // Sync immediately so a pre-existing selection (none at startup) shows.
    applyFrame(lastFrame_);
}

void MeasurementOverlay::onSelectionChanged() {
    ASSERT_THREAD(this);
    qCDebug(cMeas).noquote()
        << "selection changed | n=" << (selection_ ? static_cast<int>(selection_->count()) : 0);
    applyFrame(lastFrame_);
}

void MeasurementOverlay::setFrame(int t) {
    ASSERT_THREAD(this);
    lastFrame_ = t;
    applyFrame(t);
}

void MeasurementOverlay::setVisible(bool on) {
    ASSERT_THREAD(this);
    visible_ = on;
    applyFrame(lastFrame_);
}

void MeasurementOverlay::applyFrame(int t) {
    if (!actors_[0] || !protein_ || !conformation_)
        return;
    if (t < 0 || static_cast<std::size_t>(t) >= conformation_->frameCount())
        return;

    const std::size_t n = selection_ ? selection_->count() : 0;

    std::array<model::Vec3, kMaxSpheres> pos;
    bool                                 allInRange = true;

    for (std::size_t s = 0; s < kMaxSpheres; ++s) {
        bool show = visible_ && s < n;
        if (show) {
            const std::size_t a = selection_->atoms()[s];
            if (a >= protein_->atomCount()) {
                show       = false;
                allInRange = false;
            } else {
                const model::Vec3 p =
                    conformation_->atomPosition(static_cast<std::size_t>(t), a);
                spheres_[s]->SetCenter(p.x(), p.y(), p.z());
                pos[s] = p;
            }
        }
        actors_[s]->SetVisibility(show ? 1 : 0);
    }

    // Connecting polyline: shown when the selection defines at least a distance
    // (>= 2 atoms), every atom is in range, and the overlay is visible. Rebuilt
    // from the slot-ordered positions each frame so it holds through rotation.
    const bool showLine = visible_ && allInRange && n >= 2;
    if (showLine) {
        linePoints_->Reset();
        for (std::size_t i = 0; i < n; ++i)
            linePoints_->InsertNextPoint(pos[i].x(), pos[i].y(), pos[i].z());

        std::array<vtkIdType, kMaxSpheres> ids;
        for (std::size_t i = 0; i < n; ++i)
            ids[i] = static_cast<vtkIdType>(i);
        lineCells_->Reset();
        lineCells_->InsertNextCell(static_cast<vtkIdType>(n), ids.data());

        linePoints_->Modified();
        lineCells_->Modified();
        lineData_->Modified();
    }
    lineActor_->SetVisibility(showLine ? 1 : 0);
}

}  // namespace h5reader::app
