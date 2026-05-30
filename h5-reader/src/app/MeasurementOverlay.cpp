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

// User-mode marker palette — pastel-shifted Okabe-Ito (50% mix with white).
// Local to the overlay because the dashboard dock's slot swatches keep the
// saturated AtomSelection::SlotColorRgb values for on-screen legibility on
// a UI panel; the 3-D scene wants the lighter, more translucent pastels
// against the molecule. Hue order matches AtomSelection's slot order, so
// slot 0 here pairs with the orange swatch in the dock, etc.
//
// Distinct by design from kInstrumentRgb below: bright, opaque, 1.5 Å
// markers are the harness/debug preset (POST /selection/instrument and the
// Instrument toolbar toggle); these pastel translucent ones are normal
// user-facing selection display.
constexpr double kPastelSlotRgb[4][3] = {
    {0.951, 0.812, 0.500},  // slot 0 — pale peach   (orange + 50% white)
    {0.669, 0.853, 0.957},  // slot 1 — pale sky     (sky blue + 50% white)
    {0.500, 0.810, 0.726},  // slot 2 — pale mint    (bluish green + 50% white)
    {0.900, 0.738, 0.828},  // slot 3 — pale rose    (reddish purple + 50% white)
};

// Connecting polyline: a neutral, near-white line that reads against any slot
// colour without competing with the spheres. Width is in screen pixels.
constexpr double kLineRgb[3]  = {0.92, 0.92, 0.92};
constexpr double kLineWidth   = 2.5;
constexpr double kLineOpacity = 0.90;

// Instrument mode — CPK-distinct marker palette. Chosen 2026-05-30 (memory
// VIEWPORT_OBSERVATIONS_2026-05-30.md §5b): hues that fall outside every CPK
// element colour so a hue threshold isolates the marker against any rendered
// scene. Bright + opaque + 1.5 Å so the harness blob detector finds it
// reliably; kept distinct from kPastelSlotRgb so users can tell at a glance
// when debug/harness mode is active. RGB in 0..1.
constexpr double kInstrumentRgb[4][3] = {
    {1.000, 0.000, 1.000},  // slot 0 — pure magenta       (#FF00FF)
    {0.000, 1.000, 0.498},  // slot 1 — spring green       (#00FF7F)
    {1.000, 0.078, 0.576},  // slot 2 — deep pink          (#FF1493)
    {0.616, 0.000, 1.000},  // slot 3 — vivid violet       (#9D00FF)
};
constexpr double kInstrumentOpacity = 1.0;
constexpr double kInstrumentRadiusA = 1.5;
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

        actors_[s] = vtkSmartPointer<vtkActor>::New();
        actors_[s]->SetMapper(mapper);
        actors_[s]->GetProperty()->SetColor(
            kPastelSlotRgb[s][0], kPastelSlotRgb[s][1], kPastelSlotRgb[s][2]);
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

void MeasurementOverlay::setInstrumentMode(bool on, bool focusOnly) {
    ASSERT_THREAD(this);
    if (instrumentMode_ == on && focusOnly_ == focusOnly)
        return;
    instrumentMode_ = on;
    // focusOnly is only meaningful while instrument mode is on; clear it
    // when turning instrument mode off so a future on-without-focusOnly
    // call doesn't inherit a stale flag.
    focusOnly_ = on ? focusOnly : false;
    applyInstrumentModeToActors();
    // Reposition the visible spheres so the new radius/colour shows at the
    // current frame without waiting for the next setFrame. applyFrame()
    // also honours focusOnly_ in its visibility branch. The caller (REST
    // handler) follows with scene_->requestRender() to flush, matching the
    // setVisible() flow at ReaderMainWindow.cpp:591.
    applyFrame(lastFrame_);
    qCInfo(cMeas).noquote() << "instrument mode" << (on ? "ON" : "OFF")
                            << "focus_only=" << focusOnly_;
}

void MeasurementOverlay::applyInstrumentModeToActors() {
    // Walks the four pre-built sphere actors and switches their colour,
    // opacity, and radius wholesale. Reversible: when off, restore the
    // build-time Okabe-Ito palette + kSphereRadiusA + kOpacity (hard-coded
    // here rather than refactoring the anonymous-namespace constants — the
    // refactor is later viewport-overlay-layer work).
    //
    // focus-only variant: when both instrument and focusOnly are on, every
    // sphere gets the slot-0 magenta colour so the focus marker reads as a
    // single magenta blob regardless of which slot the focus lives in.
    // Visibility (only the focus slot rendered) is decided in applyFrame().
    for (std::size_t s = 0; s < kMaxSpheres; ++s) {
        if (!spheres_[s] || !actors_[s])
            continue;
        if (instrumentMode_) {
            if (focusOnly_) {
                actors_[s]->GetProperty()->SetColor(
                    kInstrumentRgb[0][0], kInstrumentRgb[0][1], kInstrumentRgb[0][2]);
            } else {
                actors_[s]->GetProperty()->SetColor(
                    kInstrumentRgb[s][0], kInstrumentRgb[s][1], kInstrumentRgb[s][2]);
            }
            actors_[s]->GetProperty()->SetOpacity(kInstrumentOpacity);
            spheres_[s]->SetRadius(kInstrumentRadiusA);
        } else {
            actors_[s]->GetProperty()->SetColor(
                kPastelSlotRgb[s][0], kPastelSlotRgb[s][1], kPastelSlotRgb[s][2]);
            actors_[s]->GetProperty()->SetOpacity(kOpacity);
            spheres_[s]->SetRadius(kSphereRadiusA);
        }
    }
}

void MeasurementOverlay::applyFrame(int t) {
    if (!actors_[0] || !protein_ || !conformation_)
        return;
    if (t < 0 || static_cast<std::size_t>(t) >= conformation_->frameCount())
        return;

    const std::size_t n = selection_ ? selection_->count() : 0;

    // Focus-only suppression: when the harness preset is active, only the
    // sphere at the focus slot is rendered. The slot of the focus atom is
    // selection_->slotOf(focusAtom); slots without focus do not render
    // regardless of how many atoms the selection holds.
    int focusSlot = -1;
    if (focusOnly_ && instrumentMode_ && selection_ && selection_->hasFocus()) {
        focusSlot = selection_->slotOf(selection_->focus());
    }

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
        // Focus-only override: only render the focus slot's sphere when
        // the harness preset is active. If there is no focus slot
        // (selection empty / no focus), nothing renders.
        if (focusOnly_ && instrumentMode_) {
            show = show && (static_cast<int>(s) == focusSlot);
        }
        actors_[s]->SetVisibility(show ? 1 : 0);
    }

    // Connecting polyline: shown when the selection defines at least a distance
    // (>= 2 atoms), every atom is in range, and the overlay is visible. Rebuilt
    // from the slot-ordered positions each frame so it holds through rotation.
    // Suppressed in focus-only instrument mode so the harness sees only one
    // magenta blob with no spurious connecting geometry near it.
    const bool showLine = visible_ && allInRange && n >= 2
                          && !(focusOnly_ && instrumentMode_);
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
