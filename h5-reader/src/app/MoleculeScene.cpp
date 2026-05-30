#include "MoleculeScene.h"

#include "CameraComposer.h"
#include "MeasurementOverlay.h"
#include "QtBackboneRibbonOverlay.h"
#include "QtBFieldStreamOverlay.h"
#include "QtFieldGridOverlay.h"
#include "QtRingPolygonOverlay.h"
#include "QtSelectionOverlay.h"
#include "QuietTrackballStyle.h"
#include "SceneRevealOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/Types.h"

#include <QElapsedTimer>
#include <QLoggingCategory>
#include <QMetaObject>
#include <QPointer>
#include <QVTKOpenGLNativeWidget.h>

#include <vtkActorCollection.h>

#include <cstdio>

#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>

#include <algorithm>
#include <cmath>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cScene, "h5reader.scene")

const char* RenderSourceName(MoleculeScene::RenderSource s) {
    using RS = MoleculeScene::RenderSource;
    switch (s) {
        case RS::Timer:       return "timer";
        case RS::CameraInput: return "camera-input";
        case RS::Picker:      return "picker";
        case RS::Rest:        return "rest";
        case RS::Overlay:     return "overlay";
        case RS::Reveal:      return "reveal";
        case RS::Shutdown:    return "shutdown";
        case RS::External:    return "external";
    }
    return "?";
}

// VTK bond-order encoding. vtkMolecule represents bond order as an
// unsigned short (1, 2, 3, …). We map our typed BondOrder enum onto
// those integers. Aromatic rendered as double (visually); peptide as
// single (the partial-double is handled upstream in our topology).
unsigned short VtkBondOrderFor(model::BondOrder o) {
    using model::BondOrder;
    switch (o) {
        case BondOrder::Single:   return 1;
        case BondOrder::Double:   return 2;
        case BondOrder::Triple:   return 3;
        case BondOrder::Aromatic: return 2;   // display as double
        case BondOrder::Peptide:  // display as single
        case BondOrder::Unknown:
            return 1;
    }
    return 1;
}

bool SameRevealBinding(const model::SignalBinding& a, const model::SignalBinding& b) {
    return a == b;
}
}  // namespace

MoleculeScene::MoleculeScene(QVTKOpenGLNativeWidget* vtkWidget,
                              vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
                              QObject* parent)
    : QObject(parent),
      vtkWidget_(vtkWidget),
      renderWindow_(std::move(renderWindow))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("MoleculeScene"));

    // ----- Renderer setup ------------------------------------------------
    //
    // Two-layer composition per spec §5: main scene at layer 0 + an
    // overlay layer at layer 1 sharing the same active camera. Markers
    // paint into the overlay layer so they remain visible regardless of
    // depth occlusion (the harness needs this for blob analysis).
    //
    // The existing nmr-viewer (ui/src/MainWindow.cpp:217-231) uses
    // FXAA + NO depth peeling after hitting translucency artifacts on
    // AMD hardware. We keep that choice for now because it's proven on
    // this codebase's target hardware. If a translucent overlay shows
    // sorting artifacts, flip to depth peeling per the qt6-cpp skill's
    // references/3d-vtk.md and log the change.

    renderer_ = vtkSmartPointer<vtkRenderer>::New();
    renderer_->SetLayer(0);
    renderer_->SetBackground(1.0, 1.0, 1.0);
    renderer_->SetUseFXAA(true);
    renderer_->SetUseDepthPeeling(0);

    overlayRenderer_ = vtkSmartPointer<vtkRenderer>::New();
    overlayRenderer_->SetLayer(1);
    // Shared camera with main renderer; markers always project to the
    // same screen-space coordinates as the atoms they mark.
    overlayRenderer_->SetActiveCamera(renderer_->GetActiveCamera());
    overlayRenderer_->SetBackground(0.0, 0.0, 0.0);
    overlayRenderer_->SetBackgroundAlpha(0.0);
    overlayRenderer_->SetUseFXAA(false);
    overlayRenderer_->SetUseDepthPeeling(0);
    overlayRenderer_->InteractiveOff();

    renderWindow_->SetNumberOfLayers(2);
    renderWindow_->AddRenderer(renderer_);
    renderWindow_->AddRenderer(overlayRenderer_);
    renderWindow_->SetAlphaBitPlanes(1);
    renderWindow_->SetMultiSamples(0);   // MSAA off — incompatible with
                                         // translucency; FXAA handles AA.

    if (auto* iren = renderWindow_->GetInteractor()) {
        // Quiet trackball: paint chain needs an interactor style (per
        // spec §2.7), but stock trackball mutates camera and schedules
        // renders behind our back. The subclass overrides every
        // manipulator to a no-op. EnableRender stays ON so the
        // adapter's iren->Render path still fires when paintGL drains.
        vtkNew<QuietTrackballStyle> style;
        style->AutoAdjustCameraClippingRangeOff();
        iren->SetInteractorStyle(style);
    }

    // EndEvent observer — per spec §2.8. Logs one structured line per
    // render with source + frame + ms + mode for the harness to
    // correlate render triggers against their cause. Renamed from the
    // earlier "render <ms>" form to include the source tag.
    auto endEventCb = vtkSmartPointer<vtkCallbackCommand>::New();
    endEventCb->SetClientData(this);
    endEventCb->SetCallback(
        [](vtkObject* /*caller*/, unsigned long, void* clientData, void*) {
            auto* self = static_cast<MoleculeScene*>(clientData);
            if (!self || !self->renderer_) return;
            const double ms = self->renderer_->GetLastRenderTimeInSeconds() * 1000.0;
            const char* src = RenderSourceName(self->lastRenderSource_);
            const char* modeName = "?";
            if (self->composer_) modeName = NameFor(self->composer_->mode().kind);
            qCInfo(cScene).noquote()
                << "render | source=" << src
                << "| frame=" << self->currentFrame_
                << "| ms=" << QString::number(ms, 'f', 1)
                << "| mode=" << modeName;
        });
    renderWindow_->AddObserver(vtkCommand::EndEvent, endEventCb);

    qCInfo(cScene).noquote()
        << "Renderer initialised: 2 layers (main FXAA + overlay), depth peeling OFF,"
        << "AlphaBitPlanes=1, MSAA=0, style=QuietTrackballStyle";
}

MoleculeScene::~MoleculeScene() {
    // VTK smart pointers clean up themselves. We just drop the references.
}

void MoleculeScene::Build(const model::QtProtein& protein,
                          model::Conformation&    conformation) {
    ASSERT_THREAD(this);

    if (protein_ == &protein && conformation_ == &conformation && molecule_) {
        return;   // already built with these inputs
    }

    QElapsedTimer timer;
    timer.start();

    // Fresh build — remove any prior actor.
    if (actor_) {
        renderer_->RemoveActor(actor_);
        actor_ = nullptr;
    }

    protein_      = &protein;
    conformation_ = &conformation;
    currentFrame_ = -1;

    molecule_ = vtkSmartPointer<vtkMolecule>::New();

    // Atoms — positions come from frame 0 so the first render is
    // consistent before any setFrame call. atomPosition is the shared
    // seam: the H5 for a trajectory, the snapshot's Pos column for a pose.
    for (size_t i = 0; i < protein.atomCount(); ++i) {
        const auto& atom = protein.atom(i);
        const model::Vec3 pos = conformation.atomPosition(0, i);
        const unsigned short z = static_cast<unsigned short>(
            model::AtomicNumberForElement(atom.element));
        molecule_->AppendAtom(z, pos.x(), pos.y(), pos.z());
    }

    // Bonds — connectivity is static across the trajectory.
    for (size_t i = 0; i < protein.bondCount(); ++i) {
        const auto& bond = protein.bond(i);
        molecule_->AppendBond(
            static_cast<vtkIdType>(bond.atomIndexA),
            static_cast<vtkIdType>(bond.atomIndexB),
            VtkBondOrderFor(bond.order));
    }

    // Mapper — GPU imposters for ball-and-stick, scales with molecule size.
    mapper_ = vtkSmartPointer<vtkOpenGLMoleculeMapper>::New();
    mapper_->SetInputData(molecule_);
    mapper_->UseBallAndStickSettings();

    actor_ = vtkSmartPointer<vtkActor>::New();
    actor_->SetMapper(mapper_);
    renderer_->AddActor(actor_);

    currentFrame_ = 0;

    // Camera composer — owns absolute per-frame camera writes per spec
    // §2.3. Constructed after the renderer is wired so it can hold the
    // smart pointer; protein + conformation pointers are non-owning.
    // Starts in CameraMode::Free (agent default per the implementation
    // prompt §4-c).
    if (!composer_) {
        composer_ = new CameraComposer(renderer_, &protein, &conformation, this);
    }
    composer_->setMode(FreeMode(), DefaultPolicy(), 0);

    // Overlays. Added AFTER the molecule actor so they render on top.
    // QObject parent = this — MoleculeScene's destruction destroys them.
    if (!ribbon_) {
        ribbon_ = new QtBackboneRibbonOverlay(renderer_, this);
    }
    ribbon_->Build(protein, conformation);

    if (!ringPolygons_) {
        ringPolygons_ = new QtRingPolygonOverlay(renderer_, renderWindow_, this);
    }
    ringPolygons_->Build(protein, conformation);

    if (!fieldGrid_) {
        fieldGrid_ = new QtFieldGridOverlay(renderer_, renderWindow_, this);
    }
    fieldGrid_->Build(protein, conformation);

    if (!bfieldStream_) {
        bfieldStream_ = new QtBFieldStreamOverlay(renderer_, renderWindow_, this);
    }
    bfieldStream_->Build(protein, conformation);

    if (!selection_) {
        selection_ = new QtSelectionOverlay(renderer_, renderWindow_, this);
    }
    selection_->Build(protein, conformation);

    // Markers go on the overlay-layer renderer (spec §5.2). The overlay
    // takes the main renderer for symmetry but only adds actors to the
    // overlay layer — the harness's marker-blob analysis needs the
    // marker to be findable regardless of depth occlusion by the
    // molecule's imposter spheres.
    if (!measurement_) {
        measurement_ = new MeasurementOverlay(overlayRenderer_, renderWindow_, this);
    }
    measurement_->Build(protein, conformation);

    if (!reveal_) {
        reveal_ = new SceneRevealOverlay(overlayRenderer_, this);
    }
    reveal_->Build(protein, conformation);

    qCInfo(cScene).noquote()
        << "Built molecule + overlays |"
        << "atoms=" << molecule_->GetNumberOfAtoms()
        << "| bonds=" << molecule_->GetNumberOfBonds()
        << "| rings=" << protein.ringCount()
        << "| build=" << timer.elapsed() << "ms";
}

void MoleculeScene::ResetCamera() {
    ASSERT_THREAD(this);
    if (!renderer_) return;
    renderer_->ResetCamera();
    // After the framing reset, sync the composer to the new camera
    // distance so a subsequent setMode(Atom/Bond/Plane) captures the
    // right zoom level. Composer Free mode is a no-op; this keeps the
    // initial view stable.
    if (composer_) {
        composer_->setMode(FreeMode(), DefaultPolicy(), 0);
    }
    requestRender(RenderSource::External);
}

void MoleculeScene::requestRender(RenderSource src) {
    ASSERT_THREAD(this);
    lastRenderSource_ = src;
    if (renderPending_ || !renderWindow_) return;
    renderPending_ = true;
    QPointer<MoleculeScene> self(this);
    QMetaObject::invokeMethod(this, [self]() {
        if (!self || !self->renderWindow_) return;
        self->renderPending_ = false;
        // widget->update() alone only schedules a Qt paint; paint then blits
        // the stale FBO because QVTKRenderWindowAdapter::paint() (lines 241-266)
        // only calls iren->Render() when DoVTKRenderInPaintGL is true, and that
        // flag is set by vtkRenderWindow::Render() → Frame() → adapter::frame().
        // So we have to call into VTK to actually re-render the FBO; the adapter
        // then schedules the widget paint that blits it.
        if (auto* iren = self->renderWindow_->GetInteractor()) {
            iren->Render();
        } else {
            self->renderWindow_->Render();
        }
    }, Qt::QueuedConnection);
}

void MoleculeScene::requestRender() {
    requestRender(RenderSource::External);
}

void MoleculeScene::refreshCurrentFrame() {
    ASSERT_THREAD(this);
    const int t = currentFrame_;
    // Bypass the same-frame guard by clearing currentFrame_ first.
    // setFrame then runs every update path (atom positions, each
    // overlay, camera write, render) for the current frame.
    currentFrame_ = -1;
    setFrame(t);
}

void MoleculeScene::PushAtomPositions(int t, double bounds[6]) {
    const size_t st = static_cast<size_t>(t);
    const size_t N  = protein_->atomCount();
    bounds[0] = bounds[2] = bounds[4] = +1e30;
    bounds[1] = bounds[3] = bounds[5] = -1e30;
    for (size_t i = 0; i < N; ++i) {
        const model::Vec3 p = conformation_->atomPosition(st, i);
        if (p.x() < bounds[0]) bounds[0] = p.x();
        if (p.x() > bounds[1]) bounds[1] = p.x();
        if (p.y() < bounds[2]) bounds[2] = p.y();
        if (p.y() > bounds[3]) bounds[3] = p.y();
        if (p.z() < bounds[4]) bounds[4] = p.z();
        if (p.z() > bounds[5]) bounds[5] = p.z();
        molecule_->SetAtomPosition(
            static_cast<vtkIdType>(i), p.x(), p.y(), p.z());
    }
}

void MoleculeScene::setFrame(int t) {
    ASSERT_THREAD(this);
    if (!molecule_ || !protein_ || !conformation_) return;
    if (t == currentFrame_) return;
    if (t < 0 || static_cast<size_t>(t) >= conformation_->frameCount()) return;

    QElapsedTimer timer;
    timer.start();

    // 1. Position push — one pass through atoms updates molecule
    //    positions and accumulates per-frame bounds. Per
    //    feedback_vtk_bounds_cache, we compute bounds ourselves because
    //    vtkMolecule::GetBounds() / vtkActor::GetBounds() cache from the
    //    mapper's input on first query and don't invalidate on
    //    SetAtomPosition + Modified().
    double bounds[6];
    PushAtomPositions(t, bounds);

    // 2. Modified bumps. vtkMolecule::SetAtomPosition calls
    //    Points->SetPoint(...) and Modified() on the molecule, but
    //    Points->SetPoint does NOT bump the points array's own MTime
    //    explicitly (per vtkMolecule.cxx:184-197). Some VBO-gate paths
    //    in vtkOpenGLPolyDataMapper consult the input data's Points
    //    MTime through the trivial-producer chain inside
    //    vtkOpenGLMoleculeMapper; with the molecule's MTime up but the
    //    points' MTime stale, the gate can miss when other conditions
    //    on its OR chain fall out. PROBE per spec §6.1: explicitly bump
    //    the points' MTime. If this clears the end-of-trajectory
    //    atom-render drop, the missed points-MTime check is the cause
    //    and this stays as the settled fix. If the drop persists, this
    //    line is removed and the investigation moves to
    //    notes/RESIDUAL_RENDER_DROP.md.
    molecule_->Modified();
    if (auto* points = molecule_->GetPoints())
        points->Modified();

    // 3. Move currentFrame_ BEFORE the camera write and render schedule
    //    so the EndEvent observer (Stage 8) reads the correct frame.
    currentFrame_ = t;

    // 4. Camera composer writes absolute camera state for frame t.
    //    Free mode is a no-op; non-Free modes apply the per-frame fit.
    //    Returns false on degenerate input (collinear atoms, missing
    //    indices); we log and continue with stale camera state in that
    //    case rather than teleporting.
    if (composer_) {
        const bool wrote = composer_->write(static_cast<std::size_t>(t));
        if (!wrote && composer_->mode().kind != CameraMode::Kind::Free) {
            qCWarning(cScene).noquote()
                << "camera composer reported degenerate mode at frame" << t
                << "; keeping previous camera state";
        }
    }

    // 5. Fan to overlays — each updates its own backing data.
    if (ribbon_)       ribbon_->setFrame(t);
    if (ringPolygons_) ringPolygons_->setFrame(t);
    if (fieldGrid_)    fieldGrid_->setFrame(t);
    if (bfieldStream_) bfieldStream_->setFrame(t);
    if (selection_)    selection_->setFrame(t);
    if (measurement_)  measurement_->setFrame(t);
    if (reveal_)       reveal_->setFrame(t);

    // 6. Resync near/far clipping planes from THIS FRAME's actual atom
    //    bounds (computed above), not from vtkActor::GetBounds() which
    //    stays pinned at frame-0 values. Pad each axis by 5 Å so
    //    overlays extending past the molecule (ring polygons, butterfly
    //    isosurfaces, streamlines) also stay inside the frustum.
    constexpr double pad = 5.0;
    double padded[6] = {
        bounds[0] - pad, bounds[1] + pad,
        bounds[2] - pad, bounds[3] + pad,
        bounds[4] - pad, bounds[5] + pad,
    };
    renderer_->ResetCameraClippingRange(padded);

    // 7. Schedule one render via the Qt paint chain.
    requestRender(RenderSource::Timer);

    // 8. Per-frame timing at DEBUG. Every 50 frames a diagnostic
    //    snapshot (RSS, actor count, mol bounds, visibility) ALSO at
    //    DEBUG — kept around because it caught the VTK bounds-cache
    //    bug; raise to qCInfo temporarily if a similar progressive-
    //    rendering issue recurs. See feedback_vtk_bounds_cache memory
    //    for the story.
    qCDebug(cScene).noquote()
        << "scene | frame=" << t
        << "| atoms=" << protein_->atomCount()
        << "| bounds=[" << bounds[0] << "," << bounds[1]
        << "][" << bounds[2] << "," << bounds[3]
        << "][" << bounds[4] << "," << bounds[5] << "]"
        << "| dt_ms=" << timer.elapsed();
    if (t % 50 == 0) {
        long rssKb = 0;
        if (FILE* f = std::fopen("/proc/self/statm", "r")) {
            long pages = 0;
            if (std::fscanf(f, "%ld %ld", &pages, &pages) >= 1) {
                rssKb = pages * 4;   // statm col 2 is resident pages, 4 KB each
            }
            std::fclose(f);
        }
        const int nActors = renderer_->GetActors()->GetNumberOfItems();
        const int molVis  = actor_->GetVisibility();
        qCDebug(cScene).noquote()
            << "snapshot @ frame" << t
            << "| rss=" << rssKb << "KB"
            << "| actors=" << nActors
            << "| mol vis=" << molVis
            << "| atom bounds=[" << bounds[0] << "," << bounds[1]
            << "][" << bounds[2] << "," << bounds[3]
            << "][" << bounds[4] << "," << bounds[5] << "]";
    }
}

bool MoleculeScene::lockCameraToSelectionPlane(const std::vector<std::size_t>& atoms) {
    ASSERT_THREAD(this);
    if (!composer_ || atoms.size() != 3) return false;
    const std::size_t frame = currentFrame_ >= 0
        ? static_cast<std::size_t>(currentFrame_) : 0;
    composer_->setMode(PlaneMode(atoms[0], atoms[1], atoms[2]),
                        DefaultPolicy(), frame);
    if (composer_->mode().kind != CameraMode::Kind::Plane) {
        // Composer rejected the lock (degenerate geometry / out-of-range
        // atoms). Mirror the original return path.
        emit cameraPlaneLockChanged(false);
        return false;
    }
    // Apply the per-frame fit immediately so the camera state is
    // updated on this tick (the original lockCameraToSelectionPlane
    // wrote the camera in the enable path; consumers that GET the
    // camera right after enable should see the locked state, not a
    // stale free-camera state).
    (void)composer_->write(frame);
    // Refresh clipping from current bounds — the camera moved, so the
    // stale clipping range may now cut into the molecule.
    if (renderer_) {
        double bounds[6];
        if (protein_) {
            PushAtomPositions(static_cast<int>(frame), bounds);
            constexpr double pad = 5.0;
            double padded[6] = {
                bounds[0] - pad, bounds[1] + pad,
                bounds[2] - pad, bounds[3] + pad,
                bounds[4] - pad, bounds[5] + pad,
            };
            renderer_->ResetCameraClippingRange(padded);
        }
    }
    requestRender(RenderSource::External);
    qCInfo(cScene).noquote()
        << "camera plane lock enabled | atoms="
        << atoms[0] << atoms[1] << atoms[2];
    emit cameraPlaneLockChanged(true);
    return true;
}

void MoleculeScene::clearCameraPlaneLock() {
    ASSERT_THREAD(this);
    if (!composer_) return;
    if (composer_->mode().kind != CameraMode::Kind::Plane) return;
    const std::size_t frame = currentFrame_ >= 0
        ? static_cast<std::size_t>(currentFrame_) : 0;
    composer_->setMode(FreeMode(), DefaultPolicy(), frame);
    requestRender(RenderSource::External);
    qCInfo(cScene).noquote() << "camera plane lock disabled";
    emit cameraPlaneLockChanged(false);
}

bool MoleculeScene::isCameraPlaneLocked() const {
    return composer_ && composer_->mode().kind == CameraMode::Kind::Plane;
}

std::vector<std::size_t> MoleculeScene::cameraPlaneLockAtoms() const {
    if (!composer_ || composer_->mode().kind != CameraMode::Kind::Plane)
        return {};
    return composer_->mode().atoms;
}

void MoleculeScene::revealBinding(const model::SignalBinding& binding) {
    ASSERT_THREAD(this);
    if (!reveal_ || !protein_ || !conformation_)
        return;

    if (activeRevealBinding_ && SameRevealBinding(*activeRevealBinding_, binding)) {
        clearReveal();
        return;
    }

    const int frame = currentFrame_ >= 0 ? currentFrame_ : 0;
    reveal_->reveal(binding, frame);
    if (reveal_->isActive()) {
        activeRevealBinding_ = binding;
        focusCameraOnReveal(binding, reveal_->activeAtoms(), frame);
    } else {
        activeRevealBinding_.reset();
    }
    requestRender(RenderSource::Reveal);
}

void MoleculeScene::clearReveal() {
    ASSERT_THREAD(this);
    if (!reveal_)
        return;
    activeRevealBinding_.reset();
    reveal_->clear();
    requestRender(RenderSource::Reveal);
}

void MoleculeScene::focusCameraOnReveal(const model::SignalBinding& binding,
                                        const std::vector<std::size_t>& atoms,
                                        int frame) {
    if (!composer_ || !protein_ || !conformation_ || atoms.empty()) return;
    if (frame < 0 || static_cast<std::size_t>(frame) >= conformation_->frameCount())
        return;

    // Per spec §3.2: dihedral reveal lifts to sustained
    // CameraMode::Dihedral + DownAxis policy so the sight-down view
    // persists frame-to-frame rather than being a one-shot. Other
    // anchors get CameraMode::Subset with the active atoms (focal at
    // centroid; sight inherited from current camera).
    const auto* tuple = std::get_if<model::AtomTupleAnchor>(&binding.anchor);
    const bool canSightDown = tuple && tuple->atoms.size() >= 4
        && tuple->atoms[0] < protein_->atomCount()
        && tuple->atoms[1] < protein_->atomCount()
        && tuple->atoms[2] < protein_->atomCount()
        && tuple->atoms[3] < protein_->atomCount();

    if (canSightDown) {
        composer_->setMode(
            DihedralMode(tuple->atoms[0], tuple->atoms[1],
                          tuple->atoms[2], tuple->atoms[3]),
            DownAxisPolicy(tuple->atoms[1], tuple->atoms[2]),
            static_cast<std::size_t>(frame));
    } else if (atoms.size() >= 3) {
        composer_->setMode(SubsetMode(atoms), DefaultPolicy(),
                            static_cast<std::size_t>(frame));
    } else if (atoms.size() == 1) {
        composer_->setMode(AtomMode(atoms[0]), DefaultPolicy(),
                            static_cast<std::size_t>(frame));
    } else if (atoms.size() == 2) {
        composer_->setMode(BondMode(atoms[0], atoms[1]), DefaultPolicy(),
                            static_cast<std::size_t>(frame));
    }
    // No requestRender — revealBinding caller already issues one.
}

}  // namespace h5reader::app
