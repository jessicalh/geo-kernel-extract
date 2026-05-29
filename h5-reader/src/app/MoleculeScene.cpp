#include "MoleculeScene.h"

#include "QtBackboneRibbonOverlay.h"
#include "QtBFieldStreamOverlay.h"
#include "QtFieldGridOverlay.h"
#include "QtRingPolygonOverlay.h"
#include "QtSelectionOverlay.h"
#include "MeasurementOverlay.h"
#include "SceneRevealOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/Types.h"

#include <QElapsedTimer>
#include <QLoggingCategory>

#include <vtkActorCollection.h>

#include <cstdio>

#include <vtkCallbackCommand.h>
#include <vtkCamera.h>
#include <vtkCommand.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNew.h>
#include <vtkProperty.h>
#include <vtkRenderWindowInteractor.h>

#include <algorithm>
#include <cmath>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cScene, "h5reader.scene")

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

MoleculeScene::MoleculeScene(vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
                             QObject* parent)
    : QObject(parent),
      renderWindow_(std::move(renderWindow))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("MoleculeScene"));

    // ----- Renderer setup ------------------------------------------------
    //
    // The existing nmr-viewer (ui/src/MainWindow.cpp:217-231) uses
    // FXAA + NO depth peeling after hitting translucency artifacts on
    // AMD hardware. We keep that choice for now because it's proven on
    // this codebase's target hardware. If a translucent overlay shows
    // sorting artifacts, flip to depth peeling per the qt6-cpp skill's
    // references/3d-vtk.md and log the change.

    renderer_ = vtkSmartPointer<vtkRenderer>::New();
    renderer_->SetBackground(1.0, 1.0, 1.0);
    renderer_->SetUseFXAA(true);
    renderer_->SetUseDepthPeeling(0);

    renderWindow_->AddRenderer(renderer_);
    renderWindow_->SetAlphaBitPlanes(1);
    renderWindow_->SetMultiSamples(0);   // MSAA off — incompatible with
                                         // translucency; FXAA handles AA.

    if (auto* iren = renderWindow_->GetInteractor()) {
        vtkNew<vtkInteractorStyleTrackballCamera> style;
        iren->SetInteractorStyle(style);
    }

    // Render-time observer: vtkRenderer::GetLastRenderTimeInSeconds()
    // reports CPU wall-time spent in the last Render() call (includes
    // any lazy filter Update() that runs during render, plus GPU dispatch
    // overhead). Logged after every EndEvent of the render window so the
    // log stream gives BS-overlay-pipe-time alongside total render time.
    // Diagnostic during Windows-vs-Linux perf investigation; the
    // callback is cheap (one double read + one log line) and can be
    // left in.
    auto renderTimeCb = vtkSmartPointer<vtkCallbackCommand>::New();
    renderTimeCb->SetClientData(this);
    renderTimeCb->SetCallback(
        [](vtkObject* /*caller*/, unsigned long, void* clientData, void*) {
            auto* self = static_cast<MoleculeScene*>(clientData);
            if (!self || !self->renderer_) return;
            const double ms = self->renderer_->GetLastRenderTimeInSeconds() * 1000.0;
            qCInfo(cScene).noquote() << "render"
                                      << QString::number(ms, 'f', 1) << "ms";
        });
    renderWindow_->AddObserver(vtkCommand::EndEvent, renderTimeCb);

    qCInfo(cScene).noquote()
        << "Renderer initialised: FXAA on, depth peeling OFF, AlphaBitPlanes=1,"
        << "MSAA=0, style=vtkInteractorStyleTrackballCamera";
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
    if (cameraPlaneLock_) {
        cameraPlaneLock_.reset();
        emit cameraPlaneLockChanged(false);
    }

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

    if (!measurement_) {
        measurement_ = new MeasurementOverlay(renderer_, renderWindow_, this);
    }
    measurement_->Build(protein, conformation);

    if (!reveal_) {
        reveal_ = new SceneRevealOverlay(renderer_, this);
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

    // Seed the camera-follow baseline at frame 0. Every subsequent
    // setFrame() will shift the camera by the centroid delta.
    if (conformation_ && conformation_->frameCount() > 0) {
        lastCentroid_     = ComputeCentroid(0);
        haveLastCentroid_ = true;
    }

    renderWindow_->Render();
}

void MoleculeScene::requestRender() {
    ASSERT_THREAD(this);
    if (renderWindow_) renderWindow_->Render();
}

void MoleculeScene::refreshCurrentFrame() {
    ASSERT_THREAD(this);
    const int t = currentFrame_;
    // Bypass the same-frame guard by clearing currentFrame_ first.
    // setFrame then runs every update path (atom positions, each
    // overlay, camera follow, render) for the current frame.
    currentFrame_ = -1;
    setFrame(t);
}

model::Vec3 MoleculeScene::ComputeCentroid(size_t tIndex) const {
    if (!protein_ || !conformation_) return model::Vec3::Zero();
    const size_t N = protein_->atomCount();
    if (N == 0) return model::Vec3::Zero();
    model::Vec3 sum = model::Vec3::Zero();
    for (size_t i = 0; i < N; ++i) sum += conformation_->atomPosition(tIndex, i);
    return sum / static_cast<double>(N);
}

std::optional<math::PlaneFrame> MoleculeScene::computePlaneFrame(
    const std::vector<std::size_t>& atoms,
    std::size_t frame) const {
    if (!protein_ || !conformation_ || atoms.size() != 3)
        return std::nullopt;
    if (frame >= conformation_->frameCount())
        return std::nullopt;
    for (std::size_t atom : atoms) {
        if (atom >= protein_->atomCount())
            return std::nullopt;
    }
    return math::computePlaneFrame({
        conformation_->atomPosition(frame, atoms[0]),
        conformation_->atomPosition(frame, atoms[1]),
        conformation_->atomPosition(frame, atoms[2]),
    });
}

bool MoleculeScene::lockCameraToSelectionPlane(const std::vector<std::size_t>& atoms) {
    ASSERT_THREAD(this);
    if (!renderer_ || !conformation_ || atoms.size() != 3) {
        clearCameraPlaneLock();
        return false;
    }

    const std::size_t frame = currentFrame_ >= 0 ? static_cast<std::size_t>(currentFrame_) : 0;
    const auto basis = computePlaneFrame(atoms, frame);
    auto* camera = renderer_->GetActiveCamera();
    if (!basis || !camera) {
        clearCameraPlaneLock();
        return false;
    }

    double posRaw[3];
    double upRaw[3];
    double directionRaw[3];
    camera->GetPosition(posRaw);
    camera->GetViewUp(upRaw);
    camera->GetDirectionOfProjection(directionRaw);
    const model::Vec3 position(posRaw[0], posRaw[1], posRaw[2]);
    model::Vec3 viewUp(upRaw[0], upRaw[1], upRaw[2]);
    model::Vec3 viewDirection(directionRaw[0], directionRaw[1], directionRaw[2]);
    if (viewDirection.norm() < 1e-6)
        viewDirection = basis->origin - position;
    if (viewDirection.norm() < 1e-6)
        viewDirection = basis->z;
    viewDirection.normalize();

    CameraPlaneLock lock;
    lock.atoms = atoms;
    lock.normalSign = viewDirection.dot(basis->z) < 0.0 ? -1.0 : 1.0;

    const model::Vec3 lockedNormal = basis->z * lock.normalSign;
    viewUp -= viewUp.dot(lockedNormal) * lockedNormal;
    if (viewUp.norm() < 1e-6)
        viewUp = basis->y;
    viewUp.normalize();

    lock.localViewUp = model::Vec3(viewUp.dot(basis->x),
                                   viewUp.dot(basis->y),
                                   0.0);
    if (lock.localViewUp.norm() < 1e-6)
        lock.localViewUp = model::Vec3(0.0, 1.0, 0.0);
    lock.distance = std::max(1.0, (position - basis->origin).norm());
    cameraPlaneLock_ = std::move(lock);

    if (!applyCameraPlaneLock(frame)) {
        clearCameraPlaneLock();
        return false;
    }

    renderer_->ResetCameraClippingRange();
    if (renderWindow_)
        renderWindow_->Render();

    qCInfo(cScene).noquote() << "camera plane lock enabled | atoms="
                              << atoms[0] << atoms[1] << atoms[2];
    emit cameraPlaneLockChanged(true);
    return true;
}

void MoleculeScene::clearCameraPlaneLock() {
    ASSERT_THREAD(this);
    if (!cameraPlaneLock_)
        return;
    cameraPlaneLock_.reset();
    qCInfo(cScene).noquote() << "camera plane lock disabled";
    emit cameraPlaneLockChanged(false);
}

bool MoleculeScene::isCameraPlaneLocked() const {
    return cameraPlaneLock_.has_value();
}

std::vector<std::size_t> MoleculeScene::cameraPlaneLockAtoms() const {
    return cameraPlaneLock_ ? cameraPlaneLock_->atoms : std::vector<std::size_t>{};
}

bool MoleculeScene::applyCameraPlaneLock(std::size_t frame) {
    if (!cameraPlaneLock_ || !renderer_)
        return false;

    const auto basis = computePlaneFrame(cameraPlaneLock_->atoms, frame);
    auto* camera = renderer_->GetActiveCamera();
    if (!basis || !camera)
        return false;

    auto vectorToWorld = [&basis](const model::Vec3& localVector) {
        return basis->x * localVector.x()
               + basis->y * localVector.y()
               + basis->z * localVector.z();
    };

    // Per-frame normal-sign continuity guard. The basis is rebuilt each
    // frame from the natural cross product (b-a)×(c-a); that direction
    // can flip sign across a near-degenerate configuration (ring flip,
    // third atom crossing the line through the first two). Without the
    // guard the camera would teleport to the other side of the plane.
    // Compare the candidate view direction against the previous frame's
    // direction and, if they disagree, flip the stored normalSign so
    // subsequent frames stay continuous.
    model::Vec3 viewDirection = basis->z * (cameraPlaneLock_->normalSign < 0.0 ? -1.0 : 1.0);
    if (viewDirection.norm() < 1e-6)
        return false;
    viewDirection.normalize();

    if (cameraPlaneLock_->lastDirection
        && viewDirection.dot(*cameraPlaneLock_->lastDirection) < 0.0) {
        cameraPlaneLock_->normalSign *= -1.0;
        viewDirection = -viewDirection;
    }
    cameraPlaneLock_->lastDirection = viewDirection;

    model::Vec3 viewUp = vectorToWorld(cameraPlaneLock_->localViewUp);
    viewUp -= viewUp.dot(viewDirection) * viewDirection;
    if (viewUp.norm() < 1e-6) {
        viewUp = basis->y - basis->y.dot(viewDirection) * viewDirection;
    }
    if (viewUp.norm() < 1e-6) {
        viewUp = basis->x - basis->x.dot(viewDirection) * viewDirection;
    }
    if (viewUp.norm() < 1e-6)
        return false;
    viewUp.normalize();

    const double distance = std::max(1.0, cameraPlaneLock_->distance);
    const model::Vec3 position = basis->origin - viewDirection * distance;
    camera->SetFocalPoint(basis->origin.x(), basis->origin.y(), basis->origin.z());
    camera->SetPosition(position.x(), position.y(), position.z());
    camera->SetViewUp(viewUp.x(), viewUp.y(), viewUp.z());
    camera->OrthogonalizeViewUp();
    return true;
}

void MoleculeScene::setFrame(int t) {
    ASSERT_THREAD(this);
    if (!molecule_ || !protein_ || !conformation_) return;
    if (t == currentFrame_) return;
    if (t < 0 || static_cast<size_t>(t) >= conformation_->frameCount()) return;

    QElapsedTimer timer;
    timer.start();

    model::Conformation* conf = conformation_;
    const size_t st = static_cast<size_t>(t);
    const size_t N = protein_->atomCount();

    // Update atom positions AND accumulate centroid + bounds in one
    // pass so we don't touch frame.position(i) more than once.
    //
    // Why we compute bounds ourselves: vtkMolecule::GetBounds() /
    // vtkActor::GetBounds() cache from the mapper's input on first
    // query and don't invalidate on SetAtomPosition + Modified().
    // Observed: across frames 0..600 on 1B1V_4292 the actor bounds
    // stay pinned at their frame-0 values (diagnostic snapshot log).
    // This makes renderer_->ResetCameraClippingRange() compute near/far
    // from stale geometry, so as the protein diffuses across the MD
    // simulation box over 25 ns the actual atoms drift outside the
    // clip range and disappear progressively. Computing bounds here
    // from the authoritative per-frame positions + passing them
    // explicitly to ResetCameraClippingRange(double[6]) bypasses the
    // cache entirely.
    model::Vec3 sum = model::Vec3::Zero();
    double bounds[6] = { +1e30, -1e30, +1e30, -1e30, +1e30, -1e30 };
    for (size_t i = 0; i < N; ++i) {
        const model::Vec3 p = conf->atomPosition(st, i);
        sum += p;
        if (p.x() < bounds[0]) bounds[0] = p.x();
        if (p.x() > bounds[1]) bounds[1] = p.x();
        if (p.y() < bounds[2]) bounds[2] = p.y();
        if (p.y() > bounds[3]) bounds[3] = p.y();
        if (p.z() < bounds[4]) bounds[4] = p.z();
        if (p.z() > bounds[5]) bounds[5] = p.z();
        molecule_->SetAtomPosition(
            static_cast<vtkIdType>(i), p.x(), p.y(), p.z());
    }
    model::Vec3 centroid = model::Vec3::Zero();
    if (N > 0) centroid = sum / static_cast<double>(N);

    bool planeLockApplied = false;
    if (cameraPlaneLock_) {
        planeLockApplied = applyCameraPlaneLock(st);
        if (!planeLockApplied) {
            qCWarning(cScene).noquote()
                << "camera plane lock dropped: selected atoms no longer define a stable plane";
            clearCameraPlaneLock();
        }
    }

    // Camera-follow: translate focal point and camera position by the
    // delta between this frame's centroid and the previous one unless
    // the explicit three-atom plane lock owns the camera this frame.
    if (!planeLockApplied && haveLastCentroid_) {
        const model::Vec3 delta = centroid - lastCentroid_;
        if (delta.norm() > 0.0) {
            auto* camera = renderer_->GetActiveCamera();
            double fp[3]; camera->GetFocalPoint(fp);
            double pos[3]; camera->GetPosition(pos);
            camera->SetFocalPoint(fp[0] + delta.x(),
                                  fp[1] + delta.y(),
                                  fp[2] + delta.z());
            camera->SetPosition(pos[0] + delta.x(),
                                pos[1] + delta.y(),
                                pos[2] + delta.z());
            // Clipping-plane resync is owned by the explicit-bounds
            // call below (after position + overlay updates). The old
            // zero-arg ResetCameraClippingRange() lived here and relied
            // on vtkActor::GetBounds, whose cache does not invalidate on
            // SetAtomPosition + Modified — see feedback_vtk_bounds_cache.
            // One call with real per-frame bounds is the clean version.
        }
    }
    lastCentroid_     = centroid;
    haveLastCentroid_ = true;

    molecule_->Modified();

    // Explicitly mark the composite mapper modified.
    //
    // DIAGNOSTIC PROBE, not a settled fix. vtkOpenGLMoleculeMapper builds
    // internal sphere- and cylinder-imposter mappers from the vtkMolecule
    // it was given via SetInputData. In principle, molecule_->Modified()
    // should propagate through the mapper's input chain and force the
    // internal mappers to re-upload their VBOs. In practice we have one
    // residual intermittent end-of-trajectory atom-render drop (overlays
    // render, spheres do not) that could be explained by the composite
    // chain occasionally missing the re-upload.
    //
    // If this line eliminates the drop entirely, the Modified() chain
    // through the composite mapper is broken under some condition — the
    // band-aid works, but "why does it need a band-aid?" is the real
    // question and should be answered before viva. Likely culprits:
    //   (a) a VTK bug in vtkOpenGLMoleculeMapper's Update-on-input path,
    //   (b) an ordering issue between SetAtomPosition and Render on
    //       the GUI thread when a GPU context stall coincides with a
    //       pending overlay update,
    //   (c) our own omission of Modified() on a different upstream
    //       object (e.g., vtkPoints inside the molecule).
    //
    // If the drop persists after this change, the probe is noise —
    // remove and pursue hypotheses C / D in notes/RESIDUAL_RENDER_DROP.md.
    //
    // Either way, leave the UDP log running on long playback; the
    // per-50-frame snapshot already flags bounds/actors/visibility.
    mapper_->Modified();

    // Propagate to overlays BEFORE the render, so the single Render()
    // below picks up everyone's modified data in one pass.
    if (ribbon_)       ribbon_->setFrame(t);
    if (ringPolygons_) ringPolygons_->setFrame(t);
    if (fieldGrid_)    fieldGrid_->setFrame(t);
    if (bfieldStream_) bfieldStream_->setFrame(t);
    if (selection_)    selection_->setFrame(t);
    if (measurement_)  measurement_->setFrame(t);
    if (reveal_)       reveal_->setFrame(t);

    // Resync near/far clipping planes from THIS FRAME's actual atom
    // bounds (computed above), not from vtkActor::GetBounds() which
    // stays pinned at frame-0 values. Pad each axis by 5 Å so
    // overlays extending past the molecule (ring polygons, butterfly
    // isosurfaces, streamlines) also stay inside the frustum.
    if (N > 0) {
        constexpr double pad = 5.0;
        double padded[6] = {
            bounds[0] - pad, bounds[1] + pad,
            bounds[2] - pad, bounds[3] + pad,
            bounds[4] - pad, bounds[5] + pad,
        };
        renderer_->ResetCameraClippingRange(padded);
    }

    renderWindow_->Render();

    currentFrame_ = t;

    // Per-frame timing at DEBUG. Every 50 frames a diagnostic snapshot
    // (RSS, actor count, mol bounds, visibility) ALSO at DEBUG
    // — kept around because it caught the VTK bounds-cache bug; raise
    // to qCInfo temporarily if a similar progressive-rendering issue
    // recurs. See feedback_vtk_bounds_cache memory for the story.
    //
    // The bounds reported are the per-frame atom-position bounds
    // computed earlier in this function, NOT actor_->GetBounds() —
    // the actor's bounds cache is pinned to frame 0 (that's the bug
    // this snapshot exists to catch). Using the live values lets the
    // snapshot show actual frame-to-frame motion.
    qCDebug(cScene).noquote()
        << "frame" << t << "applied |" << timer.elapsed() << "ms";
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
    requestRender();
}

void MoleculeScene::clearReveal() {
    ASSERT_THREAD(this);
    if (!reveal_)
        return;
    activeRevealBinding_.reset();
    reveal_->clear();
    requestRender();
}

void MoleculeScene::focusCameraOnReveal(const model::SignalBinding& binding,
                                        const std::vector<std::size_t>& atoms,
                                        int frame) {
    if (!renderer_ || !protein_ || !conformation_ || atoms.empty())
        return;
    if (frame < 0 || static_cast<std::size_t>(frame) >= conformation_->frameCount())
        return;

    auto* camera = renderer_->GetActiveCamera();
    if (!camera)
        return;

    auto atomPosition = [this, frame](std::size_t atom) {
        return conformation_->atomPosition(static_cast<std::size_t>(frame), atom);
    };

    model::Vec3 center = model::Vec3::Zero();
    for (std::size_t atom : atoms)
        center += atomPosition(atom);
    center /= static_cast<double>(atoms.size());

    double oldFpRaw[3];
    double oldPosRaw[3];
    camera->GetFocalPoint(oldFpRaw);
    camera->GetPosition(oldPosRaw);
    const model::Vec3 oldFp(oldFpRaw[0], oldFpRaw[1], oldFpRaw[2]);
    const model::Vec3 oldPos(oldPosRaw[0], oldPosRaw[1], oldPosRaw[2]);
    const double distance = std::max(12.0, (oldPos - oldFp).norm());

    const auto* tuple = std::get_if<model::AtomTupleAnchor>(&binding.anchor);
    const bool canLookDownDihedral = tuple && tuple->atoms.size() >= 4
                                      && tuple->atoms[1] < protein_->atomCount()
                                      && tuple->atoms[2] < protein_->atomCount();

    if (canLookDownDihedral) {
        const model::Vec3 b = atomPosition(tuple->atoms[1]);
        const model::Vec3 c = atomPosition(tuple->atoms[2]);
        model::Vec3 axis = c - b;
        if (axis.norm() > 1e-9) {
            axis.normalize();
            center = (b + c) * 0.5;

            model::Vec3 currentView = oldFp - oldPos;
            if (currentView.norm() > 1e-9)
                currentView.normalize();
            const model::Vec3 viewDir = currentView.dot(axis) >= currentView.dot(-axis) ? axis : -axis;

            camera->SetFocalPoint(center.x(), center.y(), center.z());
            const model::Vec3 pos = center - viewDir * distance;
            camera->SetPosition(pos.x(), pos.y(), pos.z());

            double oldUpRaw[3];
            camera->GetViewUp(oldUpRaw);
            model::Vec3 up(oldUpRaw[0], oldUpRaw[1], oldUpRaw[2]);
            up -= up.dot(viewDir) * viewDir;
            if (up.norm() < 1e-6) {
                up = model::Vec3(0.0, 0.0, 1.0);
                up -= up.dot(viewDir) * viewDir;
            }
            if (up.norm() < 1e-6) {
                up = model::Vec3(0.0, 1.0, 0.0);
                up -= up.dot(viewDir) * viewDir;
            }
            if (up.norm() > 1e-6) {
                up.normalize();
                camera->SetViewUp(up.x(), up.y(), up.z());
            }
            camera->OrthogonalizeViewUp();
            renderer_->ResetCameraClippingRange();
            return;
        }
    }

    const model::Vec3 delta = center - oldFp;
    camera->SetFocalPoint(center.x(), center.y(), center.z());
    camera->SetPosition(oldPos.x() + delta.x(),
                        oldPos.y() + delta.y(),
                        oldPos.z() + delta.z());
    renderer_->ResetCameraClippingRange();
}

}  // namespace h5reader::app
