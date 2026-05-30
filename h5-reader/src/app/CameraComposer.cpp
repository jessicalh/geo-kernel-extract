#include "CameraComposer.h"

#include "FitTargetMath.h"
#include "PlaneFrameMath.h"

#include "../diagnostics/ObjectCensus.h"

#include <QLoggingCategory>

#include <vtkCamera.h>

#include <algorithm>
#include <array>
#include <cmath>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cComposer, "h5reader.composer")

// Rotate a vector v around a unit axis by angle (rad), Rodrigues form.
// Used for azimuth / elevation / roll deltas applied on top of the fit.
model::Vec3 RodriguesRotate(const model::Vec3& v,
                             const model::Vec3& axis,
                             double angleRad) {
    const double c = std::cos(angleRad);
    const double s = std::sin(angleRad);
    return v * c + axis.cross(v) * s + axis * axis.dot(v) * (1.0 - c);
}
}  // namespace

CameraComposer::CameraComposer(vtkSmartPointer<vtkRenderer>  renderer,
                                const model::QtProtein*       protein,
                                model::Conformation*          conformation,
                                QObject*                      parent)
    : QObject(parent),
      renderer_(std::move(renderer)),
      protein_(protein),
      conformation_(conformation) {
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("CameraComposer"));
}

CameraComposer::~CameraComposer() = default;

void CameraComposer::setMode(CameraMode mode,
                              OrientationPolicy policy,
                              std::size_t currentFrame) {
    ASSERT_THREAD(this);

    // Reset accumulated user deltas — each lock acquisition is a fresh
    // start (agent decision per the implementation prompt §4-b).
    accumAzimuthRad_   = 0.0;
    accumElevationRad_ = 0.0;
    accumRollRad_      = 0.0;
    accumPan_          = model::Vec3::Zero();
    accumDolly_        = 1.0;

    // Capture the new mode + policy. captureInitialState uses currentFrame
    // and the new mode to compute things like distance and subset
    // reference positions.
    mode_   = std::move(mode);
    policy_ = policy;

    planeLocalViewUp_  = model::Vec3::Zero();
    planeNormalSign_   = 1.0;
    planeLastDirection_.reset();
    subsetReference_.clear();

    if (!validateAtomsForCurrentMode()) {
        qCWarning(cComposer).noquote() << "atoms out of range for mode" << NameFor(mode_.kind)
                                        << "; falling back to free";
        mode_ = FreeMode();
        emit modeChanged();
        return;
    }
    if (!conformation_ || conformation_->frameCount() == 0) {
        emit modeChanged();
        return;
    }
    if (currentFrame >= conformation_->frameCount())
        currentFrame = 0;

    captureInitialState(currentFrame);

    qCInfo(cComposer).noquote() << "camera mode" << NameFor(mode_.kind)
                                 << "| policy=" << NameFor(policy_.kind)
                                 << "| atoms=" << static_cast<int>(mode_.atoms.size())
                                 << "| ref_frame=" << static_cast<qlonglong>(currentFrame);
    emit modeChanged();
}

void CameraComposer::captureInitialState(std::size_t referenceFrame) {
    auto* camera = renderer_ ? renderer_->GetActiveCamera() : nullptr;
    if (!camera) return;

    double posRaw[3];
    double fpRaw[3];
    double upRaw[3];
    camera->GetPosition(posRaw);
    camera->GetFocalPoint(fpRaw);
    camera->GetViewUp(upRaw);
    const model::Vec3 pos(posRaw[0], posRaw[1], posRaw[2]);
    const model::Vec3 fp(fpRaw[0], fpRaw[1], fpRaw[2]);
    const model::Vec3 up(upRaw[0], upRaw[1], upRaw[2]);

    // Distance from focal to position — used by Atom/Bond/Plane modes to
    // hold the camera at the same zoom across the lock acquisition.
    distance_ = std::max(1.0, (pos - fp).norm());

    if (mode_.kind == CameraMode::Kind::Plane && mode_.atoms.size() == 3) {
        // Re-derive plane basis + plane-local viewUp + initial sign from
        // the captured camera state. Matches the original
        // lockCameraToSelectionPlane logic at MoleculeScene.cpp:282-343.
        const std::array<model::Vec3, 3> positions{{
            conformation_->atomPosition(referenceFrame, mode_.atoms[0]),
            conformation_->atomPosition(referenceFrame, mode_.atoms[1]),
            conformation_->atomPosition(referenceFrame, mode_.atoms[2]),
        }};
        const auto basis = math::computePlaneFrame(positions);
        if (!basis) return;

        // View direction: prefer the camera's current sight; fall back to
        // (origin - position) if the camera has no useful direction.
        double dirRaw[3];
        camera->GetDirectionOfProjection(dirRaw);
        model::Vec3 viewDir(dirRaw[0], dirRaw[1], dirRaw[2]);
        if (viewDir.norm() < 1e-6)
            viewDir = basis->origin - pos;
        if (viewDir.norm() < 1e-6)
            viewDir = basis->z;
        viewDir.normalize();

        planeNormalSign_ = viewDir.dot(basis->z) < 0.0 ? -1.0 : 1.0;
        const model::Vec3 lockedNormal = basis->z * planeNormalSign_;

        model::Vec3 viewUp = up - up.dot(lockedNormal) * lockedNormal;
        if (viewUp.norm() < 1e-6)
            viewUp = basis->y;
        viewUp.normalize();

        planeLocalViewUp_ = model::Vec3(viewUp.dot(basis->x),
                                         viewUp.dot(basis->y),
                                         0.0);
        if (planeLocalViewUp_.norm() < 1e-6)
            planeLocalViewUp_ = model::Vec3(0.0, 1.0, 0.0);
    } else if (mode_.kind == CameraMode::Kind::Subset && mode_.atoms.size() >= 3) {
        subsetReference_.reserve(mode_.atoms.size());
        for (std::size_t a : mode_.atoms) {
            if (a >= protein_->atomCount()) continue;
            subsetReference_.push_back(conformation_->atomPosition(referenceFrame, a));
        }
    }
}

bool CameraComposer::validateAtomsForCurrentMode() const {
    if (mode_.kind == CameraMode::Kind::Free) return true;
    if (!protein_) return false;
    const std::size_t atomCount = protein_->atomCount();
    const std::size_t expected = ExpectedAtomCount(mode_.kind);
    if (expected != static_cast<std::size_t>(-1) && mode_.atoms.size() != expected)
        return false;
    if (mode_.kind == CameraMode::Kind::Subset && mode_.atoms.size() < 3)
        return false;
    for (std::size_t a : mode_.atoms) {
        if (a >= atomCount) return false;
    }
    return true;
}

std::vector<model::Vec3> CameraComposer::readAtomPositions(std::size_t t) const {
    std::vector<model::Vec3> out;
    if (!conformation_ || t >= conformation_->frameCount())
        return out;
    out.reserve(mode_.atoms.size());
    for (std::size_t a : mode_.atoms) {
        if (a >= (protein_ ? protein_->atomCount() : 0)) return {};
        out.push_back(conformation_->atomPosition(t, a));
    }
    return out;
}

bool CameraComposer::write(std::size_t t) {
    ASSERT_THREAD(this);
    if (!renderer_ || !renderer_->GetActiveCamera()) return false;
    if (!conformation_ || t >= conformation_->frameCount()) return false;

    switch (mode_.kind) {
        case CameraMode::Kind::Free:     return writeFree(t);
        case CameraMode::Kind::Atom:     return writeAtom(t);
        case CameraMode::Kind::Bond:     return writeBond(t);
        case CameraMode::Kind::Dihedral: return writeDihedral(t);
        case CameraMode::Kind::Plane:    return writePlane(t);
        case CameraMode::Kind::Subset:   return writeSubset(t);
    }
    return false;
}

bool CameraComposer::writeFree(std::size_t /*t*/) {
    // Free camera state is owned by accumulated user gestures. On every
    // tick we re-write the camera to (captured initial) + (delta) so
    // that the camera stays put without the trackball interfering.
    auto* camera = renderer_->GetActiveCamera();
    if (!camera) return false;
    double fpRaw[3];
    double posRaw[3];
    double upRaw[3];
    camera->GetFocalPoint(fpRaw);
    camera->GetPosition(posRaw);
    camera->GetViewUp(upRaw);
    model::Vec3 fp(fpRaw[0], fpRaw[1], fpRaw[2]);
    model::Vec3 pos(posRaw[0], posRaw[1], posRaw[2]);
    model::Vec3 up(upRaw[0], upRaw[1], upRaw[2]);
    // For Free mode, the deltas were already applied at applyGesture
    // time directly to the camera. Per-frame re-application would
    // double-count them. Free mode is a strict no-op here; the per-frame
    // path is left intact for the molecule's intrinsic motion.
    (void)fp; (void)pos; (void)up;
    return true;
}

bool CameraComposer::writeAtom(std::size_t t) {
    if (mode_.atoms.size() != 1) return false;
    const auto positions = readAtomPositions(t);
    if (positions.size() != 1) return false;

    const std::array<model::Vec3, 1> arr{{positions[0]}};
    auto anchor = math::ComputeAtomAnchor(arr);
    if (!anchor) return false;

    auto* camera = renderer_->GetActiveCamera();
    if (!camera) return false;
    double posRaw[3];
    double upRaw[3];
    camera->GetPosition(posRaw);
    camera->GetViewUp(upRaw);
    const model::Vec3 pos(posRaw[0], posRaw[1], posRaw[2]);
    const model::Vec3 up(upRaw[0], upRaw[1], upRaw[2]);

    // Inherit the view direction from the current camera; the focal
    // moves to the atom, the position is at the same distance along the
    // same direction. This is the no-orientation-override default for
    // Atom mode (per spec §2.3.3).
    double fpRaw[3];
    camera->GetFocalPoint(fpRaw);
    const model::Vec3 oldFp(fpRaw[0], fpRaw[1], fpRaw[2]);
    model::Vec3 viewDir = oldFp - pos;
    if (viewDir.norm() < 1e-9)
        viewDir = model::Vec3(0.0, 0.0, -1.0);
    else
        viewDir.normalize();

    const model::Vec3 newFp  = anchor->focal;
    const model::Vec3 newPos = newFp - viewDir * distance_;

    writeCameraComposed(newFp, newPos, up);
    return true;
}

bool CameraComposer::writeBond(std::size_t t) {
    if (mode_.atoms.size() != 2) return false;
    const auto positions = readAtomPositions(t);
    if (positions.size() != 2) return false;

    const std::array<model::Vec3, 2> arr{{positions[0], positions[1]}};
    auto anchor = math::ComputeBondAnchor(arr);
    if (!anchor || !anchor->axis) return false;

    auto* camera = renderer_->GetActiveCamera();
    if (!camera) return false;
    double posRaw[3];
    double upRaw[3];
    camera->GetPosition(posRaw);
    camera->GetViewUp(upRaw);
    const model::Vec3 pos(posRaw[0], posRaw[1], posRaw[2]);
    model::Vec3 up(upRaw[0], upRaw[1], upRaw[2]);

    // Default orientation for Bond: ViewUp parallel to the bond axis
    // (sight perpendicular to the bond). Resolve to the natural pair via
    // OrientationPolicy::Default; explicit Free leaves the gesture's up.
    const model::Vec3 bondAxis = *anchor->axis;
    const model::Vec3 fp       = anchor->focal;

    if (policy_.kind == OrientationPolicy::Kind::Default ||
        policy_.kind == OrientationPolicy::Kind::PerpendicularToBond) {
        // View direction perpendicular to the bond axis; inherit horizontal
        // direction from the current camera by projecting (oldFp - pos)
        // perpendicular to the bond axis. ViewUp = bond axis (so the bond
        // is vertical in screen space, perp to sight).
        double fpRaw[3];
        camera->GetFocalPoint(fpRaw);
        const model::Vec3 oldFp(fpRaw[0], fpRaw[1], fpRaw[2]);
        model::Vec3 sightCandidate = oldFp - pos;
        const double along = sightCandidate.dot(bondAxis);
        sightCandidate -= along * bondAxis;
        if (sightCandidate.norm() < 1e-9) {
            // Fallback: pick any vector perpendicular to bondAxis.
            model::Vec3 ref(1.0, 0.0, 0.0);
            if (std::abs(bondAxis.dot(ref)) > 0.99)
                ref = model::Vec3(0.0, 1.0, 0.0);
            sightCandidate = ref - ref.dot(bondAxis) * bondAxis;
        }
        sightCandidate.normalize();

        const model::Vec3 newPos = fp - sightCandidate * distance_;
        const model::Vec3 newUp  = bondAxis;
        writeCameraComposed(fp, newPos, newUp);
        return true;
    }

    // Free / DownAxis variants: inherit the sight direction.
    double fpRaw[3];
    camera->GetFocalPoint(fpRaw);
    const model::Vec3 oldFp(fpRaw[0], fpRaw[1], fpRaw[2]);
    model::Vec3 sightDir = oldFp - pos;
    if (sightDir.norm() < 1e-9)
        sightDir = model::Vec3(0.0, 0.0, -1.0);
    else
        sightDir.normalize();
    const model::Vec3 newPos = fp - sightDir * distance_;
    writeCameraComposed(fp, newPos, up);
    return true;
}

bool CameraComposer::writeDihedral(std::size_t t) {
    if (mode_.atoms.size() != 4) return false;
    const auto positions = readAtomPositions(t);
    if (positions.size() != 4) return false;

    const std::array<model::Vec3, 4> arr{{positions[0], positions[1],
                                           positions[2], positions[3]}};
    auto anchor = math::ComputeDihedralAnchor(arr);
    if (!anchor || !anchor->axis) return false;

    auto* camera = renderer_->GetActiveCamera();
    if (!camera) return false;
    double posRaw[3];
    double upRaw[3];
    camera->GetPosition(posRaw);
    camera->GetViewUp(upRaw);
    const model::Vec3 pos(posRaw[0], posRaw[1], posRaw[2]);
    const model::Vec3 oldUp(upRaw[0], upRaw[1], upRaw[2]);

    // Default policy for Dihedral: sight DOWN the (atoms[1], atoms[2])
    // axis (Newman projection). DownAxis policy picks a different pair.
    model::Vec3 axisVec;
    if (policy_.kind == OrientationPolicy::Kind::DownAxis) {
        const std::size_t a = policy_.axisAtoms[0];
        const std::size_t b = policy_.axisAtoms[1];
        if (a >= protein_->atomCount() || b >= protein_->atomCount())
            return false;
        const model::Vec3 p0 = conformation_->atomPosition(t, a);
        const model::Vec3 p1 = conformation_->atomPosition(t, b);
        axisVec = p1 - p0;
        if (axisVec.norm() < 1e-9) return false;
        axisVec.normalize();
    } else {
        axisVec = *anchor->axis;  // (c - b).normalized()
    }

    // Sign continuity: keep sight direction on the same side as last frame.
    model::Vec3 sightDir = axisVec;
    // Inherit sign from the current camera's view direction so the first
    // write doesn't flip — match dot-product against (oldFp - pos).
    double oldFpRaw[3];
    camera->GetFocalPoint(oldFpRaw);
    const model::Vec3 oldFp(oldFpRaw[0], oldFpRaw[1], oldFpRaw[2]);
    const model::Vec3 currentView = (oldFp - pos).norm() > 1e-9
        ? (oldFp - pos).normalized()
        : model::Vec3(0.0, 0.0, -1.0);
    if (currentView.dot(sightDir) < 0.0)
        sightDir = -sightDir;

    const model::Vec3 fp = anchor->focal;

    // ViewUp: prefer dihedral's natural up (from a-b leg projected perp
    // to axis); fall back to current up orthogonalised against sight.
    model::Vec3 up;
    if (anchor->viewUp) {
        up = *anchor->viewUp;
    } else {
        if (auto orthog = math::OrthogonalizeViewUp(sightDir, oldUp))
            up = *orthog;
        else
            up = model::Vec3(0.0, 1.0, 0.0);
    }
    if (auto orthog = math::OrthogonalizeViewUp(sightDir, up))
        up = *orthog;

    const model::Vec3 newPos = fp - sightDir * distance_;
    writeCameraComposed(fp, newPos, up);
    return true;
}

bool CameraComposer::writePlane(std::size_t t) {
    if (mode_.atoms.size() != 3) return false;
    const auto positions = readAtomPositions(t);
    if (positions.size() != 3) return false;

    const std::array<model::Vec3, 3> arr{{positions[0], positions[1], positions[2]}};
    auto anchor = math::ComputePlaneAnchor(arr);
    if (!anchor || !anchor->frame || !anchor->axis) return false;
    const auto& basis = *anchor->frame;

    // Per-frame normal-sign continuity guard (lifted verbatim from the
    // original applyCameraPlaneLock at MoleculeScene.cpp:382-394). The
    // basis is rebuilt each frame from (b-a)x(c-a); that direction can
    // flip sign across near-degenerate configurations (ring flip,
    // third atom crossing the line through the first two). Without the
    // guard the camera teleports to the other side of the plane.
    model::Vec3 viewDir = basis.z * (planeNormalSign_ < 0.0 ? -1.0 : 1.0);
    if (viewDir.norm() < 1e-9) return false;
    viewDir.normalize();
    if (planeLastDirection_
        && viewDir.dot(*planeLastDirection_) < 0.0) {
        planeNormalSign_ *= -1.0;
        viewDir = -viewDir;
    }
    planeLastDirection_ = viewDir;

    // Reconstruct ViewUp from plane-local coordinates captured at setMode.
    auto vectorToWorld = [&basis](const model::Vec3& local) {
        return basis.x * local.x() + basis.y * local.y() + basis.z * local.z();
    };
    model::Vec3 viewUp = vectorToWorld(planeLocalViewUp_);
    viewUp -= viewUp.dot(viewDir) * viewDir;
    if (viewUp.norm() < 1e-6) {
        viewUp = basis.y - basis.y.dot(viewDir) * viewDir;
    }
    if (viewUp.norm() < 1e-6) {
        viewUp = basis.x - basis.x.dot(viewDir) * viewDir;
    }
    if (viewUp.norm() < 1e-6) return false;
    viewUp.normalize();

    const model::Vec3 fp  = basis.origin;
    const model::Vec3 pos = fp - viewDir * distance_;
    writeCameraComposed(fp, pos, viewUp);
    return true;
}

bool CameraComposer::writeSubset(std::size_t t) {
    if (mode_.atoms.size() < 3 || subsetReference_.size() < 3) return false;

    std::vector<model::Vec3> current;
    current.reserve(mode_.atoms.size());
    for (std::size_t a : mode_.atoms) {
        if (a >= (protein_ ? protein_->atomCount() : 0)) return false;
        current.push_back(conformation_->atomPosition(t, a));
    }
    const std::size_t n = std::min(current.size(), subsetReference_.size());
    if (n < 3) return false;
    std::vector<model::Vec3> ref(subsetReference_.begin(),
                                  subsetReference_.begin() + n);
    current.resize(n);

    // Kabsch: R, T such that R * current[i] + T approximates ref[i].
    auto transform = math::ComputeSubsetTransform(current, ref);
    if (!transform) return false;

    // The camera should follow the stabilised local frame: the focal
    // moves to the centroid of the current subset (so the visible
    // protein stays centred); the view direction + viewUp are inherited
    // from the user's current camera state (Free orientation default),
    // rotated by the Kabsch transform's inverse so the camera "sees"
    // the protein from the same relative angle each frame.
    auto* camera = renderer_->GetActiveCamera();
    double posRaw[3];
    double upRaw[3];
    double fpRaw[3];
    camera->GetPosition(posRaw);
    camera->GetViewUp(upRaw);
    camera->GetFocalPoint(fpRaw);
    const model::Vec3 pos(posRaw[0], posRaw[1], posRaw[2]);
    const model::Vec3 up(upRaw[0], upRaw[1], upRaw[2]);
    const model::Vec3 fp(fpRaw[0], fpRaw[1], fpRaw[2]);

    // Subset centroid at current frame = mean(current). Move the focal
    // there, keep the same relative camera-to-focal vector rotated to
    // follow the molecule's orientation.
    model::Vec3 currentCentroid = model::Vec3::Zero();
    for (const auto& p : current) currentCentroid += p;
    currentCentroid /= static_cast<double>(current.size());

    // For Subset with no orientation override, simply translate the
    // focal to the current subset centroid and the position by the same
    // amount, preserving the gesture-driven sight direction. This is
    // the centroid-follow done absolutely, not as a delta from frame
    // to frame.
    const model::Vec3 delta = currentCentroid - fp;
    const model::Vec3 newFp  = currentCentroid;
    const model::Vec3 newPos = pos + delta;
    writeCameraComposed(newFp, newPos, up);
    return true;
}

void CameraComposer::writeCameraComposed(const model::Vec3& focal,
                                          const model::Vec3& position,
                                          const model::Vec3& viewUp) {
    model::Vec3 f = focal;
    model::Vec3 p = position;
    model::Vec3 u = viewUp;
    applyDeltaToState(f, p, u);

    auto* camera = renderer_->GetActiveCamera();
    if (!camera) return;
    camera->SetFocalPoint(f.x(), f.y(), f.z());
    camera->SetPosition(p.x(), p.y(), p.z());
    if (u.norm() > 1e-9) {
        u.normalize();
        camera->SetViewUp(u.x(), u.y(), u.z());
    }
    camera->OrthogonalizeViewUp();
}

void CameraComposer::applyDeltaToState(model::Vec3& focal,
                                        model::Vec3& position,
                                        model::Vec3& viewUp) const {
    // Pan: shift focal and position by accumulated world-space pan.
    focal    += accumPan_;
    position += accumPan_;

    // Build a local frame from the current state.
    model::Vec3 sight = focal - position;
    if (sight.norm() < 1e-9) return;
    sight.normalize();
    model::Vec3 right = sight.cross(viewUp);
    if (right.norm() < 1e-9) {
        // viewUp parallel to sight; choose an arbitrary perpendicular.
        model::Vec3 ref(0.0, 1.0, 0.0);
        if (std::abs(sight.dot(ref)) > 0.99) ref = model::Vec3(1.0, 0.0, 0.0);
        right = sight.cross(ref).normalized();
    } else {
        right.normalize();
    }
    model::Vec3 up = right.cross(sight).normalized();

    // Apply azimuth around current up, elevation around right, roll around sight.
    const model::Vec3 axisAzi   = up;
    const model::Vec3 axisElev  = right;
    const model::Vec3 axisRoll  = sight;

    auto rotateAround = [](const model::Vec3& focal,
                            const model::Vec3& point,
                            const model::Vec3& axis,
                            double angleRad) -> model::Vec3 {
        if (std::abs(angleRad) < 1e-12) return point;
        const model::Vec3 rel = point - focal;
        const model::Vec3 rotated = RodriguesRotate(rel, axis, angleRad);
        return focal + rotated;
    };
    position = rotateAround(focal, position, axisAzi,  accumAzimuthRad_);
    viewUp   = RodriguesRotate(up, axisAzi,  accumAzimuthRad_);
    position = rotateAround(focal, position, axisElev, accumElevationRad_);
    viewUp   = RodriguesRotate(viewUp, axisElev, accumElevationRad_);
    viewUp   = RodriguesRotate(viewUp, axisRoll, accumRollRad_);

    // Dolly: scale (position - focal) by 1/accumDolly_ (>1 zooms in).
    if (std::abs(accumDolly_) > 1e-9 && std::abs(accumDolly_ - 1.0) > 1e-12) {
        const model::Vec3 toCam = position - focal;
        position = focal + toCam / accumDolly_;
    }
}

void CameraComposer::applyGesture(const CameraGesture& g) {
    ASSERT_THREAD(this);
    switch (g.kind) {
        case CameraGesture::Kind::Azimuth:
            accumAzimuthRad_ += g.deltaRadians;
            break;
        case CameraGesture::Kind::Elevation:
            accumElevationRad_ += g.deltaRadians;
            break;
        case CameraGesture::Kind::Roll:
            accumRollRad_ += g.deltaRadians;
            break;
        case CameraGesture::Kind::Pan:
            // Pan is screen-space pixel delta; convert to world-space
            // using the camera's parallel scale (height of the viewport
            // in world units). The right vector and the up vector are
            // re-derived per call from the live camera, so accumulated
            // gestures always pan relative to the current view.
            if (auto* camera = renderer_ ? renderer_->GetActiveCamera() : nullptr) {
                double posRaw[3];
                double fpRaw[3];
                double upRaw[3];
                camera->GetPosition(posRaw);
                camera->GetFocalPoint(fpRaw);
                camera->GetViewUp(upRaw);
                const model::Vec3 pos(posRaw[0], posRaw[1], posRaw[2]);
                const model::Vec3 fp(fpRaw[0], fpRaw[1], fpRaw[2]);
                const model::Vec3 up(upRaw[0], upRaw[1], upRaw[2]);
                model::Vec3 sight = fp - pos;
                if (sight.norm() > 1e-9) sight.normalize();
                model::Vec3 right = sight.cross(up);
                if (right.norm() > 1e-9) right.normalize();
                // Convert px deltas to world units via the camera height.
                const double height = std::max(1.0, (pos - fp).norm()) * 0.001;
                accumPan_ += -right * (g.dxScreenPx * height) + up * (g.dyScreenPx * height);
            }
            break;
        case CameraGesture::Kind::Dolly:
            accumDolly_ *= g.dollyFactor;
            if (accumDolly_ < 1e-3) accumDolly_ = 1e-3;
            if (accumDolly_ > 1e3)  accumDolly_ = 1e3;
            break;
    }

    // For Free mode there is no per-frame fit to compose with; apply the
    // delta to the camera directly so the visible response is immediate
    // and not deferred to the next setFrame.
    if (mode_.kind == CameraMode::Kind::Free) {
        auto* camera = renderer_ ? renderer_->GetActiveCamera() : nullptr;
        if (!camera) return;
        double fpRaw[3];
        double posRaw[3];
        double upRaw[3];
        camera->GetFocalPoint(fpRaw);
        camera->GetPosition(posRaw);
        camera->GetViewUp(upRaw);
        model::Vec3 fp(fpRaw[0], fpRaw[1], fpRaw[2]);
        model::Vec3 pos(posRaw[0], posRaw[1], posRaw[2]);
        model::Vec3 up(upRaw[0], upRaw[1], upRaw[2]);
        // Reset deltas so we don't double-apply on the next setFrame.
        const double az  = accumAzimuthRad_;
        const double el  = accumElevationRad_;
        const double rl  = accumRollRad_;
        const model::Vec3 pan = accumPan_;
        const double dly = accumDolly_;
        accumAzimuthRad_   = 0.0;
        accumElevationRad_ = 0.0;
        accumRollRad_      = 0.0;
        accumPan_          = model::Vec3::Zero();
        accumDolly_        = 1.0;

        // Apply the just-extracted delta directly.
        fp  += pan;
        pos += pan;
        model::Vec3 sight = fp - pos;
        if (sight.norm() < 1e-9) return;
        sight.normalize();
        model::Vec3 right = sight.cross(up);
        if (right.norm() < 1e-9) {
            model::Vec3 ref(0.0, 1.0, 0.0);
            if (std::abs(sight.dot(ref)) > 0.99) ref = model::Vec3(1.0, 0.0, 0.0);
            right = sight.cross(ref).normalized();
        } else {
            right.normalize();
        }
        model::Vec3 upLocal = right.cross(sight).normalized();
        auto rotateAround = [](const model::Vec3& focal,
                                const model::Vec3& point,
                                const model::Vec3& axis,
                                double angleRad) -> model::Vec3 {
            if (std::abs(angleRad) < 1e-12) return point;
            const model::Vec3 rel = point - focal;
            return focal + RodriguesRotate(rel, axis, angleRad);
        };
        pos = rotateAround(fp, pos, upLocal, az);
        up  = RodriguesRotate(up, upLocal, az);
        pos = rotateAround(fp, pos, right,   el);
        up  = RodriguesRotate(up, right,   el);
        up  = RodriguesRotate(up, sight,   rl);
        if (std::abs(dly - 1.0) > 1e-12) {
            const model::Vec3 toCam = pos - fp;
            pos = fp + toCam / dly;
        }
        camera->SetFocalPoint(fp.x(), fp.y(), fp.z());
        camera->SetPosition(pos.x(), pos.y(), pos.z());
        if (up.norm() > 1e-9) {
            up.normalize();
            camera->SetViewUp(up.x(), up.y(), up.z());
        }
        camera->OrthogonalizeViewUp();
    }
    qCDebug(cComposer).noquote() << "gesture | kind=" << static_cast<int>(g.kind)
                                  << "| accumAz=" << accumAzimuthRad_
                                  << "| accumEl=" << accumElevationRad_;
}

}  // namespace h5reader::app
