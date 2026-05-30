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
    // Dihedral sign-continuity reset (Codex finding #1) — same lifecycle
    // as planeLastDirection_; first write after setMode picks the natural
    // axis direction, subsequent writes flip the sign only if the axis
    // crosses through perpendicular to the stored reference.
    dihedralLastDirection_.reset();
    // Atom/Bond reference captures (Codex finding #2) — zero out so
    // captureInitialState's per-mode arm sets the right values; without
    // this reset, switching Atom -> Bond would inherit the Atom mode's
    // captured sight as the bond's initial reference.
    atomReferenceSight_    = model::Vec3::Zero();
    atomReferenceUp_       = model::Vec3::Zero();
    atomReferenceCamRel_   = model::Vec3::Zero();
    bondReferenceSight_    = model::Vec3::Zero();
    bondReferenceUp_       = model::Vec3::Zero();
    bondReferenceCamRel_   = model::Vec3::Zero();
    bondReferenceMidpoint_ = model::Vec3::Zero();
    subsetReference_.clear();
    subsetReferenceCamRel_   = model::Vec3::Zero();
    subsetReferenceUp_       = model::Vec3::Zero();
    subsetReferenceCentroid_ = model::Vec3::Zero();

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

    if (mode_.kind == CameraMode::Kind::Atom && mode_.atoms.size() == 1) {
        // Atom mode reference capture (Codex finding #2). The prior
        // implementation derived each frame's sight from the live camera
        // — which already contained accumulated user gestures — so the
        // gesture re-applied on top of itself frame after frame
        // (visible drift even with no further input). Now we capture
        // sight/up/cam-relative once here and writeAtom composes
        // gestures on top of the captured pose exactly once per frame.
        const std::size_t a = mode_.atoms[0];
        if (a < protein_->atomCount()) {
            const model::Vec3 atomPos =
                conformation_->atomPosition(referenceFrame, a);
            atomReferenceCamRel_ = pos - atomPos;
            model::Vec3 sight = fp - pos;
            if (sight.norm() > 1e-9) {
                sight.normalize();
                atomReferenceSight_ = sight;
                // Use safeViewUp so the captured up is guaranteed
                // perpendicular even when the live camera presents a
                // degenerate up (Codex finding #3).
                atomReferenceUp_ = math::safeViewUp(sight, up);
            } else {
                atomReferenceSight_ = model::Vec3(0.0, 0.0, -1.0);
                atomReferenceUp_    = math::safeViewUp(atomReferenceSight_, up);
            }
        }
    } else if (mode_.kind == CameraMode::Kind::Bond && mode_.atoms.size() == 2) {
        // Bond mode reference capture (Codex finding #2). Mirror of the
        // Atom arm above; the bond's anchor is the midpoint, so we
        // cache that for the per-frame fallback alongside the
        // sight/up/cam-relative triple.
        const std::size_t a = mode_.atoms[0];
        const std::size_t b = mode_.atoms[1];
        if (a < protein_->atomCount() && b < protein_->atomCount()) {
            const model::Vec3 pa = conformation_->atomPosition(referenceFrame, a);
            const model::Vec3 pb = conformation_->atomPosition(referenceFrame, b);
            bondReferenceMidpoint_ = 0.5 * (pa + pb);
            bondReferenceCamRel_   = pos - bondReferenceMidpoint_;
            model::Vec3 sight = fp - pos;
            if (sight.norm() > 1e-9) {
                sight.normalize();
                bondReferenceSight_ = sight;
                bondReferenceUp_    = math::safeViewUp(sight, up);
            } else {
                bondReferenceSight_ = model::Vec3(0.0, 0.0, -1.0);
                bondReferenceUp_    = math::safeViewUp(bondReferenceSight_, up);
            }
        }
    } else if (mode_.kind == CameraMode::Kind::Plane && mode_.atoms.size() == 3) {
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
        // Capture the camera's pose relative to the subset's centroid at
        // lock acquisition. Each frame the Kabsch fit gives R such that
        // R * current[i] + T ≈ reference[i]; rotating the captured
        // camera-relative vector by R^T moves the camera into the
        // current frame's body axes. This is the rotation half of
        // writeSubset (subsetReferenceCamRel_ + subsetReferenceUp_ +
        // subsetReferenceCentroid_ are the persistent reference state).
        if (subsetReference_.size() >= 3) {
            subsetReferenceCentroid_ = model::Vec3::Zero();
            for (const auto& p : subsetReference_) subsetReferenceCentroid_ += p;
            subsetReferenceCentroid_ /= static_cast<double>(subsetReference_.size());
            subsetReferenceCamRel_ = pos - subsetReferenceCentroid_;
            // Orthogonalise the captured up against the captured sight
            // (camera-to-focal direction). Falls back to canonical up if
            // the captured camera has no usable sight direction.
            model::Vec3 sight = fp - pos;
            if (sight.norm() > 1e-9) {
                sight.normalize();
                if (auto orthog = math::OrthogonalizeViewUp(sight, up))
                    subsetReferenceUp_ = *orthog;
                else
                    subsetReferenceUp_ = model::Vec3(0.0, 1.0, 0.0);
            } else {
                subsetReferenceUp_ = up.norm() > 1e-9 ? up.normalized()
                                                       : model::Vec3(0.0, 1.0, 0.0);
            }
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

    // Codex finding #2: derive the per-frame natural pose from the
    // captured reference at lock acquisition, NOT from the live camera.
    // The live camera already contains any composed user gestures from
    // the previous frame; deriving sight from it would re-apply those
    // gestures on top of themselves each tick, producing drift even
    // with no further input. The captured reference gives the same
    // natural pose every frame and writeCameraComposed applies the
    // accumulated gesture delta on top exactly once.
    const model::Vec3 newFp  = anchor->focal;
    const model::Vec3 newPos = newFp + atomReferenceCamRel_;
    writeCameraComposed(newFp, newPos, atomReferenceUp_);
    return true;
}

bool CameraComposer::writeBond(std::size_t t) {
    if (mode_.atoms.size() != 2) return false;
    const auto positions = readAtomPositions(t);
    if (positions.size() != 2) return false;

    const std::array<model::Vec3, 2> arr{{positions[0], positions[1]}};
    auto anchor = math::ComputeBondAnchor(arr);
    if (!anchor || !anchor->axis) return false;

    const model::Vec3 bondAxis = *anchor->axis;
    const model::Vec3 fp       = anchor->focal;

    // Codex finding #2: derive each frame's natural camera pose from
    // the captured reference at lock acquisition, NOT from the live
    // camera. The captured cam-relative offset rides on top of the
    // current midpoint; writeCameraComposed applies the accumulated
    // gesture delta on top exactly once. This kills the gesture
    // double-apply drift the prior live-camera-inherits implementation
    // had.
    //
    // Policy still drives the natural view-up:
    //   Default / PerpendicularToBond — ViewUp = bondAxis (the bond is
    //     vertical in screen space; sight is perpendicular to it). The
    //     captured cam-relative offset already encodes a perpendicular
    //     sight from the reference frame and we re-anchor at the
    //     current midpoint each frame.
    //   Free / DownAxis — use the captured reference up; sight is
    //     inherited from the captured cam-relative offset (which the
    //     user oriented when they acquired the lock).
    if (policy_.kind == OrientationPolicy::Kind::Default ||
        policy_.kind == OrientationPolicy::Kind::PerpendicularToBond) {
        // Build a sight that's perpendicular to the bond axis. Project the
        // captured cam-relative direction perpendicular to the current
        // frame's bond axis (the bond axis rotates frame-to-frame with
        // the molecule; the captured sight may not be exactly perp to
        // the new axis). If the projection collapses (captured cam-rel
        // parallel to the new bond axis), pick any perpendicular via
        // safeViewUp.
        model::Vec3 sightCandidate = bondReferenceCamRel_
            - bondReferenceCamRel_.dot(bondAxis) * bondAxis;
        if (sightCandidate.norm() < 1e-9) {
            sightCandidate = math::safeViewUp(bondAxis, bondReferenceCamRel_);
        } else {
            sightCandidate.normalize();
        }

        // Distance: use the captured cam-relative magnitude so the bond
        // lock holds the captured zoom even after the bond axis rotates.
        const double dist = std::max(1e-6, bondReferenceCamRel_.norm());
        const model::Vec3 newPos = fp + sightCandidate * dist;
        // ViewUp = bondAxis per the perp-to-bond convention. The bond
        // is vertical in screen space; safeViewUp guards against the
        // degenerate (axis parallel to sight) case that
        // OrthogonalizeViewUp would otherwise have to handle.
        const model::Vec3 newUp = math::safeViewUp(
            /*sight=*/(fp - newPos).normalized(), bondAxis);
        writeCameraComposed(fp, newPos, newUp);
        return true;
    }

    // Free / DownAxis: use the captured reference cam-relative offset
    // and up. Camera follows the bond's midpoint at the captured offset
    // (preserving distance + orientation captured at lock acquisition).
    const model::Vec3 newPos = fp + bondReferenceCamRel_;
    writeCameraComposed(fp, newPos, bondReferenceUp_);
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
    double upRaw[3];
    camera->GetViewUp(upRaw);
    const model::Vec3 oldUp(upRaw[0], upRaw[1], upRaw[2]);

    // Default policy for Dihedral: sight DOWN the (atoms[1], atoms[2])
    // axis (Newman projection). DownAxis policy picks a different pair.
    // The override changes the axis, which means the natural viewUp
    // computed against the default axis no longer applies — viewUp must
    // be recomputed perpendicular to the override axis.
    model::Vec3 axisVec;
    const bool axisOverridden =
        policy_.kind == OrientationPolicy::Kind::DownAxis;
    if (axisOverridden) {
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

    // Codex finding #1: sign continuity for the sight axis is an
    // EXPLICIT state guard (dihedralLastDirection_), not an implicit
    // feedback loop through the live camera. The old code dotted axisVec
    // against the camera's current view direction — which already
    // contained both the prior frame's output AND any accumulated user
    // gesture. That's a feedback loop: a user rotation that pushed the
    // view across the dot-product sign boundary would flip the next
    // frame's axis, producing a 180° camera teleport.
    //
    // Now the guard uses the POST-flip axis we wrote last frame
    // (mirroring writePlane / planeLastDirection_). On the first write
    // after setMode (dihedralLastDirection_ == nullopt), we use the
    // natural axis direction and STORE it. Subsequent frames flip only
    // if axisVec dots negative against the stored reference, then we
    // update the reference with the new post-flip axis.
    model::Vec3 sightDir = axisVec;
    if (dihedralLastDirection_
        && sightDir.dot(*dihedralLastDirection_) < 0.0) {
        sightDir = -sightDir;
    }
    dihedralLastDirection_ = sightDir;

    const model::Vec3 fp = anchor->focal;

    // ViewUp: when the axis is overridden, the dihedral's natural viewUp
    // (computed against (c-b)) doesn't apply. Recompute from the flanking
    // dihedral atom (a) projected perpendicular to the override axis;
    // fall back via safeViewUp using the live camera's current up as
    // the preferred direction. Without the override, prefer the
    // anchor's natural up.
    //
    // Codex finding #3: every fallback now routes through safeViewUp,
    // which is guaranteed perpendicular to sight even when the
    // preferred candidate is parallel to sight (the old (0,1,0)
    // fallback failed silently when sight aligned with world Y; VTK's
    // OrthogonalizeViewUp can't recover from a sight-parallel up).
    model::Vec3 up;
    if (axisOverridden) {
        // Project (a - midpoint) perpendicular to the override axis to
        // get a stable up that follows the molecule's geometry rather
        // than the gesture's last orientation. If (a - midpoint) is
        // parallel to the override axis, fall back via safeViewUp.
        const model::Vec3& aPos = positions[0];
        model::Vec3 candidate = aPos - fp;
        candidate -= candidate.dot(sightDir) * sightDir;
        if (candidate.norm() > 1e-9) {
            up = candidate.normalized();
        } else {
            up = math::safeViewUp(sightDir, oldUp);
        }
    } else if (anchor->viewUp) {
        up = *anchor->viewUp;
    } else {
        up = math::safeViewUp(sightDir, oldUp);
    }
    // Final guard: orthogonalise up against sight. safeViewUp guarantees
    // a non-parallel result; this is belt-and-suspenders against any
    // numerical drift between the projection above and the final write.
    up = math::safeViewUp(sightDir, up);

    const model::Vec3 newPos = fp - sightDir * distance_;
    writeCameraComposed(fp, newPos, up);
    return true;
}

bool CameraComposer::writePlane(std::size_t t) {
    if (mode_.atoms.size() != 3) return false;
    // Plane mode's only meaningful orientation policies are Default and
    // PerpendicularToPlane — both reduce to "sight along the plane
    // normal". Other policies don't apply (Free / PerpToBond / DownAxis
    // would override the normal-sign continuity that's the whole point
    // of writePlane). Reject loud rather than silently honour them.
    if (policy_.kind != OrientationPolicy::Kind::Default
        && policy_.kind != OrientationPolicy::Kind::PerpendicularToPlane) {
        qCWarning(cComposer).noquote()
            << "writePlane: policy" << NameFor(policy_.kind)
            << "not applicable to plane mode; rejecting frame";
        return false;
    }
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
    model::Vec3 viewUpCandidate = vectorToWorld(planeLocalViewUp_);
    // Codex finding #3: route through safeViewUp so the fallback chain
    // is deterministic and guaranteed-non-degenerate. The prior ad-hoc
    // sequence (try plane Y, then plane X, then bail) failed silently
    // when the plane axes themselves were parallel to sight (degenerate
    // captured plane); safeViewUp's world-axis fallback chain always
    // finds a perpendicular.
    const model::Vec3 viewUp = math::safeViewUp(viewDir, viewUpCandidate);

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
    // The body-to-current rotation is R^T (so that R^T * v_in_reference
    // = v_in_current_frame). The camera was captured in the reference
    // frame; rotating its relative pose by R^T moves it into the
    // current frame so the molecule appears stationary while the camera
    // follows its rotation.
    //
    // Codex finding #4 handles rank-deficient inputs inside
    // ComputeSubsetTransform — nullopt return means the fit can't be
    // trusted; freezing the frame (return false) keeps the camera at
    // its last-good state.
    auto transform = math::ComputeSubsetTransform(current, ref);
    if (!transform) {
        qCWarning(cComposer).noquote()
            << "writeSubset | frame=" << static_cast<qlonglong>(t)
            << "| Kabsch returned nullopt (rank-deficient or det invalid); freezing";
        return false;
    }

    // Codex finding #5: belt-and-suspenders validation before applying
    // R^T. ComputeSubsetTransform's own guards should already guarantee
    // R^T * R ~ I and det(R) = +1, but a defensive check at the use
    // site catches propagation bugs (e.g. if a future refactor adds a
    // post-processing step). If either guard fires we freeze the frame
    // and log; with the upstream guards in place these should never
    // fire in normal use.
    constexpr double kOrthoTol = 1e-6;
    constexpr double kDetTol   = 1e-6;
    const model::Mat3 RtR_minus_I = transform->R.transpose() * transform->R
                                     - model::Mat3::Identity();
    if (RtR_minus_I.norm() > kOrthoTol) {
        qCWarning(cComposer).noquote()
            << "writeSubset | frame=" << static_cast<qlonglong>(t)
            << "| R^T*R - I Frobenius norm=" << RtR_minus_I.norm()
            << "exceeds" << kOrthoTol << "; freezing";
        return false;
    }
    if (std::abs(transform->R.determinant() - 1.0) > kDetTol) {
        qCWarning(cComposer).noquote()
            << "writeSubset | frame=" << static_cast<qlonglong>(t)
            << "| det(R)=" << transform->R.determinant()
            << "not ~+1; freezing";
        return false;
    }
    const model::Mat3 Rinv = transform->R.transpose();

    // Subset centroid at current frame = mean(current). Focal lands
    // there; position and view-up are the captured reference state
    // rotated by R^T into the molecule's current orientation.
    model::Vec3 currentCentroid = model::Vec3::Zero();
    for (const auto& p : current) currentCentroid += p;
    currentCentroid /= static_cast<double>(current.size());

    const model::Vec3 newCamRel = Rinv * subsetReferenceCamRel_;
    const model::Vec3 newUp     = Rinv * subsetReferenceUp_;
    const model::Vec3 newFp     = currentCentroid;
    const model::Vec3 newPos    = currentCentroid + newCamRel;

    // Orthogonalise the rotated up against the rotated sight; both came
    // from a single rigid rotation of orthogonal reference vectors, so
    // this is a guard against accumulated floating-point drift across
    // many frames. safeViewUp (Codex finding #3) guarantees a non-
    // degenerate result even if newUp drifts parallel to sight.
    model::Vec3 sight = newFp - newPos;
    if (sight.norm() < 1e-9) return false;
    sight.normalize();
    const model::Vec3 orthoUp = math::safeViewUp(sight, newUp);

    writeCameraComposed(newFp, newPos, orthoUp);
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
            // Pan is a screen-space pixel delta; convert to world-space
            // by dividing the world-units-per-pixel implied by the
            // current camera + viewport. For a perspective camera, the
            // world height seen at the focal distance is
            //     2 * distance * tan(view_angle/2)
            // so one pixel = that height / viewport_height. For a
            // parallel camera, ParallelScale already IS half the world
            // height of the viewport in world units. Both cases derive
            // from live camera + renderer state so pans scale correctly
            // on proteins of any size.
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

                const int* sizePx = renderer_->GetSize();
                const double viewportHeightPx =
                    sizePx && sizePx[1] > 0 ? static_cast<double>(sizePx[1]) : 600.0;
                double worldPerPx = 1.0;
                if (camera->GetParallelProjection()) {
                    // ParallelScale = half the world-height of the viewport.
                    worldPerPx = (2.0 * camera->GetParallelScale()) / viewportHeightPx;
                } else {
                    const double distance = std::max(1e-6, (pos - fp).norm());
                    const double halfAngleRad =
                        camera->GetViewAngle() * 0.5 * M_PI / 180.0;
                    const double worldHeight = 2.0 * distance * std::tan(halfAngleRad);
                    worldPerPx = worldHeight / viewportHeightPx;
                }
                accumPan_ += -right * (g.dxScreenPx * worldPerPx)
                              + up    * (g.dyScreenPx * worldPerPx);
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
