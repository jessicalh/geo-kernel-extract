#include "SceneRevealOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/QtTrajectoryH5.h"
#include "../model/QtBondVectorBuffers.h"
#include "../model/TrajectoryConformation.h"
#include "TensorGlyphMath.h"

#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <cstdint>
#include <optional>
#include <utility>

namespace h5reader::app {

namespace {
constexpr std::size_t kMaxRevealAtoms = 96;
constexpr double kSphereRadiusA = 1.05;
constexpr double kSphereOpacity = 0.45;
constexpr double kLineWidth = 3.0;
constexpr double kLineOpacity = 0.92;
constexpr double kRevealRgb[3] = {0.0, 0.72, 0.78};
constexpr double kLineRgb[3] = {0.78, 1.0, 0.96};
// Distinct hue from the sphere/line highlights so the ellipsoid reads
// as "the orientation tensor" rather than another sphere. Soft amber.
constexpr double kTensorRgb[3] = {0.95, 0.74, 0.32};
constexpr double kTensorOpacity = 0.40;
// Scale factor on the eigendecomposition radii — orientation tensors
// have eigenvalues that sum to 1 (PSD trace-normalised), so sqrt
// values land near ~0.3-0.7. A 2.0× scale gives an ellipsoid roughly
// the size of a backbone bond — same visual budget as the highlight
// spheres.
constexpr double kTensorScale = 2.0;

// Resolve a (residue, kind) bond-vector anchor to its {tail, head}
// atom pair by walking the TR identity tables that own it. iRED is
// checked first (NH-only), then Reorient (NH + Cα-Hα + C=O). Returns
// nullopt when no TR knows the (residue, kind) tuple.
//
// Centralises the lookup so both `atomsForBinding` (atom-sphere
// highlights) and `reveal` (line endpoints) consume one source of
// truth — when a third BondVector-keyed TR lands, both call sites
// pick it up by extending this helper, not by editing both branches.
std::optional<std::pair<std::int32_t, std::int32_t>>
lookupBondVector(const h5reader::io::QtTrajectoryH5* h5,
                 std::size_t residue,
                 std::uint8_t kind) {
    if (!h5)
        return std::nullopt;
    auto fromTable = [&](const h5reader::model::QtBondVectorTable* table)
        -> std::optional<std::pair<std::int32_t, std::int32_t>> {
        if (!table || table->n_vectors == 0)
            return std::nullopt;
        const auto row = table->rowFor(residue, kind);
        if (!row)
            return std::nullopt;
        return std::make_pair(table->tail_atom[*row], table->head_atom[*row]);
    };
    const auto* ired = h5->iredOrderParameters();
    if (auto pair = fromTable(ired ? &ired->identity : nullptr))
        return pair;
    const auto* rd = h5->reorientationalDynamics();
    return fromTable(rd ? &rd->identity : nullptr);
}
}  // namespace

SceneRevealOverlay::SceneRevealOverlay(vtkSmartPointer<vtkRenderer> renderer,
                                       QObject* parent)
    : QObject(parent),
      renderer_(std::move(renderer))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("SceneRevealOverlay"));
}

SceneRevealOverlay::~SceneRevealOverlay()
{
    if (!renderer_)
        return;
    for (auto& actor : sphereActors_) {
        if (actor)
            renderer_->RemoveActor(actor);
    }
    if (lineActor_)
        renderer_->RemoveActor(lineActor_);
    if (tensorActor_)
        renderer_->RemoveActor(tensorActor_);
}

void SceneRevealOverlay::Build(const model::QtProtein& protein,
                               model::Conformation& conformation)
{
    ASSERT_THREAD(this);
    if (protein_ == &protein && conformation_ == &conformation && lineActor_)
        return;

    for (auto& actor : sphereActors_) {
        if (actor)
            renderer_->RemoveActor(actor);
    }
    spheres_.clear();
    sphereActors_.clear();

    if (lineActor_) {
        renderer_->RemoveActor(lineActor_);
        lineActor_ = nullptr;
    }

    protein_ = &protein;
    conformation_ = &conformation;
    activeAtoms_.clear();
    lineAtoms_.clear();
    active_ = false;
    lastFrame_ = 0;

    linePoints_ = vtkSmartPointer<vtkPoints>::New();
    lineCells_ = vtkSmartPointer<vtkCellArray>::New();
    lineData_ = vtkSmartPointer<vtkPolyData>::New();
    lineData_->SetPoints(linePoints_);
    lineData_->SetLines(lineCells_);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(lineData_);

    lineActor_ = vtkSmartPointer<vtkActor>::New();
    lineActor_->SetMapper(mapper);
    lineActor_->GetProperty()->SetColor(kLineRgb[0], kLineRgb[1], kLineRgb[2]);
    lineActor_->GetProperty()->SetLineWidth(kLineWidth);
    lineActor_->GetProperty()->SetOpacity(kLineOpacity);
    lineActor_->SetVisibility(0);
    lineActor_->PickableOff();
    renderer_->AddActor(lineActor_);

    // L-3a (2026-05-29): tensor glyph pipeline. vtkSphereSource → per-
    // frame vtkTransform (scale = eigendecomposition radii; rotation =
    // eigenvector frame) → vtkTransformPolyDataFilter → mapper →
    // actor. Mat3 update happens in applyTensorFrame; we set up the
    // pipeline once here so per-frame ticks are cheap.
    if (tensorActor_) {
        renderer_->RemoveActor(tensorActor_);
        tensorActor_ = nullptr;
    }
    tensorSphere_ = vtkSmartPointer<vtkSphereSource>::New();
    tensorSphere_->SetRadius(1.0);
    tensorSphere_->SetPhiResolution(24);
    tensorSphere_->SetThetaResolution(24);

    tensorTransform_ = vtkSmartPointer<vtkTransform>::New();
    tensorTransform_->Identity();

    tensorFilter_ = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    tensorFilter_->SetInputConnection(tensorSphere_->GetOutputPort());
    tensorFilter_->SetTransform(tensorTransform_);

    auto tensorMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    tensorMapper->SetInputConnection(tensorFilter_->GetOutputPort());

    tensorActor_ = vtkSmartPointer<vtkActor>::New();
    tensorActor_->SetMapper(tensorMapper);
    tensorActor_->GetProperty()->SetColor(kTensorRgb[0], kTensorRgb[1], kTensorRgb[2]);
    tensorActor_->GetProperty()->SetOpacity(kTensorOpacity);
    tensorActor_->SetVisibility(0);
    tensorActor_->PickableOff();
    renderer_->AddActor(tensorActor_);
    tensorActive_ = false;
}

void SceneRevealOverlay::ensureSphereCount(std::size_t count)
{
    ASSERT_THREAD(this);
    while (spheres_.size() < count) {
        auto sphere = vtkSmartPointer<vtkSphereSource>::New();
        sphere->SetRadius(kSphereRadiusA);
        sphere->SetPhiResolution(18);
        sphere->SetThetaResolution(18);

        auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputConnection(sphere->GetOutputPort());

        auto actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(kRevealRgb[0], kRevealRgb[1], kRevealRgb[2]);
        actor->GetProperty()->SetOpacity(kSphereOpacity);
        actor->SetVisibility(0);
        actor->PickableOff();
        renderer_->AddActor(actor);

        spheres_.push_back(sphere);
        sphereActors_.push_back(actor);
    }
}

std::vector<std::size_t> SceneRevealOverlay::atomsForBinding(const model::SignalBinding& binding) const
{
    std::vector<std::size_t> atoms;
    if (!protein_)
        return atoms;

    auto addAtom = [&](std::size_t atom) {
        if (atom >= protein_->atomCount())
            return;
        if (atoms.size() >= kMaxRevealAtoms)
            return;
        if (std::find(atoms.begin(), atoms.end(), atom) == atoms.end())
            atoms.push_back(atom);
    };
    auto addBond = [&](std::size_t bondIndex) {
        if (bondIndex >= protein_->bondCount())
            return;
        const auto& bond = protein_->bond(bondIndex);
        if (bond.atomIndexA >= 0)
            addAtom(static_cast<std::size_t>(bond.atomIndexA));
        if (bond.atomIndexB >= 0)
            addAtom(static_cast<std::size_t>(bond.atomIndexB));
    };
    auto addRing = [&](std::size_t ringIndex) {
        if (ringIndex >= protein_->ringCount())
            return;
        const auto& ring = protein_->ring(ringIndex);
        for (int32_t atom : ring.atomIndices) {
            if (atom >= 0)
                addAtom(static_cast<std::size_t>(atom));
        }
    };
    auto addTypedRing = [&](model::QtRingAxis axis, std::size_t ringIndex) {
        const auto absolute = protein_->topology().absoluteRingIndex(axis, ringIndex);
        if (absolute)
            addRing(*absolute);
    };
    auto addRingMembership = [&](std::size_t membershipIndex) {
        if (membershipIndex >= protein_->ringMembershipCount())
            return;
        const auto& membership = protein_->topology().ringMembershipAt(membershipIndex);
        if (membership.ringId >= 0)
            addRing(static_cast<std::size_t>(membership.ringId));
        if (membership.atomIndex >= 0)
            addAtom(static_cast<std::size_t>(membership.atomIndex));
    };

    const model::SignalAnchor& anchor = binding.anchor;
    if (const auto* atom = std::get_if<model::AtomAnchor>(&anchor)) {
        addAtom(atom->atom);
    } else if (const auto* residueAnchor = std::get_if<model::ResidueAnchor>(&anchor)) {
        if (residueAnchor->residue < protein_->residueCount()) {
            const auto& residue = protein_->residue(residueAnchor->residue);
            for (int32_t atom : residue.atomIndices) {
                if (atom >= 0)
                    addAtom(static_cast<std::size_t>(atom));
            }
        }
    } else if (const auto* tuple = std::get_if<model::AtomTupleAnchor>(&anchor)) {
        for (std::size_t atom : tuple->atoms)
            addAtom(atom);
    } else if (const auto* bond = std::get_if<model::BondAnchor>(&anchor)) {
        addBond(bond->bond);
    } else if (const auto* ring = std::get_if<model::RingAnchor>(&anchor)) {
        addRing(ring->ring);
    } else if (const auto* ring = std::get_if<model::AromaticRingAnchor>(&anchor)) {
        addTypedRing(model::QtRingAxis::AromaticRing, ring->ring);
    } else if (const auto* ring = std::get_if<model::SaturatedRingAnchor>(&anchor)) {
        addTypedRing(model::QtRingAxis::SaturatedRing, ring->ring);
    } else if (const auto* membership = std::get_if<model::RingMembershipAnchor>(&anchor)) {
        addRingMembership(membership->membership);
    } else if (const auto* vec = std::get_if<model::BondVectorAnchor>(&anchor)) {
        const auto* trajectory = conformation_ ? conformation_->asTrajectory() : nullptr;
        if (const auto pair = lookupBondVector(trajectory ? trajectory->h5() : nullptr,
                                               vec->residue, vec->kind)) {
            if (pair->first >= 0)
                addAtom(static_cast<std::size_t>(pair->first));
            if (pair->second >= 0)
                addAtom(static_cast<std::size_t>(pair->second));
        }
    }
    return atoms;
}

void SceneRevealOverlay::reveal(const model::SignalBinding& binding, int frame)
{
    ASSERT_THREAD(this);
    activeAtoms_ = atomsForBinding(binding);
    lineAtoms_.clear();
    auto setLineFromBond = [this](std::size_t bondIndex) {
        if (!protein_ || bondIndex >= protein_->bondCount())
            return;
        const auto& bond = protein_->bond(bondIndex);
        if (bond.atomIndexA >= 0)
            lineAtoms_.push_back(static_cast<std::size_t>(bond.atomIndexA));
        if (bond.atomIndexB >= 0)
            lineAtoms_.push_back(static_cast<std::size_t>(bond.atomIndexB));
    };
    auto setLineFromRing = [this](std::size_t ringIndex) {
        if (!protein_ || ringIndex >= protein_->ringCount())
            return;
        const auto& ring = protein_->ring(ringIndex);
        for (int32_t atom : ring.atomIndices) {
            if (atom >= 0)
                lineAtoms_.push_back(static_cast<std::size_t>(atom));
        }
        if (lineAtoms_.size() > 2)
            lineAtoms_.push_back(lineAtoms_.front());
    };
    const model::SignalAnchor& anchor = binding.anchor;
    if (const auto* tuple = std::get_if<model::AtomTupleAnchor>(&anchor)) {
        for (std::size_t atom : tuple->atoms) {
            if (protein_ && atom < protein_->atomCount())
                lineAtoms_.push_back(atom);
        }
    } else if (const auto* bond = std::get_if<model::BondAnchor>(&anchor)) {
        setLineFromBond(bond->bond);
    } else if (const auto* ring = std::get_if<model::RingAnchor>(&anchor)) {
        setLineFromRing(ring->ring);
    } else if (const auto* ring = std::get_if<model::AromaticRingAnchor>(&anchor)) {
        const auto absolute = protein_ ? protein_->topology().absoluteRingIndex(model::QtRingAxis::AromaticRing, ring->ring)
                                       : std::nullopt;
        if (absolute)
            setLineFromRing(*absolute);
    } else if (const auto* ring = std::get_if<model::SaturatedRingAnchor>(&anchor)) {
        const auto absolute = protein_ ? protein_->topology().absoluteRingIndex(model::QtRingAxis::SaturatedRing, ring->ring)
                                       : std::nullopt;
        if (absolute)
            setLineFromRing(*absolute);
    } else if (const auto* membership = std::get_if<model::RingMembershipAnchor>(&anchor)) {
        if (protein_ && membership->membership < protein_->ringMembershipCount()) {
            const auto& row = protein_->topology().ringMembershipAt(membership->membership);
            if (row.ringId >= 0)
                setLineFromRing(static_cast<std::size_t>(row.ringId));
        }
    } else if (const auto* vec = std::get_if<model::BondVectorAnchor>(&anchor)) {
        // Reuse the same helper that atomsForBinding uses, so the
        // sphere highlights and the line endpoints always agree on
        // which row the (residue, kind) anchor resolves to.
        const auto* trajectory = conformation_ ? conformation_->asTrajectory() : nullptr;
        if (const auto pair = lookupBondVector(trajectory ? trajectory->h5() : nullptr,
                                               vec->residue, vec->kind)) {
            if (pair->first >= 0)
                lineAtoms_.push_back(static_cast<std::size_t>(pair->first));
            if (pair->second >= 0)
                lineAtoms_.push_back(static_cast<std::size_t>(pair->second));
        }
    }
    if (activeAtoms_.empty()) {
        clear();
        return;
    }
    active_ = true;
    ensureSphereCount(activeAtoms_.size());
    applyFrame(frame);
}

void SceneRevealOverlay::clear()
{
    ASSERT_THREAD(this);
    active_ = false;
    activeAtoms_.clear();
    lineAtoms_.clear();
    for (auto& actor : sphereActors_) {
        if (actor)
            actor->SetVisibility(0);
    }
    if (lineActor_)
        lineActor_->SetVisibility(0);
}

void SceneRevealOverlay::setFrame(int t)
{
    ASSERT_THREAD(this);
    lastFrame_ = t;
    applyFrame(t);
    applyTensorFrame(t);
}

void SceneRevealOverlay::revealTensor(std::size_t tailAtom,
                                      std::size_t headAtom,
                                      const std::array<double, 9>& tensor,
                                      int frame)
{
    ASSERT_THREAD(this);
    if (!protein_ || tailAtom >= protein_->atomCount() || headAtom >= protein_->atomCount()) {
        clearTensor();
        return;
    }
    tensorTail_ = tailAtom;
    tensorHead_ = headAtom;
    tensorData_ = tensor;
    tensorActive_ = true;
    applyTensorFrame(frame);
}

void SceneRevealOverlay::clearTensor()
{
    ASSERT_THREAD(this);
    tensorActive_ = false;
    if (tensorActor_)
        tensorActor_->SetVisibility(0);
}

void SceneRevealOverlay::applyTensorFrame(int t)
{
    if (!tensorActor_ || !tensorTransform_)
        return;
    if (!tensorActive_ || !protein_ || !conformation_) {
        tensorActor_->SetVisibility(0);
        return;
    }
    if (t < 0 || static_cast<std::size_t>(t) >= conformation_->frameCount()) {
        tensorActor_->SetVisibility(0);
        return;
    }
    // Body→lab transform per frame. The bond_orientation_tensor is
    // accumulated in the body frame (producer-side: trajectory rigidly
    // aligned to frame 0). The eigenvectors decomposeSymmetric3x3
    // returns are body-frame; writing them into the VTK world
    // transform unchanged would produce a glyph whose orientation is
    // wrong once the molecule rotates relative to frame 0 (Codex
    // NOW-4, 2026-05-29). composeEllipsoidTransform takes the current
    // lab-frame bond direction and rotates the body-frame
    // eigenvectors to align the primary axis with the current bond,
    // approximating the per-frame body→lab rotation without needing
    // the full Kabsch matrix.
    const model::Vec3 tail =
        conformation_->atomPosition(static_cast<std::size_t>(t), tensorTail_);
    const model::Vec3 head =
        conformation_->atomPosition(static_cast<std::size_t>(t), tensorHead_);
    const model::Vec3 mid = (tail + head) * 0.5;
    const model::Vec3 bondVec = head - tail;
    if (bondVec.norm() < 1e-9) {
        tensorActor_->SetVisibility(0);
        return;
    }
    const model::Vec3 bondDir = bondVec.normalized();

    const auto eig = math::decomposeSymmetric3x3(tensorData_);
    const auto M = math::composeEllipsoidTransform(eig, bondDir, mid, kTensorScale);

    vtkNew<vtkMatrix4x4> vtkM;
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
            vtkM->SetElement(r, c, M[r * 4 + c]);
    tensorTransform_->SetMatrix(vtkM);
    tensorTransform_->Modified();
    tensorActor_->SetVisibility(1);
}

void SceneRevealOverlay::applyFrame(int t)
{
    if (!protein_ || !conformation_ || !lineActor_)
        return;
    if (t < 0 || static_cast<std::size_t>(t) >= conformation_->frameCount())
        return;

    const std::size_t visibleCount = active_ ? activeAtoms_.size() : 0;
    ensureSphereCount(visibleCount);
    for (std::size_t i = 0; i < sphereActors_.size(); ++i) {
        const bool show = i < visibleCount;
        if (show) {
            const model::Vec3 p =
                conformation_->atomPosition(static_cast<std::size_t>(t), activeAtoms_[i]);
            spheres_[i]->SetCenter(p.x(), p.y(), p.z());
        }
        sphereActors_[i]->SetVisibility(show ? 1 : 0);
    }

    bool showLine = active_ && lineAtoms_.size() >= 2;
    if (showLine) {
        linePoints_->Reset();
        for (std::size_t atom : lineAtoms_) {
            if (atom >= protein_->atomCount()) {
                showLine = false;
                break;
            }
            const model::Vec3 p =
                conformation_->atomPosition(static_cast<std::size_t>(t), atom);
            linePoints_->InsertNextPoint(p.x(), p.y(), p.z());
        }
        if (showLine) {
            std::vector<vtkIdType> ids(lineAtoms_.size());
            for (std::size_t i = 0; i < ids.size(); ++i)
                ids[i] = static_cast<vtkIdType>(i);
            lineCells_->Reset();
            lineCells_->InsertNextCell(static_cast<vtkIdType>(ids.size()), ids.data());
            linePoints_->Modified();
            lineCells_->Modified();
            lineData_->Modified();
        }
    }
    lineActor_->SetVisibility(showLine ? 1 : 0);
}

}  // namespace h5reader::app
