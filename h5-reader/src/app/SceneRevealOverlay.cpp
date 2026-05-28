#include "SceneRevealOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

#include <algorithm>
#include <cstdint>
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

    switch (binding.anchorKind) {
        case model::SignalAnchorKind::Atom:
            if (binding.atom)
                addAtom(*binding.atom);
            break;
        case model::SignalAnchorKind::Residue:
            if (binding.residue && *binding.residue < protein_->residueCount()) {
                const auto& residue = protein_->residue(*binding.residue);
                for (int32_t atom : residue.atomIndices) {
                    if (atom >= 0)
                        addAtom(static_cast<std::size_t>(atom));
                }
            }
            break;
        case model::SignalAnchorKind::AtomTuple:
            for (std::size_t atom : binding.atomTuple)
                addAtom(atom);
            break;
        case model::SignalAnchorKind::None:
            break;
    }
    return atoms;
}

void SceneRevealOverlay::reveal(const model::SignalBinding& binding, int frame)
{
    ASSERT_THREAD(this);
    activeAtoms_ = atomsForBinding(binding);
    lineAtoms_.clear();
    if (binding.anchorKind == model::SignalAnchorKind::AtomTuple) {
        for (std::size_t atom : binding.atomTuple) {
            if (protein_ && atom < protein_->atomCount())
                lineAtoms_.push_back(atom);
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
