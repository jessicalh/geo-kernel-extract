#include "QtSelectionOverlay.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include "../model/QtFrame.h"

#include <QLoggingCategory>

#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

namespace h5reader::app {

namespace {
Q_LOGGING_CATEGORY(cSel, "h5reader.overlay.selection")

constexpr double kSphereRadiusA = 1.0;    // same as library viewer
constexpr double kOpacity       = 0.30;
}

QtSelectionOverlay::QtSelectionOverlay(
    vtkSmartPointer<vtkRenderer>                  renderer,
    vtkSmartPointer<vtkGenericOpenGLRenderWindow> renderWindow,
    QObject* parent)
    : QObject(parent),
      renderer_(std::move(renderer)),
      renderWindow_(std::move(renderWindow))
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("QtSelectionOverlay"));
}

QtSelectionOverlay::~QtSelectionOverlay() {
    if (actor_) renderer_->RemoveActor(actor_);
}

void QtSelectionOverlay::Build(const model::QtProtein&      protein,
                                const model::QtConformation& conformation) {
    ASSERT_THREAD(this);
    if (protein_ == &protein && conformation_ == &conformation && actor_)
        return;

    if (actor_) renderer_->RemoveActor(actor_);
    actor_ = nullptr;

    protein_      = &protein;
    conformation_ = &conformation;
    hasSelection_ = false;

    sphere_ = vtkSmartPointer<vtkSphereSource>::New();
    sphere_->SetRadius(kSphereRadiusA);
    sphere_->SetPhiResolution(16);
    sphere_->SetThetaResolution(16);

    auto mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(sphere_->GetOutputPort());

    actor_ = vtkSmartPointer<vtkActor>::New();
    actor_->SetMapper(mapper);
    actor_->GetProperty()->SetColor(1.0, 1.0, 0.0);
    actor_->GetProperty()->SetOpacity(kOpacity);
    actor_->SetVisibility(0);
    renderer_->AddActor(actor_);
}

void QtSelectionOverlay::setPickedAtom(std::size_t atomIdx) {
    ASSERT_THREAD(this);
    if (!protein_ || atomIdx >= protein_->atomCount()) return;
    pickedAtom_   = atomIdx;
    hasSelection_ = true;
    // Position comes from the playback controller via setFrame; here we
    // just mark the selection as present. Caller is expected to call
    // setFrame(currentFrame) immediately after for the sphere to appear.
    qCInfo(cSel).noquote() << "picked atom" << atomIdx;
}

void QtSelectionOverlay::clearSelection() {
    ASSERT_THREAD(this);
    hasSelection_ = false;
    if (actor_) actor_->SetVisibility(0);
}

void QtSelectionOverlay::setFrame(int t) {
    ASSERT_THREAD(this);
    if (!hasSelection_ || !protein_ || !conformation_) {
        if (actor_) actor_->SetVisibility(0);
        return;
    }
    applyCurrentPosition(t);
}

void QtSelectionOverlay::applyCurrentPosition(int t) {
    if (t < 0 || static_cast<size_t>(t) >= conformation_->frameCount()) return;
    if (pickedAtom_ >= protein_->atomCount()) {
        actor_->SetVisibility(0);
        return;
    }
    const auto& frame = conformation_->frame(static_cast<size_t>(t));
    const model::Vec3 p = frame.position(pickedAtom_);
    sphere_->SetCenter(p.x(), p.y(), p.z());
    actor_->SetVisibility(1);
}

}  // namespace h5reader::app
