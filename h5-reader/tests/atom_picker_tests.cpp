#include <QtTest>
#include <QSurfaceFormat>
#include <QVTKOpenGLNativeWidget.h>

#include "model/QtAtom.h"
#include "model/QtAtomNames.h"
#include "model/QtBond.h"
#include "model/QtResidue.h"
#include "model/QtResidueNames.h"
#include "model/QtRing.h"
#include "model/QtRingMembership.h"
#include "model/QtTopology.h"

#define private public
#include "model/QtProtein.h"
#undef private

#include "app/MeasurementOverlay.h"
#include "app/MoleculeScene.h"
#include "app/QtAtomPicker.h"
#include "model/AtomSelection.h"
#include "model/QtConformationSnapshot.h"
#include "model/SingleConformation.h"

#include <vtkActorCollection.h>
#include <vtkCamera.h>
#include <vtkMapper.h>

using namespace h5reader;

namespace {

std::unique_ptr<model::QtProtein> twoAtomProtein() {
    auto protein = std::make_unique<model::QtProtein>();
    protein->atoms_.resize(2);
    protein->atomNames_.resize(2);
    protein->topology_ = std::make_unique<model::QtTopology>(
        2, std::vector<model::QtBond>{},
        std::vector<std::unique_ptr<model::QtRing>>{},
        std::vector<model::QtRingMembership>{}, 0, 0);
    return protein;
}

QPoint widgetPosition(vtkRenderer* renderer, const QVTKOpenGLNativeWidget& widget,
                      const model::Vec3& position) {
    renderer->SetWorldPoint(position.x(), position.y(), position.z(), 1.0);
    renderer->WorldToDisplay();
    const double* display = renderer->GetDisplayPoint();
    const double scale = widget.devicePixelRatioF();
    return {qRound(display[0] / scale), qRound(widget.height() - display[1] / scale)};
}

}  // namespace

class AtomPickerTests final : public QObject {
    Q_OBJECT

private slots:
    void filteredPickingAndRotatedMarker() {
        auto protein = twoAtomProtein();
        auto snapshot = std::make_shared<model::QtConformationSnapshot>(protein.get(), 0, 0.0);
        auto& positions = snapshot->mutableColumn(io::FieldKind::Pos);
        positions.present = true;
        positions.rows = 2;
        positions.cols = 3;
        positions.data = {0.0, 0.0, 0.0, 0.0, 0.0, -4.0};
        model::SingleConformation conformation(protein.get(), snapshot);
        QVTKOpenGLNativeWidget widget;
        widget.resize(640, 480);
        auto window = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
        widget.setRenderWindow(window);
        app::MoleculeScene scene(&widget, window);
        scene.Build(*protein, conformation);
        model::AtomSelection selection(protein.get());
        auto* marker = scene.measurementOverlay();
        marker->setSelection(&selection);
        QObject::connect(&selection, &model::AtomSelection::changed,
                         marker, &app::MeasurementOverlay::onSelectionChanged);
        app::QtAtomPicker picker(&widget, &scene, &conformation, nullptr);
        QObject::connect(&picker, &app::QtAtomPicker::atomPicked,
                         &selection, &model::AtomSelection::applyPick);
        QSignalSpy picks(&picker, &app::QtAtomPicker::atomPicked);
        widget.show();
        QVERIFY(QTest::qWaitForWindowExposed(&widget));

        auto* camera = scene.Renderer()->GetActiveCamera();
        camera->SetPosition(0.0, 0.0, 20.0);
        camera->SetFocalPoint(0.0, 0.0, -4.0);
        camera->SetViewUp(0.0, 1.0, 0.0);
        scene.setAtomFilter({1});
        window->Render();
        const auto atomPosition = conformation.atomPosition(0, 1);
        auto point = widgetPosition(scene.Renderer(), widget, atomPosition);

        // The hidden atom is first in storage and lies on the same click ray.
        QCOMPARE(picker.atomAt(point.x(), point.y()), std::optional<std::size_t>{1});
        QTest::mouseDClick(&widget, Qt::LeftButton, Qt::NoModifier, point);
        QCOMPARE(picks.count(), 1);
        QCOMPARE(picks.last().at(0).value<std::size_t>(), std::size_t{1});

        auto* actors = scene.OverlayRenderer()->GetActors();
        actors->InitTraversal();
        auto* markerActor = actors->GetNextActor();
        QVERIFY(markerActor != nullptr);
        auto* sphere = vtkSphereSource::SafeDownCast(markerActor->GetMapper()->GetInputAlgorithm());
        QVERIFY(sphere != nullptr);
        QCOMPARE(scene.OverlayRenderer()->GetActiveCamera(), camera);
        for (double angle : {37.0, 53.0, 90.0}) {
            camera->Azimuth(angle);
            window->Render();
            point = widgetPosition(scene.Renderer(), widget, atomPosition);
            QCOMPARE(picker.atomAt(point.x(), point.y()), std::optional<std::size_t>{1});
            QTest::mouseDClick(&widget, Qt::LeftButton, Qt::NoModifier, point);
            QCOMPARE(picks.last().at(0).value<std::size_t>(), std::size_t{1});
            QVERIFY(markerActor->GetVisibility());
            const model::Vec3 markerPosition(sphere->GetCenter());
            QVERIFY((markerPosition - atomPosition).norm() < 1e-12);
            QCOMPARE(widgetPosition(scene.OverlayRenderer(), widget, markerPosition), point);
        }

        scene.clearAtomFilter();
        camera->SetPosition(0.0, 0.0, 20.0);
        camera->SetViewUp(0.0, 1.0, 0.0);
        window->Render();
        point = widgetPosition(scene.Renderer(), widget, conformation.atomPosition(0, 0));
        QCOMPARE(picker.atomAt(point.x(), point.y()), std::optional<std::size_t>{0});
    }
};

int main(int argc, char** argv) {
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication application(argc, argv);
    AtomPickerTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "atom_picker_tests.moc"
