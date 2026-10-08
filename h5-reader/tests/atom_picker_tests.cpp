#include <QtTest>
#include <QSurfaceFormat>
#include <QTreeWidget>
#include <QScrollBar>
#include <QSplitter>
#include <QLabel>
#include <QToolButton>
#include <QDir>
#include <QContextMenuEvent>
#include <QHeaderView>
#include <QStyleOptionViewItem>
#include <QPointer>
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
#include "app/CameraComposer.h"
#include "app/CameraInputFilter.h"
#include "app/MoleculeScene.h"
#include "app/QtAtomPicker.h"
#include "app/QtAtomInspectorDock.h"
#include "app/ReaderMainWindow.h"
#include "app/AtomInspectorGlossary.h"
#include "app/MetricGlossaryPopup.h"
#include "app/QtPlaybackController.h"
#include "app/SelectionContextWidget.h"
#include "app/TensorGlyphActor.h"
#include "app/CsaTensorOverlay.h"
#include "model/AtomSelection.h"
#include "model/ConformationGeometry.h"
#include "model/QtConformationSnapshot.h"
#include "model/SingleConformation.h"
#include "model/TrajectorySignalCatalog.h"
#include "model/TrajectoryConformation.h"
#include "model/TransformedConformation.h"

#include <vtkActorCollection.h>
#include <vtkCamera.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkMapper.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>

using namespace h5reader;

namespace {

std::unique_ptr<model::QtProtein> smallProtein(
    std::size_t atomCount = 2, std::vector<model::QtBond> bonds = {}) {
    auto protein = std::make_unique<model::QtProtein>();
    protein->atoms_.resize(atomCount);
    protein->atomNames_.resize(atomCount);
    for (auto& atom : protein->atoms_) atom.element = model::Element::C;
    protein->topology_ = std::make_unique<model::QtTopology>(
        atomCount, std::move(bonds),
        std::vector<std::unique_ptr<model::QtRing>>{},
        std::vector<model::QtRingMembership>{}, 0, 0);
    return protein;
}

QPoint widgetPosition(vtkRenderer* renderer, const QVTKOpenGLNativeWidget& widget,
                      const model::Vec3& position) {
    renderer->SetWorldPoint(position.x(), position.y(), position.z(), 1.0);
    renderer->WorldToDisplay();
    const double* display = renderer->GetDisplayPoint();
    const double scale = widget.effectiveDevicePixelRatio();
    return {qRound(display[0] / scale),
            qRound((widget.renderWindow()->GetSize()[1] - 1 - display[1]) / scale)};
}

}  // namespace

class AtomPickerTests final : public QObject {
    Q_OBJECT

private slots:
    void tensorPalettesMatchInspector() {
        auto protein = smallProtein();
        auto snapshot = std::make_shared<model::QtConformationSnapshot>(protein.get(), 0, 0.0);
        model::SingleConformation conformation(protein.get(), snapshot);
        app::QtAtomInspectorDock inspector;
        inspector.setContext(protein.get(), &conformation);
        inspector.setPickedAtom(1);
        app::OrientationTensorInfo orientation;
        inspector.setOrientationTensor(1, orientation);
        auto* tree = inspector.findChild<QTreeWidget*>(QStringLiteral("tensorFields"));
        QVERIFY(tree);
        for (const auto& source : {QStringLiteral("ORCA DFT"), QStringLiteral("Predicted")}) {
            app::CsaTensorInfo shielding;
            shielding.sourceLabel = source;
            inspector.setCsaTensor(1, shielding);
            auto* principalValues = tree->topLevelItem(0)->child(2);
            QCOMPARE(principalValues->text(0), QStringLiteral("Principal values"));
            for (int i = 0; i < 3; ++i) {
                const auto& warm = app::kShieldingTensorColours[i];
                const auto& cool = app::kOrientationTensorColours[i];
                const QColor shieldingColour = QColor::fromRgbF(warm[0], warm[1], warm[2]);
                const QColor orientationColour = QColor::fromRgbF(cool[0], cool[1], cool[2]);
                QCOMPARE(principalValues->child(i)->icon(0).pixmap(14, 14).toImage().pixelColor(6, 6).rgb(),
                         shieldingColour.rgb());
                auto* axis = tree->topLevelItem(1)->child(3 + i);
                QCOMPARE(axis->text(0), QStringLiteral("lambda_%1").arg(i + 1));
                QCOMPARE(axis->icon(0).pixmap(14, 14).toImage().pixelColor(6, 6).rgb(),
                         orientationColour.rgb());
                QVERIFY(shieldingColour != orientationColour);
            }
        }
        const auto checkActors = [](vtkRenderer* renderer, const app::TensorAxisColours& colours) {
            auto* actors = renderer->GetActors();
            actors->InitTraversal();
            int arrow = 0;
            while (auto* actor = actors->GetNextActor()) {
                if (!actor->GetVisibility()) continue;
                QVERIFY(arrow < 6);
                const auto& rgb = colours[arrow / 2];
                for (int channel = 0; channel < 3; ++channel)
                    QCOMPARE(actor->GetProperty()->GetColor()[channel], rgb[channel]);
                ++arrow;
            }
            QCOMPARE(arrow, 6);
        };
        auto renderer = vtkSmartPointer<vtkRenderer>::New();
        {
            app::CsaTensorOverlay shielding(renderer);
            model::Mat3 matrix = model::Mat3::Zero();
            matrix.diagonal() << 1.0, 2.0, 4.0;
            shielding.show(model::Vec3::Zero(), model::ComputeCsaShape(matrix));
            checkActors(renderer, app::kShieldingTensorColours);
        }
        {
            app::TensorGlyphActor orientationGlyph(renderer);
            app::TensorGlyphActor::Style style;
            style.axisColours = app::kOrientationTensorColours;
            orientationGlyph.show(model::Vec3::Zero(), {0.8, 0.15, 0.05}, model::Mat3::Identity(), 1.0 / 3, 1.0, style);
            checkActors(renderer, app::kOrientationTensorColours);
            orientationGlyph.show(model::Vec3::Zero(), {1.0, 2.0, 4.0}, model::Mat3::Identity(), 7.0 / 3);
            checkActors(renderer, app::kDefaultTensorColours);
        }

        app::CsaTensorInfo pending;
        pending.sourceLabel = QStringLiteral("Predicted");
        pending.status = QStringLiteral("Calculating...");
        inspector.setCsaTensor(1, pending);
        QCOMPARE(tree->topLevelItem(0)->text(1), pending.status);
        QCOMPARE(tree->topLevelItem(0)->childCount(), 2); // Atom and source, no placeholder numerical values.
        QVERIFY(tree->topLevelItem(0)->flags().testFlag(Qt::ItemIsUserCheckable));
    }

    void inspectorGlossaryOnRightClick() {
        auto protein = smallProtein();
        auto snapshot = std::make_shared<model::QtConformationSnapshot>(protein.get(), 0, 0.0);
        model::SingleConformation conformation(protein.get(), snapshot);
        app::QtAtomInspectorDock inspector;
        inspector.resize(520, 650);
        inspector.setContext(protein.get(), &conformation);
        inspector.setPickedAtom(1);
        app::CsaTensorInfo shielding;
        shielding.sourceLabel = QStringLiteral("ORCA DFT");
        inspector.setCsaTensor(1, shielding);
        inspector.show();
        QVERIFY(QTest::qWaitForWindowExposed(&inspector));
        auto* fields = inspector.findChild<QTreeWidget*>(QStringLiteral("atomFields"));
        auto* tensors = inspector.findChild<QTreeWidget*>(QStringLiteral("tensorFields"));
        QVERIFY(fields && tensors);
        const auto selectionBefore = inspector.dumpTree();
        const auto rightClick = [](QTreeWidget* tree, QTreeWidgetItem* item, int column) {
            tree->scrollToItem(item);
            const QRect row = tree->visualItemRect(item);
            const QPoint point(tree->header()->sectionViewportPosition(column) + 12, row.center().y());
            QContextMenuEvent event(QContextMenuEvent::Mouse, point, tree->viewport()->mapToGlobal(point));
            QApplication::sendEvent(tree->viewport(), &event);
        };

        auto* element = fields->topLevelItem(0)->child(0)->child(0);
        QCOMPARE(element->text(0), QStringLiteral("Element"));
        for (int column : {0, 1}) {
            rightClick(fields, element, column);
            QPointer<QWidget> popup = QApplication::activePopupWidget();
            QVERIFY(popup);
            QCOMPARE(popup->objectName(), QStringLiteral("MetricGlossaryPopup"));
            QCOMPARE(popup->findChild<QLabel*>(QStringLiteral("glossaryTitle"))->text(), QStringLiteral("Element"));
            QCOMPARE(popup->findChild<QLabel*>(QStringLiteral("glossaryMeaning"))->text(), QStringLiteral("Chemical element of the atom."));
            QCOMPARE(inspector.dumpTree(), selectionBefore);
            QTest::keyClick(popup, Qt::Key_Escape);
            QTRY_VERIFY(popup.isNull());
        }

        rightClick(tensors, tensors->topLevelItem(0)->child(1), 0);
        QPointer<QWidget> popup = QApplication::activePopupWidget();
        QVERIFY(popup);
        QCOMPARE(popup->findChild<QLabel*>(QStringLiteral("glossaryTitle"))->text(), QStringLiteral("sigma_iso"));
        QVERIFY(popup->findChild<QLabel*>(QStringLiteral("glossaryOrigin"))->text().contains(QStringLiteral("ORCA DFT")));
        QCOMPARE(inspector.dumpTree(), selectionBefore);
        const QString calculation = popup->findChild<QLabel*>(QStringLiteral("glossaryCalculation"))->text();
        inspector.setCsaTensor(1, shielding);
        inspector.clearSelection();
        QVERIFY(popup);
        QCOMPARE(popup->findChild<QLabel*>(QStringLiteral("glossaryCalculation"))->text(), calculation);
        const QString shots = qEnvironmentVariable("H5READER_UI_TEST_SHOTS");
        if (!shots.isEmpty())
            QVERIFY(popup->grab().save(QDir(shots).filePath(QStringLiteral("right-click-glossary.png"))));
        QTest::keyClick(popup, Qt::Key_Escape);
        QTRY_VERIFY(popup.isNull());

        QContextMenuEvent blank(QContextMenuEvent::Mouse, QPoint(20, 150), tensors->viewport()->mapToGlobal(QPoint(20, 150)));
        QApplication::sendEvent(tensors->viewport(), &blank);
        QVERIFY(!QApplication::activePopupWidget());

        model::SignalDescriptor descriptor;
        descriptor.family = QStringLiteral("sasa");
        descriptor.conceptKey = QStringLiteral("atom_sasa");
        descriptor.label = QStringLiteral("SASA");
        app::ShowMetricGlossaryPopup(descriptor, inspector.mapToGlobal(QPoint(30, 30)), &inspector);
        popup = QApplication::activePopupWidget();
        QVERIFY(popup);
        QCOMPARE(popup->findChild<QLabel*>(QStringLiteral("glossaryTitle"))->text(), descriptor.label);
        QVERIFY(popup->findChild<QLabel*>(QStringLiteral("glossaryMeaning"))->text().contains(QStringLiteral("solvent")));
        QTest::keyClick(popup, Qt::Key_Escape);
        QTRY_VERIFY(popup.isNull());
    }

    void inspectorGlossaryKeepsTensorSourcesDistinct() {
        using app::AtomInspectorGlossary;
        const QStringList tensorFields = {
            QStringLiteral("T0 signed iso"), QStringLiteral("|T2| anisotropy"),
            QStringLiteral("|T1| antisymmetric"), QStringLiteral("span"),
            QStringLiteral("eta"), QStringLiteral("skew"), QStringLiteral("sigma_11"),
            QStringLiteral("T1 antisym vector"), QStringLiteral("m=-2")};
        const QStringList sources = {
            QStringLiteral("bs_shielding"), QStringLiteral("hm_shielding"),
            QStringLiteral("coulomb_shielding"), QStringLiteral("mopac_coulomb_shielding"),
            QStringLiteral("mopac_mc_shielding"), QStringLiteral("\u03c3 total"),
            QStringLiteral("\u03c3 diamagnetic"), QStringLiteral("\u03c3 paramagnetic"),
            QStringLiteral("\u0394\u03c3 total")};
        for (const auto& source : sources) {
            QVERIFY2(AtomInspectorGlossary({source}, {}).has_value(), qPrintable(source));
            for (const auto& field : tensorFields) {
                const auto help = AtomInspectorGlossary({source, field}, {});
                QVERIFY2(help.has_value(), qPrintable(source + QLatin1Char('/') + field));
                QVERIFY(!help->meaning.isEmpty() && !help->calculation.isEmpty() && !help->origin.isEmpty());
            }
        }
        const auto dia = AtomInspectorGlossary({QStringLiteral("\u03c3 diamagnetic"), QStringLiteral("T0")}, {});
        const auto para = AtomInspectorGlossary({QStringLiteral("\u03c3 paramagnetic"), QStringLiteral("T0")}, {});
        QVERIFY(dia->meaning.contains(QStringLiteral("diamagnetic")));
        QVERIFY(para->meaning.contains(QStringLiteral("paramagnetic")));
        const auto efg = AtomInspectorGlossary({QStringLiteral("aimnet2_efg"), QStringLiteral("PAS / EFG convention"), QStringLiteral("eta")}, {});
        QVERIFY(efg);
        QVERIFY(efg->calculation.contains(QStringLiteral("(Vyy - Vxx) / Vzz")));
        const auto coulomb = AtomInspectorGlossary({QStringLiteral("coulomb_shielding")}, {});
        QVERIFY(coulomb->calculation.contains(QStringLiteral("second spatial derivatives")));
        const auto bondScale = AtomInspectorGlossary({QStringLiteral("Bond orientation tensor"), QStringLiteral("Glyph size")}, {});
        QVERIFY(bondScale);
        QVERIFY(bondScale->calculation.contains(QStringLiteral("absolute deviations")));
        QVERIFY(bondScale->calculation.contains(QStringLiteral("from their mean")));
        QVERIFY(!AtomInspectorGlossary({QStringLiteral("No tensor for this atom")}, {}));
    }

    void inspectorConnectsTensorNumbersToGlyphs() {
        auto protein = smallProtein();
        auto snapshot = std::make_shared<model::QtConformationSnapshot>(protein.get(), 0, 0.0);
        model::SingleConformation conformation(protein.get(), snapshot);
        app::QtAtomInspectorDock inspector;
        inspector.resize(520, 450);
        inspector.setContext(protein.get(), &conformation);
        inspector.setPickedAtom(0);
        app::CsaTensorInfo shielding;
        shielding.sourceLabel = QStringLiteral("ORCA DFT");
        shielding.sigma11 = -10.0;
        shielding.sigma22 = 20.0;
        shielding.sigma33 = 80.0;
        shielding.sigmaIso = 30.0;
        shielding.span = 90.0;
        inspector.setCsaTensor(0, shielding);
        app::OrientationTensorInfo orientation;
        orientation.bond = QStringLiteral("A:ALA7 N-H");
        orientation.lambda1 = 0.8;
        orientation.lambda2 = 0.15;
        orientation.lambda3 = 0.05;
        orientation.s2 = 0.4975;
        inspector.setOrientationTensor(0, orientation);
        inspector.setTensorVisibility(true, true);
        inspector.show();
        QVERIFY(QTest::qWaitForWindowExposed(&inspector));

        auto* tree = inspector.findChild<QTreeWidget*>(QStringLiteral("atomFields"));
        auto* tensors = inspector.findChild<QTreeWidget*>(QStringLiteral("tensorFields"));
        auto* splitter = inspector.findChild<QSplitter*>(QStringLiteral("atomInfoSplitter"));
        QVERIFY(tree);
        QVERIFY(tensors && splitter);
        QCOMPARE(splitter->count(), 2);
        QCOMPARE(splitter->orientation(), Qt::Vertical);
        QVERIFY(!splitter->childrenCollapsible());
        QVERIFY(!tensors->isWindow());
        splitter->setSizes({180, 240});
        const auto sizes = splitter->sizes();
        const auto child = [](QTreeWidgetItem* parent, const QString& text) -> QTreeWidgetItem* {
            for (int i = 0; i < parent->childCount(); ++i)
                if (parent->child(i)->text(0) == text) return parent->child(i);
            return nullptr;
        };
        auto* root = tree->topLevelItem(0);
        QCOMPARE(root->child(0)->text(0), QStringLiteral("Identity"));
        QVERIFY(!child(root, QStringLiteral("Shielding tensor (ORCA DFT)")));
        auto* csa = tensors->topLevelItem(0);
        auto* bond = tensors->topLevelItem(1);
        QCOMPARE(csa->text(0), QStringLiteral("Shielding tensor (ORCA DFT)"));
        QCOMPARE(bond->text(0), QStringLiteral("Bond orientation tensor"));
        QCOMPARE(csa->text(1), QStringLiteral("Shown"));
        QVERIFY(csa->font(0).bold());
        auto* axes = child(csa, QStringLiteral("Principal values"));
        QVERIFY(axes && axes->isExpanded());
        QCOMPARE(axes->text(1), QStringLiteral("ppm (from mean)"));
        QCOMPARE(axes->child(0)->text(1), QStringLiteral("-10 (-40)"));
        QCOMPARE(axes->child(1)->text(1), QStringLiteral("20 (-10)"));
        QCOMPARE(axes->child(2)->text(1), QStringLiteral("80 (+50)"));
        QVERIFY(!axes->child(0)->icon(0).isNull());
        auto* details = child(csa, QStringLiteral("Details"));
        QVERIFY(details && !details->isExpanded());
        QCOMPARE(child(bond, QStringLiteral("Average over"))->text(1), QStringLiteral("Trajectory"));
        QCOMPARE(child(bond, QStringLiteral("Bond"))->text(1), orientation.bond);
        QCOMPARE(child(bond, QStringLiteral("Main axis follows"))->text(1), QStringLiteral("Current bond"));

        tree->verticalScrollBar()->setValue(5);
        tensors->verticalScrollBar()->setValue(2);
        const int tensorScroll = tensors->verticalScrollBar()->value();
        QVERIFY(tensorScroll > 0);
        const int scroll = tree->verticalScrollBar()->value();
        QVERIFY(scroll > 0);
        inspector.setTensorVisibility(false, true);
        QVERIFY(csa->text(1).isEmpty());
        QVERIFY(!csa->font(0).bold());
        QCOMPARE(bond->text(1), QStringLiteral("Shown"));
        QCOMPARE(tree->verticalScrollBar()->value(), scroll);
        QCOMPARE(tensors->verticalScrollBar()->value(), tensorScroll);
        inspector.setCsaTensor(0, shielding);
        QCOMPARE(tree->verticalScrollBar()->value(), scroll);
        QCOMPARE(tensors->verticalScrollBar()->value(), tensorScroll);
        QCOMPARE(splitter->sizes(), sizes);
        inspector.clearCsaTensor();
        inspector.setTensorVisibility(false, false);
        inspector.clearSelection();
        inspector.setTensorVisibility(true, true);
        QCOMPARE(tree->topLevelItemCount(), 1);
        QCOMPARE(tree->topLevelItem(0)->childCount(), 0);
        QCOMPARE(tensors->topLevelItemCount(), 0);
    }

    void tensorCheckboxesKeepTheirChoicesAcrossRebuilds() {
        auto protein = smallProtein();
        auto snapshot = std::make_shared<model::QtConformationSnapshot>(protein.get(), 0, 0.0);
        model::SingleConformation conformation(protein.get(), snapshot);
        app::QtAtomInspectorDock inspector;
        inspector.resize(520, 650);
        inspector.setContext(protein.get(), &conformation);
        inspector.setPickedAtom(0);
        QSignalSpy changes(&inspector, &app::QtAtomInspectorDock::tensorDisplayChanged);
        QSignalSpy requests(&inspector, &app::QtAtomInspectorDock::shieldingTensorRequested);
        inspector.setCsaTensor(0, {});
        inspector.setOrientationTensor(0, {});
        inspector.setTensorVisibility(true, true);
        inspector.show();
        QVERIFY(QTest::qWaitForWindowExposed(&inspector));
        auto* tree = inspector.findChild<QTreeWidget*>(QStringLiteral("tensorFields"));
        QVERIFY(tree);
        auto* shielding = tree->topLevelItem(0);
        auto* orientation = tree->topLevelItem(1);
        QCOMPARE(shielding->checkState(0), Qt::Checked);
        QCOMPARE(orientation->checkState(0), Qt::Checked);
        QCOMPARE(changes.count(), 0);

        tree->scrollToItem(shielding);
        QStyleOptionViewItem option;
        option.initFrom(tree);
        option.rect = tree->visualItemRect(shielding);
        option.features = QStyleOptionViewItem::HasCheckIndicator;
        option.checkState = Qt::Checked;
        const QRect checkbox = tree->style()->subElementRect(
            QStyle::SE_ItemViewItemCheckIndicator, &option, tree);
        QTest::mouseClick(tree->viewport(), Qt::LeftButton, Qt::NoModifier, checkbox.center());
        QVERIFY(!inspector.shieldingTensorEnabled());
        QVERIFY(inspector.orientationTensorEnabled());
        QCOMPARE(changes.count(), 1);
        inspector.setTensorVisibility(false, true);
        QCOMPARE(changes.count(), 1);

        tree->setCurrentItem(orientation);
        QTest::keyClick(tree, Qt::Key_Space);
        QVERIFY(!inspector.orientationTensorEnabled());
        QCOMPARE(changes.count(), 2);
        inspector.setTensorVisibility(false, false);

        inspector.setFrame(1);
        inspector.setCsaTensor(0, {});
        inspector.setOrientationTensor(0, {});
        QCOMPARE(tree->topLevelItem(0)->checkState(0), Qt::Unchecked);
        QCOMPARE(tree->topLevelItem(1)->checkState(0), Qt::Unchecked);
        QCOMPARE(changes.count(), 2);
        inspector.clearSelection();
        inspector.setPickedAtom(1);
        inspector.setOrientationTensor(1, {});
        QCOMPARE(tree->topLevelItemCount(), 1);
        QCOMPARE(tree->topLevelItem(0)->checkState(0), Qt::Unchecked);
        QCOMPARE(changes.count(), 2);

        tree->setCurrentItem(tree->topLevelItem(0));
        QTest::keyClick(tree, Qt::Key_Space);
        QVERIFY(inspector.orientationTensorEnabled());
        QVERIFY(!inspector.shieldingTensorEnabled());
        QCOMPARE(changes.count(), 3);
        inspector.setCsaTensor(1, {});
        QCOMPARE(tree->topLevelItem(0)->checkState(0), Qt::Unchecked);
        QCOMPARE(tree->topLevelItem(1)->checkState(0), Qt::Checked);
        QCOMPARE(changes.count(), 3);
        QCOMPARE(requests.count(), 0);
        inspector.setTensorDisplayEnabled(true, true);
        QCOMPARE(tree->topLevelItem(0)->checkState(0), Qt::Checked);
        QCOMPARE(changes.count(), 4);
        QCOMPARE(requests.count(), 1);
    }

    void tensorCheckboxesInWindowPreserveTreeItems() {
        const QString fixture = qEnvironmentVariable("H5READER_REST_FIXTURE");
        if (fixture.isEmpty())
            QSKIP("Set H5READER_REST_FIXTURE to a trajectory with ORCA and bond tensors.");
        app::ReaderMainWindow window;
        QVERIFY2(window.loadRunPath(fixture, false), qPrintable(window.lastLoadError()));
        window.show();
        QVERIFY(QTest::qWaitForWindowExposed(&window));

        const auto* trajectory = window.transformedConformation()->asTrajectory();
        QVERIFY(trajectory);
        const auto* dynamics = trajectory->h5()->reorientationalDynamics();
        QVERIFY(dynamics && dynamics->identity.n_vectors > 0);
        auto* selection = window.findChild<model::AtomSelection*>();
        auto* inspector = window.findChild<app::QtAtomInspectorDock*>();
        auto* scene = window.findChild<app::MoleculeScene*>();
        QVERIFY(selection && inspector && scene);
        selection->applyPick(std::size_t(dynamics->identity.tail_atom[0]), Qt::NoModifier);
        QTRY_VERIFY_WITH_TIMEOUT(scene->csaOverlay()->isActive(), 30000);
        auto* tree = inspector->findChild<QTreeWidget*>(QStringLiteral("tensorFields"));
        QVERIFY(tree);
        QCOMPARE(tree->topLevelItemCount(), 2);
        QSignalSpy resets(tree->model(), &QAbstractItemModel::modelReset);

        for (int row = 0; row < tree->topLevelItemCount(); ++row) {
            auto* group = tree->topLevelItem(row);
            auto* details = group->child(group->childCount() - 1);
            QCOMPARE(details->text(0), QStringLiteral("Details"));
            details->setExpanded(true);
            tree->setCurrentItem(group);
            const QPersistentModelIndex heading(tree->currentIndex());
            for (const auto state : {Qt::Unchecked, Qt::Checked}) {
                QTest::keyClick(tree, Qt::Key_Space);
                bool inputHandled = false;
                QMetaObject::invokeMethod(&window, [&] { inputHandled = true; },
                                          Qt::QueuedConnection);
                QTRY_VERIFY(inputHandled);
                QCOMPARE(resets.count(), 0);
                QVERIFY(heading.isValid());
                QCOMPARE(tree->currentIndex(), QModelIndex(heading));
                QCOMPARE(group->checkState(0), state);
                QVERIFY(details->isExpanded());
            }
        }
        window.shutdown();
    }

    void selectionContextKeepsHighlightAndSelectionSeparate() {
        auto protein = smallProtein(4);
        protein->atomNames_[0].iupac = QStringLiteral("CA");
        protein->atomNames_[1].iupac = QStringLiteral("CB");
        protein->atomNames_[2].iupac = QStringLiteral("CG");
        protein->atomNames_[3].iupac = QStringLiteral("CD1");
        auto snapshot = std::make_shared<model::QtConformationSnapshot>(protein.get(), 0, 0.0);
        auto& positions = snapshot->mutableColumn(io::FieldKind::Pos);
        positions.present = true;
        positions.rows = 4;
        positions.cols = 3;
        positions.data = {0.0, 0.0, 0.0, 0.0, 0.0, -4.0,
                          1.0, 0.0, -4.0, 1.0, 1.0, -4.0};
        model::SingleConformation raw(protein.get(), snapshot);
        model::TransformedConformation conformation(&raw);
        QVTKOpenGLNativeWidget widget;
        widget.resize(640, 480);
        auto window = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
        widget.setRenderWindow(window);
        app::MoleculeScene scene(&widget, window);
        scene.Build(*protein, conformation);
        model::AtomSelection selection(protein.get());
        app::QtPlaybackController playback(1);
        model::TrajectorySignalCatalog catalog;
        app::QtAtomInspectorDock inspector;
        inspector.setContext(protein.get(), &conformation);
        app::SelectionContextWidget context(*protein, selection, conformation,
                                          scene, playback, catalog, &inspector);
        inspector.setSelectionContextWidget(&context);
        inspector.resize(360, 560);
        inspector.show();
        widget.show();
        QVERIFY(QTest::qWaitForWindowExposed(&widget));
        QVERIFY(QTest::qWaitForWindowExposed(&inspector));
        QVERIFY(context.isVisible());
        QVERIFY(!context.isWindow());

        selection.bulkSet({0, 1});
        QVERIFY(context.isVisible());
        QCOMPARE(context.stateJson()["meaning"].toString(),
                 QStringLiteral("Atoms are not bonded."));
        model::SignalBinding binding;
        binding.anchor = model::AtomAnchor{1};
        scene.revealBinding(binding);
        QVERIFY(context.stateJson()["highlight"].toString().contains(QStringLiteral("CB")));
        const auto mode = scene.cameraComposer()->mode();
        const auto atoms = selection.atoms();

        QTest::keyClick(&context, Qt::Key_Escape);
        QVERIFY(context.isVisible());
        QVERIFY(scene.activeRevealBinding());
        QVERIFY(scene.cameraComposer()->mode() == mode);
        QCOMPARE(selection.atoms(), atoms);
        playback.setFrame(0);
        QVERIFY(context.isVisible());
        QVERIFY(!context.triggerAction(QStringLiteral("dismiss")));
        context.setInteractionEnabled(false);
        QVERIFY(!context.triggerAction(QStringLiteral("clear_highlight")));
        QVERIFY(scene.activeRevealBinding());
        context.setInteractionEnabled(true);
        QVERIFY(context.triggerAction(QStringLiteral("clear_highlight")));
        QVERIFY(!scene.activeRevealBinding());
        QCOMPARE(selection.atoms(), atoms);
        QVERIFY(scene.cameraComposer()->mode() == mode);
        QVERIFY(context.stateJson()["highlight"].toString().isEmpty());
        QVERIFY(context.triggerAction(QStringLiteral("release_camera")));
        QCOMPARE(selection.atoms(), atoms);
        QCOMPARE(scene.cameraComposer()->mode().kind, app::CameraMode::Kind::Free);
        QVERIFY(context.triggerAction(QStringLiteral("clear_selection")));
        QVERIFY(selection.empty());
        QVERIFY(context.isVisible());
        QVERIFY(!context.isWindow());

        selection.bulkSet({0, 1, 2, 3});
        inspector.setPickedAtom(3);
        scene.revealBinding(binding);
        app::CsaTensorInfo shielding;
        shielding.sourceLabel = QStringLiteral("ORCA DFT");
        inspector.setCsaTensor(3, shielding);
        inspector.setTensorVisibility(true, false);
        for (const QSize size : {QSize(360, 560), QSize(300, 480)}) {
            inspector.resize(size);
            QCoreApplication::processEvents();
            QCOMPARE(inspector.size(), size);
            QVERIFY(context.isVisible());
            QCOMPARE(context.stateJson()["measurement"].toObject()["kind"].toString(), QStringLiteral("Dihedral"));
            QVERIFY(!context.stateJson()["highlight"].toString().isEmpty());
            QVERIFY(context.stateJson()["camera"].toString() != QStringLiteral("Camera: free"));
            for (auto* label : context.findChildren<QLabel*>()) {
                if (label->isVisible())
                    QVERIFY(context.rect().contains(QRect(label->mapTo(&context, QPoint()), label->size())));
            }
            for (auto* button : context.findChildren<QToolButton*>()) {
                if (button->isVisible())
                    QVERIFY(context.rect().contains(QRect(button->mapTo(&context, QPoint()), button->size())));
            }
            const QString shots = qEnvironmentVariable("H5READER_UI_TEST_SHOTS");
            if (!shots.isEmpty())
                QVERIFY(inspector.grab().save(QDir(shots).filePath(QStringLiteral("split-context-%1.png").arg(size.width()))));
        }
    }

    void measurementGeometry() {
        using model::Vec3;
        QCOMPARE(model::Distance(Vec3(0, 0, 0), Vec3(3, 4, 0)), 5.0);
        QVERIFY(std::abs(model::AngleDegrees(Vec3(1, 0, 0), Vec3(0, 0, 0),
                                             Vec3(0.5, std::sqrt(3.0) / 2.0, 0)) - 60.0) < 1e-12);
        const Vec3 a(0, 1, 0), b(0, 0, 0), c(1, 0, 0), d(1, 0, 1);
        QCOMPARE(model::DihedralDegrees(a, b, c, d), -90.0);
        QCOMPARE(model::DihedralDegrees(a, b, c, Vec3(1, 0, -1)), 90.0);
        QVERIFY(std::isnan(model::AngleDegrees(a, a, c)));
        QVERIFY(std::isnan(model::DihedralDegrees(a, b, b, d)));
    }

    void hardwarePicking_data() {
        QTest::addColumn<double>("pixelRatio");
        QTest::addColumn<bool>("parallelCamera");
        for (double ratio : {1.0, 1.5, 2.0}) {
            for (bool parallel : {false, true}) {
                const QByteArray name = QByteArray::number(ratio) + (parallel ? "-parallel" : "-perspective");
                QTest::newRow(name.constData()) << ratio << parallel;
            }
        }
    }

    void hardwarePicking() {
        QFETCH(double, pixelRatio);
        QFETCH(bool, parallelCamera);
        model::QtBond bond;
        bond.atomIndexA = 2;
        bond.atomIndexB = 3;
        bond.order = model::BondOrder::Single;
        auto protein = smallProtein(4, {bond});
        auto snapshot = std::make_shared<model::QtConformationSnapshot>(protein.get(), 0, 0.0);
        auto& positions = snapshot->mutableColumn(io::FieldKind::Pos);
        positions.present = true;
        positions.rows = 4;
        positions.cols = 3;
        positions.data = {0.25, 0.0, 0.0, 0.0, 0.0, -3.0,
                         -1.0, -2.0, 0.0, 1.0, -2.0, 0.0};
        model::SingleConformation conformation(protein.get(), snapshot);
        QVTKOpenGLNativeWidget widget;
        widget.resize(640, 480);
        widget.setCustomDevicePixelRatio(pixelRatio);
        auto window = vtkSmartPointer<vtkGenericOpenGLRenderWindow>::New();
        widget.setRenderWindow(window);
        app::MoleculeScene scene(&widget, window);
        scene.Build(*protein, conformation);
        app::QtAtomPicker picker(&widget, &scene);
        app::CameraInputFilter cameraInput(&widget, &scene, scene.cameraComposer());
        QSignalSpy picks(&picker, &app::QtAtomPicker::atomPicked);
        int emptyClicks = 0;
        connect(&cameraInput, &app::CameraInputFilter::viewportClicked, &scene, [&](QPointF point) {
            const auto hit = scene.pickAt(point);
            if (hit && hit->empty()) ++emptyClicks;
        });
        widget.show();
        QVERIFY(QTest::qWaitForWindowExposed(&widget));
        auto* camera = scene.Renderer()->GetActiveCamera();
        camera->SetPosition(0.0, 0.0, 20.0);
        camera->SetFocalPoint(0.0, 0.0, 0.0);
        camera->SetViewUp(0.0, 1.0, 0.0);
        camera->SetParallelProjection(parallelCamera);
        camera->SetParallelScale(4.0);
        scene.syncCameraClippingRange();
        window->Render();
        const QPoint overlap = widgetPosition(scene.Renderer(), widget, model::Vec3(0, 0, 0));
        const auto foreground = scene.pickAt(overlap);
        QVERIFY(foreground);
        QCOMPARE(foreground->atom, std::optional<std::size_t>{0});
        QVERIFY(!foreground->bond && !foreground->empty());
        QTest::mouseDClick(&widget, Qt::LeftButton, Qt::ShiftModifier, overlap);
        QCOMPARE(picks.count(), 1);
        QCOMPARE(picks.last().at(0).value<std::size_t>(), std::size_t{0});
        QCOMPARE(picks.last().at(1).value<Qt::KeyboardModifiers>(), Qt::KeyboardModifiers(Qt::ShiftModifier));
        QTest::mouseDClick(&widget, Qt::RightButton, Qt::NoModifier, overlap);
        QCOMPARE(picks.count(), 1);

        scene.setAtomFilter({1, 0, 3, 2});
        QCOMPARE(scene.pickAt(overlap)->atom, std::optional<std::size_t>{0});
        scene.setAtomFilter({1});
        QCOMPARE(scene.pickAt(overlap)->atom, std::optional<std::size_t>{1});
        scene.clearAtomFilter();
        const QPoint bondPoint = widgetPosition(scene.Renderer(), widget, model::Vec3(0, -2, 0));
        const auto bondHit = scene.pickAt(bondPoint);
        QVERIFY(bondHit && bondHit->bond && !bondHit->atom && !bondHit->empty());
        QTest::mouseDClick(&widget, Qt::LeftButton, Qt::NoModifier, bondPoint);
        QCOMPARE(picks.count(), 1);
        QTest::mouseClick(&widget, Qt::LeftButton, Qt::NoModifier, bondPoint);
        QCOMPARE(emptyClicks, 0);
        QTest::mouseClick(&widget, Qt::LeftButton, Qt::NoModifier, QPoint(5, 5));
        QCOMPARE(emptyClicks, 1);
        QVERIFY(!scene.pickAt(QPointF(-1, 0)));
        QVERIFY(!scene.pickAt(QPointF(widget.width(), 0)));

        auto style = scene.moleculeStyle();
        style.renderAtoms = false;
        scene.applyMoleculeStyle(style);
        QVERIFY(scene.pickAt(overlap)->empty());
        style.renderAtoms = true;
        style.atomicRadiusScaleFactor = 0.12f;
        scene.applyMoleculeStyle(style);
        QCOMPARE(scene.pickAt(overlap)->atom, std::optional<std::size_t>{1});
        style.atomicRadiusScaleFactor = 0.3f;
        scene.applyMoleculeStyle(style);

        vtkNew<vtkSphereSource> decoration;
        decoration->SetCenter(0.0, 0.0, 2.0);
        decoration->SetRadius(0.6);
        vtkNew<vtkPolyDataMapper> decorationMapper;
        decorationMapper->SetInputConnection(decoration->GetOutputPort());
        vtkNew<vtkActor> mainDecoration;
        mainDecoration->SetMapper(decorationMapper);
        mainDecoration->PickableOff();
        scene.Renderer()->AddActor(mainDecoration);
        vtkNew<vtkActor> layerDecoration;
        layerDecoration->SetMapper(decorationMapper);
        scene.OverlayRenderer()->AddActor(layerDecoration);
        window->Render();
        const QImage before = widget.grabFramebuffer();
        QSignalSpy rendered(&scene, &app::MoleculeScene::renderCompleted);
        struct PickFrames { int hidden = 0; int displayed = 0; } pickFrames;
        vtkNew<vtkCallbackCommand> frameObserver;
        frameObserver->SetClientData(&pickFrames);
        frameObserver->SetCallback([](vtkObject* caller, unsigned long, void* data, void*) {
            auto* window = vtkRenderWindow::SafeDownCast(caller);
            auto& frames = *static_cast<PickFrames*>(data);
            // QVTKRenderWindowAdapter::frame ignores exactly this condition.
            if (window->GetDoubleBuffer() && !window->GetSwapBuffers()) ++frames.hidden;
            else ++frames.displayed;
        });
        const auto observerTag = window->AddObserver(vtkCommand::WindowFrameEvent, frameObserver);
        for (int i = 0; i < 5; ++i)
            QCOMPARE(scene.pickAt(overlap)->atom, std::optional<std::size_t>{0});
        window->RemoveObserver(observerTag);
        QVERIFY(pickFrames.hidden > 0);
        QCOMPARE(pickFrames.displayed, 0);
        QCOMPARE(rendered.count(), 0);
        QCOMPARE(mainDecoration->GetVisibility(), 1);
        QCOMPARE(mainDecoration->GetPickable(), 0);
        QCOMPARE(layerDecoration->GetVisibility(), 1);
        QVERIFY(scene.OverlayRenderer()->GetDraw());
        QTRY_VERIFY(rendered.count() > 0);
        QCOMPARE(widget.grabFramebuffer(), before);
        mainDecoration->VisibilityOff();
        scene.OverlayRenderer()->DrawOff();
        QCOMPARE(scene.pickAt(overlap)->atom, std::optional<std::size_t>{0});
        QCOMPARE(mainDecoration->GetVisibility(), 0);
        QVERIFY(!scene.OverlayRenderer()->GetDraw());

        positions.data[2] = -6.0;
        scene.refreshCurrentFrame();
        QCOMPARE(scene.pickAt(overlap)->atom, std::optional<std::size_t>{1});
    }

    void filteredPickingAndRotatedMarker() {
        auto protein = smallProtein();
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
        app::QtAtomPicker picker(&widget, &scene);
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
        QCOMPARE(scene.pickAt(point)->atom, std::optional<std::size_t>{1});
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
            QCOMPARE(scene.pickAt(point)->atom, std::optional<std::size_t>{1});
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
        QCOMPARE(scene.pickAt(point)->atom, std::optional<std::size_t>{0});

        selection.bulkSet({0, 1});
        window->Render();
        vtkActor* connector = nullptr;
        actors->InitTraversal();
        while (auto* actor = actors->GetNextActor()) {
            auto* data = vtkPolyData::SafeDownCast(actor->GetMapper()->GetInput());
            if (data && data->GetNumberOfLines() == 1) {
                connector = actor;
                QCOMPARE(data->GetNumberOfPoints(), vtkIdType{2});
                for (vtkIdType i = 0; i < 2; ++i) {
                    const model::Vec3 endpoint(data->GetPoint(i));
                    QVERIFY((endpoint - conformation.atomPosition(0, std::size_t(i))).norm() < 1e-12);
                }
                break;
            }
        }
        QVERIFY(connector);
        QVERIFY(connector->GetVisibility());
        const double* color = connector->GetProperty()->GetColor();
        for (int i = 0; i < 3; ++i)
            QVERIFY(color[i] < 0.5);
    }
};

int main(int argc, char** argv) {
    QSurfaceFormat::setDefaultFormat(QVTKOpenGLNativeWidget::defaultFormat());
    QApplication application(argc, argv);
    AtomPickerTests tests;
    return QTest::qExec(&tests, argc, argv);
}

#include "atom_picker_tests.moc"
