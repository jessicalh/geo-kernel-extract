// h5reader_app_tests — app-controller robustness tests.

#include "app/DashboardDisplayController.h"
#include "model/Conformation.h"
#include "model/DashboardPanelModel.h"
#include "model/DashboardSignalModel.h"
#include "model/SignalTimeSeries.h"
#include "model/TrajectorySignalCatalog.h"

#include <QtTest>

#include <cstddef>
#include <memory>
#include <vector>

using namespace h5reader;

namespace {

class CountingConformation final : public model::Conformation {
public:
    explicit CountingConformation(std::size_t frames)
        : model::Conformation(nullptr),
          frames_(frames) {}

    std::size_t frameCount() const override { return frames_; }
    double timePicoseconds(std::size_t frame) const override {
        return static_cast<double>(frame);
    }
    model::Vec3 atomPosition(std::size_t, std::size_t) const override {
        return model::Vec3::Zero();
    }

    void resetCounts() {
        snapshotRequests = 0;
        requestedFrames.clear();
    }

    int snapshotRequests = 0;
    std::vector<std::size_t> requestedFrames;

protected:
    std::shared_ptr<const model::QtConformationSnapshot> loadSnapshot(std::size_t frame) override {
        ++snapshotRequests;
        requestedFrames.push_back(frame);
        return nullptr;
    }

private:
    std::size_t frames_ = 0;
};

}  // namespace

class DashboardControllerTests : public QObject {
    Q_OBJECT

private slots:
    void scrubDefersFrameSnapshotRequestsUntilRelease();
    void stripHistorySurvivesRebuildByLegacyModeId();
    void replacingPendingSampleRecomputesValidityAndRange();
    void f003TensorBindingTracksActivePanelReference();
};

void DashboardControllerTests::scrubDefersFrameSnapshotRequestsUntilRelease() {
    CountingConformation conformation(1000);
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    app::DashboardDisplayController controller;

    const model::SignalDescriptor* descriptor =
        catalog.findDescriptor(QStringLiteral("npy:bs_total_B"));
    QVERIFY(descriptor != nullptr);
    signalModel.addSignal(*descriptor,
                          model::AtomAnchor{0},
                          QString(),
                          {QStringLiteral("strip.vector.component")},
                          false,
                          QStringLiteral("Frame-local magnetic field"));

    controller.setContext(nullptr, &conformation);
    controller.setSignalModels(&catalog, &signalModel);
    conformation.resetCounts();

    controller.setScrubActive(true);
    controller.setFrame(750);
    QCOMPARE(conformation.snapshotRequests, 0);
    QVERIFY(conformation.requestedFrames.empty());

    controller.setScrubActive(false);
    QCOMPARE(conformation.snapshotRequests, 1);
    QCOMPARE(conformation.requestedFrames.size(), std::size_t{1});
    QCOMPARE(conformation.requestedFrames.front(), std::size_t{750});
}

void DashboardControllerTests::stripHistorySurvivesRebuildByLegacyModeId() {
    CountingConformation conformation(1000);
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    app::DashboardDisplayController controller;

    const model::SignalDescriptor* descriptor =
        catalog.findDescriptor(QStringLiteral("npy:bs_total_B"));
    QVERIFY(descriptor != nullptr);
    signalModel.addSignal(*descriptor,
                          model::AtomAnchor{0},
                          QString(),
                          {QStringLiteral("strip.vector.component")},
                          false,
                          QStringLiteral("Frame-local magnetic field"));

    controller.setContext(nullptr, &conformation);
    controller.setSignalModels(&catalog, &signalModel);
    controller.setFrame(3);
    const app::DashboardSmokeSummary before = controller.smokeSummary();
    QCOMPARE(before.seriesSparseness.size(), 3);
    QCOMPARE(before.seriesSparseness.front().samples, 4);

    controller.rebuild();
    const app::DashboardSmokeSummary after = controller.smokeSummary();
    QCOMPARE(after.seriesSparseness.size(), before.seriesSparseness.size());
    for (int i = 0; i < after.seriesSparseness.size(); ++i) {
        QCOMPARE(after.seriesSparseness.at(i).displayModeId,
                 before.seriesSparseness.at(i).displayModeId);
        QCOMPARE(after.seriesSparseness.at(i).channelId,
                 before.seriesSparseness.at(i).channelId);
        QCOMPARE(after.seriesSparseness.at(i).samples,
                 before.seriesSparseness.at(i).samples);
    }
}

void DashboardControllerTests::replacingPendingSampleRecomputesValidityAndRange() {
    model::SignalBuffer buffer;
    buffer.append(model::FrameSignalSample::Valid(1.0));
    buffer.append(model::FrameSignalSample::Valid(5.0));
    buffer.append(model::FrameSignalSample::Gap(model::GapReason::Pending));

    QVERIFY(buffer.channel.hasRange);
    QCOMPARE(buffer.channel.yMin, 1.0);
    QCOMPARE(buffer.channel.yMax, 5.0);

    buffer.replace(1, model::FrameSignalSample::Gap(model::GapReason::FrameSourceAbsent));
    QVERIFY(!buffer.isValidAt(1));
    QCOMPARE(buffer.channel.yMin, 1.0);
    QCOMPARE(buffer.channel.yMax, 1.0);

    buffer.replace(2, model::FrameSignalSample::Valid(-2.0));
    QVERIFY(buffer.isValidAt(2));
    QCOMPARE(buffer.channel.yMin, -2.0);
    QCOMPARE(buffer.channel.yMax, 1.0);
    QCOMPARE(buffer.statuses[2], model::SampleStatus::Valid);
    QCOMPARE(buffer.gapReasons[2], model::GapReason::None);
}

void DashboardControllerTests::f003TensorBindingTracksActivePanelReference() {
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    model::DashboardPanelModel panelModel;
    app::DashboardDisplayController controller;

    controller.setPanelModel(&panelModel);
    controller.setSignalModels(&catalog, &signalModel);

    const model::SignalDescriptor* descriptor =
        catalog.findDescriptor(QStringLiteral("ml:experimental_shielding_t2"));
    QVERIFY(descriptor != nullptr);

    QSignalSpy bindingSpy(
        &controller,
        &app::DashboardDisplayController::sceneTensorBindingChanged);
    QVERIFY(bindingSpy.isValid());

    const QUuid signalId =
        signalModel.addSignal(*descriptor,
                              model::AtomAnchor{16},
                              QString(),
                              {QStringLiteral("static.tensor")});
    QVERIFY(!signalId.isNull());
    QCOMPARE(bindingSpy.count(), 0);

    const model::DashboardDisplayRef ref{
        signalId,
        QStringLiteral("static.tensor"),
        QStringLiteral("panel")};
    QVERIFY(panelModel.addDisplayRef(panelModel.activePanelId(), ref));
    QCOMPARE(bindingSpy.count(), 1);
    QCOMPARE(bindingSpy.at(0).at(0).toString(),
             QStringLiteral("ml:experimental_shielding_t2"));
    QCOMPARE(bindingSpy.at(0).at(1).toLongLong(), qint64{16});

    QVERIFY(panelModel.removeDisplayRef(panelModel.activePanelId(), ref));
    QCOMPARE(bindingSpy.count(), 2);
    QVERIFY(bindingSpy.at(1).at(0).toString().isEmpty());
    QCOMPARE(bindingSpy.at(1).at(1).toLongLong(), qint64{-1});
}

QTEST_GUILESS_MAIN(DashboardControllerTests)

#include "dashboard_controller_tests.moc"
