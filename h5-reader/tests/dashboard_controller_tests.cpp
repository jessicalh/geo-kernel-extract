// h5reader_app_tests — app-controller robustness tests.

#include "app/DashboardDisplayController.h"
#include "model/Conformation.h"
#include "model/QtConformationSnapshot.h"
#include "model/DashboardPanelModel.h"
#include "model/DashboardSignalModel.h"
#include "model/SignalTimeSeries.h"
#include "model/TrajectorySignalCatalog.h"

#include <QtTest>

#include <cstddef>
#include <cmath>
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

class FieldConformation final : public model::Conformation {
public:
    explicit FieldConformation(double firstValue = 1.0)
        : model::Conformation(nullptr), firstValue_(firstValue) {}
    std::size_t frameCount() const override { return 4; }
    double timePicoseconds(std::size_t frame) const override { return double(frame); }
    model::Vec3 atomPosition(std::size_t, std::size_t) const override {
        return model::Vec3::Zero();
    }

protected:
    std::shared_ptr<const model::QtConformationSnapshot> loadSnapshot(std::size_t frame) override {
        auto snapshot = std::make_shared<model::QtConformationSnapshot>(nullptr, frame, double(frame));
        auto& column = snapshot->mutableColumn(io::FieldKind::BSTotalB);
        column.present = true;
        column.rows = 2;
        column.cols = 3;
        column.data = {firstValue_ + double(frame), 0.0, 0.0,
                       firstValue_ + 99.0 + double(frame), 0.0, 0.0};
        return snapshot;
    }

private:
    double firstValue_;
};

}  // namespace

class DashboardControllerTests : public QObject {
    Q_OBJECT

private slots:
    void scrubDefersFrameSnapshotRequestsUntilRelease();
    void stripHistorySurvivesRebuildByModeId();
    void replacingPendingSampleRecomputesValidityAndRange();
    void f003TensorBindingTracksActivePanelReference();
    void retargetingStripRecomputesHistory();
    void newConformationStartsAtFrameZero();
    void revisitingScrubbedFramesFillsPendingSamples();
};

void DashboardControllerTests::revisitingScrubbedFramesFillsPendingSamples() {
    FieldConformation conformation;
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    app::DashboardDisplayController controller;
    const auto* descriptor = catalog.findDescriptor(QStringLiteral("npy:bs_total_B"));
    QVERIFY(descriptor);
    signalModel.addSignal(*descriptor, model::AtomAnchor{0}, QString(),
                          {QStringLiteral("strip.vector.component")});
    controller.setContext(nullptr, &conformation);
    controller.setSignalModels(&catalog, &signalModel);

    controller.setScrubActive(true);
    controller.setFrame(3);
    controller.setScrubActive(false);
    const auto& values = controller.stripTracks()[0].buffer->values;
    QCOMPARE(values.size(), std::size_t{4});
    QCOMPARE(values[0], 1.0);
    QVERIFY(std::isnan(values[1]));
    QVERIFY(std::isnan(values[2]));
    QCOMPARE(values[3], 4.0);

    controller.setFrame(1);
    QCOMPARE(values[1], 2.0);
    QVERIFY(std::isnan(values[2]));
    controller.setFrame(2);
    QCOMPARE(values, (std::vector<double>{1.0, 2.0, 3.0, 4.0}));
}

void DashboardControllerTests::newConformationStartsAtFrameZero() {
    FieldConformation first, second(20.0);
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    app::DashboardDisplayController controller;
    const auto* descriptor = catalog.findDescriptor(QStringLiteral("npy:bs_total_B"));
    QVERIFY(descriptor);
    signalModel.addSignal(*descriptor, model::AtomAnchor{0}, QString(),
                          {QStringLiteral("strip.vector.component")});
    controller.setContext(nullptr, &first);
    controller.setSignalModels(&catalog, &signalModel);
    controller.setFrame(2);
    QCOMPARE(controller.stripTracks()[0].buffer->values, (std::vector<double>{1.0, 2.0, 3.0}));

    controller.setContext(nullptr, &second);
    QCOMPARE(controller.stripTracks()[0].buffer->values, (std::vector<double>{20.0}));
}

void DashboardControllerTests::retargetingStripRecomputesHistory() {
    FieldConformation conformation;
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    app::DashboardDisplayController controller;
    const auto* descriptor = catalog.findDescriptor(QStringLiteral("npy:bs_total_B"));
    QVERIFY(descriptor);
    const QUuid id = signalModel.addSignal(*descriptor, model::AtomAnchor{0}, QString(),
                                      {QStringLiteral("strip.vector.component")});
    controller.setContext(nullptr, &conformation);
    controller.setSignalModels(&catalog, &signalModel);
    controller.setFrame(2);
    auto tracks = controller.stripTracks();
    QCOMPARE(tracks.size(), 3);
    QCOMPARE(tracks[0].buffer->values, (std::vector<double>{1.0, 2.0, 3.0}));

    auto binding = signalModel.signalById(id)->binding;
    binding.anchor = model::AtomAnchor{1};
    QVERIFY(signalModel.updateBinding(id, binding));
    tracks = controller.stripTracks();
    QCOMPARE(tracks.size(), 3);
    QCOMPARE(tracks[0].buffer->values, (std::vector<double>{100.0, 101.0, 102.0}));
}

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

void DashboardControllerTests::stripHistorySurvivesRebuildByModeId() {
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
