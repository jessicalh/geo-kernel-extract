// h5reader_app_tests — app-controller robustness tests.

#include "app/DashboardDisplayController.h"
#include "model/Conformation.h"
#include "model/QtConformationSnapshot.h"
#include "model/DashboardPanelModel.h"
#include "model/DashboardSignalModel.h"
#include "model/DftShieldingStore.h"
#include "model/SignalTimeSeries.h"
#include "model/TrajectorySignalCatalog.h"

#include <QTemporaryDir>
#include <QtTest>

#include <cstddef>
#include <cmath>
#include <memory>
#include <optional>
#include <stdexcept>
#include <vector>

using namespace h5reader;

namespace {

template <typename Provider>
bool waitUntilIdle(Provider& provider) {
    QSignalSpy idle(&provider, &Provider::becameIdle);
    return !provider.isBusy() || idle.wait();
}

class CountingConformation final : public model::Conformation {
public:
    explicit CountingConformation(std::size_t frames, bool failReads = false)
        : model::Conformation(nullptr),
          frames_(frames), failReads_(failReads) {}

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

    mutable int snapshotRequests = 0;
    mutable std::vector<std::size_t> requestedFrames;

protected:
    SnapshotReader snapshotReader(std::size_t frame) const override {
        ++snapshotRequests;
        requestedFrames.push_back(frame);
        return [snapshot = std::shared_ptr<const model::QtConformationSnapshot>{},
                fail = failReads_] {
            if (fail)
                throw std::runtime_error("Test snapshot read failure");
            return snapshot;
        };
    }

private:
    std::size_t frames_ = 0;
    bool failReads_ = false;
};

class FieldConformation final : public model::Conformation {
public:
    explicit FieldConformation(double firstValue = 1.0, std::size_t originalStride = 1,
                               std::size_t frames = 4)
        : model::Conformation(nullptr), firstValue_(firstValue), originalStride_(originalStride), frames_(frames) {}
    std::size_t frameCount() const override { return frames_; }
    double timePicoseconds(std::size_t frame) const override { return double(frame); }
    model::Vec3 atomPosition(std::size_t, std::size_t) const override {
        return model::Vec3::Zero();
    }
    std::size_t originalFrameIndex(std::size_t frame) const override {
        return frame * originalStride_;
    }
    std::optional<std::size_t> frameRowForOriginalIndex(std::size_t original) const override {
        if (original % originalStride_ != 0 || original / originalStride_ >= frameCount())
            return std::nullopt;
        return original / originalStride_;
    }

    mutable std::vector<std::size_t> requestedFrames;

protected:
    SnapshotReader snapshotReader(std::size_t frame) const override {
        requestedFrames.push_back(frame);
        auto snapshot = std::make_shared<model::QtConformationSnapshot>(nullptr, frame, double(frame));
        auto& column = snapshot->mutableColumn(io::FieldKind::BSTotalB);
        column.present = true;
        column.rows = 2;
        column.cols = 3;
        column.data = {firstValue_ + double(frame), 0.0, 0.0,
                       firstValue_ + 99.0 + double(frame), 0.0, 0.0};
        return [snapshot] { return snapshot; };
    }

private:
    double firstValue_;
    std::size_t originalStride_;
    std::size_t frames_;
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
    void snapshotCompletionFillsPendingSamples();
    void everyPlaybackFrameReachesTheStrip();
    void snapshotCompletionResolvesAbsentInput_data();
    void snapshotCompletionResolvesAbsentInput();
    void providerCompletionsStayIndependent();
    void replacingContextDisconnectsOldCompletions();
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
    QVERIFY(waitUntilIdle(conformation));

    controller.setScrubActive(true);
    controller.setFrame(3);
    controller.setScrubActive(false);
    QVERIFY(waitUntilIdle(conformation));
    const auto& values = controller.stripTracks()[0].buffer->values;
    QCOMPARE(values.size(), std::size_t{4});
    QCOMPARE(values[0], 1.0);
    QVERIFY(std::isnan(values[1]));
    QVERIFY(std::isnan(values[2]));
    QCOMPARE(values[3], 4.0);

    controller.setFrame(1);
    QVERIFY(waitUntilIdle(conformation));
    QCOMPARE(values[1], 2.0);
    QVERIFY(std::isnan(values[2]));
    controller.setFrame(2);
    QVERIFY(waitUntilIdle(conformation));
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
    QVERIFY(waitUntilIdle(first));
    QCOMPARE(controller.stripTracks()[0].buffer->values, (std::vector<double>{1.0, 2.0, 3.0}));

    controller.setContext(nullptr, &second);
    QVERIFY(waitUntilIdle(second));
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
    QVERIFY(waitUntilIdle(conformation));
    auto tracks = controller.stripTracks();
    QCOMPARE(tracks.size(), 3);
    QCOMPARE(tracks[0].buffer->values, (std::vector<double>{1.0, 2.0, 3.0}));

    auto binding = signalModel.signalById(id)->binding;
    binding.anchor = model::AtomAnchor{1};
    QVERIFY(signalModel.updateBinding(id, binding));
    QVERIFY(waitUntilIdle(conformation));
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
    QVERIFY(waitUntilIdle(conformation));
    conformation.resetCounts();

    controller.setScrubActive(true);
    controller.setFrame(750);
    QCOMPARE(conformation.snapshotRequests, 0);
    QVERIFY(conformation.requestedFrames.empty());

    controller.setScrubActive(false);
    QCOMPARE(conformation.snapshotRequests, 1);
    QCOMPARE(conformation.requestedFrames.size(), std::size_t{1});
    QCOMPARE(conformation.requestedFrames.front(), std::size_t{750});
    QVERIFY(waitUntilIdle(conformation));
    QCOMPARE(controller.smokeSummary().pendingGapSamples, 3LL * 749);
    QCOMPARE(controller.smokeSummary().frameSourceAbsentGapSamples, 6LL);
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
    QVERIFY(waitUntilIdle(conformation));
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

void DashboardControllerTests::everyPlaybackFrameReachesTheStrip() {
    FieldConformation conformation(1.0, 15, 100);
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    app::DashboardDisplayController controller;
    const auto* descriptor = catalog.findDescriptor(QStringLiteral("npy:bs_total_B"));
    QVERIFY(descriptor);
    signalModel.addSignal(*descriptor, model::AtomAnchor{0}, QString(),
                          {QStringLiteral("strip.vector.component")});
    controller.setContext(nullptr, &conformation);
    controller.setSignalModels(&catalog, &signalModel);
    for (int frame = 0; frame < 100; ++frame)
        controller.setFrame(frame);
    QTRY_VERIFY_WITH_TIMEOUT(!conformation.isBusy(), 15000);
    QCOMPARE(conformation.requestedFrames.size(), std::size_t(100));
    QCOMPARE(controller.smokeSummary().pendingGapSamples, 0LL);
    QCOMPARE(controller.smokeSummary().validSamples, 300LL);
    const auto tracks = controller.stripTracks();
    QCOMPARE(tracks[0].buffer->values.size(), std::size_t(100));
    for (std::size_t frame = 0; frame < 100; ++frame) {
        QCOMPARE(conformation.requestedFrames[frame], frame);
        QCOMPARE(tracks[0].buffer->values[frame], 1.0 + double(frame));
    }
}

void DashboardControllerTests::snapshotCompletionFillsPendingSamples() {
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
    controller.setFrame(3);
    controller.rebuild();
    controller.setFrame(3);

    QCOMPARE(controller.smokeSummary().pendingGapSamples, 12LL);
    QSignalSpy changed(&controller, &app::DashboardDisplayController::stripTracksChanged);
    QVERIFY(waitUntilIdle(conformation));
    QCOMPARE(changed.count(), 4);
    QCOMPARE(conformation.requestedFrames, (std::vector<std::size_t>{0, 1, 2, 3}));
    QCOMPARE(controller.smokeSummary().pendingGapSamples, 0LL);
    QCOMPARE(controller.smokeSummary().validSamples, 12LL);
    const auto tracks = controller.stripTracks();
    QCOMPARE(tracks[0].buffer->values, (std::vector<double>{1.0, 2.0, 3.0, 4.0}));
    QCOMPARE(tracks[0].buffer->yMin, 1.0);
    QCOMPARE(tracks[0].buffer->yMax, 4.0);
    QVERIFY(!conformation.snapshot(0));
    QVERIFY(conformation.snapshot(3));

    // An unrelated consumer loading an old frame must not erase strip history.
    conformation.requestSnapshotAsync(0);
    QVERIFY(waitUntilIdle(conformation));
    QCOMPARE(changed.count(), 4);
    QCOMPARE(tracks[0].buffer->values, (std::vector<double>{1.0, 2.0, 3.0, 4.0}));
}

void DashboardControllerTests::snapshotCompletionResolvesAbsentInput_data() {
    QTest::addColumn<bool>("failReads");
    QTest::newRow("absent") << false;
    QTest::newRow("failed") << true;
}

void DashboardControllerTests::snapshotCompletionResolvesAbsentInput() {
    QFETCH(bool, failReads);
    CountingConformation conformation(4, failReads);
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    app::DashboardDisplayController controller;
    const auto* descriptor = catalog.findDescriptor(QStringLiteral("npy:bs_total_B"));
    QVERIFY(descriptor);
    signalModel.addSignal(*descriptor, model::AtomAnchor{0}, QString(),
                          {QStringLiteral("strip.vector.component")});
    controller.setContext(nullptr, &conformation);
    controller.setSignalModels(&catalog, &signalModel);
    controller.setFrame(3);
    QCOMPARE(controller.smokeSummary().pendingGapSamples, 12LL);

    QVERIFY(waitUntilIdle(conformation));
    const auto summary = controller.smokeSummary();
    QCOMPARE(summary.pendingGapSamples, 0LL);
    QCOMPARE(summary.frameSourceAbsentGapSamples, 12LL);
    QCOMPARE(summary.validSamples, 0LL);
    QCOMPARE(conformation.requestedFrames, (std::vector<std::size_t>{0, 1, 2, 3}));
}

void DashboardControllerTests::providerCompletionsStayIndependent() {
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    io::DftFrame failedJob;
    failedJob.frame_index = 20;
    failedJob.meta_json_abspath = directory.filePath(QStringLiteral("missing-meta.json"));
    model::DftShieldingStore dftStore(nullptr, {failedJob});
    FieldConformation conformation(1.0, 10);
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    app::DashboardDisplayController controller;
    const auto* field = catalog.findDescriptor(QStringLiteral("npy:bs_total_B"));
    const auto* dft = catalog.findDescriptor(QStringLiteral("orca_dft:total"));
    const auto* ml = catalog.findDescriptor(QStringLiteral("ml:experimental_shielding_t2"));
    QVERIFY(field);
    QVERIFY(dft);
    QVERIFY(ml);
    signalModel.addSignal(*field, model::AtomAnchor{0}, QString(),
                          {QStringLiteral("strip.vector.component")});
    signalModel.addSignal(*dft, model::AtomAnchor{0}, QString(),
                          {QStringLiteral("strip.tensor.T0")});
    signalModel.addSignal(*ml, model::AtomAnchor{0}, QString(),
                          {QStringLiteral("strip.tensor.T2")});
    controller.setContext(nullptr, &conformation);
    controller.setSignalModels(&catalog, &signalModel);
    // Attaching the provider replaces earlier SourceAbsent samples.
    controller.setDftStore(&dftStore);
    controller.setFrame(3);

    const auto pending = controller.smokeSummary();
    QCOMPARE(pending.pendingGapSamples, 14LL); // NPY: 12; DFT: original frames 20, 30.
    QCOMPARE(pending.frameSourceAbsentGapSamples, 2LL); // Inline absent DFT 0, 10.
    QCOMPARE(pending.sourceAbsentGapSamples, 4LL); // No ML provider.
    QVERIFY(waitUntilIdle(conformation));
    QVERIFY(waitUntilIdle(dftStore));
    const auto resolved = controller.smokeSummary();
    QCOMPARE(resolved.pendingGapSamples, 0LL);
    QCOMPARE(resolved.validSamples, 12LL);
    QCOMPARE(resolved.orcaDftFrameSourceAbsentGapSamples, 4LL);
    QCOMPARE(resolved.sourceAbsentGapSamples, 4LL);
    QVERIFY(dftStore.hasFailedFrame(20));
    QCOMPARE(controller.stripTracks()[0].buffer->values,
             (std::vector<double>{1.0, 2.0, 3.0, 4.0}));
}

void DashboardControllerTests::replacingContextDisconnectsOldCompletions() {
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
    QVERIFY(waitUntilIdle(first));

    controller.setContext(nullptr, &second);
    QVERIFY(waitUntilIdle(second));
    controller.setScrubActive(true);
    controller.setFrame(3);
    controller.setScrubActive(false);
    QVERIFY(waitUntilIdle(second));
    QCOMPARE(controller.smokeSummary().pendingGapSamples, 6LL);

    QSignalSpy changed(&controller, &app::DashboardDisplayController::stripTracksChanged);
    first.requestSnapshotAsync(1);
    QVERIFY(waitUntilIdle(first));
    QCOMPARE(changed.count(), 0);
    QCOMPARE(controller.smokeSummary().pendingGapSamples, 6LL);
    QCOMPARE(controller.stripTracks()[0].buffer->values[0], 20.0);
    QCOMPARE(controller.stripTracks()[0].buffer->values[3], 23.0);
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
