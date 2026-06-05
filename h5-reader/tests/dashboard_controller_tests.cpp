// h5reader_app_tests — app-controller robustness tests.

#include "app/DashboardDisplayController.h"
#include "model/Conformation.h"
#include "model/DashboardSignalModel.h"
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
};

void DashboardControllerTests::scrubDefersFrameSnapshotRequestsUntilRelease() {
    CountingConformation conformation(1000);
    model::TrajectorySignalCatalog catalog;
    model::DashboardSignalModel signalModel;
    app::DashboardDisplayController controller;

    const model::SignalDescriptor* descriptor =
        catalog.findDescriptor(QStringLiteral("npy:pos"));
    QVERIFY(descriptor != nullptr);
    signalModel.addSignal(*descriptor,
                          model::AtomAnchor{0},
                          QString(),
                          {QStringLiteral("strip.vector.component")},
                          false,
                          QStringLiteral("Snapshot positions"));

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
        catalog.findDescriptor(QStringLiteral("npy:pos"));
    QVERIFY(descriptor != nullptr);
    signalModel.addSignal(*descriptor,
                          model::AtomAnchor{0},
                          QString(),
                          {QStringLiteral("strip.vector.component")},
                          false,
                          QStringLiteral("Snapshot positions"));

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

QTEST_GUILESS_MAIN(DashboardControllerTests)

#include "dashboard_controller_tests.moc"
