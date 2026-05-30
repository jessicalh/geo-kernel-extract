// h5reader_trajectory_frame_map_tests — direct tests for sampled-row and
// original-frame-index mapping helpers.

#include "io/TrajectoryFrameMap.h"

#include "highfive/H5DataSpace.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5Group.hpp"

#include <QDir>
#include <QTemporaryDir>
#include <QtTest>

#include <cstdint>
#include <string>
#include <vector>

using h5reader::io::QtTrajectoryH5;
using h5reader::io::TrajectoryFrameMap;

namespace {

void writeMinimalTrajectory(const QString& path, const std::vector<std::uint64_t>& originalIndices) {
    const std::size_t frameCount = originalIndices.size();
    HighFive::File file(path.toStdString(), HighFive::File::Overwrite);
    file.createAttribute("n_atoms", std::uint64_t{1});
    file.createAttribute("protein_id", std::string("unit"));

    auto atoms = file.createGroup("/atoms");
    atoms.createDataSet("element", std::vector<std::int32_t>{1});
    atoms.createDataSet("residue_index", std::vector<std::uint64_t>{0});
    atoms.createDataSet("pdb_atom_name", std::vector<std::string>{"H"});

    std::vector<double> frameTimes(frameCount, 0.0);
    for (std::size_t frame = 0; frame < frameCount; ++frame)
        frameTimes[frame] = static_cast<double>(frame);

    auto frames = file.createGroup("/trajectory/frames");
    frames.createDataSet("time_ps", frameTimes);
    frames.createDataSet("original_index", originalIndices);

    auto positions = file.createGroup("/trajectory/positions");
    std::vector<double> xyz(frameCount * 3, 0.0);
    for (std::size_t frame = 0; frame < frameCount; ++frame)
        xyz[frame * 3] = static_cast<double>(frame);
    auto ds = positions.createDataSet<double>("xyz", HighFive::DataSpace({1, frameCount, 3}));
    ds.write_raw(xyz.data());
    positions.createDataSet("frame_times", frameTimes);
    positions.createDataSet("frame_indices", originalIndices);
}

}  // namespace

class TrajectoryFrameMapTests : public QObject {
    Q_OBJECT

private slots:
    void emptyPerFrameNpysDirReturnsEmptyRows();
    void scanSampledRowsReturnsSortedH5Rows();
    void originalIndexUsesFrameMapAndIdentityFallback();
};

void TrajectoryFrameMapTests::emptyPerFrameNpysDirReturnsEmptyRows() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    const QString h5Path = dir.filePath(QStringLiteral("trajectory.h5"));
    writeMinimalTrajectory(h5Path, {20, 10, 30});
    const QtTrajectoryH5 h5(h5Path);

    const auto rows = TrajectoryFrameMap::ScanSampledRows(QString(), h5);

    QVERIFY(rows.empty());
}

void TrajectoryFrameMapTests::scanSampledRowsReturnsSortedH5Rows() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    const QString h5Path = dir.filePath(QStringLiteral("trajectory.h5"));
    writeMinimalTrajectory(h5Path, {20, 10, 30, 40});
    const QtTrajectoryH5 h5(h5Path);

    QDir root(dir.path());
    QVERIFY(root.mkpath(QStringLiteral("npys/frame_000030")));
    QVERIFY(root.mkpath(QStringLiteral("npys/frame_000020")));
    QVERIFY(root.mkpath(QStringLiteral("npys/frame_000010")));
    QVERIFY(root.mkpath(QStringLiteral("npys/frame_999999")));
    QVERIFY(root.mkpath(QStringLiteral("npys/frame_bad")));

    const auto rows = TrajectoryFrameMap::ScanSampledRows(dir.filePath(QStringLiteral("npys")), h5);

    QCOMPARE(rows, (std::vector<std::size_t>{0, 1, 2}));
}

void TrajectoryFrameMapTests::originalIndexUsesFrameMapAndIdentityFallback() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    const QString h5Path = dir.filePath(QStringLiteral("trajectory.h5"));
    writeMinimalTrajectory(h5Path, {20, 10, 30});
    const QtTrajectoryH5 h5(h5Path);

    QCOMPARE(TrajectoryFrameMap::OriginalIndex(0, h5), std::size_t{20});
    QCOMPARE(TrajectoryFrameMap::OriginalIndex(1, h5), std::size_t{10});
    QCOMPARE(TrajectoryFrameMap::OriginalIndex(99, h5), std::size_t{99});
}

QTEST_GUILESS_MAIN(TrajectoryFrameMapTests)

#include "trajectory_frame_map_tests.moc"
