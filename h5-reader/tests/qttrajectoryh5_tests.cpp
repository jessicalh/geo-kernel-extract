// h5reader_io_tests — QtTrajectoryH5 boundary robustness tests.

#include "io/QtTrajectoryH5.h"

#include "highfive/H5DataSpace.hpp"
#include "highfive/H5File.hpp"
#include "highfive/H5Group.hpp"

#include <QTemporaryDir>
#include <QtTest>

#include <cstdint>
#include <string>
#include <vector>

using h5reader::io::QtTrajectoryH5;

namespace {

void writeMinimalTrajectory(const QString& path, bool includePositions = true) {
    HighFive::File file(path.toStdString(), HighFive::File::Overwrite);
    file.createAttribute("n_atoms", std::uint64_t{1});
    file.createAttribute("protein_id", std::string("unit"));

    auto atoms = file.createGroup("/atoms");
    atoms.createDataSet("element", std::vector<std::int32_t>{1});
    atoms.createDataSet("residue_index", std::vector<std::uint64_t>{0});
    atoms.createDataSet("pdb_atom_name", std::vector<std::string>{"H"});

    auto frames = file.createGroup("/trajectory/frames");
    frames.createDataSet("time_ps", std::vector<double>{0.0, 1.0});
    frames.createDataSet("original_index", std::vector<std::uint64_t>{0, 1});

    if (includePositions) {
        auto positions = file.createGroup("/trajectory/positions");
        const std::vector<double> xyz{
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
        };
        auto ds = positions.createDataSet<double>("xyz", HighFive::DataSpace({1, 2, 3}));
        ds.write_raw(xyz.data());
        positions.createDataSet("frame_times", std::vector<double>{0.0, 1.0});
        positions.createDataSet("frame_indices", std::vector<std::uint64_t>{0, 1});
    }

}

void addOptionalShieldingGroupWithMalformedFrameMeta(const QString& path) {
    HighFive::File file(path.toStdString(), HighFive::File::ReadWrite);
    auto grp = file.createGroup("/trajectory/bs_shielding_time_series");
    const std::vector<double> xyz(1 * 2 * 9, 0.0);
    auto ds = grp.createDataSet<double>("xyz", HighFive::DataSpace({1, 2, 9}));
    ds.write_raw(xyz.data());
    grp.createDataSet("frame_times", std::vector<std::string>{"bad", "bad"});
}

}  // namespace

class QtTrajectoryH5Tests : public QObject {
    Q_OBJECT

private slots:
    void optionalMalformedReaderDoesNotAbortLoad();
    void missingPositionsIsHardLoadError();
};

void QtTrajectoryH5Tests::optionalMalformedReaderDoesNotAbortLoad() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    const QString path = dir.filePath(QStringLiteral("trajectory.h5"));
    writeMinimalTrajectory(path);
    addOptionalShieldingGroupWithMalformedFrameMeta(path);

    QtTrajectoryH5 h5(path);

    QCOMPARE(h5.atomCount(), std::size_t{1});
    QCOMPARE(h5.frameCount(), std::size_t{2});
    QVERIFY(h5.positions() != nullptr);
    QVERIFY(h5.bsShielding() == nullptr);
}

void QtTrajectoryH5Tests::missingPositionsIsHardLoadError() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    const QString path = dir.filePath(QStringLiteral("trajectory.h5"));
    writeMinimalTrajectory(path, /*includePositions=*/false);

    try {
        const QtTrajectoryH5 h5(path);
        Q_UNUSED(h5);
        QFAIL("QtTrajectoryH5 accepted a trajectory without /trajectory/positions");
    } catch (const std::runtime_error& e) {
        const QString message = QString::fromUtf8(e.what());
        QVERIFY2(message.contains(QStringLiteral("/trajectory/positions")),
                 qPrintable(message));
        QVERIFY2(message.contains(QStringLiteral("positions buffer is null")),
                 qPrintable(message));
    }
}

QTEST_GUILESS_MAIN(QtTrajectoryH5Tests)

#include "qttrajectoryh5_tests.moc"
