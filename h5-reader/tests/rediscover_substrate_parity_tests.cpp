#include "rediscover/SubstrateParity.h"
#include "rediscover/SphericalBasis.h"

#include <QByteArray>
#include <QDir>
#include <QSaveFile>
#include <QTemporaryDir>
#include <QtTest>

#include <array>
#include <vector>

using h5reader::model::Mat3;
using h5reader::rediscover::AuditAllAtomEquivariantSidecars;
using h5reader::rediscover::AuditPerAtomSubstrateSidecars;
using h5reader::rediscover::DecomposeLibrary;
using h5reader::rediscover::ProducerTargetLookup;

namespace {

bool writeText(const QString& path, const QString& text) {
    QSaveFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    f.write(text.toUtf8());
    return f.commit();
}

template <typename T>
bool writeNpy(const QString& path, const std::vector<std::size_t>& shape,
              const std::vector<T>& data, const QByteArray& descr) {
    std::size_t n = 1;
    for (std::size_t dim : shape) n *= dim;
    if (n != data.size()) return false;

    QByteArray header;
    header += "{'descr': '";
    header += descr;
    header += "', 'fortran_order': False, 'shape': (";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        if (i) header += ", ";
        header += QByteArray::number(static_cast<qulonglong>(shape[i]));
    }
    if (shape.size() == 1) header += ",";
    header += "), }";
    constexpr int kPreambleBytes = 10;
    const int pad = (16 - ((kPreambleBytes + header.size() + 1) % 16)) % 16;
    header += QByteArray(pad, ' ');
    header += '\n';

    QByteArray prefix;
    prefix.append("\x93NUMPY", 6);
    prefix.append(char(1));
    prefix.append(char(0));
    const quint16 headerLen = static_cast<quint16>(header.size());
    prefix.append(char(headerLen & 0xff));
    prefix.append(char((headerLen >> 8) & 0xff));

    QSaveFile f(path);
    if (!f.open(QIODevice::WriteOnly)) return false;
    if (f.write(prefix) != prefix.size()) return false;
    if (f.write(header) != header.size()) return false;
    const qsizetype bytes = static_cast<qsizetype>(data.size() * sizeof(T));
    if (bytes > 0 && f.write(reinterpret_cast<const char*>(data.data()), bytes) != bytes)
        return false;
    return f.commit();
}

Mat3 targetTensor() {
    Mat3 m;
    m << 10.0, 1.0, 2.0,
         3.0, 20.0, 4.0,
         5.0, 6.0, 30.0;
    return m;
}

ProducerTargetLookup fixedLookup(const Mat3& m) {
    return [m](std::size_t, std::size_t, QString* reason) -> std::optional<Mat3> {
        if (reason) reason->clear();
        return m;
    };
}

void writePerAtomFixture(const QString& dir, bool badRowContract = false) {
    const Mat3 raw = targetTensor();
    const auto decomp = DecomposeLibrary(raw);
    const std::size_t rowId = badRowContract ? 5 : 1;
    QVERIFY(writeText(QStringLiteral("%1/per_atom_substrate_rows.csv").arg(dir),
                      QStringLiteral("row_id,atom_index,h5_row,frame_slot,original_frame_index,dft_present\n"
                                     "0,0,7,0,77,1\n"
                                     "%1,1,7,0,77,1\n").arg(rowId)));
    QVERIFY(writeNpy<double>(QStringLiteral("%1/per_atom_substrate_target_T0.npy").arg(dir),
                             {2}, {decomp.T0, decomp.T0}, QByteArray("<f8")));
    QVERIFY(writeNpy<double>(QStringLiteral("%1/per_atom_substrate_target_T1.npy").arg(dir),
                             {2, 3},
                             {decomp.T1[0], decomp.T1[1], decomp.T1[2],
                              decomp.T1[0], decomp.T1[1], decomp.T1[2]},
                             QByteArray("<f8")));
    QVERIFY(writeNpy<double>(QStringLiteral("%1/per_atom_substrate_target_T2.npy").arg(dir),
                             {2, 5},
                             {decomp.T2[0], decomp.T2[1], decomp.T2[2], decomp.T2[3], decomp.T2[4],
                              decomp.T2[0], decomp.T2[1], decomp.T2[2], decomp.T2[3], decomp.T2[4]},
                             QByteArray("<f8")));
}

void writeAllAtomFixture(const QString& dir) {
    const Mat3 raw = targetTensor();
    const auto decomp = DecomposeLibrary(raw);
    QVERIFY(writeText(QStringLiteral("%1/all_atom_equivariant_targets.csv").arg(dir),
                      QStringLiteral("row_id,atom_index,h5_row,original_index,dft_present\n"
                                     "0,0,7,77,1\n")));
    QVERIFY(writeNpy<double>(QStringLiteral("%1/all_atom_equivariant_target_sigma_iso.npy").arg(dir),
                             {1}, {decomp.T0}, QByteArray("<f8")));
    QVERIFY(writeNpy<double>(QStringLiteral("%1/all_atom_equivariant_target_T2.npy").arg(dir),
                             {1, 5},
                             {decomp.T2[0], decomp.T2[1], decomp.T2[2], decomp.T2[3], decomp.T2[4]},
                             QByteArray("<f8")));
    std::vector<double> flatRaw;
    flatRaw.reserve(9);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) flatRaw.push_back(raw(i, j));
    QVERIFY(writeNpy<double>(QStringLiteral("%1/all_atom_equivariant_target_raw.npy").arg(dir),
                             {1, 9}, flatRaw, QByteArray("<f8")));
}

}  // namespace

class RediscoverSubstrateParityTests : public QObject {
    Q_OBJECT

private slots:
    void perAtomParityAcceptsProducerEquivalentTargets();
    void perAtomParityRejectsRowContractDrift();
    void allAtomParityChecksRawAndT2Targets();
};

void RediscoverSubstrateParityTests::perAtomParityAcceptsProducerEquivalentTargets() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    writePerAtomFixture(dir.path());

    const auto stats = AuditPerAtomSubstrateSidecars(dir.path(), 2, fixedLookup(targetTensor()));
    QVERIFY2(stats.ok(), qPrintable(stats.errors.join(QLatin1Char('\n'))));
    QCOMPARE(stats.rows_seen, std::size_t(2));
    QCOMPARE(stats.rows_checked, std::size_t(2));
    QCOMPARE(stats.target_T0_checked, std::size_t(2));
    QCOMPARE(stats.target_T1_checked, std::size_t(6));
    QCOMPARE(stats.target_T2_checked, std::size_t(10));
}

void RediscoverSubstrateParityTests::perAtomParityRejectsRowContractDrift() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    writePerAtomFixture(dir.path(), true);

    const auto stats = AuditPerAtomSubstrateSidecars(dir.path(), 2, fixedLookup(targetTensor()));
    QVERIFY(!stats.ok());
    QVERIFY(stats.errors.join(QLatin1Char('\n')).contains(QStringLiteral("row contract mismatch")));
}

void RediscoverSubstrateParityTests::allAtomParityChecksRawAndT2Targets() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    writeAllAtomFixture(dir.path());

    const auto stats = AuditAllAtomEquivariantSidecars(dir.path(), fixedLookup(targetTensor()));
    QVERIFY2(stats.ok(), qPrintable(stats.errors.join(QLatin1Char('\n'))));
    QCOMPARE(stats.rows_seen, std::size_t(1));
    QCOMPARE(stats.rows_checked, std::size_t(1));
    QCOMPARE(stats.target_T0_checked, std::size_t(1));
    QCOMPARE(stats.target_T2_checked, std::size_t(5));
    QCOMPARE(stats.target_raw_components_checked, std::size_t(9));
}

QTEST_GUILESS_MAIN(RediscoverSubstrateParityTests)

#include "rediscover_substrate_parity_tests.moc"
