// rediscover_buckingham_tests — pin the APBS E-field per_atom_feature carrier.
//
// Load-bearing checks: local projection is the local z component, |E| is
// rotation-invariant, T1 is emitted only as unverified audit payload, and T2
// keeps the library √6 fixture.

#include "io/QtNpyReader.h"
#include "rediscover/BuckinghamEfieldSink.h"
#include "rediscover/LocalFrameBasis.h"
#include "rediscover/SphericalBasis.h"

#include <QFile>
#include <QTemporaryDir>
#include <QtTest>

#include <cmath>

using h5reader::io::QtNpyReader;
using h5reader::model::Mat3;
using h5reader::model::SphericalTensor;
using h5reader::model::Vec3;
using h5reader::rediscover::BuckinghamEfieldRow;
using h5reader::rediscover::BuckinghamEfieldSink;
using h5reader::rediscover::DecomposeLibrary;
using h5reader::rediscover::FiniteVec3;
using h5reader::rediscover::LocalFrame;

namespace {

constexpr double kTol = 1e-9;

LocalFrame fixedFrame() {
    LocalFrame f;
    f.x = Vec3::UnitX();
    f.z = Vec3::UnitY();
    f.y = f.z.cross(f.x).normalized();
    f.is_valid = true;
    return f;
}

Mat3 sqrt6Tensor() {
    Mat3 m = Mat3::Zero();
    m(0, 0) = -1.0;
    m(1, 1) = -1.0;
    m(2, 2) = 2.0;
    return m;
}

Mat3 asymmetricTensor() {
    Mat3 m;
    m << 10.0, 1.0, -2.0,
         4.0, 11.0, 3.0,
        -5.0, 7.0, 12.0;
    return m;
}

BuckinghamEfieldRow rowWithFieldAndTarget(const Vec3& fieldLocal,
                                          const SphericalTensor& target,
                                          const Mat3& totalRaw) {
    BuckinghamEfieldRow row;
    row.atom_index = 10;
    row.residue_index = 1;
    row.residue_number = 2;
    row.amino_acid = 3;
    row.element = 2;
    row.atom_name = QStringLiteral("N");
    row.frame_variant = 4;
    row.frame_valid = true;
    row.h5_row = 5;
    row.original_index = 7;
    row.time_ps = 1.25;
    row.dft_present = true;
    row.apbs_efield_present = true;
    row.efield_local = fieldLocal;
    row.e_proj = fieldLocal.z();
    row.e_mag = fieldLocal.norm();
    row.efield_units = QStringLiteral("V/Angstrom");
    row.dft_total_raw = totalRaw;
    row.dft_total_decomp = target;
    return row;
}

}  // namespace

class RediscoverBuckinghamTests : public QObject {
    Q_OBJECT

private slots:
    void projectionIsLocalZAndMagnitudeInvariant();
    void targetBasisKeepsSqrt6Fixture();
    void sinkWritesScalarCsvAndAuditSidecars();
};

void RediscoverBuckinghamTests::projectionIsLocalZAndMagnitudeInvariant() {
    const LocalFrame frame = fixedFrame();
    const Vec3 eLab(1.0, 2.0, -3.0);
    const Vec3 eLocal = frame.ToLocal(eLab);

    QVERIFY(FiniteVec3(eLocal));
    QVERIFY(std::abs(eLocal.z() - eLab.dot(frame.z)) < kTol);
    QVERIFY(std::abs(eLocal.z() - 2.0) < kTol);
    QVERIFY(std::abs(eLocal.norm() - eLab.norm()) < kTol);
}

void RediscoverBuckinghamTests::targetBasisKeepsSqrt6Fixture() {
    const SphericalTensor st = DecomposeLibrary(sqrt6Tensor());
    QVERIFY(std::abs(st.T0) < kTol);
    QVERIFY(std::abs(st.T2[2] - std::sqrt(6.0)) < kTol);
}

void RediscoverBuckinghamTests::sinkWritesScalarCsvAndAuditSidecars() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());

    BuckinghamEfieldSink sink(dir.path(), QStringLiteral("buckingham_efield"));
    QVERIFY(sink.Ok());

    const Mat3 total = asymmetricTensor();
    const SphericalTensor target = DecomposeLibrary(total);
    sink.Write(rowWithFieldAndTarget(Vec3(1.0, 2.0, 3.0), target, total));
    sink.Write(rowWithFieldAndTarget(Vec3(-1.0, 0.5, 4.0), target, total));

    QVERIFY(sink.Commit());
    QCOMPARE(sink.rowsWritten(), static_cast<std::size_t>(2));

    const QtNpyReader::WidenedArray fieldNpy = QtNpyReader::ReadArrayWidened(
        dir.filePath(QStringLiteral("buckingham_efield_feature_field_local.npy")));
    QVERIFY2(fieldNpy.ok, qPrintable(fieldNpy.error));
    QCOMPARE(fieldNpy.rows, static_cast<std::size_t>(2));
    QCOMPARE(fieldNpy.cols, static_cast<std::size_t>(3));
    QVERIFY(std::abs(fieldNpy.data[2] - 3.0) < kTol);
    QVERIFY(std::abs(fieldNpy.data[5] - 4.0) < kTol);

    const QtNpyReader::WidenedArray t1Npy = QtNpyReader::ReadArrayWidened(
        dir.filePath(QStringLiteral("buckingham_efield_target_T1_unverified.npy")));
    QVERIFY2(t1Npy.ok, qPrintable(t1Npy.error));
    QCOMPARE(t1Npy.rows, static_cast<std::size_t>(2));
    QCOMPARE(t1Npy.cols, static_cast<std::size_t>(3));
    for (double v : t1Npy.data) QVERIFY(std::isfinite(v));

    const QtNpyReader::WidenedArray t2Npy = QtNpyReader::ReadArrayWidened(
        dir.filePath(QStringLiteral("buckingham_efield_target_T2.npy")));
    QVERIFY2(t2Npy.ok, qPrintable(t2Npy.error));
    QCOMPARE(t2Npy.rows, static_cast<std::size_t>(2));
    QCOMPARE(t2Npy.cols, static_cast<std::size_t>(5));
    for (double v : t2Npy.data) QVERIFY(std::isfinite(v));

    QFile csv(dir.filePath(QStringLiteral("buckingham_efield_aggregated.csv")));
    QVERIFY(csv.open(QIODevice::ReadOnly | QIODevice::Text));
    const QByteArray payload = csv.readAll();
    QVERIFY(payload.contains("E_proj,E_mag"));
    QVERIFY(payload.contains("unverified_emit_only"));
    const QList<QByteArray> lines = payload.split('\n');
    int nonEmpty = 0;
    for (const QByteArray& line : lines)
        if (!line.trimmed().isEmpty()) ++nonEmpty;
    QCOMPARE(nonEmpty, 3);  // header + two rows
}

QTEST_GUILESS_MAIN(RediscoverBuckinghamTests)

#include "rediscover_buckingham_tests.moc"
