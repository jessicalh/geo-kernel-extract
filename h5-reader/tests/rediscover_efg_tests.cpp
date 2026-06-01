// rediscover_efg_tests — pin the APBS EFG per_atom_feature carrier.
//
// The basis fixture is the load-bearing check: emitted EFG T2 uses the same
// library [xy, yz, zz, xz, xx-yy] isometric basis as the DFT target sidecar.

#include "io/QtNpyReader.h"
#include "rediscover/EfgFeatureSink.h"
#include "rediscover/SphericalBasis.h"

#include <QFile>
#include <QTemporaryDir>
#include <QtTest>

#include <cmath>
#include <vector>

using h5reader::io::QtNpyReader;
using h5reader::model::Mat3;
using h5reader::rediscover::DecomposeEfgLibraryT2;
using h5reader::rediscover::DecomposeLibrary;
using h5reader::rediscover::EfgFeatureRow;
using h5reader::rediscover::EfgFeatureSink;
using h5reader::rediscover::FiniteT2;
using h5reader::rediscover::ReconstructLibraryT2;
using h5reader::rediscover::T2Magnitude;

namespace {

constexpr double kTol = 1e-9;

void assertArrayNear(const std::array<double, 5>& a, const std::array<double, 5>& b) {
    for (std::size_t i = 0; i < 5; ++i)
        QVERIFY(std::abs(a[i] - b[i]) < kTol);
}

void assertMatNear(const Mat3& a, const Mat3& b, double tol = kTol) {
    QVERIFY((a - b).norm() < tol);
}

Mat3 fixedRotation() {
    const h5reader::model::Vec3 axis = h5reader::model::Vec3(1.0, -2.0, 0.5).normalized();
    const double angle = 0.9;
    Mat3 k;
    k << 0.0, -axis.z(), axis.y(),
         axis.z(), 0.0, -axis.x(),
        -axis.y(), axis.x(), 0.0;
    return Mat3::Identity() + std::sin(angle) * k
           + (1.0 - std::cos(angle)) * (k * k);
}

Mat3 knownEfgTensor() {
    Mat3 m;
    m << 1.2, 0.4, -0.7,
         0.4, -2.5, 1.1,
        -0.7, 1.1, 1.3;
    return m;  // symmetric, trace 0
}

EfgFeatureRow rowWithT2(int idx, const std::array<double, 5>& feature,
                        const std::array<double, 5>& target) {
    EfgFeatureRow row;
    row.atom_index = idx;
    row.residue_index = 0;
    row.residue_number = 1;
    row.amino_acid = 0;
    row.element = 1;
    row.atom_name = QStringLiteral("H");
    row.frame_variant = 1;
    row.frame_valid = true;
    row.frame_z = h5reader::model::Vec3::UnitZ();
    row.frame_x = h5reader::model::Vec3::UnitX();
    row.frame_y = h5reader::model::Vec3::UnitY();
    row.h5_row = idx;
    row.original_index = idx * 2;
    row.time_ps = static_cast<double>(idx);
    row.dft_present = true;
    row.apbs_efg_present = true;
    row.efg_feature_T2 = feature;
    row.dft_target_T2 = target;
    row.efg_feature_lab_T2 = feature;
    row.dft_target_lab_T2 = target;
    row.efg_units = QStringLiteral("V/Angstrom^2");
    return row;
}

}  // namespace

class RediscoverEfgTests : public QObject {
    Q_OBJECT

private slots:
    void basisFixtureMatchesTargetBasis();
    void rotationPreservesMagnitude();
    void sinkWritesFiniteRowAlignedSidecars();
};

void RediscoverEfgTests::basisFixtureMatchesTargetBasis() {
    const Mat3 efg = knownEfgTensor();
    const std::array<double, 5> emitted = DecomposeEfgLibraryT2(efg);
    const std::array<double, 5> targetBasis = DecomposeLibrary(efg).T2;

    assertArrayNear(emitted, targetBasis);
    assertMatNear(ReconstructLibraryT2(emitted), efg);
    QVERIFY(FiniteT2(emitted));
}

void RediscoverEfgTests::rotationPreservesMagnitude() {
    const Mat3 efg = knownEfgTensor();
    const Mat3 r = fixedRotation();
    const Mat3 rotated = r * efg * r.transpose();

    const std::array<double, 5> baseT2 = DecomposeEfgLibraryT2(efg);
    const std::array<double, 5> rotatedT2 = DecomposeEfgLibraryT2(rotated);

    QVERIFY(std::abs(T2Magnitude(rotatedT2) - T2Magnitude(baseT2)) < kTol);
    assertMatNear(ReconstructLibraryT2(rotatedT2), rotated);
}

void RediscoverEfgTests::sinkWritesFiniteRowAlignedSidecars() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());

    EfgFeatureSink sink(dir.path(), QStringLiteral("efg"));
    QVERIFY(sink.Ok());

    const std::array<double, 5> feature = DecomposeEfgLibraryT2(knownEfgTensor());
    const std::array<double, 5> target = DecomposeLibrary(Mat3::Identity() * 10.0 + knownEfgTensor()).T2;
    const int dftPresentAtoms = 3;
    for (int i = 0; i < dftPresentAtoms; ++i)
        sink.Write(rowWithT2(i, feature, target));

    QVERIFY(sink.Commit());
    QCOMPARE(sink.rowsWritten(), static_cast<std::size_t>(dftPresentAtoms));

    const QtNpyReader::WidenedArray featureNpy =
        QtNpyReader::ReadArrayWidened(dir.filePath(QStringLiteral("efg_feature_T2.npy")));
    QVERIFY2(featureNpy.ok, qPrintable(featureNpy.error));
    QCOMPARE(featureNpy.rows, static_cast<std::size_t>(dftPresentAtoms));
    QCOMPARE(featureNpy.cols, static_cast<std::size_t>(5));
    for (double v : featureNpy.data) QVERIFY(std::isfinite(v));

    const QtNpyReader::WidenedArray targetNpy =
        QtNpyReader::ReadArrayWidened(dir.filePath(QStringLiteral("efg_target_T2.npy")));
    QVERIFY2(targetNpy.ok, qPrintable(targetNpy.error));
    QCOMPARE(targetNpy.rows, static_cast<std::size_t>(dftPresentAtoms));
    QCOMPARE(targetNpy.cols, static_cast<std::size_t>(5));
    for (double v : targetNpy.data) QVERIFY(std::isfinite(v));

    const QtNpyReader::WidenedArray featureLabNpy =
        QtNpyReader::ReadArrayWidened(dir.filePath(QStringLiteral("efg_feature_lab_T2.npy")));
    QVERIFY2(featureLabNpy.ok, qPrintable(featureLabNpy.error));
    QCOMPARE(featureLabNpy.rows, static_cast<std::size_t>(dftPresentAtoms));
    QCOMPARE(featureLabNpy.cols, static_cast<std::size_t>(5));
    for (double v : featureLabNpy.data) QVERIFY(std::isfinite(v));

    const QtNpyReader::WidenedArray targetLabNpy =
        QtNpyReader::ReadArrayWidened(dir.filePath(QStringLiteral("efg_target_lab_T2.npy")));
    QVERIFY2(targetLabNpy.ok, qPrintable(targetLabNpy.error));
    QCOMPARE(targetLabNpy.rows, static_cast<std::size_t>(dftPresentAtoms));
    QCOMPARE(targetLabNpy.cols, static_cast<std::size_t>(5));
    for (double v : targetLabNpy.data) QVERIFY(std::isfinite(v));

    QFile csv(dir.filePath(QStringLiteral("efg_aggregated.csv")));
    QVERIFY(csv.open(QIODevice::ReadOnly | QIODevice::Text));
    const QList<QByteArray> lines = csv.readAll().split('\n');
    int nonEmpty = 0;
    for (const QByteArray& line : lines)
        if (!line.trimmed().isEmpty()) ++nonEmpty;
    QCOMPARE(nonEmpty, dftPresentAtoms + 1);  // header + DFT-present rows
}

QTEST_GUILESS_MAIN(RediscoverEfgTests)

#include "rediscover_efg_tests.moc"
