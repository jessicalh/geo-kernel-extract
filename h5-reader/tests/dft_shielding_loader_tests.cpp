// h5reader_dft_shielding_loader_tests — direct tests for the static
// DFT job loader. Post-2026-05-31 SIMPLIFY: the loader takes an
// explicit meta.json path (from `.LGS`'s `dft.frames[].meta_json`),
// not a jobs-dir + frame index pair.

#include "io/DftShieldingLoader.h"

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

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QTemporaryDir>
#include <QtTest>

#include <cmath>
#include <cstddef>
#include <memory>

using h5reader::io::DftShieldingLoader;
using h5reader::model::QtProtein;

namespace {

std::unique_ptr<QtProtein> makeProtein(std::size_t atomCount) {
    auto protein = std::make_unique<QtProtein>();
    protein->atoms_.resize(atomCount);
    return protein;
}

bool writeFile(const QString& path, const QByteArray& bytes) {
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate))
        return false;
    return file.write(bytes) == bytes.size();
}

QString matrix(double diagonal) {
    return QStringLiteral(
               "    %1 0.0000 0.0000\n"
               "    0.0000 %1 0.0000\n"
               "    0.0000 0.0000 %1\n")
        .arg(diagonal, 0, 'f', 4);
}

QString atomBlock(int index, const QString& element, double dia, double para) {
    return QStringLiteral(
               "Nucleus %1%2\n"
               "Diamagnetic contribution to the shielding tensor\n"
               "%3"
               "Paramagnetic contribution to the shielding tensor\n"
               "%4"
               "Total shielding tensor\n"
               "%5")
        .arg(index)
        .arg(element, matrix(dia), matrix(para), matrix(dia + para));
}

QString orcaOut(const QString& atomBlocks) {
    return QStringLiteral(
               "ORCA unit-test output\n"
               "CHEMICAL SHIELDINGS\n"
               "%1")
        .arg(atomBlocks);
}

// Builds a single DFT job dir + meta.json + primary .out, returning
// the absolute meta.json path that LoadAndValidate consumes.
QString writeJob(const QString& root, std::size_t originalIndex, const QString& outText) {
    const QString jobId = QStringLiteral("orca_f%1_t0").arg(originalIndex, 6, 10, QLatin1Char('0'));
    QDir parent(root);
    if (!parent.mkpath(jobId)) return {};
    const QString jobDir = parent.absoluteFilePath(jobId);
    const QString metaPath = QStringLiteral("%1/%2_meta.json").arg(jobDir, jobId);
    if (!writeFile(metaPath, QByteArrayLiteral("{\"files\":{\"out_primary\":\"primary.out\"}}\n")))
        return {};
    if (!writeFile(QStringLiteral("%1/primary.out").arg(jobDir), outText.toUtf8()))
        return {};
    return metaPath;
}

}  // namespace

class DftShieldingLoaderTests : public QObject {
    Q_OBJECT

private slots:
    void loadAndValidateHappyPath();
    void missingMetaReturnsNull();
    void parserHoleReturnsNull();
};

void DftShieldingLoaderTests::loadAndValidateHappyPath() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    const auto protein = makeProtein(2);
    const QString out = orcaOut(atomBlock(0, QStringLiteral("H"), 10.0, 1.0)
                                + atomBlock(1, QStringLiteral("C"), 20.0, 2.0));
    const QString metaPath = writeJob(dir.path(), 7, out);
    QVERIFY(!metaPath.isEmpty());

    const auto frame = DftShieldingLoader::LoadAndValidate(metaPath, protein.get());

    QVERIFY(frame != nullptr);
    QCOMPARE(frame->atoms.size(), std::size_t{2});
    QCOMPARE(frame->valid, true);
    QVERIFY(frame->atoms[0].element == h5reader::model::Element::H);
    QVERIFY(frame->atoms[1].element == h5reader::model::Element::C);
    QVERIFY(std::abs(frame->atoms[0].total.T0 - (frame->atoms[0].dia.T0 + frame->atoms[0].para.T0)) < 1e-9);
}

void DftShieldingLoaderTests::missingMetaReturnsNull() {
    // No file at the supplied path → null.
    const auto protein = makeProtein(1);
    const auto frame = DftShieldingLoader::LoadAndValidate(
        QStringLiteral("/tmp/this-meta-does-not-exist.json"), protein.get());
    QVERIFY(frame == nullptr);
}

void DftShieldingLoaderTests::parserHoleReturnsNull() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());
    const auto protein = makeProtein(3);
    const QString out = orcaOut(atomBlock(0, QStringLiteral("H"), 10.0, 1.0)
                                + atomBlock(2, QStringLiteral("C"), 20.0, 2.0));
    const QString metaPath = writeJob(dir.path(), 13, out);
    QVERIFY(!metaPath.isEmpty());

    const auto frame = DftShieldingLoader::LoadAndValidate(metaPath, protein.get());

    QVERIFY(frame == nullptr);
}

QTEST_GUILESS_MAIN(DftShieldingLoaderTests)

#include "dft_shielding_loader_tests.moc"
