#include "../src/rediscover/ReaderOutputCatalog.h"

#include "../src/rediscover/Aimnet2FeatureSink.h"
#include "../src/rediscover/AllAtomEquivariantSink.h"
#include "../src/rediscover/BroadBackboneSink.h"
#include "../src/rediscover/BuckinghamEfieldSink.h"
#include "../src/rediscover/EfgFeatureSink.h"
#include "../src/rediscover/RecordSink.h"

#include <QSet>
#include <QString>
#include <QStringList>
#include <QTemporaryDir>
#include <QtTest/QtTest>

#include <string_view>

using namespace h5reader::rediscover;

namespace {

QString qstr(std::string_view value) {
    return QString::fromLatin1(value.data(), static_cast<qsizetype>(value.size()));
}

QString stripNpy(QString fileName) {
    if (fileName.endsWith(QStringLiteral(".npy"))) fileName.chop(4);
    return fileName;
}

QSet<QString> stemsFromFiles(const QStringList& files) {
    QSet<QString> stems;
    for (const QString& file : files) stems.insert(stripNpy(file));
    return stems;
}

QSet<QString> declaredFor(ReaderOutputProducer producer) {
    QSet<QString> stems;
    for (const ReaderOutputSpec& spec : kReaderOutputCatalog) {
        if (spec.produced_by == producer) stems.insert(qstr(spec.stem));
    }
    return stems;
}

QStringList sortedValues(const QSet<QString>& values) {
    QStringList out;
    for (const QString& value : values) out << value;
    out.sort();
    return out;
}

QString diffMessage(const QSet<QString>& declared, const QSet<QString>& emitted) {
    const QSet<QString> stale = declared - emitted;
    const QSet<QString> missing = emitted - declared;
    QStringList parts;
    if (!stale.isEmpty()) {
        const QStringList values = sortedValues(stale);
        parts << QStringLiteral("declared_not_emitted=[%1]").arg(values.join(QStringLiteral(", ")));
    }
    if (!missing.isEmpty()) {
        const QStringList values = sortedValues(missing);
        parts << QStringLiteral("emitted_not_declared=[%1]").arg(values.join(QStringLiteral(", ")));
    }
    return parts.join(QStringLiteral("; "));
}

void expectCatalogMatches(ReaderOutputProducer producer, const QStringList& emittedFiles) {
    const QSet<QString> declared = declaredFor(producer);
    const QSet<QString> emitted = stemsFromFiles(emittedFiles);
    if (declared != emitted) {
        QFAIL(qPrintable(diffMessage(declared, emitted)));
    }
}

FeatureSchema recordSchema(const QString& caseName, bool includeBareKernel) {
    FeatureSchema schema;
    schema.caseName = caseName;
    schema.includeBareKernel = includeBareKernel;
    return schema;
}

QStringList perAtomWriterSidecars() {
    QStringList files;
    for (const std::string_view file : kPerAtomSubstrateAlwaysSidecars) files << qstr(file);
    files << qstr(kPerAtomSubstrateEmbeddingSidecar);
    return files;
}

}  // namespace

class ReaderOutputCatalogTests : public QObject {
    Q_OBJECT

private slots:
    void entriesAreCompleteAndUnique() {
        QCOMPARE(kReaderOutputCatalog.size(), std::size_t{73});

        QSet<QString> seen;
        for (const ReaderOutputSpec& spec : kReaderOutputCatalog) {
            QVERIFY2(!spec.stem.empty(), "catalog stem must not be empty");
            QVERIFY2(!spec.purpose.empty(), qPrintable(qstr(spec.stem) + QStringLiteral(" has empty purpose")));
            QVERIFY2(!spec.mechanism.empty(),
                     qPrintable(qstr(spec.stem) + QStringLiteral(" has empty mechanism")));
            QVERIFY2(!spec.units.empty(), qPrintable(qstr(spec.stem) + QStringLiteral(" has empty units")));
            QVERIFY2(!qstr(spec.stem).endsWith(QStringLiteral(".npy")),
                     qPrintable(qstr(spec.stem) + QStringLiteral(" should be a stem, not a filename")));
            QVERIFY2(ReaderOutputProducerName(spec.produced_by) != std::string_view("unknown"),
                     qPrintable(qstr(spec.stem) + QStringLiteral(" has unknown producer")));
            QVERIFY2(ReaderOutputTrailingElements(spec) > 0,
                     qPrintable(qstr(spec.stem) + QStringLiteral(" has zero trailing elements")));

            const QString stem = qstr(spec.stem);
            QVERIFY2(!seen.contains(stem), qPrintable(stem + QStringLiteral(" is duplicated")));
            seen.insert(stem);
            QCOMPARE(FindReaderOutputByStem(spec.stem), &spec);
        }
    }

    void declarationsMatchWriterSidecars() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());

        RecordSink ring(dir.path(), recordSchema(QStringLiteral("ring_current"), true));
        expectCatalogMatches(ReaderOutputProducer::RingCurrent, ring.sidecarFiles());

        RecordSink mc(dir.path(), recordSchema(QStringLiteral("mcconnell"), true));
        expectCatalogMatches(ReaderOutputProducer::McConnell, mc.sidecarFiles());

        RecordSink charge(dir.path(), recordSchema(QStringLiteral("charge_dipole"), false));
        expectCatalogMatches(ReaderOutputProducer::ChargeDipole, charge.sidecarFiles());

        BroadBackboneSink broad(dir.path(), QStringLiteral("broad_backbone"));
        expectCatalogMatches(ReaderOutputProducer::BroadBackbone, broad.sidecarFiles());

        EfgFeatureSink efg(dir.path(), QStringLiteral("efg"));
        expectCatalogMatches(ReaderOutputProducer::Efg, efg.sidecarFiles());

        BuckinghamEfieldSink buckingham(dir.path(), QStringLiteral("buckingham_efield"));
        expectCatalogMatches(ReaderOutputProducer::BuckinghamEfield, buckingham.sidecarFiles());

        AllAtomEquivariantSink allAtom(dir.path(), QStringLiteral("all_atom_equivariant"));
        expectCatalogMatches(ReaderOutputProducer::AllAtomEquivariant, allAtom.sidecarFiles());

        Aimnet2FeatureSink aimnet2(dir.path(), QStringLiteral("aimnet2_features"));
        expectCatalogMatches(ReaderOutputProducer::Aimnet2Features, aimnet2.sidecarFiles());

        expectCatalogMatches(ReaderOutputProducer::PerAtomSubstrate, perAtomWriterSidecars());
    }

    void writerConstantsDriveWidths() {
        const ReaderOutputSpec* classical =
            FindReaderOutputByStem("per_atom_substrate_features_classical");
        QVERIFY(classical);
        QCOMPARE(classical->extent0, kPerAtomClassicalCols);

        const ReaderOutputSpec* conditioning =
            FindReaderOutputByStem("per_atom_substrate_features_conditioning");
        QVERIFY(conditioning);
        QCOMPARE(conditioning->extent0, kPerAtomConditioningCols);

        const ReaderOutputSpec* methodPaths =
            FindReaderOutputByStem("per_atom_substrate_features_method_paths");
        QVERIFY(methodPaths);
        QCOMPARE(methodPaths->extent0, kPerAtomMethodPathCols);

        const ReaderOutputSpec* hbond =
            FindReaderOutputByStem("per_atom_substrate_features_hbond_conditioning");
        QVERIFY(hbond);
        QCOMPARE(hbond->extent0, kPerAtomHbondConditioningCols);

        const ReaderOutputSpec* embedding =
            FindReaderOutputByStem("per_atom_substrate_aimnet2_embedding");
        QVERIFY(embedding);
        QCOMPARE(embedding->extent0, kPerAtomEmbeddingDims);

        const ReaderOutputSpec* raw = FindReaderOutputByStem("all_atom_equivariant_target_raw");
        QVERIFY(raw);
        QCOMPARE(raw->rank, ReaderOutputRank::Rank3);
        QCOMPARE(raw->extent0, std::size_t{3});
        QCOMPARE(raw->extent1, std::size_t{3});
    }
};

QTEST_MAIN(ReaderOutputCatalogTests)
#include "reader_output_catalog_tests.moc"
