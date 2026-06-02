// rediscover_all_atom_equivariant_tests -- contract test for the corrected
// all-atom equivariant sink. The important invariant is frame discipline:
// source vectors and target tensors are emitted in the molecular/lab frame, not
// a per-atom local frame.

#include "io/QtNpyReader.h"
#include "rediscover/AllAtomEquivariantSink.h"

#include <QtTest>

#include <array>

using h5reader::model::Mat3;
using h5reader::model::Vec3;
using h5reader::rediscover::AllAtomEquivariantSink;
using h5reader::rediscover::AllAtomEquivariantSourceRecord;
using h5reader::rediscover::AllAtomEquivariantTargetRecord;

class RediscoverAllAtomEquivariantTests : public QObject {
    Q_OBJECT

private slots:
    void sinkWritesLabFrameTargetAndSourcePayloads();
};

void RediscoverAllAtomEquivariantTests::sinkWritesLabFrameTargetAndSourcePayloads() {
    QTemporaryDir dir;
    QVERIFY(dir.isValid());

    AllAtomEquivariantSink sink(dir.path(), QStringLiteral("all_atom_equivariant"), 4);
    QVERIFY(sink.Ok());

    AllAtomEquivariantTargetRecord target;
    target.atom_index = 7;
    target.residue_index = 2;
    target.residue_number = 14;
    target.amino_acid = 3;
    target.element = 1;
    target.atom_name = QStringLiteral("HA");
    target.h5_row = 11;
    target.original_index = 101;
    target.time_ps = 22.5;
    target.target.present = true;
    target.target.total_raw = (Mat3() << 1.0, 0.2, 0.3,
                                          0.2, 2.0, 0.4,
                                          0.3, 0.4, 3.0).finished();
    target.target.total_decomp.T0 = 2.0;
    target.target.total_decomp.T2 = {0.1, 0.2, 0.3, 0.4, 0.5};
    target.apbs_efield_present = true;
    target.apbs_efield_lab = Vec3(1.0, 2.0, 3.0);
    target.apbs_efg_present = true;
    target.apbs_efg_T2 = {1.0, 2.0, 3.0, 4.0, 5.0};
    target.aimnet2_charge_present = true;
    target.aimnet2_charge = -0.12;
    target.aimnet2_crg_present = true;
    target.aimnet2_crg_scalar = 9.0;
    target.aimnet2_crg_lab = Vec3(4.0, 5.0, 6.0);
    const std::array<float, 4> emb = {1.0f, 2.0f, 3.0f, 4.0f};
    target.aimnet2_embedding_present = true;
    target.aimnet2_embedding = emb.data();
    target.aimnet2_embedding_dims = emb.size();

    const int64_t rowId = sink.WriteTarget(target);
    QVERIFY(rowId == 0);

    AllAtomEquivariantSourceRecord source;
    source.target_atom_index = target.atom_index;
    source.target_residue_index = target.residue_index;
    source.h5_row = target.h5_row;
    source.original_index = target.original_index;
    source.time_ps = target.time_ps;
    source.mechanism = QStringLiteral("ring");
    source.source_kind = QStringLiteral("ring_center");
    source.category = QStringLiteral("PheBenzene");
    source.category_ord = 0;
    source.source_id = 3;
    source.disp = Vec3(1.0, 0.0, 0.0);
    source.r = 1.0;
    source.inv_r3 = 1.0;
    source.orientation_a = Vec3(0.0, 0.0, 1.0);
    source.source_value = -12.0;
    source.source_units = QStringLiteral("nA/T");
    source.ring_index = 3;
    source.ring_type_index = 0;
    sink.WriteSource(rowId, source);

    QVERIFY(sink.Commit());
    QCOMPARE(sink.targetRowsWritten(), std::size_t(1));
    QCOMPARE(sink.sourceRowsWritten(), std::size_t(1));

    QFile targets(QStringLiteral("%1/all_atom_equivariant_targets.csv").arg(dir.path()));
    QVERIFY(targets.open(QIODevice::ReadOnly | QIODevice::Text));
    const QString targetsCsv = QString::fromUtf8(targets.readAll());
    QVERIFY(targetsCsv.contains(QStringLiteral("molecular_lab_h5_orca_aligned")));
    QVERIFY(!targetsCsv.contains(QStringLiteral("target_local")));

    QFile sources(QStringLiteral("%1/all_atom_equivariant_sources.csv").arg(dir.path()));
    QVERIFY(sources.open(QIODevice::ReadOnly | QIODevice::Text));
    const QString sourcesCsv = QString::fromUtf8(sources.readAll());
    QVERIFY(sourcesCsv.contains(QStringLiteral("orientation_a_x")));
    QVERIFY(sourcesCsv.contains(QStringLiteral("ring_center")));
    QVERIFY(!sourcesCsv.contains(QStringLiteral("disp_local")));

    const auto t2 = h5reader::io::QtNpyReader::ReadArrayWidened(
        QStringLiteral("%1/all_atom_equivariant_target_T2.npy").arg(dir.path()));
    QVERIFY2(t2.ok, qPrintable(t2.error));
    QCOMPARE(t2.rows, std::size_t(1));
    QCOMPARE(t2.cols, std::size_t(5));
    QVERIFY(std::abs(t2.data[4] - 0.5) < 1e-12);

    const auto embedding = h5reader::io::QtNpyReader::ReadArrayWidened(
        QStringLiteral("%1/all_atom_equivariant_aimnet2_embedding.npy").arg(dir.path()));
    QVERIFY2(embedding.ok, qPrintable(embedding.error));
    QCOMPARE(embedding.rows, std::size_t(1));
    QCOMPARE(embedding.cols, std::size_t(4));
    QVERIFY(std::abs(embedding.data[3] - 4.0) < 1e-12);
}

QTEST_MAIN(RediscoverAllAtomEquivariantTests)
#include "rediscover_all_atom_equivariant_tests.moc"
