// h5reader_calcset_manifest_tests — QtTest binary for
// CalcsetManifest, the `.LGS` (Lowly Graduate Student) loader that
// supersedes ReaderInputManifest.
//
// Covers:
//   * directory resolution (zero / one / multiple .LGS files)
//   * schema_version strict equality
//   * malformed JSON
//   * kind dispatch (trajectory / single_pose / mutant_pair)
//   * required vs optional field handling
//   * declared paths must exist on disk (manifest doesn't lie)
//   * optional reference_pdb (present + absent)
//   * dft block: frames[] enumeration, frame_index validation
//   * mutant_pair: nested .LGS path validation
//
// Fixtures are built per-test using QTemporaryDir.

#include "io/CalcsetManifest.h"

#include <QDir>
#include <QFile>
#include <QObject>
#include <QTemporaryDir>
#include <QTextStream>
#include <QtTest>

using h5reader::io::CalcsetManifest;

namespace {

bool writeFile(const QString& path, const QString& content) {
    QFileInfo fi(path);
    if (!QDir().mkpath(fi.absolutePath())) return false;
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text)) return false;
    QTextStream out(&f);
    out << content;
    return true;
}

bool makeDir(const QString& path) {
    return QDir().mkpath(path);
}

// Minimum trajectory .LGS body builder. Caller appends an optional
// dft block / closing brace.
QString trajectoryLgs(const QString& extraField = QString()) {
    return QStringLiteral(R"({
        "schema_version": 1,
        "kind": "trajectory",
        "dataset_id": "test-traj",
        "protein_id": "X",
        "human_name": "Test trajectory",
        "trajectory": {
            "md_dir": "md",
            "topology_top": "topol.top",
            "extraction_dir": "extract",
            "trajectory_h5": "extract/trajectory.h5",
            "extraction_manifest": "extract/extraction_manifest.json",
            "frame_dt_ps": 10.0,
            "frame_index_basis": "trr_frame_index"%1
        }
    })").arg(extraField);
}

// Sets up the file/dir scaffolding a minimum trajectory .LGS expects
// to exist on disk.
bool scaffoldTrajectoryLayout(const QString& root) {
    if (!makeDir(root + "/md")) return false;
    if (!makeDir(root + "/extract")) return false;
    if (!writeFile(root + "/topol.top", QString())) return false;
    if (!writeFile(root + "/extract/trajectory.h5", QString())) return false;
    if (!writeFile(root + "/extract/extraction_manifest.json",
                   QStringLiteral("{\"protein_id\":\"X\"}"))) return false;
    return true;
}

}  // namespace

class CalcsetManifestTests : public QObject {
    Q_OBJECT

private slots:
    void testNoLgsInDirectory();
    void testMultipleLgsInDirectory();
    void testLgsFileDirectly();
    void testMissingFile();
    void testMalformedJson();
    void testWrongSchemaVersion();
    void testWrongKind();
    void testMissingDatasetId();

    void testTrajectoryHappyPath();
    void testTrajectoryMissingRequiredPath();
    void testTrajectoryOptionalReferencePdbPresent();
    void testTrajectoryOptionalReferencePdbAbsent();
    void testTrajectoryOptionalReferencePdbDeclaredButMissing();

    void testSinglePoseHappyPath();
    void testSinglePoseWrongPoseKind();

    void testMutantPairNestedLgs();

    void testDftBlockHappyPath();
    void testDftBlockNegativeFrameIndex();
    void testDftBlockMissingMetaJsonFile();

    void testMetadataPassthrough();
};

// --- discovery / file resolution ------------------------------------

void CalcsetManifestTests::testNoLgsInDirectory() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("no .LGS")));
}

void CalcsetManifestTests::testMultipleLgsInDirectory() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(writeFile(tmp.path() + "/a.LGS", QStringLiteral("{}")));
    QVERIFY(writeFile(tmp.path() + "/b.LGS", QStringLiteral("{}")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("multiple .LGS")));
}

void CalcsetManifestTests::testLgsFileDirectly() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    const QString lgs = tmp.path() + "/test.LGS";
    QVERIFY(writeFile(lgs, trajectoryLgs()));
    QString err;
    auto m = CalcsetManifest::Load(lgs, &err);
    QVERIFY2(m.has_value(), qPrintable(err));
    QCOMPARE(m->kind, CalcsetManifest::Kind::Trajectory);
}

void CalcsetManifestTests::testMissingFile() {
    QString err;
    auto m = CalcsetManifest::Load(QStringLiteral("/nonexistent/path"), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("does not exist")));
}

void CalcsetManifestTests::testMalformedJson() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(writeFile(tmp.path() + "/x.LGS", QStringLiteral("this is not json")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("parse error")));
}

void CalcsetManifestTests::testWrongSchemaVersion() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(writeFile(tmp.path() + "/x.LGS",
        QStringLiteral("{\"schema_version\":42,\"kind\":\"trajectory\","
                       "\"dataset_id\":\"x\",\"protein_id\":\"x\",\"human_name\":\"x\"}")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("schema_version=42 unsupported")));
}

void CalcsetManifestTests::testWrongKind() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(writeFile(tmp.path() + "/x.LGS",
        QStringLiteral("{\"schema_version\":1,\"kind\":\"hologram\","
                       "\"dataset_id\":\"x\",\"protein_id\":\"x\",\"human_name\":\"x\"}")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("kind=hologram unknown")));
}

void CalcsetManifestTests::testMissingDatasetId() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(writeFile(tmp.path() + "/x.LGS",
        QStringLiteral("{\"schema_version\":1,\"kind\":\"trajectory\","
                       "\"protein_id\":\"x\",\"human_name\":\"x\"}")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("dataset_id")));
}

// --- trajectory -----------------------------------------------------

void CalcsetManifestTests::testTrajectoryHappyPath() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    QVERIFY(writeFile(tmp.path() + "/test.LGS", trajectoryLgs()));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY2(m.has_value(), qPrintable(err));
    QCOMPARE(m->schema_version, 1);
    QCOMPARE(m->kind, CalcsetManifest::Kind::Trajectory);
    QCOMPARE(m->dataset_id, QStringLiteral("test-traj"));
    QCOMPARE(m->protein_id, QStringLiteral("X"));
    QCOMPARE(m->human_name, QStringLiteral("Test trajectory"));
    QVERIFY(m->trajectory.has_value());
    QVERIFY(m->trajectory->trajectory_h5_abspath.endsWith(
        QStringLiteral("/extract/trajectory.h5")));
    QVERIFY(m->trajectory->extraction_dir_abspath.endsWith(QStringLiteral("/extract")));
    QVERIFY(m->trajectory->extraction_manifest_abspath.endsWith(
        QStringLiteral("/extract/extraction_manifest.json")));
    QVERIFY(m->trajectory->md_dir_abspath.endsWith(QStringLiteral("/md")));
    QVERIFY(m->trajectory->topology_top_abspath.endsWith(QStringLiteral("/topol.top")));
    QCOMPARE(m->trajectory->frame_dt_ps, 10.0);
    QCOMPARE(m->trajectory->frame_index_basis, QStringLiteral("trr_frame_index"));
    QVERIFY(m->trajectory->reference_pdb_abspath.isEmpty());
    QVERIFY(!m->dft.has_value());
    QVERIFY(!m->single_pose.has_value());
    QVERIFY(!m->mutant_pair.has_value());
}

void CalcsetManifestTests::testTrajectoryMissingRequiredPath() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    // Delete the trajectory.h5 to make the manifest lie.
    QFile::remove(tmp.path() + "/extract/trajectory.h5");
    QVERIFY(writeFile(tmp.path() + "/test.LGS", trajectoryLgs()));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("trajectory_h5")));
    QVERIFY(err.contains(QStringLiteral("does not exist")));
}

void CalcsetManifestTests::testTrajectoryOptionalReferencePdbPresent() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    QVERIFY(writeFile(tmp.path() + "/extract/reference.pdb", QStringLiteral("REMARK\n")));
    QVERIFY(writeFile(tmp.path() + "/test.LGS",
        trajectoryLgs(QStringLiteral(",\n        \"reference_pdb\": \"extract/reference.pdb\""))));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY2(m.has_value(), qPrintable(err));
    QVERIFY(m->trajectory->reference_pdb_abspath.endsWith(
        QStringLiteral("/extract/reference.pdb")));
}

void CalcsetManifestTests::testTrajectoryOptionalReferencePdbAbsent() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    QVERIFY(writeFile(tmp.path() + "/test.LGS", trajectoryLgs()));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY2(m.has_value(), qPrintable(err));
    QVERIFY(m->trajectory->reference_pdb_abspath.isEmpty());
}

void CalcsetManifestTests::testTrajectoryOptionalReferencePdbDeclaredButMissing() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    // Declare reference_pdb but don't create it.
    QVERIFY(writeFile(tmp.path() + "/test.LGS",
        trajectoryLgs(QStringLiteral(",\n        \"reference_pdb\": \"extract/missing.pdb\""))));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("reference_pdb")));
    QVERIFY(err.contains(QStringLiteral("does not exist")));
}

// --- single pose ----------------------------------------------------

void CalcsetManifestTests::testSinglePoseHappyPath() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/pose"));
    QVERIFY(writeFile(tmp.path() + "/pose/extraction_manifest.json",
                      QStringLiteral("{\"protein_id\":\"X\"}")));
    QVERIFY(writeFile(tmp.path() + "/test.LGS",
        QStringLiteral(R"({
            "schema_version": 1,
            "kind": "single_pose",
            "dataset_id": "ds",
            "protein_id": "X",
            "human_name": "n",
            "single_pose": {
                "pose_kind": "orca",
                "pose_dir": "pose",
                "extraction_manifest": "pose/extraction_manifest.json"
            }
        })")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY2(m.has_value(), qPrintable(err));
    QCOMPARE(m->kind, CalcsetManifest::Kind::SinglePose);
    QVERIFY(m->single_pose.has_value());
    QCOMPARE(m->single_pose->pose_kind, CalcsetManifest::PoseKind::Orca);
    QVERIFY(m->single_pose->pose_dir_abspath.endsWith(QStringLiteral("/pose")));
}

void CalcsetManifestTests::testSinglePoseWrongPoseKind() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/pose"));
    QVERIFY(writeFile(tmp.path() + "/pose/extraction_manifest.json", QStringLiteral("{}")));
    QVERIFY(writeFile(tmp.path() + "/test.LGS",
        QStringLiteral(R"({
            "schema_version": 1,
            "kind": "single_pose",
            "dataset_id": "ds",
            "protein_id": "X",
            "human_name": "n",
            "single_pose": {
                "pose_kind": "xyz",
                "pose_dir": "pose",
                "extraction_manifest": "pose/extraction_manifest.json"
            }
        })")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("pose_kind=xyz unknown")));
}

// --- mutant pair ----------------------------------------------------

void CalcsetManifestTests::testMutantPairNestedLgs() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/wt"));
    QVERIFY(makeDir(tmp.path() + "/ala"));
    QVERIFY(writeFile(tmp.path() + "/wt/wt.LGS", QStringLiteral("{}")));
    QVERIFY(writeFile(tmp.path() + "/ala/ala.LGS", QStringLiteral("{}")));
    QVERIFY(writeFile(tmp.path() + "/test.LGS",
        QStringLiteral(R"({
            "schema_version": 1,
            "kind": "mutant_pair",
            "dataset_id": "ds",
            "protein_id": "X",
            "human_name": "n",
            "mutant_pair": {
                "wt_lgs": "wt/wt.LGS",
                "ala_lgs": "ala/ala.LGS"
            }
        })")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY2(m.has_value(), qPrintable(err));
    QCOMPARE(m->kind, CalcsetManifest::Kind::MutantPair);
    QVERIFY(m->mutant_pair.has_value());
    QVERIFY(m->mutant_pair->wt_lgs_abspath.endsWith(QStringLiteral("/wt/wt.LGS")));
    QVERIFY(m->mutant_pair->ala_lgs_abspath.endsWith(QStringLiteral("/ala/ala.LGS")));
}

// --- dft block ------------------------------------------------------

void CalcsetManifestTests::testDftBlockHappyPath() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    QVERIFY(makeDir(tmp.path() + "/dft/jobs/j0"));
    QVERIFY(makeDir(tmp.path() + "/dft/jobs/j2"));
    QVERIFY(writeFile(tmp.path() + "/dft/jobs/j0/m.json",
                      QStringLiteral("{\"frame_index\":0,\"frame_ps\":0.0}")));
    QVERIFY(writeFile(tmp.path() + "/dft/jobs/j2/m.json",
                      QStringLiteral("{\"frame_index\":2,\"frame_ps\":20.0}")));
    QVERIFY(writeFile(tmp.path() + "/test.LGS",
        QStringLiteral(R"({
            "schema_version": 1,
            "kind": "trajectory",
            "dataset_id": "ds",
            "protein_id": "X",
            "human_name": "n",
            "trajectory": {
                "md_dir": "md",
                "topology_top": "topol.top",
                "extraction_dir": "extract",
                "trajectory_h5": "extract/trajectory.h5",
                "extraction_manifest": "extract/extraction_manifest.json",
                "frame_dt_ps": 10.0,
                "frame_index_basis": "trr_frame_index"
            },
            "dft": {
                "method": "r2SCAN",
                "campaign_target_frames": 100,
                "frame_stride": {"first": 0, "last": 1000, "step": 2},
                "frames": [
                    {"frame_index": 0, "meta_json": "dft/jobs/j0/m.json"},
                    {"frame_index": 2, "meta_json": "dft/jobs/j2/m.json"}
                ]
            }
        })")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY2(m.has_value(), qPrintable(err));
    QVERIFY(m->dft.has_value());
    QCOMPARE(m->dft->method, QStringLiteral("r2SCAN"));
    QCOMPARE(m->dft->campaign_target_frames, 100);
    QCOMPARE(m->dft->frame_stride.first, 0);
    QCOMPARE(m->dft->frame_stride.last, 1000);
    QCOMPARE(m->dft->frame_stride.step, 2);
    QCOMPARE(m->dft->frames.size(), std::size_t{2});
    QCOMPARE(m->dft->frames[0].frame_index, 0);
    QCOMPARE(m->dft->frames[1].frame_index, 2);
    QVERIFY(m->dft->frames[0].meta_json_abspath.endsWith(QStringLiteral("/j0/m.json")));
    // Lazy meta load.
    QString metaErr;
    QVERIFY2(m->dft->frames[1].LoadMeta(&metaErr), qPrintable(metaErr));
    QCOMPARE(m->dft->frames[1].framePs(), 20.0);
}

void CalcsetManifestTests::testDftBlockNegativeFrameIndex() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    QVERIFY(makeDir(tmp.path() + "/dft/jobs/j"));
    QVERIFY(writeFile(tmp.path() + "/dft/jobs/j/m.json", QStringLiteral("{}")));
    QVERIFY(writeFile(tmp.path() + "/test.LGS",
        QStringLiteral(R"({
            "schema_version": 1,
            "kind": "trajectory",
            "dataset_id": "ds",
            "protein_id": "X",
            "human_name": "n",
            "trajectory": {
                "md_dir": "md",
                "topology_top": "topol.top",
                "extraction_dir": "extract",
                "trajectory_h5": "extract/trajectory.h5",
                "extraction_manifest": "extract/extraction_manifest.json",
                "frame_dt_ps": 10.0,
                "frame_index_basis": "trr_frame_index"
            },
            "dft": {
                "method": "r2SCAN",
                "campaign_target_frames": 0,
                "frame_stride": {"first": 0, "last": 0, "step": 1},
                "frames": [
                    {"frame_index": -1, "meta_json": "dft/jobs/j/m.json"}
                ]
            }
        })")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("frame_index=-1")));
}

void CalcsetManifestTests::testDftBlockMissingMetaJsonFile() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    // Declare meta_json but DON'T create it.
    QVERIFY(writeFile(tmp.path() + "/test.LGS",
        QStringLiteral(R"({
            "schema_version": 1,
            "kind": "trajectory",
            "dataset_id": "ds",
            "protein_id": "X",
            "human_name": "n",
            "trajectory": {
                "md_dir": "md",
                "topology_top": "topol.top",
                "extraction_dir": "extract",
                "trajectory_h5": "extract/trajectory.h5",
                "extraction_manifest": "extract/extraction_manifest.json",
                "frame_dt_ps": 10.0,
                "frame_index_basis": "trr_frame_index"
            },
            "dft": {
                "method": "r2SCAN",
                "campaign_target_frames": 0,
                "frame_stride": {"first": 0, "last": 0, "step": 1},
                "frames": [
                    {"frame_index": 0, "meta_json": "dft/jobs/missing/m.json"}
                ]
            }
        })")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY(!m.has_value());
    QVERIFY(err.contains(QStringLiteral("meta_json")));
    QVERIFY(err.contains(QStringLiteral("does not exist")));
}

// --- metadata -------------------------------------------------------

void CalcsetManifestTests::testMetadataPassthrough() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(scaffoldTrajectoryLayout(tmp.path()));
    QVERIFY(writeFile(tmp.path() + "/test.LGS",
        QStringLiteral(R"({
            "schema_version": 1,
            "kind": "trajectory",
            "dataset_id": "ds",
            "protein_id": "X",
            "human_name": "n",
            "trajectory": {
                "md_dir": "md",
                "topology_top": "topol.top",
                "extraction_dir": "extract",
                "trajectory_h5": "extract/trajectory.h5",
                "extraction_manifest": "extract/extraction_manifest.json",
                "frame_dt_ps": 10.0,
                "frame_index_basis": "trr_frame_index"
            },
            "metadata": {
                "generated_at_utc": "2026-05-31T22:00:00Z",
                "lgs_writer": "lgs-tools 0.1.0",
                "producer_extractor_version": "0.2.0"
            }
        })")));
    QString err;
    auto m = CalcsetManifest::Load(tmp.path(), &err);
    QVERIFY2(m.has_value(), qPrintable(err));
    QCOMPARE(m->generated_at_utc, QStringLiteral("2026-05-31T22:00:00Z"));
    QCOMPARE(m->lgs_writer, QStringLiteral("lgs-tools 0.1.0"));
    QCOMPARE(m->producer_extractor_version, QStringLiteral("0.2.0"));
}

QTEST_GUILESS_MAIN(CalcsetManifestTests)
#include "calcset_manifest_tests.moc"
