// h5reader_reader_input_manifest_tests — QtTest binary for
// ReaderInputManifest, the TOML-based input-directory descriptor.
//
// Covers:
//   * ExistsAt() — positive + negative
//   * missing file       -> ok=false with "not found" in error
//   * malformed TOML     -> ok=false with parse error in error
//   * wrong schema_version, wrong run_kind, missing required key
//   * trajectory + single_pose + mutant_pair happy paths
//   * required path missing on disk -> ok=false
//   * optional per_frame_npys_dir present + absent
//   * optional [dft] present (with missing dir → error; with valid dir → populated)
//   * primaryPoseDir / alternatePoseDir helpers
//
// Fixtures are built per-test using QTemporaryDir — each test owns its
// own dir tree so tests can run in any order and one failing won't leak
// into the next.

#include "io/ReaderInputManifest.h"

#include <QDir>
#include <QFile>
#include <QObject>
#include <QTemporaryDir>
#include <QTextStream>
#include <QtTest>

using h5reader::io::ReaderInputManifest;

namespace {

// Write a file with `content` at `path`. Creates parent dirs as needed.
// Returns false on any I/O failure.
bool writeFile(const QString& path, const QString& content) {
    QFileInfo fi(path);
    if (!QDir().mkpath(fi.absolutePath()))
        return false;
    QFile f(path);
    if (!f.open(QIODevice::WriteOnly | QIODevice::Text))
        return false;
    QTextStream out(&f);
    out << content;
    return true;
}

// Create an empty directory at `path` (including parents).
bool makeDir(const QString& path) {
    return QDir().mkpath(path);
}

}  // namespace

class ReaderInputManifestTests : public QObject {
    Q_OBJECT

private slots:
    void testExistsAt();
    void testMissingFile();
    void testMalformedToml();
    void testWrongSchemaVersion();
    void testWrongRunKind();
    void testMissingProteinId();
    void testTrajectoryHappyPath();
    void testTrajectoryMissingRequiredPath();
    void testTrajectoryOptionalPerFrameAbsent();
    void testTrajectoryOptionalPerFrameMissingOnDisk();
    void testSinglePoseHappyPath();
    void testMutantPairHappyPath();
    void testMutantPairHelpers();
    void testDftSectionPresent();
    void testDftSectionMissingDir();
    void testReferencePdbOptional();
    void testMetadataPassthrough();
};

void ReaderInputManifestTests::testExistsAt() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(!ReaderInputManifest::ExistsAt(tmp.path()));
    writeFile(tmp.path() + "/h5reader_manifest.toml", "schema_version = 1\n");
    QVERIFY(ReaderInputManifest::ExistsAt(tmp.path()));
}

void ReaderInputManifestTests::testMissingFile() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY(!m.ok);
    QVERIFY(m.error.contains(QStringLiteral("not found")));
}

void ReaderInputManifestTests::testMalformedToml() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    writeFile(tmp.path() + "/h5reader_manifest.toml",
              "schema_version = oops\n");   // not an int
    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY(!m.ok);
    // Either toml++ flagged it as parse error, or our reader said the
    // key is missing because the cast to int64_t failed. Either way ok=false.
    QVERIFY(!m.error.isEmpty());
}

void ReaderInputManifestTests::testWrongSchemaVersion() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    writeFile(tmp.path() + "/h5reader_manifest.toml",
              "schema_version = 42\n"
              "run_kind = \"trajectory\"\n"
              "protein_id = \"X\"\n");
    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY(!m.ok);
    QVERIFY(m.error.contains(QStringLiteral("schema_version=42 unsupported")));
}

void ReaderInputManifestTests::testWrongRunKind() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    writeFile(tmp.path() + "/h5reader_manifest.toml",
              "schema_version = 1\n"
              "run_kind = \"hologram\"\n"
              "protein_id = \"X\"\n");
    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY(!m.ok);
    QVERIFY(m.error.contains(QStringLiteral("run_kind=hologram unknown")));
}

void ReaderInputManifestTests::testMissingProteinId() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    writeFile(tmp.path() + "/h5reader_manifest.toml",
              "schema_version = 1\n"
              "run_kind = \"trajectory\"\n");
    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY(!m.ok);
    QVERIFY(m.error.contains(QStringLiteral("protein_id")));
}

void ReaderInputManifestTests::testTrajectoryHappyPath() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    // Build a layout with trajectory.h5 + sidecar dir + per_frame dir.
    QVERIFY(makeDir(tmp.path() + "/extract"));
    QVERIFY(makeDir(tmp.path() + "/extract/npys"));
    QVERIFY(writeFile(tmp.path() + "/extract/trajectory.h5", QString()));

    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"trajectory\"\n"
        "protein_id     = \"FAKE\"\n"
        "human_name     = \"Fake test trajectory\"\n"
        "\n"
        "[trajectory]\n"
        "trajectory_h5        = \"extract/trajectory.h5\"\n"
        "topology_sidecar_dir = \"extract\"\n"
        "per_frame_npys_dir   = \"extract/npys\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY2(m.ok, qPrintable(m.error));
    QCOMPARE(m.schemaVersion, 1);
    QCOMPARE(m.runKind, ReaderInputManifest::RunKind::Trajectory);
    QCOMPARE(m.proteinId, QStringLiteral("FAKE"));
    QCOMPARE(m.humanName, QStringLiteral("Fake test trajectory"));
    QVERIFY(m.trajectoryH5.endsWith(QStringLiteral("/extract/trajectory.h5")));
    QVERIFY(m.topologySidecarDir.endsWith(QStringLiteral("/extract")));
    QVERIFY(m.perFrameNpysDir.endsWith(QStringLiteral("/extract/npys")));
    QVERIFY(m.dftJobsDir.isEmpty());
    QVERIFY(m.referencePdb.isEmpty());
}

void ReaderInputManifestTests::testTrajectoryMissingRequiredPath() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/extract"));
    // Deliberately don't create trajectory.h5
    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"trajectory\"\n"
        "protein_id     = \"X\"\n"
        "[trajectory]\n"
        "trajectory_h5        = \"extract/trajectory.h5\"\n"
        "topology_sidecar_dir = \"extract\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY(!m.ok);
    QVERIFY(m.error.contains(QStringLiteral("trajectory_h5")));
    QVERIFY(m.error.contains(QStringLiteral("does not exist")));
}

void ReaderInputManifestTests::testTrajectoryOptionalPerFrameAbsent() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/extract"));
    QVERIFY(writeFile(tmp.path() + "/extract/trajectory.h5", QString()));

    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"trajectory\"\n"
        "protein_id     = \"X\"\n"
        "[trajectory]\n"
        "trajectory_h5        = \"extract/trajectory.h5\"\n"
        "topology_sidecar_dir = \"extract\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY2(m.ok, qPrintable(m.error));
    QVERIFY(m.perFrameNpysDir.isEmpty());
}

void ReaderInputManifestTests::testTrajectoryOptionalPerFrameMissingOnDisk() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/extract"));
    QVERIFY(writeFile(tmp.path() + "/extract/trajectory.h5", QString()));
    // Declare per_frame_npys_dir but DON'T create it on disk

    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"trajectory\"\n"
        "protein_id     = \"X\"\n"
        "[trajectory]\n"
        "trajectory_h5        = \"extract/trajectory.h5\"\n"
        "topology_sidecar_dir = \"extract\"\n"
        "per_frame_npys_dir   = \"extract/npys_does_not_exist\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY(!m.ok);
    QVERIFY(m.error.contains(QStringLiteral("per_frame_npys_dir")));
}

void ReaderInputManifestTests::testSinglePoseHappyPath() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/pose"));

    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"single_pose\"\n"
        "protein_id     = \"P\"\n"
        "[single_pose]\n"
        "pose_dir = \"pose\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY2(m.ok, qPrintable(m.error));
    QCOMPARE(m.runKind, ReaderInputManifest::RunKind::SinglePose);
    QVERIFY(m.poseDir.endsWith(QStringLiteral("/pose")));
    QVERIFY(m.trajectoryH5.isEmpty());
    QVERIFY(m.topologySidecarDir.isEmpty());
}

void ReaderInputManifestTests::testMutantPairHappyPath() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/wt"));
    QVERIFY(makeDir(tmp.path() + "/ala"));

    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"mutant_pair\"\n"
        "protein_id     = \"M\"\n"
        "[mutant_pair]\n"
        "wt_pose_dir  = \"wt\"\n"
        "ala_pose_dir = \"ala\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY2(m.ok, qPrintable(m.error));
    QCOMPARE(m.runKind, ReaderInputManifest::RunKind::MutantPair);
    QVERIFY(m.wtPoseDir.endsWith(QStringLiteral("/wt")));
    QVERIFY(m.alaPoseDir.endsWith(QStringLiteral("/ala")));
}

void ReaderInputManifestTests::testMutantPairHelpers() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/wt"));
    QVERIFY(makeDir(tmp.path() + "/ala"));
    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"mutant_pair\"\n"
        "protein_id     = \"M\"\n"
        "[mutant_pair]\n"
        "wt_pose_dir  = \"wt\"\n"
        "ala_pose_dir = \"ala\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY(m.ok);
    QCOMPARE(m.primaryPoseDir(), m.wtPoseDir);
    QCOMPARE(m.alternatePoseDir(), m.alaPoseDir);

    // Non-mutant: helpers return empty.
    ReaderInputManifest single;
    single.runKind = ReaderInputManifest::RunKind::SinglePose;
    QVERIFY(single.primaryPoseDir().isEmpty());
    QVERIFY(single.alternatePoseDir().isEmpty());
}

void ReaderInputManifestTests::testDftSectionPresent() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/extract"));
    QVERIFY(writeFile(tmp.path() + "/extract/trajectory.h5", QString()));
    QVERIFY(makeDir(tmp.path() + "/dft/jobs"));

    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"trajectory\"\n"
        "protein_id     = \"X\"\n"
        "[trajectory]\n"
        "trajectory_h5        = \"extract/trajectory.h5\"\n"
        "topology_sidecar_dir = \"extract\"\n"
        "[dft]\n"
        "jobs_dir = \"dft/jobs\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY2(m.ok, qPrintable(m.error));
    QVERIFY(m.dftJobsDir.endsWith(QStringLiteral("/dft/jobs")));
}

void ReaderInputManifestTests::testDftSectionMissingDir() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/extract"));
    QVERIFY(writeFile(tmp.path() + "/extract/trajectory.h5", QString()));
    // Declare dft/jobs but don't create it

    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"trajectory\"\n"
        "protein_id     = \"X\"\n"
        "[trajectory]\n"
        "trajectory_h5        = \"extract/trajectory.h5\"\n"
        "topology_sidecar_dir = \"extract\"\n"
        "[dft]\n"
        "jobs_dir = \"dft/jobs\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY(!m.ok);
    QVERIFY(m.error.contains(QStringLiteral("dft.jobs_dir")));
}

void ReaderInputManifestTests::testReferencePdbOptional() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/extract"));
    QVERIFY(writeFile(tmp.path() + "/extract/trajectory.h5", QString()));
    QVERIFY(writeFile(tmp.path() + "/extract/reference.pdb",
                      QStringLiteral("REMARK fake\nEND\n")));

    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"trajectory\"\n"
        "protein_id     = \"X\"\n"
        "reference_pdb  = \"extract/reference.pdb\"\n"
        "[trajectory]\n"
        "trajectory_h5        = \"extract/trajectory.h5\"\n"
        "topology_sidecar_dir = \"extract\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY2(m.ok, qPrintable(m.error));
    QVERIFY(m.referencePdb.endsWith(QStringLiteral("/extract/reference.pdb")));
}

void ReaderInputManifestTests::testMetadataPassthrough() {
    QTemporaryDir tmp; QVERIFY(tmp.isValid());
    QVERIFY(makeDir(tmp.path() + "/extract"));
    QVERIFY(writeFile(tmp.path() + "/extract/trajectory.h5", QString()));

    QVERIFY(writeFile(tmp.path() + "/h5reader_manifest.toml",
        "schema_version = 1\n"
        "run_kind       = \"trajectory\"\n"
        "protein_id     = \"X\"\n"
        "[trajectory]\n"
        "trajectory_h5        = \"extract/trajectory.h5\"\n"
        "topology_sidecar_dir = \"extract\"\n"
        "[metadata]\n"
        "extract_version = \"0.5.0\"\n"
        "fixture_kind    = \"calibration\"\n"));

    auto m = ReaderInputManifest::Load(tmp.path());
    QVERIFY2(m.ok, qPrintable(m.error));
    QCOMPARE(m.extractVersion, QStringLiteral("0.5.0"));
    QCOMPARE(m.fixtureKind, QStringLiteral("calibration"));
}

QTEST_GUILESS_MAIN(ReaderInputManifestTests)
#include "reader_input_manifest_tests.moc"
