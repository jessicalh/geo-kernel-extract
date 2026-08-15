#include "model/ScientificAlignment.h"

#include <Eigen/Geometry>

#include <QtTest>

#include <cmath>
#include <limits>
#include <vector>

namespace {

using h5reader::model::Mat3;
using h5reader::model::ScientificFitStatus;
using h5reader::model::Vec3;

std::vector<Vec3> transformPoints(const std::vector<Vec3>& points,
                                  const Mat3& rotation,
                                  const Vec3& translation) {
    std::vector<Vec3> transformed;
    transformed.reserve(points.size());
    for (const Vec3& point : points)
        transformed.push_back(rotation * point + translation);
    return transformed;
}

void appendFrame(h5reader::model::ScientificPositionTable* table,
                 const std::vector<Vec3>& points) {
    table->values.insert(table->values.end(), points.begin(), points.end());
    ++table->frameCount;
}

}  // namespace

class ScientificAlignmentTests final : public QObject {
    Q_OBJECT

private slots:
    void fullRankFitRecoversKnownTransform() {
        const std::vector<Vec3> current{
            {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0},
            {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0},
        };
        const Mat3 expectedRotation =
            Eigen::AngleAxisd(0.71, Vec3(1.0, 2.0, 3.0).normalized())
                .toRotationMatrix();
        const Vec3 expectedTranslation(4.0, -2.5, 0.75);
        const std::vector<Vec3> reference =
            transformPoints(current, expectedRotation, expectedTranslation);

        const auto fit = h5reader::model::ComputeScientificKabsch(
            current, reference);
        QCOMPARE(fit.status, ScientificFitStatus::Valid);
        QCOMPARE(fit.numericalRank, 3);
        QVERIFY((fit.rotation - expectedRotation).norm() < 1e-12);
        QVERIFY((fit.translation - expectedTranslation).norm() < 1e-12);
        QVERIFY(fit.rmsdAngstrom < 1e-12);
        QVERIFY(std::abs(fit.rotation.determinant() - 1.0) < 1e-12);
    }

    void planarNonCollinearFitIsAccepted() {
        const std::vector<Vec3> current{
            {-2.0, -1.0, 0.0}, {2.0, -1.0, 0.0},
            {2.0, 1.0, 0.0}, {-2.0, 1.0, 0.0},
        };
        const Mat3 expectedRotation =
            Eigen::AngleAxisd(-0.93, Vec3(0.5, 1.0, -0.25).normalized())
                .toRotationMatrix();
        const Vec3 expectedTranslation(-3.0, 1.5, 7.0);
        const std::vector<Vec3> reference =
            transformPoints(current, expectedRotation, expectedTranslation);

        const auto fit = h5reader::model::ComputeScientificKabsch(
            current, reference);
        QCOMPARE(fit.status, ScientificFitStatus::Valid);
        QCOMPARE(fit.numericalRank, 2);
        QVERIFY((fit.rotation - expectedRotation).norm() < 1e-12);
        QVERIFY((fit.translation - expectedTranslation).norm() < 1e-12);
        QVERIFY(fit.rmsdAngstrom < 1e-12);
    }

    void collinearFitIsRejected() {
        const std::vector<Vec3> current{
            {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
        };
        const std::vector<Vec3> reference{
            {1.0, 2.0, 3.0}, {1.0, 3.0, 3.0}, {1.0, 4.0, 3.0},
        };

        const auto fit = h5reader::model::ComputeScientificKabsch(
            current, reference);
        QCOMPARE(fit.status, ScientificFitStatus::RankBelowTwo);
        QCOMPARE(fit.numericalRank, 1);
    }

    void tooSmallFitIsRejected() {
        const std::vector<Vec3> points{{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}};
        const auto fit = h5reader::model::ComputeScientificKabsch(points, points);
        QCOMPARE(fit.status, ScientificFitStatus::TooFewPoints);
    }

    void nonfiniteFitIsRejected() {
        std::vector<Vec3> points{
            {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0},
        };
        points[1].x() = std::numeric_limits<double>::quiet_NaN();
        const auto fit = h5reader::model::ComputeScientificKabsch(points, points);
        QCOMPARE(fit.status, ScientificFitStatus::NonFiniteInput);
    }

    void planarTrajectoryBuildsConvergedMean() {
        const std::vector<Vec3> base{
            {-2.0, -1.0, 0.0}, {2.0, -1.0, 0.0},
            {2.0, 1.0, 0.0}, {-2.0, 1.0, 0.0}, {0.25, 0.4, 0.0},
        };
        h5reader::model::ScientificPositionTable table;
        table.atomCount = base.size();
        appendFrame(&table, base);
        appendFrame(&table, transformPoints(
            base, Eigen::AngleAxisd(0.4, Vec3::UnitX()).toRotationMatrix(),
            Vec3(5.0, -3.0, 2.0)));
        appendFrame(&table, transformPoints(
            base, Eigen::AngleAxisd(-1.1, Vec3(1.0, 1.0, 0.5).normalized())
                      .toRotationMatrix(),
            Vec3(-2.0, 8.0, -4.0)));

        const std::vector<std::size_t> atoms{0, 1, 2, 3, 4};
        const auto alignment =
            h5reader::model::BuildScientificAlignment(table, atoms);
        QVERIFY2(alignment.ok, qPrintable(alignment.error));
        QVERIFY(alignment.mean.converged);
        QVERIFY(alignment.mean.iterations >= 1);
        QCOMPARE(alignment.frameFits.size(), table.frameCount);
        for (const auto& fit : alignment.frameFits) {
            QCOMPARE(fit.status, ScientificFitStatus::Valid);
            QCOMPARE(fit.numericalRank, 2);
            QVERIFY(fit.rmsdAngstrom < 1e-11);
            QVERIFY(std::abs(fit.rotation.determinant() - 1.0) < 1e-12);
        }
    }

    void unconvergedMeanFailsLoudly() {
        h5reader::model::ScientificPositionTable table;
        table.atomCount = 4;
        appendFrame(&table, {{0.0, 0.0, 0.0}, {2.0, 0.0, 0.0},
                             {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}});
        appendFrame(&table, {{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0},
                             {0.0, 2.0, 0.0}, {0.0, 0.0, 4.0}});
        h5reader::model::ScientificAlignmentPolicy policy;
        policy.maxIterations = 1;
        policy.convergenceToleranceAngstrom = 1e-20;

        const auto alignment = h5reader::model::BuildScientificAlignment(
            table, {0, 1, 2, 3}, policy);
        QVERIFY(!alignment.ok);
        QVERIFY(!alignment.mean.converged);
        QVERIFY(alignment.error.contains(QStringLiteral("did not converge")));
    }
};

QTEST_MAIN(ScientificAlignmentTests)
#include "scientific_alignment_tests.moc"
