#include "io/QtProteinLoader.h"
#include "model/QtConformationSnapshot.h"

#include <QtTest>
#include <cstring>

using namespace h5reader::io;

class PublicationTests final : public QObject {
    Q_OBJECT
private slots:
    void publishedTrajectoryMatchesAllOriginalFrames() {
        const QString originalPath = qEnvironmentVariable("H5READER_PUBLICATION_ORIGINAL");
        const QString portablePath = qEnvironmentVariable("H5READER_PUBLICATION_PORTABLE");
        if (originalPath.isEmpty() && portablePath.isEmpty())
            QSKIP("Set H5READER_PUBLICATION_ORIGINAL and H5READER_PUBLICATION_PORTABLE to compare a published trajectory.");
        QVERIFY(!originalPath.isEmpty());
        QVERIFY(!portablePath.isEmpty());
        const auto original = QtProteinLoader::LoadRunPath(originalPath);
        const auto portable = QtProteinLoader::LoadRunPath(portablePath);
        QVERIFY2(original.ok, qPrintable(original.error));
        QVERIFY2(portable.ok, qPrintable(portable.error));
        QCOMPARE(portable.protein->atomCount(), original.protein->atomCount());
        QCOMPARE(portable.conformation->frameCount(), original.conformation->frameCount());
        for (std::size_t frame = 0; frame < original.conformation->frameCount(); ++frame) {
            QCOMPARE(portable.conformation->originalFrameIndex(frame), original.conformation->originalFrameIndex(frame));
            QCOMPARE(portable.conformation->timePicoseconds(frame), original.conformation->timePicoseconds(frame));
            for (std::size_t atom = 0; atom < original.protein->atomCount(); ++atom) {
                const auto expected = original.conformation->atomPosition(frame, atom);
                const auto actual = portable.conformation->atomPosition(frame, atom);
                QVERIFY((actual.array() == expected.array()).all());
            }
            original.conformation->requestSnapshot(frame);
            portable.conformation->requestSnapshot(frame);
            const auto expected = original.conformation->snapshot(frame);
            const auto actual = portable.conformation->snapshot(frame);
            QCOMPARE(bool(actual), bool(expected));
            if (!actual)
                continue;
            for (const auto& spec : kFieldCatalog) {
                const auto& a = actual->column(spec.kind);
                const auto& b = expected->column(spec.kind);
                const QString where =
                    QStringLiteral("frame %1, field %2").arg(frame).arg(QString::fromUtf8(spec.stem.data(), spec.stem.size()));
                QCOMPARE(a.present, b.present);
                QCOMPARE(a.rows, b.rows);
                QCOMPARE(a.cols, b.cols);
                QCOMPARE(a.data.size(), b.data.size());
                QVERIFY2(a.data.empty() || std::memcmp(a.data.data(), b.data.data(), a.data.size() * sizeof(double)) == 0,
                         qPrintable(where));
            }
            qInfo() << "Exact snapshot comparison passed: frame" << frame;
        }
    }
};

QTEST_GUILESS_MAIN(PublicationTests)
#include "publication_tests.moc"
