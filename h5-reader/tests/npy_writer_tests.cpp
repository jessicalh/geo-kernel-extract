#include "io/QtNpyReader.h"
#include "io/QtNpyWriter.h"

#include <QTemporaryDir>
#include <QtTest>

#include <limits>
#include <vector>

class NpyWriterTests final : public QObject {
    Q_OBJECT

private slots:
    void float32RoundTripPreservesShapeAndOrder() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString path = directory.filePath(QStringLiteral("matrix32.npy"));
        const std::vector<float> values{1.25f, -2.5f, 3.75f, 4.0f};

        const auto written = h5reader::io::QtNpyWriter::WriteFloat32(
            path, {2, 2}, values);
        QVERIFY2(written.ok, qPrintable(written.error));
        const auto read =
            h5reader::io::QtNpyReader::ReadNumericArrayWidened(path);
        QVERIFY2(read.ok, qPrintable(read.error));
        QCOMPARE(read.shape, (std::vector<std::size_t>{2, 2}));
        QCOMPARE(read.data, (std::vector<double>{1.25, -2.5, 3.75, 4.0}));
        QVERIFY(read.descr.find("<f4") != std::string::npos);
    }

    void float64RoundTripPreservesShapeAndOrder() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString path = directory.filePath(QStringLiteral("matrix.npy"));
        const std::vector<double> values{1.0, 2.0, 3.0, 4.0, 5.0, 6.0};

        const auto written = h5reader::io::QtNpyWriter::WriteFloat64(
            path, {2, 3}, values);
        QVERIFY2(written.ok, qPrintable(written.error));
        const auto read =
            h5reader::io::QtNpyReader::ReadNumericArrayWidened(path);
        QVERIFY2(read.ok, qPrintable(read.error));
        QCOMPARE(read.shape, (std::vector<std::size_t>{2, 3}));
        QCOMPARE(read.data, values);
        QVERIFY(read.descr.find("<f8") != std::string::npos);
    }

    void uint8RoundTrip() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString uint8Path =
            directory.filePath(QStringLiteral("status.npy"));
        const std::vector<std::uint8_t> status{0, 1, 6};

        const auto written = h5reader::io::QtNpyWriter::WriteUInt8(
            uint8Path, {status.size()}, status);
        QVERIFY2(written.ok, qPrintable(written.error));

        const auto readStatus =
            h5reader::io::QtNpyReader::ReadNumericArrayWidened(uint8Path);
        QVERIFY2(readStatus.ok, qPrintable(readStatus.error));
        QCOMPARE(readStatus.shape,
                 (std::vector<std::size_t>{status.size()}));
        QCOMPARE(readStatus.data,
                 (std::vector<double>{0.0, 1.0, 6.0}));
        QVERIFY(readStatus.descr.find("|u1") != std::string::npos);
    }

    void payloadCountMismatchIsRejected() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const auto written = h5reader::io::QtNpyWriter::WriteFloat64(
            directory.filePath(QStringLiteral("bad.npy")), {2, 3},
            std::vector<double>{1.0, 2.0});
        QVERIFY(!written.ok);
        QVERIFY(written.error.contains(QStringLiteral("payload")));
    }

    void overflowingShapeIsRejected() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const auto written = h5reader::io::QtNpyWriter::WriteFloat64(
            directory.filePath(QStringLiteral("overflow.npy")),
            {std::numeric_limits<std::size_t>::max(), 2}, {});
        QVERIFY(!written.ok);
        QVERIFY(written.error.contains(QStringLiteral("shape")));
    }
};

QTEST_MAIN(NpyWriterTests)
#include "npy_writer_tests.moc"
