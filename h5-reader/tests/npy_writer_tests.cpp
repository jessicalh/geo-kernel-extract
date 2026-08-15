#include "io/QtNpyReader.h"
#include "io/QtNpyWriter.h"

#include <QTemporaryDir>
#include <QtTest>

#include <limits>
#include <vector>

class NpyWriterTests final : public QObject {
    Q_OBJECT

private slots:
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

    void unsignedTypesRoundTrip() {
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString uint64Path =
            directory.filePath(QStringLiteral("indices.npy"));
        const QString uint8Path =
            directory.filePath(QStringLiteral("status.npy"));
        const std::vector<std::uint64_t> indices{0, 15, 1500};
        const std::vector<std::uint8_t> status{0, 1, 6};

        auto written = h5reader::io::QtNpyWriter::WriteUInt64(
            uint64Path, {indices.size()}, indices);
        QVERIFY2(written.ok, qPrintable(written.error));
        written = h5reader::io::QtNpyWriter::WriteUInt8(
            uint8Path, {status.size()}, status);
        QVERIFY2(written.ok, qPrintable(written.error));

        const auto readIndices =
            h5reader::io::QtNpyReader::ReadNumericArrayWidened(uint64Path);
        const auto readStatus =
            h5reader::io::QtNpyReader::ReadNumericArrayWidened(uint8Path);
        QVERIFY2(readIndices.ok, qPrintable(readIndices.error));
        QVERIFY2(readStatus.ok, qPrintable(readStatus.error));
        QCOMPARE(readIndices.shape,
                 (std::vector<std::size_t>{indices.size()}));
        QCOMPARE(readStatus.shape,
                 (std::vector<std::size_t>{status.size()}));
        QCOMPARE(readIndices.data,
                 (std::vector<double>{0.0, 15.0, 1500.0}));
        QCOMPARE(readStatus.data,
                 (std::vector<double>{0.0, 1.0, 6.0}));
        QVERIFY(readIndices.descr.find("<u8") != std::string::npos);
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
