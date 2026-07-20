#include "io/QtNpyReader.h"

#include <QFile>
#include <QTemporaryDir>
#include <QtTest>

#include <array>
#include <cstring>

namespace {

QByteArray makeNpy(const QByteArray& header, const std::array<float, 4>& values) {
    QByteArray bytes;
    bytes.append('\x93');
    bytes.append("NUMPY", 5);
    bytes.append('\x01');
    bytes.append('\x00');
    const auto size = static_cast<quint16>(header.size());
    bytes.append(static_cast<char>(size & 0xff));
    bytes.append(static_cast<char>((size >> 8) & 0xff));
    bytes.append(header);
    bytes.append(reinterpret_cast<const char*>(values.data()),
                 static_cast<qsizetype>(values.size() * sizeof(float)));
    return bytes;
}

QString writeNpy(QTemporaryDir& dir, const QByteArray& header) {
    const QString path = dir.filePath(QStringLiteral("array.npy"));
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly))
        return {};
    const std::array<float, 4> values{1.0F, 2.0F, 3.0F, 4.0F};
    file.write(makeNpy(header, values));
    file.close();
    return path;
}

}  // namespace

class NpyReaderTests final : public QObject {
    Q_OBJECT

private slots:
    void validRowMajorArrayLoads() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());
        const QString path = writeNpy(
            dir,
            QByteArrayLiteral("{'descr': '<f4', 'fortran_order': False, 'shape': (2, 2), }"));
        QVERIFY(!path.isEmpty());

        const auto result = h5reader::io::QtNpyReader::ReadArrayWidened(path);
        QVERIFY2(result.ok, qPrintable(result.error));
        QCOMPARE(result.rows, std::size_t{2});
        QCOMPARE(result.cols, std::size_t{2});
        QCOMPARE(result.data.size(), std::size_t{4});
        QCOMPARE(result.data[3], 4.0);
    }

    void malformedShapeTokenIsRejected() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());
        const QString path = writeNpy(
            dir,
            QByteArrayLiteral("{'descr': '<f4', 'fortran_order': False, 'shape': (2, bad, 2), }"));
        QVERIFY(!path.isEmpty());

        const auto result = h5reader::io::QtNpyReader::ReadArrayWidened(path);
        QVERIFY(!result.ok);
        QVERIFY(result.error.contains(QStringLiteral("shape")));
    }

    void doubledShapeCommaIsRejected() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());
        const QString path = writeNpy(
            dir,
            QByteArrayLiteral("{'descr': '<f4', 'fortran_order': False, 'shape': (2,, 2), }"));
        QVERIFY(!path.isEmpty());

        const auto result = h5reader::io::QtNpyReader::ReadArrayWidened(path);
        QVERIFY(!result.ok);
        QVERIFY(result.error.contains(QStringLiteral("shape")));
    }

    void missingFortranOrderIsRejected() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());
        const QString path = writeNpy(
            dir,
            QByteArrayLiteral("{'descr': '<f4', 'shape': (2, 2), }"));
        QVERIFY(!path.isEmpty());

        const auto result = h5reader::io::QtNpyReader::ReadArrayWidened(path);
        QVERIFY(!result.ok);
        QVERIFY(result.error.contains(QStringLiteral("fortran_order")));
    }

    void malformedFortranOrderIsRejected() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());
        const QString path = writeNpy(
            dir,
            QByteArrayLiteral("{'descr': '<f4', 'fortran_order': DefinitelyFalse, 'shape': (2, 2), }"));
        QVERIFY(!path.isEmpty());

        const auto result = h5reader::io::QtNpyReader::ReadArrayWidened(path);
        QVERIFY(!result.ok);
        QVERIFY(result.error.contains(QStringLiteral("fortran_order")));
    }

    void oversizedDtypeWidthIsRejectedWithoutThrowing() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());
        const QString path = writeNpy(
            dir,
            QByteArrayLiteral("{'descr': '<f999999999999', 'fortran_order': False, 'shape': (2, 2), }"));
        QVERIFY(!path.isEmpty());

        const auto result = h5reader::io::QtNpyReader::ReadArrayWidened(path);
        QVERIFY(!result.ok);
        QVERIFY(result.error.contains(QStringLiteral("dtype")));
    }

    void overflowingShapeProductIsRejected() {
        QTemporaryDir dir;
        QVERIFY(dir.isValid());
        const QString path = writeNpy(
            dir,
            QByteArrayLiteral("{'descr': '<f4', 'fortran_order': False, 'shape': (18446744073709551615, 2), }"));
        QVERIFY(!path.isEmpty());

        const auto result = h5reader::io::QtNpyReader::ReadArrayWidened(path);
        QVERIFY(!result.ok);
        QVERIFY(result.error.contains(QStringLiteral("overflows")));
    }
};

QTEST_MAIN(NpyReaderTests)
#include "npy_reader_tests.moc"
