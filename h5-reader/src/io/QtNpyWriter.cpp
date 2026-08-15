#include "QtNpyWriter.h"

#include <QByteArray>
#include <QFileInfo>
#include <QSaveFile>
#include <QSysInfo>

#include <limits>

namespace h5reader::io {

namespace {

bool ShapeElementCount(const std::vector<std::size_t>& shape,
                       std::size_t* count) {
    if (shape.empty())
        return false;
    std::size_t product = 1;
    for (std::size_t dimension : shape) {
        if (dimension != 0
            && product > std::numeric_limits<std::size_t>::max() / dimension) {
            return false;
        }
        product *= dimension;
    }
    *count = product;
    return true;
}

QByteArray ShapeLiteral(const std::vector<std::size_t>& shape) {
    QByteArray literal("(");
    for (std::size_t i = 0; i < shape.size(); ++i) {
        if (i != 0)
            literal.append(", ");
        literal.append(QByteArray::number(static_cast<qulonglong>(shape[i])));
    }
    if (shape.size() == 1)
        literal.append(',');
    literal.append(')');
    return literal;
}

QByteArray Header(const QByteArray& dtype,
                  const std::vector<std::size_t>& shape,
                  QString* error) {
    QByteArray dictionary("{'descr': '");
    dictionary.append(dtype);
    dictionary.append("', 'fortran_order': False, 'shape': ");
    dictionary.append(ShapeLiteral(shape));
    dictionary.append(", }");

    constexpr qsizetype kPrefixSize = 10;
    constexpr qsizetype kAlignment = 64;
    const qsizetype sizeWithNewline = kPrefixSize + dictionary.size() + 1;
    const qsizetype padding =
        (kAlignment - (sizeWithNewline % kAlignment)) % kAlignment;
    dictionary.append(QByteArray(padding, ' '));
    dictionary.append('\n');
    if (dictionary.size() > std::numeric_limits<std::uint16_t>::max()) {
        *error = QStringLiteral("NPY v1 header is too large");
        return {};
    }
    return dictionary;
}

}  // namespace

QtNpyWriter::Result QtNpyWriter::WriteFloat64(
    const QString& path,
    const std::vector<std::size_t>& shape,
    const std::vector<double>& values) {
    return WriteRaw(path, QByteArrayLiteral("<f8"), shape, values.data(),
                    values.size(), sizeof(double));
}

QtNpyWriter::Result QtNpyWriter::WriteUInt64(
    const QString& path,
    const std::vector<std::size_t>& shape,
    const std::vector<std::uint64_t>& values) {
    return WriteRaw(path, QByteArrayLiteral("<u8"), shape, values.data(),
                    values.size(), sizeof(std::uint64_t));
}

QtNpyWriter::Result QtNpyWriter::WriteUInt8(
    const QString& path,
    const std::vector<std::size_t>& shape,
    const std::vector<std::uint8_t>& values) {
    return WriteRaw(path, QByteArrayLiteral("|u1"), shape, values.data(),
                    values.size(), sizeof(std::uint8_t));
}

QtNpyWriter::Result QtNpyWriter::WriteRaw(
    const QString& path,
    const QByteArray& dtype,
    const std::vector<std::size_t>& shape,
    const void* values,
    std::size_t valueCount,
    std::size_t elementSize) {
    Result result;
    if (QSysInfo::ByteOrder != QSysInfo::LittleEndian && elementSize > 1) {
        result.error = QStringLiteral("NPY writer currently supports little-endian hosts only");
        return result;
    }

    std::size_t expectedCount = 0;
    if (!ShapeElementCount(shape, &expectedCount)) {
        result.error = QStringLiteral("invalid or overflowing NPY shape");
        return result;
    }
    if (expectedCount != valueCount) {
        result.error = QStringLiteral("NPY shape contains %1 values but payload contains %2")
            .arg(static_cast<qulonglong>(expectedCount))
            .arg(static_cast<qulonglong>(valueCount));
        return result;
    }
    if (valueCount != 0 && values == nullptr) {
        result.error = QStringLiteral("NPY payload pointer is null");
        return result;
    }
    if (elementSize != 0
        && valueCount > static_cast<std::size_t>(
            std::numeric_limits<qint64>::max()) / elementSize) {
        result.error = QStringLiteral("NPY payload byte count overflows qint64");
        return result;
    }

    QString headerError;
    const QByteArray header = Header(dtype, shape, &headerError);
    if (header.isEmpty()) {
        result.error = headerError;
        return result;
    }

    QByteArray prefix;
    prefix.append("\x93NUMPY", 6);
    prefix.append(char{1});
    prefix.append(char{0});
    const std::uint16_t headerSize = static_cast<std::uint16_t>(header.size());
    prefix.append(static_cast<char>(headerSize & 0xff));
    prefix.append(static_cast<char>((headerSize >> 8) & 0xff));

    QSaveFile file(path);
    if (!file.open(QIODevice::WriteOnly)) {
        result.error = QStringLiteral("could not open NPY for writing: %1")
            .arg(file.errorString());
        return result;
    }
    if (file.write(prefix) != prefix.size()
        || file.write(header) != header.size()) {
        result.error = QStringLiteral("could not write NPY header: %1")
            .arg(file.errorString());
        file.cancelWriting();
        return result;
    }

    const qint64 payloadBytes = static_cast<qint64>(valueCount * elementSize);
    if (payloadBytes > 0
        && file.write(static_cast<const char*>(values), payloadBytes) != payloadBytes) {
        result.error = QStringLiteral("could not write NPY payload: %1")
            .arg(file.errorString());
        file.cancelWriting();
        return result;
    }
    if (!file.commit()) {
        result.error = QStringLiteral("could not commit NPY: %1")
            .arg(file.errorString());
        return result;
    }

    result.fileSize = QFileInfo(path).size();
    result.ok = true;
    return result;
}

}  // namespace h5reader::io
