#pragma once

#include <QByteArray>
#include <QDir>
#include <QFileInfo>
#include <QIODevice>
#include <QSaveFile>
#include <QString>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace h5reader::rediscover {

class NpyWriter {
public:
    template <typename T>
    static bool Write(const QString& path,
                      const std::vector<std::size_t>& shape,
                      const std::vector<T>& data,
                      const QByteArray& dtypeDescr);

    static bool WriteF64(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<double>& data);
    static bool WriteF32(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<float>& data);
    static bool WriteI64(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<std::int64_t>& data);
    static bool WriteI32(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<std::int32_t>& data);
    static bool WriteI16(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<std::int16_t>& data);
    static bool WriteU8(const QString& path,
                        const std::vector<std::size_t>& shape,
                        const std::vector<std::uint8_t>& data);
};

template <typename T>
bool NpyWriter::Write(const QString& path,
                      const std::vector<std::size_t>& shape,
                      const std::vector<T>& data,
                      const QByteArray& dtypeDescr) {
    if (shape.empty()) return false;
    std::size_t expected = 1;
    for (const std::size_t dim : shape) expected *= dim;
    if (expected != data.size()) return false;

    QDir().mkpath(QFileInfo(path).dir().absolutePath());

    QByteArray header;
    header += "{'descr': '";
    header += dtypeDescr;
    header += "', 'fortran_order': False, 'shape': (";
    for (std::size_t i = 0; i < shape.size(); ++i) {
        if (i) header += ", ";
        header += QByteArray::number(static_cast<qulonglong>(shape[i]));
    }
    if (shape.size() == 1) header += ",";
    header += "), }";

    constexpr int kPreambleBytes = 10;
    const int newlineBytes = 1;
    int pad = 16 - ((kPreambleBytes + header.size() + newlineBytes) % 16);
    if (pad == 16) pad = 0;
    header += QByteArray(pad, ' ');
    header += '\n';
    if (header.size() > 65535) return false;

    QByteArray prefix;
    prefix.append("\x93NUMPY", 6);
    prefix.append(char(1));
    prefix.append(char(0));
    const quint16 headerLen = static_cast<quint16>(header.size());
    prefix.append(char(headerLen & 0xff));
    prefix.append(char((headerLen >> 8) & 0xff));

    QSaveFile f(path);
    if (!f.open(QIODevice::WriteOnly)) return false;
    if (f.write(prefix) != prefix.size()) return false;
    if (f.write(header) != header.size()) return false;
    const qsizetype payloadBytes = static_cast<qsizetype>(data.size() * sizeof(T));
    if (payloadBytes > 0
        && f.write(reinterpret_cast<const char*>(data.data()), payloadBytes) != payloadBytes)
        return false;
    return f.commit();
}

}  // namespace h5reader::rediscover
