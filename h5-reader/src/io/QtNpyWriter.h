#pragma once

#include <QByteArray>
#include <QString>

#include <cstddef>
#include <cstdint>
#include <vector>

namespace h5reader::io {

class QtNpyWriter final {
public:
    struct Result {
        bool ok = false;
        QString error;
        qint64 fileSize = 0;
    };

    static Result WriteFloat64(const QString& path,
                               const std::vector<std::size_t>& shape,
                               const std::vector<double>& values);
    static Result WriteUInt64(const QString& path,
                              const std::vector<std::size_t>& shape,
                              const std::vector<std::uint64_t>& values);
    static Result WriteUInt8(const QString& path,
                             const std::vector<std::size_t>& shape,
                             const std::vector<std::uint8_t>& values);

private:
    static Result WriteRaw(const QString& path,
                           const QByteArray& dtype,
                           const std::vector<std::size_t>& shape,
                           const void* values,
                           std::size_t valueCount,
                           std::size_t elementSize);
};

}  // namespace h5reader::io
