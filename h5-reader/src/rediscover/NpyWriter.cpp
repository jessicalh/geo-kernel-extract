#include "NpyWriter.h"

namespace h5reader::rediscover {

bool NpyWriter::WriteF64(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<double>& data) {
    return Write<double>(path, shape, data, QByteArray("<f8"));
}

bool NpyWriter::WriteF32(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<float>& data) {
    return Write<float>(path, shape, data, QByteArray("<f4"));
}

bool NpyWriter::WriteI64(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<std::int64_t>& data) {
    return Write<std::int64_t>(path, shape, data, QByteArray("<i8"));
}

bool NpyWriter::WriteI32(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<std::int32_t>& data) {
    return Write<std::int32_t>(path, shape, data, QByteArray("<i4"));
}

bool NpyWriter::WriteI16(const QString& path,
                         const std::vector<std::size_t>& shape,
                         const std::vector<std::int16_t>& data) {
    return Write<std::int16_t>(path, shape, data, QByteArray("<i2"));
}

bool NpyWriter::WriteU8(const QString& path,
                        const std::vector<std::size_t>& shape,
                        const std::vector<std::uint8_t>& data) {
    return Write<std::uint8_t>(path, shape, data, QByteArray("|u1"));
}

}  // namespace h5reader::rediscover
