#pragma once
//
// NpyWriter: write arrays to numpy .npy format.
//
// Format only. No domain knowledge. The callers (ConformationResult
// subclasses) know what data to write and what to call it.
//
// NPY v1.0: 6-byte magic, 2-byte version, 2-byte header length,
// Python dict header, padding to 64-byte alignment, raw data.
//

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdint>
#include <limits>

namespace nmr {

class NpyWriter {
public:
    // Write a 2D array of float64. cols=1 writes as 1D shape (rows,).
    static bool WriteFloat64(const std::string& path,
                             const double* data,
                             size_t rows, size_t cols = 1) {
        return Write(path, "<f8", data, rows * cols * sizeof(double),
                     cols == 1 ? std::vector<size_t>{rows}
                               : std::vector<size_t>{rows, cols});
    }

    // Write an arbitrary-rank float64 array. The caller supplies data in
    // C-contiguous order matching shape.
    static bool WriteFloat64(const std::string& path,
                             const double* data,
                             const std::vector<size_t>& shape) {
        size_t count = 0;
        if (!ElementCount(shape, count)) return false;
        return Write(path, "<f8", data, count * sizeof(double), shape);
    }

    // Write a 2D array of float32. cols=1 writes as 1D shape (rows,).
    static bool WriteFloat32(const std::string& path,
                             const float* data,
                             size_t rows, size_t cols = 1) {
        return Write(path, "<f4", data, rows * cols * sizeof(float),
                     cols == 1 ? std::vector<size_t>{rows}
                               : std::vector<size_t>{rows, cols});
    }

    // Write a 1D array of int32.
    static bool WriteInt32(const std::string& path,
                           const int32_t* data,
                           size_t count) {
        return Write(path, "<i4", data, count * sizeof(int32_t),
                     {count});
    }

    // Write a 2D array of int32. cols=1 writes the same 1D shape as the
    // historical overload.
    static bool WriteInt32(const std::string& path,
                           const int32_t* data,
                           size_t rows, size_t cols) {
        return Write(path, "<i4", data, rows * cols * sizeof(int32_t),
                     cols == 1 ? std::vector<size_t>{rows}
                               : std::vector<size_t>{rows, cols});
    }

    // Write an arbitrary-rank int32 array. The caller supplies data in
    // C-contiguous order matching shape.
    static bool WriteInt32(const std::string& path,
                           const int32_t* data,
                           const std::vector<size_t>& shape) {
        size_t count = 0;
        if (!ElementCount(shape, count)) return false;
        return Write(path, "<i4", data, count * sizeof(int32_t), shape);
    }

    // Write a 1D array of int8. Useful for boolean masks where int32
    // would be 4× the disk size for no information gain.
    static bool WriteInt8(const std::string& path,
                          const int8_t* data,
                          size_t count) {
        return Write(path, "|i1", data, count * sizeof(int8_t),
                     {count});
    }

    // Write a 2D array of int8. cols=1 writes the same 1D shape as the
    // historical overload.
    static bool WriteInt8(const std::string& path,
                          const int8_t* data,
                          size_t rows, size_t cols) {
        return Write(path, "|i1", data, rows * cols * sizeof(int8_t),
                     cols == 1 ? std::vector<size_t>{rows}
                               : std::vector<size_t>{rows, cols});
    }

    // Write unsigned-byte masks without changing their sentinel/value range.
    static bool WriteUInt8(const std::string& path,
                           const std::uint8_t* data,
                           size_t count) {
        return Write(path, "|u1", data, count * sizeof(std::uint8_t),
                     {count});
    }

    static bool WriteUInt8(const std::string& path,
                           const std::uint8_t* data,
                           size_t rows, size_t cols) {
        return Write(path, "|u1", data,
                     rows * cols * sizeof(std::uint8_t),
                     cols == 1 ? std::vector<size_t>{rows}
                               : std::vector<size_t>{rows, cols});
    }

private:
    static bool ElementCount(const std::vector<size_t>& shape,
                             size_t& count) {
        if (shape.empty()) return false;
        count = 1;
        for (size_t extent : shape) {
            if (extent != 0 &&
                count > std::numeric_limits<size_t>::max() / extent) {
                return false;
            }
            count *= extent;
        }
        return true;
    }

    static bool Write(const std::string& path,
                      const char* descr,
                      const void* data,
                      size_t data_bytes,
                      const std::vector<size_t>& shape) {

        std::ostringstream hdr;
        hdr << "{'descr': '" << descr
            << "', 'fortran_order': False, 'shape': (";
        for (size_t i = 0; i < shape.size(); ++i) {
            hdr << shape[i];
            if (i + 1 < shape.size() || shape.size() == 1) hdr << ",";
        }
        hdr << "), }";

        std::string h = hdr.str();
        size_t total = 10 + h.size() + 1;
        size_t pad = (64 - (total % 64)) % 64;
        h.append(pad, ' ');
        h += '\n';

        std::ofstream out(path, std::ios::binary);
        if (!out.is_open()) return false;

        const char magic[] = "\x93NUMPY";
        out.write(magic, 6);
        const char version[] = "\x01\x00";
        out.write(version, 2);
        uint16_t hlen = static_cast<uint16_t>(h.size());
        out.write(reinterpret_cast<const char*>(&hlen), 2);
        out.write(h.data(), static_cast<std::streamsize>(h.size()));
        out.write(reinterpret_cast<const char*>(data),
                  static_cast<std::streamsize>(data_bytes));
        return out.good();
    }
};

}  // namespace nmr
