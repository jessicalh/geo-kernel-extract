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

    // Write a 1D array of int32.
    static bool WriteInt32(const std::string& path,
                           const int32_t* data,
                           size_t count) {
        return Write(path, "<i4", data, count * sizeof(int32_t),
                     {count});
    }

private:
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
