// QtNpyReader — NPY v1.0 binary loader.
//
// Reads a numpy .npy file's header (magic + version + header dict),
// validates the dtype descriptor string matches an expected value and
// the shape matches an expected (rows,) 1-D layout, then loads the raw
// bytes into a caller-supplied typed buffer via memcpy. Used by
// QtTopologySidecar for the 5 sidecar files (residues, bonds, rings,
// ring_membership, atoms_category_info).
//
// Why hand-roll this rather than vendor a library: NPY 1.0 is a
// stable, ~80-line-of-code format (6-byte magic + 2-byte version +
// 2-byte LE header length + ASCII Python dict literal + raw bytes).
// Cross-platform via QFile (no POSIX); zero dependency surface beyond
// Qt; failure modes log via ErrorBus with the specific dtype mismatch
// for forensics. See src/NpyWriter.h for the writer-side spec.
//
// Cross-platform discipline: QFile + QByteArray (no fread, no mmap).
// All multi-byte ints assumed little-endian (Intel/AMD/ARM64 default;
// the writer emits LE explicitly per the dtype prefix).

#pragma once

#include <QString>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace h5reader::io {

class QtNpyReader {
public:
    // Result of a successful structured-record read.
    struct StructuredResult {
        bool ok = false;
        QString error;              // empty on success
        std::size_t row_count = 0;  // == out_buffer size (in records)
        std::string dtype_descr;    // verbatim from header, for diagnostics
    };

    // Read a 1-D structured NPY file into `out` (caller-supplied
    // vector). On entry `out` is resized to row_count and the raw bytes
    // are memcpy'd in. Validates:
    //   - file exists, opens
    //   - magic == "\x93NUMPY", version == 1.0
    //   - dtype descr string contains `expected_dtype_substr` (a
    //     fragment like `'descr': [('atom_index', '<i4'),`) so the
    //     caller can pin the schema. Pass empty string to skip the
    //     check (not recommended; the byte layout is the runtime
    //     contract).
    //   - sizeof(RecordT) divides the raw-data length evenly; the
    //     remainder bytes are the implicit record count.
    //   - shape is (row_count,) — 1-D only. 2-D arrays go through a
    //     different overload (not yet needed).
    //
    // Failure modes log via diagnostics::ErrorBus (Severity::Error)
    // with the specific dtype mismatch or row-count parse failure.
    template <typename RecordT>
    static StructuredResult
    ReadStructured(const QString& path, const std::string& expected_dtype_substr, std::vector<RecordT>& out);

    // Read a 1-D NPY of fixed-width scalars (e.g. <f8 doubles, <i8
    // uint64) into a raw byte buffer. The caller validates dtype via
    // ParsedHeader inspection.
    //
    // Not used by the sidecar (the structured-record overload covers
    // it); reserved for future per-trajectory-NPY reads.
    static StructuredResult
    ReadRawBytes(const QString& path, std::vector<unsigned char>& out_bytes, std::size_t& out_record_size);

private:
    // Parsed-header struct (file-scope; not exposed in the public API).
    struct ParsedHeader {
        bool ok = false;
        QString error;
        std::string descr_substring;  // verbatim header text for diagnostics
        bool fortran_order = false;
        std::vector<std::size_t> shape;
        std::size_t data_offset = 0;  // bytes from file start to raw data
    };

    // Parses the NPY header from the prefix of `bytes`. Returns
    // ok=false with an error message if anything is structurally wrong.
    static ParsedHeader ParseHeader(const QByteArray& bytes, const QString& path_for_diagnostics);
};

// ── Template implementation ────────────────────────────────────────
//
// Lives in the header for the single template parameter the call sites
// instantiate (QtNpyResidueRow / QtNpyBondRow / QtNpyRingRow /
// QtNpyRingMembershipRow / QtNpyAtomCategoryRow). If the instantiation
// set grows past five, fold the body into the .cpp via explicit
// instantiations.

}  // namespace h5reader::io


// The template implementation must be visible to every TU that
// instantiates it. The non-templated helpers live in the .cpp.

#include <QByteArray>
#include <QFile>
#include <cstring>

#include "../diagnostics/ErrorBus.h"

namespace h5reader::io {

template <typename RecordT>
QtNpyReader::StructuredResult
QtNpyReader::ReadStructured(const QString& path, const std::string& expected_dtype_substr, std::vector<RecordT>& out) {
    StructuredResult r;

    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) {
        r.error = QStringLiteral("QtNpyReader: could not open %1").arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                r.error,
                                                path);
        return r;
    }
    const QByteArray bytes = f.readAll();
    f.close();

    const ParsedHeader hdr = ParseHeader(bytes, path);
    if (!hdr.ok) {
        r.error = hdr.error;
        return r;
    }
    r.dtype_descr = hdr.descr_substring;

    // Dtype substring check.
    if (!expected_dtype_substr.empty() && hdr.descr_substring.find(expected_dtype_substr) == std::string::npos) {
        r.error = QStringLiteral("QtNpyReader: dtype descr mismatch in %1\n"
                                 "  expected substring: %2\n"
                                 "  got: %3")
                      .arg(path, QString::fromStdString(expected_dtype_substr), QString::fromStdString(hdr.descr_substring));
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                r.error,
                                                path);
        return r;
    }

    // Shape must be 1-D.
    if (hdr.shape.size() != 1) {
        r.error = QStringLiteral("QtNpyReader: expected 1-D shape, got rank=%1 in %2").arg(hdr.shape.size()).arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                r.error,
                                                path);
        return r;
    }
    const std::size_t expected_rows = hdr.shape[0];

    // Row size check against the typed struct.
    const std::size_t raw_data_len = static_cast<std::size_t>(bytes.size()) - hdr.data_offset;
    const std::size_t expected_bytes = expected_rows * sizeof(RecordT);
    if (raw_data_len != expected_bytes) {
        r.error = QStringLiteral("QtNpyReader: byte count mismatch in %1: "
                                 "raw=%2 expected=%3 (rows=%4, sizeof(RecordT)=%5)")
                      .arg(path)
                      .arg(raw_data_len)
                      .arg(expected_bytes)
                      .arg(expected_rows)
                      .arg(sizeof(RecordT));
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                r.error,
                                                path);
        return r;
    }

    out.resize(expected_rows);
    std::memcpy(out.data(), bytes.constData() + hdr.data_offset, expected_bytes);

    r.ok = true;
    r.row_count = expected_rows;
    return r;
}

}  // namespace h5reader::io
