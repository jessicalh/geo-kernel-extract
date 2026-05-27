// QtNpyReader implementation — the non-templated parts. The
// ReadStructured<RecordT> template lives in the header so callers can
// instantiate it on each NPY row struct without an explicit-
// instantiation list to maintain.

#include "QtNpyReader.h"

#include "../diagnostics/ErrorBus.h"

#include <QByteArray>
#include <QFile>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <string>

namespace h5reader::io {

// ── Local helpers ──────────────────────────────────────────────────

namespace {

// Find the (single-quoted) Python literal that follows a `'key':` token.
// Tolerates whitespace; does not unescape the value. Returns std::string::npos
// for both indices if the key is not present.
struct Located {
    std::size_t value_begin = std::string::npos;  // first char after the opening quote/bracket
    std::size_t value_end = std::string::npos;    // last char before the closing quote/bracket
};

// Extract the value of `'descr': <value>` as a verbatim slice. The
// value may be a string literal `'<<i4'` OR a list of tuples
// `[('a', '<i4'), ...]`. We capture from the first non-space after the
// colon to the matching closing character (quote for strings, bracket
// for lists). The result includes the brackets/quotes — the caller's
// substring check sees the same shape it spelled.
Located LocateDescrValue(const std::string& hdr) {
    Located out;
    const std::size_t k = hdr.find("'descr'");
    if (k == std::string::npos)
        return out;
    std::size_t colon = hdr.find(':', k);
    if (colon == std::string::npos)
        return out;
    std::size_t i = colon + 1;
    // Skip whitespace.
    while (i < hdr.size() && (hdr[i] == ' ' || hdr[i] == '\t'))
        ++i;
    if (i >= hdr.size())
        return out;

    const char open = hdr[i];
    const char close = (open == '[') ? ']' : (open == '\'' || open == '"') ? open : 0;
    if (close == 0)
        return out;

    out.value_begin = i;
    // For brackets we need to track nesting; for quotes we walk to the
    // next matching quote. NPY headers don't escape quotes inside
    // values, so a simple scan suffices.
    if (open == '[') {
        int depth = 0;
        for (std::size_t j = i; j < hdr.size(); ++j) {
            if (hdr[j] == '[')
                ++depth;
            else if (hdr[j] == ']') {
                --depth;
                if (depth == 0) {
                    out.value_end = j;
                    break;
                }
            }
        }
    } else {
        std::size_t j = i + 1;
        while (j < hdr.size() && hdr[j] != close)
            ++j;
        if (j < hdr.size())
            out.value_end = j;
    }
    return out;
}

// Extract the `'shape': (R,)` or `'shape': (R, C)` tuple as a vector
// of std::size_t. Returns empty on parse failure.
std::vector<std::size_t> ParseShape(const std::string& hdr) {
    std::vector<std::size_t> out;
    const std::size_t k = hdr.find("'shape'");
    if (k == std::string::npos)
        return out;
    std::size_t lp = hdr.find('(', k);
    if (lp == std::string::npos)
        return out;
    std::size_t rp = hdr.find(')', lp);
    if (rp == std::string::npos)
        return out;
    std::string inside = hdr.substr(lp + 1, rp - lp - 1);

    // Split by ',' and parse each token as an unsigned integer; skip
    // empties (Python tuple `(N,)` has a trailing comma).
    std::size_t pos = 0;
    while (pos < inside.size()) {
        std::size_t comma = inside.find(',', pos);
        std::string tok = (comma == std::string::npos) ? inside.substr(pos) : inside.substr(pos, comma - pos);
        // Strip whitespace.
        std::size_t lo = tok.find_first_not_of(" \t");
        std::size_t hi = tok.find_last_not_of(" \t");
        if (lo != std::string::npos && hi != std::string::npos) {
            const std::string trimmed = tok.substr(lo, hi - lo + 1);
            try {
                out.push_back(std::stoull(trimmed));
            } catch (...) {
                // ignore; will fall through to validation
            }
        }
        if (comma == std::string::npos)
            break;
        pos = comma + 1;
    }
    return out;
}

// Read the `'fortran_order': True/False` value. Defaults to False (the
// writer always emits row-major) but parse what's actually in the file.
bool ParseFortranOrder(const std::string& hdr) {
    const std::size_t k = hdr.find("'fortran_order'");
    if (k == std::string::npos)
        return false;
    const std::size_t colon = hdr.find(':', k);
    if (colon == std::string::npos)
        return false;
    // Search the rest of the header for the next True/False token.
    const std::size_t t = hdr.find("True", colon);
    const std::size_t f = hdr.find("False", colon);
    if (t == std::string::npos && f == std::string::npos)
        return false;
    if (t == std::string::npos)
        return false;
    if (f == std::string::npos)
        return true;
    return t < f;
}

}  // namespace

// ── Header parser ──────────────────────────────────────────────────

QtNpyReader::ParsedHeader QtNpyReader::ParseHeader(const QByteArray& bytes, const QString& path_for_diagnostics) {
    ParsedHeader h;

    if (bytes.size() < 10) {
        h.error = QStringLiteral("QtNpyReader: file too small to contain NPY header: %1 (size=%2)")
                      .arg(path_for_diagnostics)
                      .arg(bytes.size());
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                h.error,
                                                path_for_diagnostics);
        return h;
    }

    // Magic
    static constexpr char kMagic[6] = {static_cast<char>('\x93'), 'N', 'U', 'M', 'P', 'Y'};
    if (std::memcmp(bytes.constData(), kMagic, 6) != 0) {
        h.error = QStringLiteral("QtNpyReader: bad magic in %1").arg(path_for_diagnostics);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                h.error,
                                                path_for_diagnostics);
        return h;
    }

    // Version
    const std::uint8_t major = static_cast<std::uint8_t>(bytes[6]);
    const std::uint8_t minor = static_cast<std::uint8_t>(bytes[7]);
    if (major != 1 || minor != 0) {
        h.error = QStringLiteral("QtNpyReader: unsupported NPY version %1.%2 in %3 (only 1.0 handled)")
                      .arg(major)
                      .arg(minor)
                      .arg(path_for_diagnostics);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                h.error,
                                                path_for_diagnostics);
        return h;
    }

    // Header length (uint16, little-endian)
    const std::uint16_t header_len =
        static_cast<std::uint8_t>(bytes[8]) | (static_cast<std::uint16_t>(static_cast<std::uint8_t>(bytes[9])) << 8);

    constexpr int kPreHeader = 10;
    if (static_cast<int>(bytes.size()) < kPreHeader + header_len) {
        h.error = QStringLiteral("QtNpyReader: header truncated in %1 (need %2 bytes, got %3)")
                      .arg(path_for_diagnostics)
                      .arg(kPreHeader + header_len)
                      .arg(bytes.size());
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                h.error,
                                                path_for_diagnostics);
        return h;
    }

    const std::string hdr_str(bytes.constData() + kPreHeader, header_len);

    const Located descr = LocateDescrValue(hdr_str);
    if (descr.value_begin == std::string::npos) {
        h.error = QStringLiteral("QtNpyReader: 'descr' key missing in %1").arg(path_for_diagnostics);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                h.error,
                                                path_for_diagnostics);
        return h;
    }
    h.descr_substring = hdr_str.substr(descr.value_begin, descr.value_end - descr.value_begin + 1);

    h.fortran_order = ParseFortranOrder(hdr_str);
    if (h.fortran_order) {
        h.error = QStringLiteral("QtNpyReader: fortran_order=True not supported in %1").arg(path_for_diagnostics);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                h.error,
                                                path_for_diagnostics);
        return h;
    }

    h.shape = ParseShape(hdr_str);
    if (h.shape.empty()) {
        h.error = QStringLiteral("QtNpyReader: 'shape' parse failed in %1").arg(path_for_diagnostics);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                h.error,
                                                path_for_diagnostics);
        return h;
    }

    h.data_offset = static_cast<std::size_t>(kPreHeader) + header_len;
    h.ok = true;
    return h;
}


// ── Raw-bytes overload (unused by sidecar; reserved for raw scalar NPYs).
QtNpyReader::StructuredResult
QtNpyReader::ReadRawBytes(const QString& path, std::vector<unsigned char>& out_bytes, std::size_t& out_record_size) {
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
    if (hdr.shape.size() != 1) {
        r.error = QStringLiteral("QtNpyReader::ReadRawBytes: expected 1-D shape in %1").arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                r.error,
                                                path);
        return r;
    }
    const std::size_t raw = static_cast<std::size_t>(bytes.size()) - hdr.data_offset;
    if (hdr.shape[0] == 0 || raw % hdr.shape[0] != 0) {
        r.error = QStringLiteral("QtNpyReader::ReadRawBytes: raw bytes %1 not divisible by row count %2 in %3")
                      .arg(raw)
                      .arg(hdr.shape[0])
                      .arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"),
                                                r.error,
                                                path);
        return r;
    }
    out_record_size = raw / hdr.shape[0];
    out_bytes.assign(reinterpret_cast<const unsigned char*>(bytes.constData() + hdr.data_offset),
                     reinterpret_cast<const unsigned char*>(bytes.constData()) + bytes.size());
    r.ok = true;
    r.row_count = hdr.shape[0];
    return r;
}


// ── 2-D / 1-D numeric read with dtype widening ──────────────────────
QtNpyReader::WidenedArray QtNpyReader::ReadArrayWidened(const QString& path) {
    WidenedArray r;

    QFile f(path);
    if (!f.open(QIODevice::ReadOnly)) {
        r.error = QStringLiteral("QtNpyReader::ReadArrayWidened: could not open %1").arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"), r.error, path);
        return r;
    }
    const QByteArray bytes = f.readAll();
    f.close();

    const ParsedHeader hdr = ParseHeader(bytes, path);  // logs its own structural errors
    if (!hdr.ok) {
        r.error = hdr.error;
        return r;
    }
    r.descr = hdr.descr_substring;

    // Rank → (rows, cols). 1-D is rows x 1; the caller reshapes a Protein-axis
    // 1-D row if needed.
    std::size_t rows = 0;
    std::size_t cols = 0;
    if (hdr.shape.size() == 1) {
        rows = hdr.shape[0];
        cols = 1;
    } else if (hdr.shape.size() == 2) {
        rows = hdr.shape[0];
        cols = hdr.shape[1];
    } else {
        r.error = QStringLiteral("QtNpyReader::ReadArrayWidened: rank %1 unsupported in %2")
                      .arg(hdr.shape.size())
                      .arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"), r.error, path);
        return r;
    }

    // Reject structured dtypes (list-of-tuples descr) — those are sidecar
    // records read via ReadStructured, not numeric calculator arrays.
    if (hdr.descr_substring.find('[') != std::string::npos) {
        r.error = QStringLiteral("QtNpyReader::ReadArrayWidened: structured dtype %1 in %2")
                      .arg(QString::fromStdString(hdr.descr_substring), path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"), r.error, path);
        return r;
    }

    // Parse dtype kind + element size from the descr ('<f8' -> f,8; '|i1' ->
    // i,1). Byte order assumed little-endian (Intel/AMD/ARM64 default).
    char kind = 0;
    int esize = 0;
    for (std::size_t i = 0; i < hdr.descr_substring.size(); ++i) {
        const char c = hdr.descr_substring[i];
        if (c == 'f' || c == 'i' || c == 'u' || c == 'b') {
            kind = c;
            std::string digits;
            for (std::size_t j = i + 1;
                 j < hdr.descr_substring.size() && std::isdigit(static_cast<unsigned char>(hdr.descr_substring[j]));
                 ++j)
                digits += hdr.descr_substring[j];
            if (!digits.empty())
                esize = std::stoi(digits);
            break;
        }
    }
    const bool supported = (kind == 'f' && (esize == 4 || esize == 8)) ||
                           ((kind == 'i' || kind == 'u' || kind == 'b') &&
                            (esize == 1 || esize == 2 || esize == 4 || esize == 8));
    if (!supported) {
        r.error = QStringLiteral("QtNpyReader::ReadArrayWidened: unsupported dtype %1 in %2")
                      .arg(QString::fromStdString(hdr.descr_substring), path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"), r.error, path);
        return r;
    }

    const std::size_t n = rows * cols;
    const std::size_t raw = static_cast<std::size_t>(bytes.size()) - hdr.data_offset;
    if (raw != n * static_cast<std::size_t>(esize)) {
        r.error = QStringLiteral("QtNpyReader::ReadArrayWidened: byte count %1 != rows*cols*esize %2 in %3")
                      .arg(raw)
                      .arg(n * static_cast<std::size_t>(esize))
                      .arg(path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtNpyReader"), r.error, path);
        return r;
    }

    r.data.resize(n);
    const unsigned char* p = reinterpret_cast<const unsigned char*>(bytes.constData()) + hdr.data_offset;
    for (std::size_t i = 0; i < n; ++i) {
        const unsigned char* q = p + i * static_cast<std::size_t>(esize);
        double v = 0.0;
        if (kind == 'f' && esize == 8) {
            double t;
            std::memcpy(&t, q, 8);
            v = t;
        } else if (kind == 'f' && esize == 4) {
            float t;
            std::memcpy(&t, q, 4);
            v = static_cast<double>(t);
        } else if (kind == 'i' && esize == 8) {
            std::int64_t t;
            std::memcpy(&t, q, 8);
            v = static_cast<double>(t);
        } else if (kind == 'i' && esize == 4) {
            std::int32_t t;
            std::memcpy(&t, q, 4);
            v = static_cast<double>(t);
        } else if (kind == 'i' && esize == 2) {
            std::int16_t t;
            std::memcpy(&t, q, 2);
            v = static_cast<double>(t);
        } else if (kind == 'i' && esize == 1) {
            std::int8_t t;
            std::memcpy(&t, q, 1);
            v = static_cast<double>(t);
        } else {  // 'u' or 'b'
            if (esize == 8) {
                std::uint64_t t;
                std::memcpy(&t, q, 8);
                v = static_cast<double>(t);
            } else if (esize == 4) {
                std::uint32_t t;
                std::memcpy(&t, q, 4);
                v = static_cast<double>(t);
            } else if (esize == 2) {
                std::uint16_t t;
                std::memcpy(&t, q, 2);
                v = static_cast<double>(t);
            } else {
                std::uint8_t t;
                std::memcpy(&t, q, 1);
                v = static_cast<double>(t);
            }
        }
        r.data[i] = v;
    }

    r.rows = rows;
    r.cols = cols;
    r.ok = true;
    return r;
}

}  // namespace h5reader::io
