#pragma once
//
// BlessCompare: tolerance-aware NPY comparison for smoke / regression tests.
//
// Replaces byte-identical file comparison with three checks:
//   1. Structural agreement (dtype + shape).
//   2. Element-wise numeric agreement under (rtol, atol).
//   3. Nonzero-fraction sanity — catches silent-empty calculator output.
//
// Why: byte-identity fails on platforms where libgromacs / libtorch /
// OpenBabel ship differently (different LTO, different blas). The data
// is still scientifically equal; the bytes are not. We want a contract
// that distinguishes "the calculator drifted" from "the same number,
// different last bit." We also want to catch the failure mode where a
// calculator silently emits zeros — a byte diff would pass that.
//
// Per-array policy lives in tests/golden/blessed/bless_policy.toml.
// Layout:
//   [default]
//   rtol = 1.0e-5
//   atol = 1.0e-8
//   min_nonzero_fraction = 0.05
//
//   [arrays.<stem>]                # optional per-array overrides
//   rtol = ...
//   atol = ...
//   min_nonzero_fraction = ...
//
// The diagnostic message reports magnitudes ("max |Δ|/|val|=3.2e-7
// across 1231 elements") so a reviewer sees what the drift IS, not
// just that drift exists.
//

#include <string>

namespace nmr::test {

struct BlessPolicy {
    double rtol = 1.0e-5;
    double atol = 1.0e-8;
    double min_nonzero_fraction = 0.05;
};

enum class BlessVerdict {
    Identical,             // bit-identical files
    WithinTolerance,       // numeric drift within (rtol, atol)
    Drifted,               // drift exceeds tolerance
    ZeroOutput,            // run array suspiciously sparse vs blessed
    DtypeOrShapeMismatch,  // structural disagreement
    ReadFailed,            // I/O or parse error on either file
};

const char* VerdictName(BlessVerdict v);

struct BlessResult {
    BlessVerdict verdict = BlessVerdict::Identical;
    std::string  diagnostic;  // "" on Identical; descriptive otherwise
};

// Compare two NPY files under the given policy. Both files must exist.
// Supported dtypes: <f8, <f4, <i4, <i8 (the set NpyWriter emits).
BlessResult CompareNpy(const std::string& run_path,
                       const std::string& blessed_path,
                       const BlessPolicy& policy);

// Read policy for a specific array. Caller passes the array stem (e.g.
// "bs_shielding"); the loader consults [arrays.<stem>] if present, else
// returns the [default] policy. Missing file => default policy + warning
// once per process. Cached internally.
BlessPolicy LoadPolicy(const std::string& policy_toml_path,
                       const std::string& array_stem);

}  // namespace nmr::test
