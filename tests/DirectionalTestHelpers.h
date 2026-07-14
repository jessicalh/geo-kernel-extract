#pragma once

// Deterministic O(3) test support for the directional-output freeze gate.
//
// The helpers in this file transform calculator *inputs*.  Covariance tests
// must construct a fresh ProteinConformation and rerun the owning production
// Compute() method; applying these functions to an already-produced output is
// only suitable for forming the expected value.

#include "SolventEnvironment.h"
#include "Types.h"

#include <Eigen/Geometry>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <random>
#include <string>
#include <vector>

namespace nmr::test::directional {

struct OrthogonalTransform {
    Mat3 Q = Mat3::Identity();
    Vec3 translation = Vec3::Zero();

    double Determinant() const { return Q.determinant(); }
    bool IsImproper() const { return Determinant() < 0.0; }
};

// The seed is recorded at every call site.  A proper axis-angle rotation is
// drawn first; for an improper transform, a fixed x reflection is composed on
// the input side.  Thus Q is orthogonal to roundoff and det(Q) is exactly the
// requested sign to roundoff.
inline OrthogonalTransform SeededTransform(std::uint64_t seed,
                                           bool improper) {
    std::mt19937_64 rng(seed);
    std::normal_distribution<double> normal(0.0, 1.0);
    std::uniform_real_distribution<double> angle(-2.6, 2.6);
    std::uniform_real_distribution<double> shift(-7.0, 7.0);

    Vec3 axis(normal(rng), normal(rng), normal(rng));
    axis.normalize();
    Mat3 Q = Eigen::AngleAxisd(angle(rng), axis).toRotationMatrix();
    if (improper) {
        Mat3 reflection = Mat3::Identity();
        reflection(0, 0) = -1.0;
        Q = Q * reflection;
    }

    OrthogonalTransform out;
    out.Q = Q;
    out.translation = Vec3(shift(rng), shift(rng), shift(rng));
    return out;
}

inline Vec3 Position(const OrthogonalTransform& x, const Vec3& p) {
    return x.Q * p + x.translation;
}

inline Vec3 Polar(const OrthogonalTransform& x, const Vec3& v) {
    return x.Q * v;
}

inline Vec3 Axial(const OrthogonalTransform& x, const Vec3& v) {
    return x.Determinant() * x.Q * v;
}

inline Mat3 EvenRank2(const OrthogonalTransform& x, const Mat3& tensor) {
    return x.Q * tensor * x.Q.transpose();
}

inline std::vector<Vec3> Positions(const OrthogonalTransform& x,
                                   const std::vector<Vec3>& positions) {
    std::vector<Vec3> out;
    out.reserve(positions.size());
    for (const Vec3& p : positions) out.push_back(Position(x, p));
    return out;
}

inline SolventEnvironment Solvent(const OrthogonalTransform& x,
                                  const SolventEnvironment& source) {
    SolventEnvironment out = source;
    for (WaterMolecule& water : out.waters) {
        water.O_pos = Position(x, water.O_pos);
        water.H1_pos = Position(x, water.H1_pos);
        water.H2_pos = Position(x, water.H2_pos);
    }
    for (Ion& ion : out.ions) ion.pos = Position(x, ion.pos);
    out.water_O_positions.clear();
    out.water_O_positions.reserve(out.waters.size());
    for (const WaterMolecule& water : out.waters)
        out.water_O_positions.push_back(water.O_pos);
    return out;
}

inline Vec3 T1Vector(const SphericalTensor& tensor) {
    return Vec3(tensor.T1[0], tensor.T1[1], tensor.T1[2]);
}

inline Mat3 T2Matrix(const SphericalTensor& tensor) {
    SphericalTensor t2 = tensor;
    t2.T0 = 0.0;
    t2.T1 = {0.0, 0.0, 0.0};
    return t2.Reconstruct();
}

// Project-native-T2 round trip used by calculator tests:
// native five-vector -> Cartesian symmetric/traceless -> Q T Q^T ->
// project-native decomposition.  This deliberately does not use e3nn.
inline SphericalTensor RotateNativeT2(const OrthogonalTransform& x,
                                      const SphericalTensor& tensor) {
    return SphericalTensor::Decompose(EvenRank2(x, T2Matrix(tensor)));
}

inline double MatrixMaxAbs(const Mat3& value) {
    return value.cwiseAbs().maxCoeff();
}

inline bool Near(double actual, double expected,
                 double abs_tolerance, double rel_tolerance) {
    return std::abs(actual - expected) <=
           abs_tolerance + rel_tolerance * std::abs(expected);
}

inline bool NearVector(const Vec3& actual, const Vec3& expected,
                       double abs_tolerance, double rel_tolerance) {
    return (actual - expected).norm() <=
           abs_tolerance + rel_tolerance * expected.norm();
}

inline bool NearMatrix(const Mat3& actual, const Mat3& expected,
                       double abs_tolerance, double rel_tolerance) {
    return (actual - expected).norm() <=
           abs_tolerance + rel_tolerance * expected.norm();
}

inline double WrapPi(double angle) {
    constexpr double kPi = 3.141592653589793238462643383279502884;
    while (angle <= -kPi) angle += 2.0 * kPi;
    while (angle > kPi) angle -= 2.0 * kPi;
    return angle;
}

inline double CircularDifference(double a, double b) {
    return WrapPi(a - b);
}

inline std::string ShellQuote(const std::filesystem::path& path) {
    std::string quoted = "'";
    for (const char ch : path.string()) {
        if (ch == '\'') quoted += "'\\''";
        else quoted += ch;
    }
    quoted += "'";
    return quoted;
}

// Exercise NumPy itself at the serialized boundary.  The Python helper uses
// np.load(..., allow_pickle=False) for every named file; calculator tests keep
// their independent physics oracles in C++ and use this only as the final
// reader/format gate.
inline int RunNumpyAllowPickleFalse(
        const std::filesystem::path& python_executable,
        const std::filesystem::path& script,
        const std::vector<std::filesystem::path>& npy_paths) {
    std::string command = ShellQuote(python_executable) + " " +
                          ShellQuote(script);
    for (const auto& path : npy_paths) command += " " + ShellQuote(path);
    return std::system(command.c_str());
}

}  // namespace nmr::test::directional
