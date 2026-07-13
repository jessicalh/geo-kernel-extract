#include "LocalBackboneGeometryResult.h"

#include "NpyWriter.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace nmr {

namespace local_backbone_geometry {

namespace {

constexpr double kNearZero = 1e-12;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

Vec3 NanVec3() {
    return Vec3(kNaN, kNaN, kNaN);
}

bool Finite(const Vec3& v) {
    return v.allFinite();
}

bool HasAtom(std::size_t atom_index) {
    return atom_index != Residue::NONE;
}

}  // namespace


double BondAngleRadians(const Vec3& endpoint_a,
                        const Vec3& vertex,
                        const Vec3& endpoint_b) {
    const Vec3 u = endpoint_a - vertex;
    const Vec3 v = endpoint_b - vertex;
    if (!Finite(u) || !Finite(v)) return kNaN;
    const double un = u.norm();
    const double vn = v.norm();
    if (un <= kNearZero || vn <= kNearZero) return kNaN;
    const double cosine = std::clamp(u.dot(v) / (un * vn), -1.0, 1.0);
    return std::acos(cosine);
}


Vec3 IdealCbPosition(const Vec3& n_pos,
                     const Vec3& ca_pos,
                     const Vec3& c_pos) {
    if (!Finite(n_pos) || !Finite(ca_pos) || !Finite(c_pos))
        return NanVec3();

    const Vec3 n = n_pos - ca_pos;
    const Vec3 c = c_pos - ca_pos;
    const double nn = n.norm();
    const double cn = c.norm();
    if (nn <= kNearZero || cn <= kNearZero) return NanVec3();

    const Vec3 cross_cn = c.cross(n);
    if (cross_cn.norm() <= kNearZero * nn * cn) return NanVec3();

    return ca_pos - 0.58273431 * cross_cn
                  + 0.56802827 * n
                  - 0.54067466 * c;
}


Vec3 CbDeviationVector(const Vec3& n_pos,
                       const Vec3& ca_pos,
                       const Vec3& c_pos,
                       const Vec3& observed_cb) {
    if (!Finite(observed_cb)) return NanVec3();
    const Vec3 ideal = IdealCbPosition(n_pos, ca_pos, c_pos);
    if (!Finite(ideal)) return NanVec3();
    return observed_cb - ideal;
}


double CbDeviation(const Vec3& n_pos,
                   const Vec3& ca_pos,
                   const Vec3& c_pos,
                   const Vec3& observed_cb) {
    const Vec3 residual = CbDeviationVector(n_pos, ca_pos, c_pos, observed_cb);
    return Finite(residual) ? residual.norm() : kNaN;
}


Measurements MeasureResidue(const ProteinConformation& conf,
                            std::size_t residue_index) {
    const Protein& protein = conf.ProteinRef();
    const Residue& residue = protein.ResidueAt(residue_index);
    const auto prev_idx = protein.BackbonePredecessor(residue_index);
    const auto next_idx = protein.BackboneSuccessor(residue_index);

    Measurements m{
        kNaN, kNaN, kNaN, kNaN, kNaN, kNaN, NanVec3()
    };

    if (HasAtom(residue.N) && HasAtom(residue.CA) && HasAtom(residue.C)) {
        const Vec3& n = conf.Positions().at(residue.N);
        const Vec3& ca = conf.Positions().at(residue.CA);
        const Vec3& c = conf.Positions().at(residue.C);
        m.tau_n_ca_c = BondAngleRadians(n, ca, c);

        if (HasAtom(residue.CB)) {
            const Vec3& cb = conf.Positions().at(residue.CB);
            m.angle_n_ca_cb = BondAngleRadians(n, ca, cb);
            m.angle_cb_ca_c = BondAngleRadians(cb, ca, c);
            m.cb_residual_vector = CbDeviationVector(n, ca, c, cb);
            if (m.cb_residual_vector.allFinite())
                m.cb_deviation = m.cb_residual_vector.norm();
        }
    } else {
        // Each C-beta valence angle remains independently defined when its
        // own three typed endpoints exist.
        if (HasAtom(residue.N) && HasAtom(residue.CA) && HasAtom(residue.CB)) {
            m.angle_n_ca_cb = BondAngleRadians(
                conf.PositionAt(residue.N), conf.PositionAt(residue.CA),
                conf.PositionAt(residue.CB));
        }
        if (HasAtom(residue.CB) && HasAtom(residue.CA) && HasAtom(residue.C)) {
            m.angle_cb_ca_c = BondAngleRadians(
                conf.PositionAt(residue.CB), conf.PositionAt(residue.CA),
                conf.PositionAt(residue.C));
        }
    }

    if (prev_idx && HasAtom(residue.N) && HasAtom(residue.CA)) {
        const Residue& prev = protein.ResidueAt(*prev_idx);
        m.angle_cprev_n_ca = BondAngleRadians(
            conf.PositionAt(prev.C), conf.PositionAt(residue.N),
            conf.PositionAt(residue.CA));
    }
    if (next_idx && HasAtom(residue.CA) && HasAtom(residue.C)) {
        const Residue& next = protein.ResidueAt(*next_idx);
        m.angle_ca_c_nnext = BondAngleRadians(
            conf.PositionAt(residue.CA), conf.PositionAt(residue.C),
            conf.PositionAt(next.N));
    }

    return m;
}

}  // namespace local_backbone_geometry


std::unique_ptr<LocalBackboneGeometryResult>
LocalBackboneGeometryResult::Compute(ProteinConformation& conf) {
    auto result = std::make_unique<LocalBackboneGeometryResult>();
    const std::size_t R = conf.ProteinRef().ResidueCount();

    result->tau_n_ca_c_.reserve(R);
    result->angle_n_ca_cb_.reserve(R);
    result->angle_cb_ca_c_.reserve(R);
    result->angle_cprev_n_ca_.reserve(R);
    result->angle_ca_c_nnext_.reserve(R);
    result->cb_deviation_.reserve(R);
    result->cb_residual_vector_.reserve(R);

    result->tau_n_ca_c_valid_.reserve(R);
    result->angle_n_ca_cb_valid_.reserve(R);
    result->angle_cb_ca_c_valid_.reserve(R);
    result->angle_cprev_n_ca_valid_.reserve(R);
    result->angle_ca_c_nnext_valid_.reserve(R);
    result->cb_deviation_valid_.reserve(R);
    result->cb_residual_vector_valid_.reserve(R);

    for (std::size_t ri = 0; ri < R; ++ri) {
        const auto m = local_backbone_geometry::MeasureResidue(conf, ri);
        result->tau_n_ca_c_.push_back(m.tau_n_ca_c);
        result->angle_n_ca_cb_.push_back(m.angle_n_ca_cb);
        result->angle_cb_ca_c_.push_back(m.angle_cb_ca_c);
        result->angle_cprev_n_ca_.push_back(m.angle_cprev_n_ca);
        result->angle_ca_c_nnext_.push_back(m.angle_ca_c_nnext);
        result->cb_deviation_.push_back(m.cb_deviation);
        result->cb_residual_vector_.push_back(m.cb_residual_vector);

        result->tau_n_ca_c_valid_.push_back(std::isfinite(m.tau_n_ca_c));
        result->angle_n_ca_cb_valid_.push_back(std::isfinite(m.angle_n_ca_cb));
        result->angle_cb_ca_c_valid_.push_back(std::isfinite(m.angle_cb_ca_c));
        result->angle_cprev_n_ca_valid_.push_back(std::isfinite(m.angle_cprev_n_ca));
        result->angle_ca_c_nnext_valid_.push_back(std::isfinite(m.angle_ca_c_nnext));
        result->cb_deviation_valid_.push_back(std::isfinite(m.cb_deviation));
        result->cb_residual_vector_valid_.push_back(
            m.cb_residual_vector.allFinite());
    }

    return result;
}


int LocalBackboneGeometryResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    (void)conf;

    auto write_scalar = [&](const char* stem,
                            const std::vector<double>& values) {
        NpyWriter::WriteFloat64(output_dir + "/" + stem + ".npy",
                                values.data(), values.size());
    };
    auto write_valid = [&](const char* stem,
                           const std::vector<std::uint8_t>& values) {
        NpyWriter::WriteUInt8(output_dir + "/" + stem + ".npy",
                              values.data(), values.size());
    };

    write_scalar("tau_N_CA_C", tau_n_ca_c_);
    write_valid("tau_N_CA_C_valid", tau_n_ca_c_valid_);
    write_scalar("angle_N_CA_CB", angle_n_ca_cb_);
    write_valid("angle_N_CA_CB_valid", angle_n_ca_cb_valid_);
    write_scalar("angle_CB_CA_C", angle_cb_ca_c_);
    write_valid("angle_CB_CA_C_valid", angle_cb_ca_c_valid_);
    write_scalar("angle_Cprev_N_CA", angle_cprev_n_ca_);
    write_valid("angle_Cprev_N_CA_valid", angle_cprev_n_ca_valid_);
    write_scalar("angle_CA_C_Nnext", angle_ca_c_nnext_);
    write_valid("angle_CA_C_Nnext_valid", angle_ca_c_nnext_valid_);
    write_scalar("cb_deviation", cb_deviation_);
    write_valid("cb_deviation_valid", cb_deviation_valid_);

    std::vector<double> residual(cb_residual_vector_.size() * 3);
    for (std::size_t i = 0; i < cb_residual_vector_.size(); ++i) {
        residual[i * 3 + 0] = cb_residual_vector_[i].x();
        residual[i * 3 + 1] = cb_residual_vector_[i].y();
        residual[i * 3 + 2] = cb_residual_vector_[i].z();
    }
    NpyWriter::WriteFloat64(output_dir + "/cb_residual_vector.npy",
                            residual.data(), cb_residual_vector_.size(), 3);
    write_valid("cb_residual_vector_valid", cb_residual_vector_valid_);

    return 14;
}

}  // namespace nmr
