#include "ExtractionSupport.h"

#include "SphericalBasis.h"

#include "../model/Conformation.h"
#include "../model/DftShielding.h"
#include "../model/QtProtein.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace h5reader::rediscover {

void FillIdentity(NeighborhoodRecord& rec, const RunData& run, std::size_t atomIdx,
                  std::size_t h5Row, const QString& stratum, const LocalFrame& frame) {
    const model::QtProtein& p = *run.protein;
    const model::QtAtom& a = p.atom(atomIdx);

    rec.atom_index = a.atomIndex;
    rec.residue_index = a.residueIndex;
    rec.element = static_cast<int>(a.element);
    rec.stratum = stratum;
    rec.atom_name = p.atomLabel(atomIdx, model::NamingConvention::Bmrb);

    if (a.residueIndex >= 0 && static_cast<std::size_t>(a.residueIndex) < p.residueCount()) {
        const model::QtResidue& r = p.residue(static_cast<std::size_t>(a.residueIndex));
        rec.residue_number = r.address.residueNumber;
        rec.amino_acid = static_cast<int>(r.aminoAcid);
    }

    rec.h5_row = static_cast<int32_t>(h5Row);
    rec.original_index = static_cast<int32_t>(run.frameMap.originalIndex(h5Row));
    rec.time_ps = run.trajectory()->timePicoseconds(h5Row);

    rec.frame_z = frame.z;
    rec.frame_x = frame.x;
    rec.frame_y = frame.y;
    rec.frame_variant = static_cast<int>(frame.variant);
    rec.frame_valid = frame.is_valid;
    // frame_anchor_atom_index is set by the extraction (the aromatic-H ring
    // frame records its typed CG/CD2 anchor); left at -1 here.
}

DftTarget BuildTarget(const RunData& run, std::size_t atomIdx, std::size_t originalIndex,
                      const LocalFrame& frame) {
    DftTarget t;
    const model::DftAtomShielding* s = run.dft.AtomShielding(atomIdx, originalIndex);
    if (!s) return t;  // present stays false — honest gap

    t.present = true;
    t.total_raw = s->total_raw;
    t.dia_raw = s->dia_raw;
    t.para_raw = s->para_raw;

    // Library-basis decomposition of the RAW tensors so DFT-T2 lands in the
    // same component order as the H5 kernel-T2. T0 == σ_iso (rotation-safe).
    t.total_decomp = DecomposeLibrary(s->total_raw);
    t.dia_decomp = DecomposeLibrary(s->dia_raw);
    t.para_decomp = DecomposeLibrary(s->para_raw);

    // Total tensor in the atom's local frame (only meaningful when valid).
    if (frame.is_valid) {
        t.total_local = frame.TensorToLocal(s->total_raw);
        t.local_frame_valid = true;
    }
    return t;
}

std::vector<FeatureColumn> IdentityColumns() {
    return {
        {QStringLiteral("atom_index"), {}},
        {QStringLiteral("residue_index"), {}},
        {QStringLiteral("residue_number"), {}},
        {QStringLiteral("amino_acid_ord"), {}},
        {QStringLiteral("element_ord"), {}},
        {QStringLiteral("atom_name"), {}},
        {QStringLiteral("stratum"), {}},
        {QStringLiteral("h5_row"), {}},
        {QStringLiteral("original_index"), {}},
        {QStringLiteral("time_ps"), QStringLiteral("ps")},
        {QStringLiteral("frame_z_x"), {}}, {QStringLiteral("frame_z_y"), {}}, {QStringLiteral("frame_z_z"), {}},
        {QStringLiteral("frame_x_x"), {}}, {QStringLiteral("frame_x_y"), {}}, {QStringLiteral("frame_x_z"), {}},
        {QStringLiteral("frame_y_x"), {}}, {QStringLiteral("frame_y_y"), {}}, {QStringLiteral("frame_y_z"), {}},
        {QStringLiteral("frame_variant"), {}},
        {QStringLiteral("frame_valid"), {}},
        {QStringLiteral("frame_anchor_atom_index"), {}},
    };
}

std::vector<FeatureColumn> BareKernelColumns() {
    std::vector<FeatureColumn> c = {
        {QStringLiteral("bare_kernel_present"), {}},
        {QStringLiteral("bare_T0"), QStringLiteral("ppm")},
    };
    for (int i = 0; i < 3; ++i) c.push_back({QStringLiteral("bare_T1_%1").arg(i), QStringLiteral("ppm")});
    for (int i = 0; i < 5; ++i) c.push_back({QStringLiteral("bare_T2_%1").arg(i), QStringLiteral("ppm")});
    return c;
}

DftFrameAlignment CheckDftFrameAlignment(const RunData& run) {
    DftFrameAlignment out;
    const model::QtProtein& p = *run.protein;
    const model::Conformation& conf = *run.conformation;
    const std::size_t nAtoms = p.atomCount();
    constexpr double kRad2Deg = 57.29577951308232;

    double sumAngle = 0.0, sumRmsd = 0.0;
    for (std::size_t row : run.frameMap.dftRows()) {
        const std::size_t orig = run.frameMap.originalIndex(row);
        const model::DftShieldingFrame* fr = run.dft.Frame(orig);
        if (!fr) continue;

        // Matched atoms: ORCA-input position (the DFT-tensor frame) vs H5 position.
        std::vector<Vec3> P, Q;  // P = ORCA, Q = H5
        const std::size_t n = std::min(nAtoms, fr->atoms.size());
        P.reserve(n);
        Q.reserve(n);
        for (std::size_t i = 0; i < n; ++i) {
            const Vec3& pc = fr->atoms[i].orca_coord;
            if (pc.squaredNorm() == 0.0) continue;  // unfilled — skip
            P.push_back(pc);
            Q.push_back(conf.atomPosition(row, i));
        }
        if (P.size() < 3) continue;

        // Kabsch: optimal rotation P→Q (translation removed by centring, so PBC
        // wrapping shifts don't matter; only a real rotation survives).
        Vec3 cP = Vec3::Zero(), cQ = Vec3::Zero();
        for (const Vec3& v : P) cP += v;
        for (const Vec3& v : Q) cQ += v;
        cP /= static_cast<double>(P.size());
        cQ /= static_cast<double>(Q.size());
        Mat3 H = Mat3::Zero();
        for (std::size_t i = 0; i < P.size(); ++i) H += (P[i] - cP) * (Q[i] - cQ).transpose();
        Eigen::JacobiSVD<Mat3> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
        const Mat3 U = svd.matrixU(), V = svd.matrixV();
        const double det = (V * U.transpose()).determinant();
        Mat3 D = Mat3::Identity();
        D(2, 2) = det < 0.0 ? -1.0 : 1.0;
        const Mat3 R = V * D * U.transpose();  // rotates P onto Q

        const double cosang = std::clamp((R.trace() - 1.0) / 2.0, -1.0, 1.0);
        const double angle = std::acos(cosang) * kRad2Deg;
        double se = 0.0;
        for (std::size_t i = 0; i < P.size(); ++i)
            se += (R * (P[i] - cP) - (Q[i] - cQ)).squaredNorm();
        const double rmsd = std::sqrt(se / static_cast<double>(P.size()));

        ++out.n_frames;
        sumAngle += angle;
        sumRmsd += rmsd;
        out.max_angle_deg = std::max(out.max_angle_deg, angle);
        out.max_rmsd_A = std::max(out.max_rmsd_A, rmsd);
        out.n_atoms_used = static_cast<int>(P.size());
    }
    if (out.n_frames > 0) {
        out.mean_angle_deg = sumAngle / out.n_frames;
        out.mean_rmsd_A = sumRmsd / out.n_frames;
    }
    return out;
}

std::vector<FeatureColumn> TargetColumns() {
    std::vector<FeatureColumn> c;
    c.push_back({QStringLiteral("dft_present"), {}});
    auto mat = [&](const QString& tag) {
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                c.push_back({QStringLiteral("dft_%1_%2%3").arg(tag).arg(i).arg(j), QStringLiteral("ppm")});
    };
    mat(QStringLiteral("total_raw"));
    mat(QStringLiteral("dia_raw"));
    mat(QStringLiteral("para_raw"));
    c.push_back({QStringLiteral("dft_sigma_iso"), QStringLiteral("ppm")});  // total T0
    for (int i = 0; i < 3; ++i) c.push_back({QStringLiteral("dft_total_T1_%1").arg(i), QStringLiteral("ppm")});
    for (int i = 0; i < 5; ++i) c.push_back({QStringLiteral("dft_total_T2_%1").arg(i), QStringLiteral("ppm")});
    c.push_back({QStringLiteral("dft_dia_iso"), QStringLiteral("ppm")});
    c.push_back({QStringLiteral("dft_para_iso"), QStringLiteral("ppm")});
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            c.push_back({QStringLiteral("dft_total_local_%1%2").arg(i).arg(j), QStringLiteral("ppm")});
    c.push_back({QStringLiteral("dft_local_frame_valid"), {}});
    return c;
}

}  // namespace h5reader::rediscover
