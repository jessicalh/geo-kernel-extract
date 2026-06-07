#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "Protein.h"
#include "Ring.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"

#include <filesystem>
#include <cmath>
#include <cstdint>

namespace fs = std::filesystem;

namespace nmr {

namespace {

void PackMat3RowMajor(const Mat3& m, double* out) {
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            out[r*3 + c] = m(r, c);
        }
    }
}

}  // namespace


int ConformationResult::WriteAllFeatures(const ProteinConformation& conf,
                                          const std::string& output_dir) {
    OperationLog::Scope scope("ConformationResult::WriteAllFeatures",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " results=" + std::to_string(conf.AllResults().size()) +
        " dir=" + output_dir);

    fs::create_directories(output_dir);

    const Protein& protein = conf.ProteinRef();
    const size_t N = conf.AtomCount();
    int total = 0;

    // Identity arrays — these belong to no single result.
    {
        std::vector<double> pos(N * 3);
        std::vector<int32_t> elem(N), res_idx(N), res_type(N);
        for (size_t i = 0; i < N; ++i) {
            Vec3 p = conf.PositionAt(i);
            pos[i*3+0] = p.x(); pos[i*3+1] = p.y(); pos[i*3+2] = p.z();
            const Atom& a = protein.AtomAt(i);
            elem[i] = AtomicNumberForElement(a.element);
            res_idx[i] = static_cast<int32_t>(a.residue_index);
            res_type[i] = static_cast<int32_t>(protein.ResidueAt(a.residue_index).type);
        }
        NpyWriter::WriteFloat64(output_dir + "/pos.npy", pos.data(), N, 3);
        NpyWriter::WriteInt32(output_dir + "/element.npy", elem.data(), N);
        NpyWriter::WriteInt32(output_dir + "/residue_index.npy", res_idx.data(), N);
        NpyWriter::WriteInt32(output_dir + "/residue_type.npy", res_type.data(), N);
        total += 4;
    }

    // Per-ring contributions — sparse (atom, ring) pair array.
    // Shape (P, 58) where P = total evaluated (atom, ring) pairs.
    // Columns: [0-8] geometry, [9-17] BS G, [18-26] HM H (pure T2),
    //          [27-35] HM G (shielding), [36-44] quad, [45-53] chi,
    //          [54-55] dispersion, [56-57] azimuthal angle.
    {
        size_t P = 0;
        for (size_t i = 0; i < N; ++i)
            P += conf.AtomAt(i).ring_neighbours.size();

        if (P > 0) {
            const size_t C = 58;
            std::vector<double> data(P * C, 0.0);
            std::vector<double> bs_B(P * 3, 0.0);
            std::vector<double> bs_B_cyl(P * 3, 0.0);
            std::vector<double> hm_B(P * 3, 0.0);
            std::vector<double> ring_dir(P * 3, 0.0);
            std::vector<double> disp_tensor(P * 9, 0.0);
            std::vector<double> disp_spherical(P * 9, 0.0);
            std::vector<double> piquad_scalar(P, 0.0);
            size_t row = 0;
            for (size_t i = 0; i < N; ++i) {
                for (const auto& rn : conf.AtomAt(i).ring_neighbours) {
                    double* r = &data[row * C];
                    r[0]  = static_cast<double>(i);
                    r[1]  = static_cast<double>(rn.ring_index);
                    r[2]  = static_cast<double>(rn.ring_type);
                    r[3]  = rn.distance_to_center;
                    r[4]  = rn.rho;
                    r[5]  = rn.z;
                    r[6]  = rn.theta;

                    double cos_th = (rn.distance_to_center > 1e-12)
                        ? rn.z / rn.distance_to_center : 0.0;
                    double r3 = rn.distance_to_center * rn.distance_to_center
                              * rn.distance_to_center;
                    r[7]  = (r3 > 1e-30)
                        ? (3.0 * cos_th * cos_th - 1.0) / r3 : 0.0;
                    r[8]  = std::exp(-rn.distance_to_center / EXP_DECAY_LENGTH);

                    rn.G_spherical.PackFull9(r + 9);       // BS shielding kernel
                    rn.hm_H_spherical.PackFull9(r + 18);   // HM raw integral (pure T2)
                    rn.hm_G_spherical.PackFull9(r + 27);   // HM shielding kernel (T0+T1+T2)
                    rn.quad_spherical.PackFull9(r + 36);
                    rn.chi_spherical.PackFull9(r + 45);
                    r[54] = rn.disp_scalar;
                    r[55] = static_cast<double>(rn.disp_contacts);
                    r[56] = rn.cos_phi;
                    r[57] = rn.sin_phi;

                    bs_B[row*3+0] = rn.B_field.x();
                    bs_B[row*3+1] = rn.B_field.y();
                    bs_B[row*3+2] = rn.B_field.z();

                    bs_B_cyl[row*3+0] = rn.B_cylindrical.x();
                    bs_B_cyl[row*3+1] = rn.B_cylindrical.y();
                    bs_B_cyl[row*3+2] = rn.B_cylindrical.z();

                    hm_B[row*3+0] = rn.hm_B_field.x();
                    hm_B[row*3+1] = rn.hm_B_field.y();
                    hm_B[row*3+2] = rn.hm_B_field.z();

                    ring_dir[row*3+0] = rn.direction_to_center.x();
                    ring_dir[row*3+1] = rn.direction_to_center.y();
                    ring_dir[row*3+2] = rn.direction_to_center.z();

                    PackMat3RowMajor(rn.disp_tensor, &disp_tensor[row*9]);
                    rn.disp_spherical.PackFull9(&disp_spherical[row*9]);
                    piquad_scalar[row] = rn.quad_scalar;
                    row++;
                }
            }
            NpyWriter::WriteFloat64(output_dir + "/ring_contributions.npy",
                                    data.data(), P, C);
            NpyWriter::WriteFloat64(output_dir + "/bs_ring_B_field.npy",
                                    bs_B.data(), P, 3);
            NpyWriter::WriteFloat64(output_dir + "/bs_ring_B_cylindrical.npy",
                                    bs_B_cyl.data(), P, 3);
            NpyWriter::WriteFloat64(output_dir + "/hm_ring_B_field.npy",
                                    hm_B.data(), P, 3);
            NpyWriter::WriteFloat64(output_dir + "/ring_direction_to_center.npy",
                                    ring_dir.data(), P, 3);
            NpyWriter::WriteFloat64(output_dir + "/disp_per_ring_tensor.npy",
                                    disp_tensor.data(), P, 9);
            NpyWriter::WriteFloat64(output_dir + "/disp_per_ring_spherical.npy",
                                    disp_spherical.data(), P, 9);
            NpyWriter::WriteFloat64(output_dir + "/piquad_quad_scalar.npy",
                                    piquad_scalar.data(), P);
            total += 8;
        }
    }

    // Ring geometry reference table — shape (R, 10), one row per ring.
    {
        const size_t R = protein.RingCount();
        if (R > 0 && R <= conf.ring_geometries.size()) {
            const size_t G = 10;
            std::vector<double> data(R * G, 0.0);
            for (size_t ri = 0; ri < R; ++ri) {
                const Ring& ring = protein.RingAt(ri);
                const RingGeometry& geom = conf.ring_geometries[ri];
                double* d = &data[ri * G];
                d[0] = static_cast<double>(ri);
                d[1] = static_cast<double>(ring.type_index);
                d[2] = static_cast<double>(ring.parent_residue_index);
                d[3] = geom.center.x();
                d[4] = geom.center.y();
                d[5] = geom.center.z();
                d[6] = geom.normal.x();
                d[7] = geom.normal.y();
                d[8] = geom.normal.z();
                d[9] = geom.radius;
            }
            NpyWriter::WriteFloat64(output_dir + "/ring_geometry.npy",
                                    data.data(), R, G);
            total++;
        }
    }

    // Walk the conformation's accumulated results. Each one writes its own.
    for (const auto& [tid, result] : conf.AllResults()) {
        int n = result->WriteFeatures(conf, output_dir);
        if (n > 0) {
            OperationLog::Info(LogCalcOther, "WriteAllFeatures",
                result->Name() + " wrote " + std::to_string(n) + " arrays");
        }
        total += n;
    }

    OperationLog::Info(LogCalcOther, "WriteAllFeatures",
        "total: " + std::to_string(total) + " arrays");
    return total;
}

}  // namespace nmr
