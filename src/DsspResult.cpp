#include "DsspResult.h"
#include "Protein.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <cif++.hpp>
#include <cif++/pdb/pdb2cif.hpp>
#include <dssp.hpp>

#include <fstream>
#include <filesystem>
#include <map>
#include <cstdio>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>

namespace fs = std::filesystem;

namespace nmr {

namespace {

constexpr std::size_t kCanonicalTorsionCount = 6;
constexpr double kTorsionBondNearZero = 1e-12;
constexpr double kTorsionNormalNearZero = 1e-10;

double GuardedIupacDihedral(const Vec3& p0, const Vec3& p1,
                            const Vec3& p2, const Vec3& p3) {
    if (!p0.allFinite() || !p1.allFinite() ||
        !p2.allFinite() || !p3.allFinite()) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const Vec3 b1 = p1 - p0;
    const Vec3 b2 = p2 - p1;
    const Vec3 b3 = p3 - p2;
    const Vec3 n1 = b1.cross(b2);
    const Vec3 n2 = b2.cross(b3);
    const double b2_norm = b2.norm();
    const double n1_norm = n1.norm();
    const double n2_norm = n2.norm();
    if (b2_norm <= kTorsionBondNearZero ||
        n1_norm <= kTorsionNormalNearZero ||
        n2_norm <= kTorsionNormalNearZero) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const Vec3 m1 = n1.cross(b2 / b2_norm);
    return std::atan2(m1.dot(n2), n1.dot(n2));
}

}  // namespace

// ============================================================================
// Write a minimal PDB from our Protein + Conformation for cif++ to read.
//
// PDB LOADING BOUNDARY: string-based PDB format output. This is the tool
// boundary where typed objects are serialised to the string format that
// libdssp expects. After DSSP runs, results are mapped back to typed
// residue indices.
// ============================================================================

static std::string WriteTempPdb(const Protein& protein,
                                 const ProteinConformation& conf) {
    auto tmp_path = fs::temp_directory_path() /
        ("dssp_" + std::to_string(reinterpret_cast<uintptr_t>(&conf)) + ".pdb");

    std::ofstream out(tmp_path);
    if (!out.is_open()) return "";

    int serial = 1;
    for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
        const auto& res = protein.ResidueAt(ri);
        for (size_t ai : res.atom_indices) {
            const auto& atom = protein.AtomAt(ai);
            Vec3 pos = conf.PositionAt(ai);

            std::string elem_str;
            switch (atom.element) {
                case Element::H: elem_str = " H"; break;
                case Element::C: elem_str = " C"; break;
                case Element::N: elem_str = " N"; break;
                case Element::O: elem_str = " O"; break;
                case Element::S: elem_str = " S"; break;
                default:         elem_str = " X"; break;
            }

            // PDB fixed-width format: atom name in columns 13-16
            char atom_field[5];
            if (atom.pdb_atom_name.size() <= 3)
                snprintf(atom_field, sizeof(atom_field), " %-3s",
                         atom.pdb_atom_name.c_str());
            else
                snprintf(atom_field, sizeof(atom_field), "%-4s",
                         atom.pdb_atom_name.c_str());

            std::string res_name = ThreeLetterCodeForAminoAcid(res.type);
            char chain = res.chain_id.empty() ? 'A' : res.chain_id[0];

            char line[82];
            snprintf(line, sizeof(line),
                "ATOM  %5d %4s %3s %c%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n",
                serial++, atom_field, res_name.c_str(), chain, res.sequence_number,
                pos.x(), pos.y(), pos.z(), elem_str.c_str());
            out << line;
        }
    }
    out << "END\n";
    out.close();
    return tmp_path.string();
}


std::unique_ptr<DsspResult> DsspResult::Compute(ProteinConformation& conf) {
    const Protein& protein = conf.ProteinRef();
    auto result = std::make_unique<DsspResult>();
    result->residues_.resize(protein.ResidueCount());

    // Write temp PDB for cif++ / DSSP
    std::string tmp_path = WriteTempPdb(protein, conf);
    if (tmp_path.empty()) {
        fprintf(stderr, "DsspResult::Compute: failed to write temp PDB.\n");
        return nullptr;
    }

    // RAII guard for temp file cleanup
    struct TempGuard {
        std::string path;
        ~TempGuard() { std::error_code ec; fs::remove(path, ec); }
    } guard{tmp_path};

    // TOOL BOUNDARY: cif++ parser + libdssp may throw on unusual geometry.
    try {
        cif::file cif_file;
        {
            std::ifstream pdb_in(tmp_path);
            cif::pdb::ReadPDBFile(pdb_in, cif_file);
        }

        auto& db = cif_file.front();
        cif::mm::structure structure(db, 1, {});

        // Constructor arguments: min_poly_proline_stretch_length = 3,
        // calculateSurfaceAccessibility = true.
        dssp dssp_calc(structure, 3, true);

        // Build lookup: (chain, seq_num) -> residue index
        std::map<std::pair<std::string, int>, size_t> res_lookup;
        for (size_t ri = 0; ri < protein.ResidueCount(); ++ri) {
            const auto& res = protein.ResidueAt(ri);
            std::string chain = res.chain_id.empty() ? "A" : res.chain_id;
            res_lookup[{chain, res.sequence_number}] = ri;
        }

        // Map DSSP results to our residues
        // Phi/Psi from DSSP are in degrees; convert to radians
        for (auto& dssp_res : dssp_calc) {
            std::string chain = dssp_res.pdb_strand_id();
            int seq_id = dssp_res.pdb_seq_num();

            auto it = res_lookup.find({chain, seq_id});
            if (it == res_lookup.end()) continue;

            auto& dr = result->residues_[it->second];
            dr.observed = true;  // Codex review 2026-05-19: distinguish
                                 // mapped from default-coil for honest
                                 // downstream SS classification.
            char ss = static_cast<char>(dssp_res.type());
            dr.secondary_structure = (ss == ' ') ? 'C' : ss;

            // DSSP phi/psi are in degrees; convert to radians
            dr.phi = dssp_res.phi() * PI / 180.0;
            dr.psi = dssp_res.psi() * PI / 180.0;
            dr.sasa = dssp_res.accessibility();

            // H-bond partners
            for (int bi = 0; bi < 2; ++bi) {
                auto [acc_res, acc_energy] = dssp_res.acceptor(bi);
                if (acc_res) {
                    auto acc_it = res_lookup.find(
                        {acc_res.pdb_strand_id(), acc_res.pdb_seq_num()});
                    if (acc_it != res_lookup.end()) {
                        dr.acceptors[bi].residue_index = acc_it->second;
                        dr.acceptors[bi].energy = acc_energy;
                    }
                }

                auto [don_res, don_energy] = dssp_res.donor(bi);
                if (don_res) {
                    auto don_it = res_lookup.find(
                        {don_res.pdb_strand_id(), don_res.pdb_seq_num()});
                    if (don_it != res_lookup.end()) {
                        dr.donors[bi].residue_index = don_it->second;
                        dr.donors[bi].energy = don_energy;
                    }
                }
            }
        }
    } catch (const std::exception& e) {
        fprintf(stderr, "DsspResult::Compute: DSSP failed: %s\n", e.what());
        return nullptr;
    }

    result->PopulateCanonicalTorsions(&conf);
    return result;
}

std::unique_ptr<DsspResult> DsspResult::CreateForTesting(
        std::vector<DsspResidue> residues) {
    auto result = std::make_unique<DsspResult>();
    result->residues_ = std::move(residues);
    result->PopulateCanonicalTorsions(nullptr);
    return result;
}


void DsspResult::PopulateCanonicalTorsions(
        const ProteinConformation* conf) {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const std::size_t R = residues_.size();
    torsion_angle_.assign(R * kCanonicalTorsionCount, nan);
    torsion_sin_.assign(R * kCanonicalTorsionCount, nan);
    torsion_cos_.assign(R * kCanonicalTorsionCount, nan);
    torsion_valid_.assign(R * kCanonicalTorsionCount, 0u);

    auto store = [&](std::size_t ri, std::size_t column, double angle) {
        if (!std::isfinite(angle)) return;
        const std::size_t offset = ri * kCanonicalTorsionCount + column;
        torsion_angle_[offset] = angle;
        torsion_sin_[offset] = std::sin(angle);
        torsion_cos_[offset] = std::cos(angle);
        torsion_valid_[offset] = 1u;
    };

    for (std::size_t ri = 0; ri < R; ++ri) {
        const DsspResidue& dr = residues_[ri];
        if (dr.observed) {
            // libdssp reports the negated-IUPAC convention and uses 360
            // degrees (2*pi after conversion) for an undefined terminus.
            if (std::isfinite(dr.phi) && std::abs(dr.phi) <= PI + 1e-6)
                store(ri, 0, -dr.phi);
            if (std::isfinite(dr.psi) && std::abs(dr.psi) <= PI + 1e-6)
                store(ri, 1, -dr.psi);
        }

        if (!conf || ri >= conf->ProteinRef().ResidueCount()) continue;
        const Residue& residue = conf->ProteinRef().ResidueAt(ri);
        for (std::size_t chi = 0; chi < 4; ++chi) {
            if (!residue.chi[chi].Valid()) continue;
            const auto& a = residue.chi[chi].a;
            const double angle = GuardedIupacDihedral(
                conf->PositionAt(a[0]), conf->PositionAt(a[1]),
                conf->PositionAt(a[2]), conf->PositionAt(a[3]));
            store(ri, 2 + chi, angle);
        }
    }
}


char DsspResult::SecondaryStructure(size_t residue_index) const {
    if (residue_index >= residues_.size()) return 'C';
    return residues_[residue_index].secondary_structure;
}

double DsspResult::Phi(size_t residue_index) const {
    if (residue_index >= residues_.size()) return 0.0;
    return residues_[residue_index].phi;
}

double DsspResult::Psi(size_t residue_index) const {
    if (residue_index >= residues_.size()) return 0.0;
    return residues_[residue_index].psi;
}

double DsspResult::SASA(size_t residue_index) const {
    if (residue_index >= residues_.size()) return 0.0;
    return residues_[residue_index].sasa;
}


// ============================================================================
// WriteFeatures: 11 NPY files
//
// All per-atom, broadcast from per-residue via Protein atom→residue mapping.
//
// 1. dssp_observed.npy (N,) — residue was observed by DSSP
// 2. dssp_backbone.npy (N, 5) — -phi, -psi, sasa, ss_helix, ss_sheet
// 3. dssp_ss8.npy (N, 8) — full 8-class SS one-hot (H/G/I/E/B/T/S/C)
// 4. dssp_ppii.npy (N,) — explicit PPII flag (1/0/-1 sentinel)
// 5. dssp_hbond_energy.npy (N, 4) — H-bond energies (acc0/acc1/don0/don1)
// 6. dssp_chi.npy (N, 12) — chi1-4 cos/sin/exists (4 angles × 3 cols)
// 7-10. dssp_torsion_{angle,sin,cos,valid}.npy (R, 6)
// 11. dssp_hbond_partner_residue_index.npy (N, 4)
// ============================================================================

int DsspResult::WriteFeatures(const ProteinConformation& conf,
                               const std::string& output_dir) const {
    const Protein& protein = conf.ProteinRef();
    const size_t N = conf.AtomCount();
    const double kNaN = std::numeric_limits<double>::quiet_NaN();
    const std::size_t expected_torsion_size =
        protein.ResidueCount() * kCanonicalTorsionCount;
    if (torsion_angle_.size() != expected_torsion_size ||
        torsion_sin_.size() != expected_torsion_size ||
        torsion_cos_.size() != expected_torsion_size ||
        torsion_valid_.size() != expected_torsion_size) {
        const std::string message =
            "DsspResult::WriteFeatures canonical torsion buffer size mismatch";
        OperationLog::Error("DsspResult::WriteFeatures", message);
        throw std::runtime_error(message);
    }

    std::vector<double> data(N * 5, 0.0);
    std::vector<int8_t> observed_mask(N, 0);

    for (size_t i = 0; i < N; ++i) {
        size_t ri = protein.AtomAt(i).residue_index;
        const bool observed = (ri < residues_.size() && residues_[ri].observed);
        observed_mask[i] = observed ? 1 : 0;
        if (observed) {
            const auto& dr = residues_[ri];
            data[i*5 + 0] = -dr.phi;
            data[i*5 + 1] = -dr.psi;
            data[i*5 + 2] = dr.sasa;
            char ss = dr.secondary_structure;
            data[i*5 + 3] = (ss == 'H' || ss == 'G' || ss == 'I') ? 1.0 : 0.0;
            data[i*5 + 4] = (ss == 'E' || ss == 'B') ? 1.0 : 0.0;
        } else {
            data[i*5 + 0] = kNaN;
            data[i*5 + 1] = kNaN;
            data[i*5 + 2] = kNaN;
        }
    }

    NpyWriter::WriteInt8(output_dir + "/dssp_observed.npy", observed_mask.data(), N);
    NpyWriter::WriteFloat64(output_dir + "/dssp_backbone.npy", data.data(), N, 5);
    int files_written = 2;

    // dssp_ss8.npy — (N, 8) float64, full 8-class one-hot
    // Column order: H(alpha), G(3_10), I(pi), E(strand), B(bridge), T(turn), S(bend), C(coil)
    {
        // Map: H=0, G=1, I=2, E=3, B=4, T=5, S=6, C=7
        auto ss_col = [](char ss) -> int {
            switch (ss) {
                case 'H': return 0;
                case 'G': return 1;
                case 'I': return 2;
                case 'E': return 3;
                case 'B': return 4;
                case 'T': return 5;
                case 'S': return 6;
                default:  return 7; // C / coil
            }
        };

        std::vector<double> ss8(N * 8, 0.0);
        for (size_t i = 0; i < N; ++i) {
            size_t ri = protein.AtomAt(i).residue_index;
            if (ri < residues_.size() && residues_[ri].observed)
                ss8[i * 8 + ss_col(residues_[ri].secondary_structure)] = 1.0;
        }
        NpyWriter::WriteFloat64(output_dir + "/dssp_ss8.npy", ss8.data(), N, 8);
        files_written++;
    }

    // dssp_ppii.npy — (N,) int8. 1=PPII, 0=observed non-PPII,
    // -1=no DSSP observation for the parent residue.
    {
        std::vector<int8_t> ppii(N, -1);
        for (size_t i = 0; i < N; ++i) {
            size_t ri = protein.AtomAt(i).residue_index;
            if (ri < residues_.size() && residues_[ri].observed) {
                ppii[i] = DsspIsPpii(residues_[ri].secondary_structure)
                    ? 1 : 0;
            }
        }
        NpyWriter::WriteInt8(output_dir + "/dssp_ppii.npy", ppii.data(), N);
        files_written++;
    }

    // dssp_hbond_energy.npy — (N, 4) float64, per-atom (broadcast from residue)
    // Columns: acceptor0_energy, acceptor1_energy, donor0_energy, donor1_energy
    // Units: kcal/mol. 0.0 if no H-bond partner.
    {
        std::vector<double> hb(N * 4, kNaN);
        for (size_t i = 0; i < N; ++i) {
            size_t ri = protein.AtomAt(i).residue_index;
            if (ri < residues_.size() && residues_[ri].observed) {
                const auto& dr = residues_[ri];
                hb[i * 4 + 0] = dr.acceptors[0].energy;
                hb[i * 4 + 1] = dr.acceptors[1].energy;
                hb[i * 4 + 2] = dr.donors[0].energy;
                hb[i * 4 + 3] = dr.donors[1].energy;
            }
        }
        NpyWriter::WriteFloat64(output_dir + "/dssp_hbond_energy.npy", hb.data(), N, 4);
        files_written++;
    }

    // dssp_chi.npy — (N, 12) float64, per-atom (broadcast from residue)
    // Per chi angle (1-4): cos, sin, exists (1.0/0.0)
    // Chi angles computed from Residue::chi[k] atom indices + conformation positions.
    // Missing or degenerate chi angles keep cos=0, sin=0, exists=0.
    {
        std::vector<double> chi_data(N * 12, 0.0);
        for (size_t i = 0; i < N; ++i) {
            size_t ri = protein.AtomAt(i).residue_index;
            if (ri >= protein.ResidueCount()) continue;
            const auto& res = protein.ResidueAt(ri);

            for (int k = 0; k < 4; ++k) {
                int base = i * 12 + k * 3;
                if (res.chi[k].Valid()) {
                    Vec3 p0 = conf.PositionAt(res.chi[k].a[0]);
                    Vec3 p1 = conf.PositionAt(res.chi[k].a[1]);
                    Vec3 p2 = conf.PositionAt(res.chi[k].a[2]);
                    Vec3 p3 = conf.PositionAt(res.chi[k].a[3]);

                    // Dihedral angle via atan2
                    Vec3 b1 = p1 - p0;
                    Vec3 b2 = p2 - p1;
                    Vec3 b3 = p3 - p2;
                    Vec3 n1 = b1.cross(b2);
                    Vec3 n2 = b2.cross(b3);
                    double n1_norm = n1.norm();
                    double n2_norm = n2.norm();
                    if (n1_norm > 1e-10 && n2_norm > 1e-10) {
                        n1 /= n1_norm;
                        n2 /= n2_norm;
                        double cos_angle = n1.dot(n2);
                        // clamp for numerical safety
                        cos_angle = std::max(-1.0, std::min(1.0, cos_angle));
                        Vec3 m1 = n1.cross(b2.normalized());
                        double sin_angle = m1.dot(n2);
                        chi_data[base + 0] = cos_angle;
                        chi_data[base + 1] = sin_angle;
                        chi_data[base + 2] = 1.0; // exists
                    }
                }
                // else: cos=0, sin=0, exists=0 (default)
            }
        }
        NpyWriter::WriteFloat64(output_dir + "/dssp_chi.npy", chi_data.data(), N, 12);
        files_written++;
    }

    // Canonical residue-axis torsions.  Columns are
    // phi, psi, chi1, chi2, chi3, chi4.  Every undefined numeric entry is
    // NaN and its uint8 validity entry is zero.
    {
        const std::size_t R = protein.ResidueCount();
        NpyWriter::WriteFloat64(output_dir + "/dssp_torsion_angle.npy",
                                torsion_angle_.data(), R,
                                kCanonicalTorsionCount);
        NpyWriter::WriteFloat64(output_dir + "/dssp_torsion_sin.npy",
                                torsion_sin_.data(), R,
                                kCanonicalTorsionCount);
        NpyWriter::WriteFloat64(output_dir + "/dssp_torsion_cos.npy",
                                torsion_cos_.data(), R,
                                kCanonicalTorsionCount);
        NpyWriter::WriteUInt8(output_dir + "/dssp_torsion_valid.npy",
                              torsion_valid_.data(), R,
                              kCanonicalTorsionCount);
        files_written += 4;
    }

    // Partner rows are atom-broadcast exactly like dssp_hbond_energy.npy;
    // the four columns retain the same acceptor0/acceptor1/donor0/donor1
    // slot order.  -1 is the explicit no-partner sentinel.
    {
        std::vector<std::int32_t> partner(N * 4, -1);
        for (std::size_t i = 0; i < N; ++i) {
            const std::size_t ri = protein.AtomAt(i).residue_index;
            if (ri >= residues_.size() || !residues_[ri].observed) continue;
            const DsspResidue& dr = residues_[ri];
            const std::size_t ids[4] = {
                dr.acceptors[0].residue_index,
                dr.acceptors[1].residue_index,
                dr.donors[0].residue_index,
                dr.donors[1].residue_index,
            };
            for (std::size_t slot = 0; slot < 4; ++slot) {
                if (ids[slot] != SIZE_MAX && ids[slot] < protein.ResidueCount())
                    partner[i * 4 + slot] = static_cast<std::int32_t>(ids[slot]);
            }
        }
        NpyWriter::WriteInt32(
            output_dir + "/dssp_hbond_partner_residue_index.npy",
            partner.data(), N, 4);
        files_written++;
    }

    return files_written;
}

}  // namespace nmr
