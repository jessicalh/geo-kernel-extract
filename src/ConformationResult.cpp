#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "Protein.h"
#include "Atom.h"
#include "AtomTopology.h"
#include "Residue.h"
#include "AminoAcidType.h"
#include "Ring.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <filesystem>
#include <fstream>
#include <cmath>
#include <cstdint>

namespace fs = std::filesystem;

namespace nmr {

// Pack SphericalTensor as [T0, T1[3], T2[5]] = 9 doubles.
static void PackST(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

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
    //
    // Per-atom identity data exported here covers what post-hoc analysis
    // needs to classify results by topology rules (Markley 1998 IUPAC):
    // residue context, locant + prochirality, ring participation, NMR
    // role, partner/parent linkage, and the IUPAC atom-name strings.
    // The SDK catalog (python/nmr_extract/_catalog.py) registers each
    // file; downstream calibration and analysis code reads typed columns
    // rather than reconstructing identity from name-string heuristics.
    //
    // Required tier (always present after FinalizeConstruction + a normal
    // RunOptions pipeline including EnrichmentResult): the calibration
    // and trajectory job specs both pull EnrichmentResult, so atom_role
    // and atom_flags are part of the identity contract — not optional.
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

    // -------------------------------------------------------------------
    // Atom topology — IUPAC symbolic identity (Markley 1998).
    //
    // 14 columns per atom:
    //   0  locant                Locant enum (Backbone/Alpha/.../Special)
    //   1  branch_index          BranchIndex enum (None/One/Two/Three)
    //   2  diastereotopic_index  DiastereotopicIndex (None/Two/Three)
    //   3  prochiral_stereo      ProchiralStereo (NotProchiral/Equivalent/ProR/ProS)
    //   4  planar_stereo         PlanarStereo (None/Cis/Trans)
    //   5  pseudoatom_class      PseudoatomClass (None/MB/MG/.../QR)
    //   6  ring_position         RingPosition (NotInRing/Substituent/Member/Junction)
    //   7  polar_h_kind          PolarHKind (None/RingNH/Hydroxyl/CarboxylAcid/Sulfanyl)
    //   8  topology_stamped      1 if AtomTopology was populated, 0 if the
    //                            (residue_type, atom_name) lookup missed
    //   9  ring_count            number of aromatic rings the atom belongs to
    //                            (0/1/3 — Trp Cδ2/Cε2 are in three rings)
    //   10..13  chi_position[0..3]  position 0..3 in chi(i+1), -1 if not
    //                              participating in that chi angle
    // -------------------------------------------------------------------
    {
        const size_t TOPO_COLS = 14;
        std::vector<int32_t> data(N * TOPO_COLS, 0);
        for (size_t i = 0; i < N; ++i) {
            const Atom& a = protein.AtomAt(i);
            const AtomTopology& t = a.topology;
            int32_t* r = &data[i * TOPO_COLS];
            r[0]  = static_cast<int32_t>(t.locant);
            r[1]  = static_cast<int32_t>(t.branch_index);
            r[2]  = static_cast<int32_t>(t.diastereotopic_index);
            r[3]  = static_cast<int32_t>(t.prochiral_stereo);
            r[4]  = static_cast<int32_t>(t.planar_stereo);
            r[5]  = static_cast<int32_t>(t.pseudoatom_class);
            r[6]  = static_cast<int32_t>(t.ring_position);
            r[7]  = static_cast<int32_t>(t.polar_h_kind);
            r[8]  = t.stamped ? 1 : 0;
            r[9]  = static_cast<int32_t>(a.ring_indices.size());
            for (int c = 0; c < 4; ++c)
                r[10 + c] = static_cast<int32_t>(t.chi_position[c]);
        }
        NpyWriter::WriteInt32(output_dir + "/atom_topology.npy",
                               data.data(), N * TOPO_COLS);
        total++;
    }

    // -------------------------------------------------------------------
    // Atom relationships — partner/parent linkage.
    //
    // 2 columns per atom (int64 because SIZE_MAX → -1; atom indices fit
    // comfortably in int64):
    //   0  partner_atom_index    diastereotopic partner (HB2↔HB3 etc.);
    //                            -1 if not at a 2/3-numbered position
    //   1  parent_atom_index     for hydrogens, the bonded heavy atom;
    //                            -1 if not hydrogen
    // -------------------------------------------------------------------
    {
        std::vector<int64_t> data(N * 2, -1);
        for (size_t i = 0; i < N; ++i) {
            const Atom& a = protein.AtomAt(i);
            data[i*2 + 0] = (a.partner_atom_index == SIZE_MAX)
                ? -1 : static_cast<int64_t>(a.partner_atom_index);
            data[i*2 + 1] = (a.parent_atom_index == SIZE_MAX)
                ? -1 : static_cast<int64_t>(a.parent_atom_index);
        }
        // NpyWriter has no Int64; pack as two int32 lanes per cell.
        // Atom indices in our use case fit in int32, so emit as int32
        // with -1 sentinel and document the type in the catalog.
        std::vector<int32_t> data32(N * 2);
        for (size_t i = 0; i < N * 2; ++i)
            data32[i] = static_cast<int32_t>(data[i]);
        NpyWriter::WriteInt32(output_dir + "/atom_relationships.npy",
                               data32.data(), N * 2);
        total++;
    }

    // -------------------------------------------------------------------
    // Atom role and classification flags — from EnrichmentResult.
    //
    // atom_role.npy        (N,)    int32  AtomRole enum
    // atom_flags.npy       (N, 10) int32  packed flags:
    //   0  hybridisation         Hybridisation enum (Unassigned/sp/sp2/sp3/...)
    //   1  is_backbone           role in backbone set
    //   2  is_amide_H            role == AmideH
    //   3  is_alpha_H             role == AlphaH
    //   4  is_methyl             role == MethylH
    //   5  is_aromatic_H         role == AromaticH
    //   6  is_on_aromatic_residue residue.IsAromatic()
    //   7  is_hbond_donor        H bonded to N or O
    //   8  is_hbond_acceptor     N or O with lone pair
    //   9  parent_is_sp2          for H, parent atom hybridisation == sp2
    //
    // EnrichmentResult is required by every job spec we run (calibration
    // pipeline, trajectory pipeline). If it has not run, fields are zero
    // — that is itself diagnostic information that downstream callers
    // can detect.
    // -------------------------------------------------------------------
    {
        std::vector<int32_t> roles(N);
        const size_t FLAG_COLS = 10;
        std::vector<int32_t> flags(N * FLAG_COLS, 0);
        for (size_t i = 0; i < N; ++i) {
            const ConformationAtom& ca = conf.AtomAt(i);
            roles[i] = static_cast<int32_t>(ca.role);
            int32_t* f = &flags[i * FLAG_COLS];
            f[0] = static_cast<int32_t>(ca.hybridisation);
            f[1] = ca.is_backbone ? 1 : 0;
            f[2] = ca.is_amide_H ? 1 : 0;
            f[3] = ca.is_alpha_H ? 1 : 0;
            f[4] = ca.is_methyl ? 1 : 0;
            f[5] = ca.is_aromatic_H ? 1 : 0;
            f[6] = ca.is_on_aromatic_residue ? 1 : 0;
            f[7] = ca.is_hbond_donor ? 1 : 0;
            f[8] = ca.is_hbond_acceptor ? 1 : 0;
            f[9] = ca.parent_is_sp2 ? 1 : 0;
        }
        NpyWriter::WriteInt32(output_dir + "/atom_role.npy", roles.data(), N);
        NpyWriter::WriteInt32(output_dir + "/atom_flags.npy",
                               flags.data(), N * FLAG_COLS);
        total += 2;
    }

    // -------------------------------------------------------------------
    // Residue context projected to atoms.
    //
    // 4 columns per atom:
    //   0  prev_residue_type   AminoAcid enum, Unknown for N-terminal
    //   1  next_residue_type   AminoAcid enum, Unknown for C-terminal
    //   2  is_n_terminal       1 if first residue in chain
    //   3  is_c_terminal       1 if last residue in chain
    //
    // Projecting to per-atom (N rows) instead of per-residue avoids a
    // residue_index→residue join in downstream Python; same memory
    // footprint as the residue table for typical proteins.
    // -------------------------------------------------------------------
    {
        const size_t CTX_COLS = 4;
        std::vector<int32_t> data(N * CTX_COLS, 0);
        for (size_t i = 0; i < N; ++i) {
            const Atom& a = protein.AtomAt(i);
            const Residue& res = protein.ResidueAt(a.residue_index);
            int32_t* r = &data[i * CTX_COLS];
            r[0] = static_cast<int32_t>(res.prev_residue_type);
            r[1] = static_cast<int32_t>(res.next_residue_type);
            r[2] = res.is_n_terminal ? 1 : 0;
            r[3] = res.is_c_terminal ? 1 : 0;
        }
        NpyWriter::WriteInt32(output_dir + "/residue_context.npy",
                               data.data(), N * CTX_COLS);
        total++;
    }

    // -------------------------------------------------------------------
    // Atom names (CSV — IUPAC strings cannot fit cleanly into NPY).
    //
    // Diagnostic surface for human-readable atom identity. The integer
    // columns above are sufficient for typed analysis; the CSV is for
    // grepping logs and printing tensor entries by name. Same precedent
    // as ccd_validation/findings.csv.
    //
    // Columns: atom_index, residue_three_letter, residue_position,
    //          chain_id, atom_name.
    // -------------------------------------------------------------------
    {
        std::ofstream csv(output_dir + "/atom_name.csv");
        if (csv.is_open()) {
            csv << "atom_index,residue_three_letter,residue_position,chain_id,atom_name\n";
            for (size_t i = 0; i < N; ++i) {
                const Atom& a = protein.AtomAt(i);
                const Residue& res = protein.ResidueAt(a.residue_index);
                const char* res_code = res.AminoAcidInfo().three_letter_code;
                csv << i << ","
                    << res_code << ","
                    << res.sequence_number << ","
                    << res.chain_id << ","
                    << a.iupac_name.AsString() << "\n";
            }
            // CSV is not counted in `total` (which counts NPY arrays for
            // the pipeline-completion log line).
        } else {
            OperationLog::Error("WriteAllFeatures",
                "could not open " + output_dir + "/atom_name.csv for write");
        }
    }

    // Per-ring contributions — sparse (atom, ring) pair array.
    // Shape (P, 59) where P = total evaluated (atom, ring) pairs.
    // Columns: [0-8] geometry, [9-17] BS G, [18-26] HM H (pure T2),
    //          [27-35] HM G (shielding), [36-44] quad, [45-53] chi,
    //          [54-56] dispersion, [57-58] azimuthal angle.
    {
        size_t P = 0;
        for (size_t i = 0; i < N; ++i)
            P += conf.AtomAt(i).ring_neighbours.size();

        if (P > 0) {
            const size_t C = 59;
            std::vector<double> data(P * C, 0.0);
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
                    r[8]  = std::exp(-rn.distance_to_center / 4.0);

                    PackST(rn.G_spherical,      r + 9);   // BS shielding kernel
                    PackST(rn.hm_H_spherical,   r + 18);  // HM raw integral (pure T2)
                    PackST(rn.hm_G_spherical,   r + 27);  // HM shielding kernel (T0+T1+T2)
                    PackST(rn.quad_spherical,    r + 36);
                    PackST(rn.chi_spherical,     r + 45);
                    r[54] = rn.disp_scalar;
                    r[55] = static_cast<double>(rn.disp_contacts);
                    r[56] = rn.gaussian_density;
                    r[57] = rn.cos_phi;
                    r[58] = rn.sin_phi;
                    row++;
                }
            }
            NpyWriter::WriteFloat64(output_dir + "/ring_contributions.npy",
                                    data.data(), P, C);
            total++;
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
