#include "DsspResult.h"
#include "Protein.h"
#include "PhysicalConstants.h"
#include "NpyWriter.h"

#include <cif++.hpp>
#include <cif++/pdb/pdb2cif.hpp>
#include <dssp.hpp>

#include <fstream>
#include <sstream>
#include <filesystem>
#include <map>
#include <cstdio>
#include <cmath>

namespace fs = std::filesystem;

namespace nmr {

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

        // min_poly_proline_stretch_length = 3 (standard, Adzhubei & Sternberg 1993)
        // calculate_accessibility = true (SASA via Lee & Richards 1971, 1.4A probe)
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

    return result;
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
// WriteFeatures: dssp_backbone.npy (N, 5)
//
// Per-atom, broadcast from per-residue via Protein atom→residue mapping.
//
// Columns: phi (rad), psi (rad), sasa (A^2), ss_helix (0/1), ss_sheet (0/1)
//
// DSSP secondary structure alphabet: H(alpha-helix), G(3_10-helix),
// I(pi-helix), E(extended strand), B(isolated bridge), T(turn), S(bend),
// C(coil). For NMR purposes these collapse to helix(H/G/I), sheet(E/B),
// other(T/S/C) — the "other" case is implicit (both columns zero).
// ============================================================================

int DsspResult::WriteFeatures(const ProteinConformation& conf,
                               const std::string& output_dir) const {
    const Protein& protein = conf.ProteinRef();
    const size_t N = conf.AtomCount();

    std::vector<double> data(N * 5, 0.0);

    for (size_t i = 0; i < N; ++i) {
        size_t ri = protein.AtomAt(i).residue_index;
        if (ri < residues_.size()) {
            const auto& dr = residues_[ri];
            data[i*5 + 0] = dr.phi;
            data[i*5 + 1] = dr.psi;
            data[i*5 + 2] = dr.sasa;
            char ss = dr.secondary_structure;
            data[i*5 + 3] = (ss == 'H' || ss == 'G' || ss == 'I') ? 1.0 : 0.0;
            data[i*5 + 4] = (ss == 'E' || ss == 'B') ? 1.0 : 0.0;
        }
    }

    NpyWriter::WriteFloat64(output_dir + "/dssp_backbone.npy", data.data(), N, 5);
    return 1;
}

}  // namespace nmr
