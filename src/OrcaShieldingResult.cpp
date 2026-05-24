#include "OrcaShieldingResult.h"
#include "Protein.h"
#include "NpyWriter.h"
#include "OperationLog.h"

#include <fstream>
#include <sstream>
#include <cstdio>
#include <cmath>
#include <filesystem>

namespace fs = std::filesystem;

namespace nmr {

// ============================================================================
// PDB LOADING BOUNDARY: parse ORCA NMR shielding tensor output.
//
// Format per nucleus:
//   " Nucleus  NNNX :"              (index, element)
//   " --------------"
//   ""
//   "Diamagnetic contribution to the shielding tensor (ppm) :"
//   "   row0_col0  row0_col1  row0_col2"
//   "   row1_col0  row1_col1  row1_col2"
//   "   row2_col0  row2_col1  row2_col2"
//   ""
//   "Paramagnetic contribution to the shielding tensor (ppm):"
//   "   row0_col0  row0_col1  row0_col2"
//   "   row1_col0  row1_col1  row1_col2"
//   "   row2_col0  row2_col1  row2_col2"
//   ""
//   "Total shielding tensor (ppm):"
//   "   row0_col0  row0_col1  row0_col2"
//   "   row1_col0  row1_col1  row1_col2"
//   "   row2_col0  row2_col1  row2_col2"
//
// After parsing: Mat3 + SphericalTensor on ConformationAtom. No strings.
// ============================================================================

struct ParsedNucleus {
    int index = -1;
    Element element = Element::Unknown;
    Mat3 diamagnetic = Mat3::Zero();
    Mat3 paramagnetic = Mat3::Zero();
    Mat3 total = Mat3::Zero();
};


// Read a 3x3 matrix from three consecutive lines.
static bool ReadMatrix(std::ifstream& in, Mat3& m) {
    for (int row = 0; row < 3; ++row) {
        std::string line;
        if (!std::getline(in, line)) return false;
        double a, b, c;
        if (std::sscanf(line.c_str(), " %lf %lf %lf", &a, &b, &c) != 3)
            return false;
        m(row, 0) = a;
        m(row, 1) = b;
        m(row, 2) = c;
    }
    return true;
}


static std::vector<ParsedNucleus> ParseOrcaNmrOutput(const std::string& path) {
    std::vector<ParsedNucleus> nuclei;
    std::ifstream in(path);
    if (!in.is_open()) return nuclei;

    std::string line;

    // Scan for "CHEMICAL SHIELDINGS" section
    bool in_shielding = false;
    while (std::getline(in, line)) {
        if (line.find("CHEMICAL SHIELDINGS") != std::string::npos &&
            line.find("ppm") != std::string::npos) {
            in_shielding = true;
            break;
        }
    }
    if (!in_shielding) return nuclei;

    // Parse nucleus blocks
    while (std::getline(in, line)) {
        // PDB LOADING BOUNDARY: " Nucleus  NNNX :"
        if (line.find("Nucleus") == std::string::npos || line.find(':') == std::string::npos)
            continue;

        // Stop at the next major section
        if (line.find("CHEMICAL SHIFTS") != std::string::npos) break;
        if (line.find("TIMINGS") != std::string::npos) break;

        int idx = -1;
        char elem[4] = {};
        if (std::sscanf(line.c_str(), " Nucleus %d%3s", &idx, elem) < 1)
            continue;

        ParsedNucleus nuc;
        nuc.index = idx;
        // Strip trailing colon/space from element: "N :" -> "N"
        std::string elem_str(elem);
        auto colon = elem_str.find(':');
        if (colon != std::string::npos) elem_str = elem_str.substr(0, colon);
        while (!elem_str.empty() && elem_str.back() == ' ') elem_str.pop_back();
        nuc.element = ElementFromSymbol(elem_str);

        // Skip separator and blank lines, find "Diamagnetic"
        bool found_dia = false;
        for (int skip = 0; skip < 5 && std::getline(in, line); ++skip) {
            if (line.find("Diamagnetic") != std::string::npos) {
                found_dia = true;
                break;
            }
        }
        if (!found_dia) continue;
        if (!ReadMatrix(in, nuc.diamagnetic)) continue;

        // Find "Paramagnetic"
        bool found_para = false;
        for (int skip = 0; skip < 3 && std::getline(in, line); ++skip) {
            if (line.find("Paramagnetic") != std::string::npos) {
                found_para = true;
                break;
            }
        }
        if (!found_para) continue;
        if (!ReadMatrix(in, nuc.paramagnetic)) continue;

        // Find "Total"
        bool found_total = false;
        for (int skip = 0; skip < 3 && std::getline(in, line); ++skip) {
            if (line.find("Total shielding") != std::string::npos) {
                found_total = true;
                break;
            }
        }
        if (!found_total) continue;
        if (!ReadMatrix(in, nuc.total)) continue;

        nuclei.push_back(nuc);
    }

    return nuclei;
}


// ============================================================================
// OrcaShieldingResult::Compute
// ============================================================================

std::unique_ptr<OrcaShieldingResult> OrcaShieldingResult::Compute(
        ProteinConformation& conf,
        const std::string& nmr_out_path) {

    OperationLog::Scope scope("OrcaShieldingResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()) +
        " source=" + nmr_out_path);

    if (!fs::exists(nmr_out_path)) {
        OperationLog::Error("OrcaShieldingResult::Compute",
            "NMR output not found: " + nmr_out_path);
        return nullptr;
    }

    auto nuclei = ParseOrcaNmrOutput(nmr_out_path);
    if (nuclei.empty()) {
        OperationLog::Error("OrcaShieldingResult::Compute",
            "no shielding tensors parsed from " + nmr_out_path);
        return nullptr;
    }

    if (nuclei.size() != conf.AtomCount()) {
        OperationLog::Error("OrcaShieldingResult::Compute",
            "parsed " + std::to_string(nuclei.size()) + " nuclei but "
            "conformation has " + std::to_string(conf.AtomCount()) + " atoms. "
            "NMR output must match the protonated geometry.");
        return nullptr;
    }

    // Verify atom ordering: ORCA nucleus element must match protein atom element.
    // If these disagree, the tensors would go to the wrong atoms. Hard error.
    const Protein& protein = conf.ProteinRef();
    int mismatches = 0;
    for (size_t ai = 0; ai < nuclei.size(); ++ai) {
        Element orca_elem = nuclei[ai].element;
        Element atom_elem = protein.AtomAt(ai).element;
        if (orca_elem != Element::Unknown && orca_elem != atom_elem) {
            if (mismatches < 5) {
                OperationLog::Error("OrcaShieldingResult::Compute",
                    "element mismatch at atom " + std::to_string(ai) +
                    ": ORCA says " + SymbolForElement(orca_elem) +
                    ", protein has " + SymbolForElement(atom_elem));
            }
            mismatches++;
        }
    }
    if (mismatches > 0) {
        OperationLog::Error("OrcaShieldingResult::Compute",
            std::to_string(mismatches) + " element mismatches — "
            "atom ordering between ORCA output and protein does not match. "
            "Refusing to load: every tensor would go to the wrong atom.");
        return nullptr;
    }

    auto result = std::make_unique<OrcaShieldingResult>();
    result->conf_ = &conf;
    result->source_ = nmr_out_path;
    result->parsed_count_ = static_cast<int>(nuclei.size());

    // Store tensors on ConformationAtoms
    for (size_t ai = 0; ai < nuclei.size(); ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        const auto& nuc = nuclei[ai];

        ca.orca_shielding_total = nuc.total;
        ca.orca_shielding_total_spherical = SphericalTensor::Decompose(nuc.total);

        ca.orca_shielding_diamagnetic = nuc.diamagnetic;
        ca.orca_shielding_diamagnetic_spherical = SphericalTensor::Decompose(nuc.diamagnetic);

        ca.orca_shielding_paramagnetic = nuc.paramagnetic;
        ca.orca_shielding_paramagnetic_spherical = SphericalTensor::Decompose(nuc.paramagnetic);

        ca.has_orca_shielding = true;
    }

    OperationLog::Info(LogCharges, "OrcaShieldingResult::Compute",
        "stored " + std::to_string(nuclei.size()) + " shielding tensors from " +
        nmr_out_path);

    return result;
}

static void PackST_O(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int OrcaShieldingResult::WriteFeatures(const ProteinConformation& conf,
                                        const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    std::vector<double> total(N * 9);
    std::vector<double> dia(N * 9);
    std::vector<double> para(N * 9);

    for (size_t i = 0; i < N; ++i) {
        const auto& ca = conf.AtomAt(i);
        PackST_O(ca.orca_shielding_total_spherical, &total[i*9]);
        PackST_O(ca.orca_shielding_diamagnetic_spherical, &dia[i*9]);
        PackST_O(ca.orca_shielding_paramagnetic_spherical, &para[i*9]);
    }

    NpyWriter::WriteFloat64(output_dir + "/orca_total.npy", total.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/orca_diamagnetic.npy", dia.data(), N, 9);
    NpyWriter::WriteFloat64(output_dir + "/orca_paramagnetic.npy", para.data(), N, 9);
    return 3;
}

}  // namespace nmr
