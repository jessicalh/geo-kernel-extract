#include "MopacResult.h"
#include "Protein.h"
#include "RuntimeEnvironment.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include <fstream>
#include <sstream>
#include <filesystem>
#include <cmath>
#include <regex>
#include <algorithm>
#include <thread>
#include <cstdlib>

namespace fs = std::filesystem;

namespace nmr {

// Log channel for MOPAC computation


// ============================================================================
// Write .mop input file for MOPAC PM7+MOZYME 1SCF.
//
// Keywords: PM7 MOZYME 1SCF CHARGE=N BONDS MULLIK LET GEO-OK THREADS=8
//   PM7:     semiempirical Hamiltonian
//   MOZYME:  linear-scaling SCF (~45s for 889 atoms)
//   1SCF:    single-point, no geometry optimisation
//   CHARGE:  net formal charge
//   BONDS:   print Wiberg bond order matrix
//   MULLIK:  print Mulliken population analysis (charges + s/p pop)
//   LET:     allow unusual geometries without aborting
//   GEO-OK:  suppress geometry warnings
//   THREADS: OpenMP parallelism
// ============================================================================

static std::string WriteMopFile(const Protein& protein,
                                 const ProteinConformation& conf,
                                 int net_charge,
                                 int threads,
                                 const std::string& path) {
    std::ofstream out(path);
    if (!out.is_open()) return "Cannot open " + path + " for writing";

    out << "PM7 MOZYME 1SCF CHARGE=" << net_charge
        << " BONDS MULLIK LET GEO-OK THREADS=" << threads << "\n";

    std::string name = protein.BuildContext().pdb_source;
    if (!name.empty()) name = fs::path(name).stem().string();
    if (name.empty()) name = "protein";
    out << name << " " << conf.AtomCount() << " atoms\n";
    out << "\n";

    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        Vec3 pos = conf.PositionAt(i);
        std::string elem = SymbolForElement(protein.AtomAt(i).element);
        char line[128];
        snprintf(line, sizeof(line), "  %-2s  %14.8f 0  %14.8f 0  %14.8f 0\n",
                 elem.c_str(), pos.x(), pos.y(), pos.z());
        out << line;
    }
    out.close();
    return "";
}


// ============================================================================
// Parse MOPAC .out file for Mulliken charges, orbital populations,
// bond orders, heat of formation.
//
// This is the text-parsing path needed because run_mopac_from_input()
// does not return a mopac_properties struct — it writes the .out file.
// The parsing matches what mopac_extract.py does, preserving the same
// data and structure.
// ============================================================================

struct MopacParsed {
    std::vector<double> charges;
    std::vector<double> s_pop;
    std::vector<double> p_pop;
    std::vector<MopacBondOrder> bond_orders;
    double heat_of_formation = 0.0;
    Vec3 dipole = Vec3::Zero();
    bool success = false;
    std::string error;
};

static MopacParsed ParseMopacOutput(const std::string& out_path, size_t natoms) {
    MopacParsed result;
    result.charges.resize(natoms, 0.0);
    result.s_pop.resize(natoms, 0.0);
    result.p_pop.resize(natoms, 0.0);

    // Read the .out file using POSIX read to bypass any C++/Fortran
    // stream buffering interaction. run_mopac_from_input uses Fortran I/O
    // which may not fully flush to the C++ iostream layer.
    std::string text;
    {
        std::error_code fec;
        size_t file_size = fs::file_size(out_path, fec);
        if (fec || file_size == 0) {
            result.error = "Cannot stat " + out_path;
            return result;
        }
        text.resize(file_size);
        FILE* fp = fopen(out_path.c_str(), "rb");
        if (!fp) {
            result.error = "Cannot open " + out_path;
            return result;
        }
        size_t nread = fread(text.data(), 1, file_size, fp);
        fclose(fp);
        text.resize(nread);
    }

    // Check for normal termination
    if (text.find("JOB ENDED NORMALLY") == std::string::npos &&
        text.find("ended normally") == std::string::npos) {
        result.error = "MOPAC abnormal termination";
        return result;
    }

    // --- Heat of formation ---
    {
        std::regex hof_re(R"(FINAL HEAT OF FORMATION\s*=\s*([-\d.]+))");
        std::smatch m;
        if (std::regex_search(text, m, hof_re)) {
            result.heat_of_formation = std::stod(m[1].str());
        }
    }

    // --- NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS ---
    // This section has: index, element, charge, ?, s_pop, p_pop
    {
        auto pos = text.find("NET ATOMIC CHARGES AND DIPOLE CONTRIBUTIONS");
        if (pos != std::string::npos) {
            std::istringstream section(text.substr(pos));
            std::string line;
            // Skip header lines (title + blank + column headers)
            std::getline(section, line);  // title
            std::getline(section, line);  // blank or underline
            std::getline(section, line);  // column headers

            while (std::getline(section, line)) {
                if (line.empty() || line.find("DIPOLE") != std::string::npos) break;
                std::istringstream ss(line);
                int idx;
                std::string elem;
                double charge, dummy, sp;
                if (ss >> idx >> elem >> charge >> dummy >> sp) {
                    size_t i = static_cast<size_t>(idx - 1);
                    if (i < natoms) {
                        result.charges[i] = charge;
                        result.s_pop[i] = sp;
                        // p-Pop is absent for hydrogen (only s orbital)
                        double pp = 0.0;
                        if (ss >> pp) {
                            result.p_pop[i] = pp;
                        }
                    }
                }
            }
        }
    }

    // --- BOND ORDERS ---
    // Format: "  4  C     (3.119)    17  C 1.014     6  C 0.968 ..."
    {
        auto pos = text.find("BOND ORDERS");
        if (pos != std::string::npos) {
            std::istringstream section(text.substr(pos));
            std::string line;
            std::getline(section, line);  // "BOND ORDERS" title
            std::getline(section, line);  // blank or header

            std::regex line_re(R"(\s+(\d+)\s+\w+\s+\([\d.]+\)(.*))");
            std::regex pair_re(R"((\d+)\s+\w+\s+([\d.]+))");

            while (std::getline(section, line)) {
                if (line.empty()) continue;
                if (line.find("BOND ORDERS") != std::string::npos) continue;
                // Stop at next section
                if (line.find("CARTESIAN") != std::string::npos ||
                    line.find("ATOMIC") != std::string::npos ||
                    line.find("COMPUTATION") != std::string::npos) break;

                std::smatch lm;
                if (std::regex_match(line, lm, line_re)) {
                    int atom_i = std::stoi(lm[1].str()) - 1;
                    std::string rest = lm[2].str();

                    auto begin = std::sregex_iterator(rest.begin(), rest.end(), pair_re);
                    auto end = std::sregex_iterator();
                    for (auto it = begin; it != end; ++it) {
                        int atom_j = std::stoi((*it)[1].str()) - 1;
                        double bo = std::stod((*it)[2].str());
                        if (bo > 0.01 &&
                            atom_i >= 0 && static_cast<size_t>(atom_i) < natoms &&
                            atom_j >= 0 && static_cast<size_t>(atom_j) < natoms) {
                            MopacBondOrder mbo;
                            mbo.atom_a = static_cast<size_t>(atom_i);
                            mbo.atom_b = static_cast<size_t>(atom_j);
                            mbo.wiberg_order = bo;
                            result.bond_orders.push_back(mbo);
                        }
                    }
                }
            }
        }
    }

    // --- DIPOLE (Debye) ---
    // Table format:
    //           DIPOLE           X         Y         Z       TOTAL
    //  POINT-CHG.        ...
    //  HYBRID            ...
    //  SUM               dx        dy        dz        total
    {
        auto pos = text.find("DIPOLE           X");
        if (pos != std::string::npos) {
            std::istringstream section(text.substr(pos));
            std::string line;
            std::getline(section, line);  // header
            while (std::getline(section, line)) {
                if (line.find("SUM") != std::string::npos) {
                    std::istringstream ss(line);
                    std::string label;
                    double dx, dy, dz;
                    if (ss >> label >> dx >> dy >> dz) {
                        result.dipole = Vec3(dx, dy, dz);
                    }
                    break;
                }
            }
        }
    }

    result.success = true;
    return result;
}


// ============================================================================
// Compute: write .mop, call MOPAC in-process, parse .out
// ============================================================================

std::unique_ptr<MopacResult> MopacResult::Compute(
        ProteinConformation& conf,
        int net_charge,
        int threads) {

    const Protein& protein = conf.ProteinRef();
    const size_t natoms = conf.AtomCount();

    // Resolve thread count: 0 = auto = 3/4 of hardware concurrency.
    // This machine is dedicated to this work — give MOPAC most of the cores.
    if (threads <= 0) {
        int hw = static_cast<int>(std::thread::hardware_concurrency());
        threads = std::max(4, (hw * 3) / 4);
    }

    // Set OMP environment for MOPAC's Fortran/OpenMP internals.
    // OMP_STACKSIZE: Fortran stack overflow at ~500 atoms with the default
    // 8 MB. 2 GB is generous and safe for proteins up to ~4000 atoms.
    // OMP_NUM_THREADS: belt-and-suspenders with the THREADS keyword.
    // These are process-global but MOPAC is the only OpenMP consumer here
    // and we run one protein at a time (no concurrent mozyme_scf calls).
    setenv("OMP_STACKSIZE", "2G", 1);
    setenv("OMP_NUM_THREADS", std::to_string(threads).c_str(), 1);

    OperationLog::Scope scope("MopacResult::Compute",
        "atoms=" + std::to_string(natoms) +
        " charge=" + std::to_string(net_charge) +
        " threads=" + std::to_string(threads));

    // Generate guid-unique temp file paths
    std::string protein_name = protein.BuildContext().pdb_source;
    if (!protein_name.empty())
        protein_name = fs::path(protein_name).stem().string();
    if (protein_name.empty()) protein_name = "protein";

    std::string mop_path = RuntimeEnvironment::TempFilePath(protein_name, "mopac.mop");

    // Write .mop input
    std::string err = WriteMopFile(protein, conf, net_charge, threads, mop_path);
    if (!err.empty()) {
        OperationLog::Error("MopacResult::Compute", err);
        return nullptr;
    }

    // Call MOPAC as a subprocess. We link libmopac but use the binary
    // because run_mopac_from_input() leaves Fortran I/O buffers unflushed —
    // the .out file is incomplete when we try to read it. The subprocess
    // exits cleanly, all file handles close, the output is complete.
    // The binary stays resident in the page cache after the first call.
    const std::string& mopac_bin = RuntimeEnvironment::Mopac();
    if (mopac_bin.empty() || !fs::exists(mopac_bin)) {
        OperationLog::Error("MopacResult::Compute",
            "MOPAC binary not found (configured: " +
            (mopac_bin.empty() ? "<not set>" : mopac_bin) + ")");
        return nullptr;
    }

    OperationLog::Log(OperationLog::Level::Info, LogMopac,
        "MopacResult::Compute",
        "running PM7+MOZYME 1SCF atoms=" + std::to_string(natoms));

    std::string cmd = "OMP_STACKSIZE=2G OMP_NUM_THREADS=" +
        std::to_string(threads) + " " + mopac_bin + " " + mop_path +
        " > /dev/null 2>&1";
    int rc = std::system(cmd.c_str());

    // MOPAC writes .out alongside .mop (same stem, different extension)
    std::string out_path = mop_path.substr(0, mop_path.size() - 4) + ".out";

    std::error_code ec;

    if (rc != 0) {
        OperationLog::Error("MopacResult::Compute",
            "MOPAC returned error code " + std::to_string(rc));
        fs::remove(mop_path, ec);
        fs::remove(out_path, ec);
        return nullptr;
    }

    // Parse the .out file
    MopacParsed parsed = ParseMopacOutput(out_path, natoms);

    if (!parsed.success) {
        OperationLog::Error("MopacResult::Compute", parsed.error);
        // Leave temp files for debugging
        OperationLog::Error("MopacResult::Compute",
            "MOPAC output preserved at: " + out_path);
        return nullptr;
    }

    // Clean up temp files on success
    std::string stem = mop_path.substr(0, mop_path.size() - 4);
    fs::remove(mop_path, ec);
    fs::remove(out_path, ec);
    fs::remove(stem + ".arc", ec);
    fs::remove(stem + ".den", ec);

    // Build result
    auto result = std::make_unique<MopacResult>();
    result->charges_ = std::move(parsed.charges);
    result->s_pop_ = std::move(parsed.s_pop);
    result->p_pop_ = std::move(parsed.p_pop);
    result->heat_of_formation_ = parsed.heat_of_formation;
    result->dipole_ = parsed.dipole;
    result->bond_orders_ = std::move(parsed.bond_orders);

    // Build valencies from bond orders (sum of orders per atom = valency)
    result->valencies_.resize(natoms, 0.0);
    for (const auto& bo : result->bond_orders_) {
        result->valencies_[bo.atom_a] += bo.wiberg_order;
        // Bond orders are listed once per pair from the MOPAC output
        // (atom_i lists atom_j, but atom_j's line also lists atom_i).
        // We only add for atom_a here; atom_b's contribution comes from
        // the reversed entry elsewhere in the bond order list.
    }

    // Build O(1) bond order lookup map
    for (const auto& bo : result->bond_orders_) {
        uint64_t key = PairKey(bo.atom_a, bo.atom_b);
        // MOPAC lists each pair from both sides; keep the first (or max)
        auto it = result->bond_order_map_.find(key);
        if (it == result->bond_order_map_.end() || bo.wiberg_order > it->second) {
            result->bond_order_map_[key] = bo.wiberg_order;
        }
    }

    // Build topology bridge: parallel to protein.Bonds()
    // Map each covalent bond's atom pair to its index for reverse lookup
    std::unordered_map<uint64_t, size_t> bond_pair_to_index;
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const auto& bond = protein.BondAt(bi);
        bond_pair_to_index[PairKey(bond.atom_index_a, bond.atom_index_b)] = bi;
    }

    result->topology_bond_orders_.resize(protein.BondCount(), 0.0);
    for (const auto& [key, order] : result->bond_order_map_) {
        auto it = bond_pair_to_index.find(key);
        if (it != bond_pair_to_index.end()) {
            result->topology_bond_orders_[it->second] = order;
        }
    }

    // Store per-atom data on ConformationAtom
    for (size_t i = 0; i < natoms; ++i) {
        auto& ca = conf.MutableAtomAt(i);
        ca.mopac_charge = result->charges_[i];
        ca.mopac_s_pop = result->s_pop_[i];
        ca.mopac_p_pop = result->p_pop_[i];
        ca.mopac_valency = result->valencies_[i];
    }

    // Build per-atom bond neighbour lists on ConformationAtom
    // First, collect per atom
    std::vector<std::vector<MopacBondNeighbour>> per_atom(natoms);
    for (const auto& [key, order] : result->bond_order_map_) {
        size_t a = static_cast<size_t>(key >> 32);
        size_t b = static_cast<size_t>(key & 0xFFFFFFFF);

        // Look up topology bond index
        size_t topo_idx = SIZE_MAX;
        auto tit = bond_pair_to_index.find(key);
        if (tit != bond_pair_to_index.end()) topo_idx = tit->second;

        MopacBondNeighbour nb_a;
        nb_a.other_atom = b;
        nb_a.wiberg_order = order;
        nb_a.topology_bond_index = topo_idx;
        per_atom[a].push_back(nb_a);

        MopacBondNeighbour nb_b;
        nb_b.other_atom = a;
        nb_b.wiberg_order = order;
        nb_b.topology_bond_index = topo_idx;
        per_atom[b].push_back(nb_b);
    }

    // Sort descending by bond order and store
    for (size_t i = 0; i < natoms; ++i) {
        std::sort(per_atom[i].begin(), per_atom[i].end(),
            [](const MopacBondNeighbour& a, const MopacBondNeighbour& b) {
                return a.wiberg_order > b.wiberg_order;
            });
        conf.MutableAtomAt(i).mopac_bond_neighbours = std::move(per_atom[i]);
    }

    OperationLog::Log(OperationLog::Level::Info, LogMopac,
        "MopacResult::Compute",
        "heat=" + std::to_string(result->heat_of_formation_) +
        " kcal/mol, " + std::to_string(result->bond_order_map_.size()) +
        " bond orders");

    return result;
}


// ============================================================================
// Query methods
// ============================================================================

double MopacResult::ChargeAt(size_t i) const {
    return (i < charges_.size()) ? charges_[i] : 0.0;
}

double MopacResult::SPopAt(size_t i) const {
    return (i < s_pop_.size()) ? s_pop_[i] : 0.0;
}

double MopacResult::PPopAt(size_t i) const {
    return (i < p_pop_.size()) ? p_pop_[i] : 0.0;
}

double MopacResult::ValencyAt(size_t i) const {
    return (i < valencies_.size()) ? valencies_[i] : 0.0;
}

double MopacResult::BondOrder(size_t atom_a, size_t atom_b) const {
    auto it = bond_order_map_.find(PairKey(atom_a, atom_b));
    return (it != bond_order_map_.end()) ? it->second : 0.0;
}

double MopacResult::TopologyBondOrder(size_t bond_index) const {
    return (bond_index < topology_bond_orders_.size())
        ? topology_bond_orders_[bond_index] : 0.0;
}


// ============================================================================
// WriteFeatures: compatible with mopac_extract.py output format.
//
//   mopac_charges.npy     — (N,)   Mulliken charges
//   mopac_scalars.npy     — (N, 4) [charge, s_pop, p_pop, valency]
//   mopac_bond_orders.npy — (B, 3) [atom_i, atom_j, bond_order]
//   mopac_global.npy      — (4,)   [heat_of_formation, dipole_x, dipole_y, dipole_z]
// ============================================================================

int MopacResult::WriteFeatures(const ProteinConformation& conf,
                                const std::string& output_dir) const {
    const size_t N = conf.AtomCount();

    // mopac_charges: (N,)
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i)
            data[i] = conf.AtomAt(i).mopac_charge;
        NpyWriter::WriteFloat64(output_dir + "/mopac_charges.npy", data.data(), N);
    }

    // mopac_scalars: (N, 4) — [charge, s_pop, p_pop, valency]
    {
        std::vector<double> data(N * 4);
        for (size_t i = 0; i < N; ++i) {
            data[i*4 + 0] = conf.AtomAt(i).mopac_charge;
            data[i*4 + 1] = conf.AtomAt(i).mopac_s_pop;
            data[i*4 + 2] = conf.AtomAt(i).mopac_p_pop;
            data[i*4 + 3] = conf.AtomAt(i).mopac_valency;
        }
        NpyWriter::WriteFloat64(output_dir + "/mopac_scalars.npy", data.data(), N, 4);
    }

    // mopac_bond_orders: (B, 3) — [atom_i, atom_j, bond_order]
    // Only unique pairs (from bond_order_map_, not the raw list which has duplicates)
    {
        std::vector<double> data;
        data.reserve(bond_order_map_.size() * 3);
        for (const auto& [key, order] : bond_order_map_) {
            size_t a = static_cast<size_t>(key >> 32);
            size_t b = static_cast<size_t>(key & 0xFFFFFFFF);
            data.push_back(static_cast<double>(a));
            data.push_back(static_cast<double>(b));
            data.push_back(order);
        }
        size_t nbonds = bond_order_map_.size();
        if (nbonds > 0) {
            NpyWriter::WriteFloat64(output_dir + "/mopac_bond_orders.npy",
                                    data.data(), nbonds, 3);
        } else {
            // Write empty (0, 3) array
            NpyWriter::WriteFloat64(output_dir + "/mopac_bond_orders.npy",
                                    nullptr, 0, 3);
        }
    }

    // mopac_global: (4,) — [heat_of_formation, dipole_x, dipole_y, dipole_z]
    {
        double data[4] = { heat_of_formation_, dipole_.x(), dipole_.y(), dipole_.z() };
        NpyWriter::WriteFloat64(output_dir + "/mopac_global.npy", data, 4);
    }

    return 4;
}

}  // namespace nmr
