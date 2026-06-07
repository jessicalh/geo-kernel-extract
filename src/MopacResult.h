#pragma once
//
// MopacResult: PM7+MOZYME semiempirical calculation (MOPAC).
//
// Provides per-atom Mulliken charges, orbital populations (s, p),
// Wiberg bond orders (continuous, conformation-dependent), heat of
// formation, and dipole moment.
//
// Dependencies: none.
//
// Runs the MOPAC binary as a subprocess (RuntimeEnvironment::Mopac()):
// writes a .mop input to RuntimeEnvironment::TempFilePath, MOPAC writes
// the .out alongside it, and we parse the .out text. The in-process
// run_mopac_from_input() is deliberately not used — it leaves Fortran
// I/O buffers unflushed, so the .out is incomplete when read.
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "Types.h"
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <unordered_map>

namespace nmr {

struct MopacBondOrder {
    size_t atom_a = 0;
    size_t atom_b = 0;
    double wiberg_order = 0.0;
};

struct MopacRawArtifact {
    std::string filename;
    std::string extension;
    std::vector<unsigned char> bytes;
    bool is_text = false;
    std::uint64_t size_bytes = 0;
    std::string sha256;
};

struct MopacTextSection {
    std::string id;
    std::string label;
    std::size_t start_byte = 0;
    std::size_t end_byte = 0;
    std::size_t start_line = 0;
    std::size_t end_line = 0;
    std::vector<std::string> heading_lines;
    std::vector<std::string> raw_lines;
    std::string parser_status;
};

struct MopacGenericTable {
    std::string id;
    std::string section_id;
    std::string axis;
    std::string label;
    std::vector<std::string> column_labels;
    std::vector<std::string> units;
    std::vector<std::string> row_labels;
    std::vector<std::vector<std::string>> original_strings;
    std::vector<std::vector<double>> numeric_values;
    std::vector<std::string> parse_warnings;
};

struct MopacAuxRecord {
    std::string key;
    std::string unit;
    std::size_t declared_count = 0;
    std::size_t source_line = 0;
    std::vector<std::string> values;
    std::vector<double> numeric_values;
    std::vector<std::string> raw_lines;
};

struct MopacInputAtomLine {
    std::size_t atom_index = 0;
    std::string element;
    Vec3 position = Vec3::Zero();
    int x_flag = 0;
    int y_flag = 0;
    int z_flag = 0;
    std::string raw_line;
};

struct MopacAtomPopulation {
    std::size_t atom_index = 0;
    std::string element;
    double net_charge = std::numeric_limits<double>::quiet_NaN();
    double electron_density = std::numeric_limits<double>::quiet_NaN();
    double s_population = std::numeric_limits<double>::quiet_NaN();
    double p_population = std::numeric_limits<double>::quiet_NaN();
    double d_population = std::numeric_limits<double>::quiet_NaN();
    double f_population = std::numeric_limits<double>::quiet_NaN();
    double dipole_x = std::numeric_limits<double>::quiet_NaN();
    double dipole_y = std::numeric_limits<double>::quiet_NaN();
    double dipole_z = std::numeric_limits<double>::quiet_NaN();
    double dipole_total = std::numeric_limits<double>::quiet_NaN();
    double mopac_valency = std::numeric_limits<double>::quiet_NaN();
    double project_valency = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::string> column_labels;
    std::vector<std::string> original_cells;
};

struct MopacAOBasisFunction {
    std::size_t ao_index = 0;
    std::size_t atom_index = 0;
    std::string symmetry_type;
    double zeta = std::numeric_limits<double>::quiet_NaN();
    int principal_quantum_number = -1;
    double population = std::numeric_limits<double>::quiet_NaN();
};

struct MopacAtomicOrbitalPopulation {
    std::size_t atom_index = 0;
    std::string element;
    std::array<double, 9> populations = {
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN(),
        std::numeric_limits<double>::quiet_NaN()
    };
    std::vector<std::string> column_labels;
    std::vector<std::string> original_cells;
};

struct MopacDipoleRow {
    std::string label;
    double x = std::numeric_limits<double>::quiet_NaN();
    double y = std::numeric_limits<double>::quiet_NaN();
    double z = std::numeric_limits<double>::quiet_NaN();
    double total = std::numeric_limits<double>::quiet_NaN();
    std::string units = "Debye";
    std::string raw_line;
};

struct MopacScalarTerm {
    std::string label;
    double value = std::numeric_limits<double>::quiet_NaN();
    std::string unit;
    std::string source;
    std::string original_string;
};

struct MopacPrintedBondEntry {
    std::size_t printed_entry_index = 0;
    std::size_t row_order = 0;
    std::size_t row_atom = 0;
    std::string row_element;
    double row_valency = std::numeric_limits<double>::quiet_NaN();
    std::size_t neighbour_atom = 0;
    std::string neighbour_element;
    double order = std::numeric_limits<double>::quiet_NaN();
    std::size_t source_line = 0;
    std::string raw_line;
};

struct MopacUniqueBondOrder {
    std::size_t atom_a = 0;
    std::size_t atom_b = 0;
    double max_order = std::numeric_limits<double>::quiet_NaN();
    double mean_order = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::size_t> printed_entry_indices;
    std::size_t topology_bond_index = SIZE_MAX;
};

struct MopacTopologyBondOrderRecord {
    std::size_t bond_index = 0;
    std::size_t atom_a = 0;
    std::size_t atom_b = 0;
    double order = std::numeric_limits<double>::quiet_NaN();
    bool present = false;
    std::size_t unique_pair_index = SIZE_MAX;
    std::string absence_reason;
    std::size_t printed_entry_count = 0;
};

struct MopacMolecularOrbital {
    std::size_t orbital_index = 0;
    double energy = std::numeric_limits<double>::quiet_NaN();
    double occupation = std::numeric_limits<double>::quiet_NaN();
    double bonding_contribution = std::numeric_limits<double>::quiet_NaN();
    std::string label;
};

struct MopacMatrixBlock {
    std::string name;
    std::string unit;
    std::string storage;
    std::size_t rows = 0;
    std::size_t cols = 0;
    std::vector<double> values;
    std::vector<std::string> original_strings;
};

struct MopacRunRecord {
    std::string schema_version = "mopac-full-capture-1.0";
    std::string sidecar_format_status =
        "manifest-only; bulk numeric capture is authoritative in sibling NPY files";
    std::string mopac_binary;
    std::string temp_stem;
    std::string keyword_line;
    std::string title_line;
    std::string comment_line;
    int net_charge = 0;
    int threads = 0;
    std::size_t atom_count = 0;
    std::string method;
    std::string version;
    std::string run_date;
    std::string empirical_formula;
    std::string point_group;
    std::string termination_status;
    std::vector<std::string> warnings;
    std::vector<std::string> errors;
    std::vector<std::string> convergence_records;
    std::vector<std::string> timing_records;
    std::vector<MopacInputAtomLine> input_atoms;
    std::vector<MopacRawArtifact> artifacts;
    std::vector<MopacTextSection> sections;
    std::vector<MopacGenericTable> tables;
    std::vector<MopacAuxRecord> aux_records;
    std::map<std::string, std::size_t> aux_record_index;
    std::vector<MopacAtomPopulation> atom_populations;
    std::vector<MopacAOBasisFunction> ao_basis;
    std::vector<MopacAtomicOrbitalPopulation> atomic_orbital_populations;
    std::vector<MopacDipoleRow> dipole_rows;
    std::vector<MopacScalarTerm> scalar_terms;
    std::vector<MopacPrintedBondEntry> printed_bond_entries;
    std::vector<MopacUniqueBondOrder> unique_bond_orders;
    std::vector<MopacTopologyBondOrderRecord> topology_bond_records;
    std::vector<MopacMolecularOrbital> molecular_orbitals;
    std::vector<MopacMatrixBlock> matrix_blocks;
};

class MopacResult : public ConformationResult {
public:
    std::string Name() const override { return "MopacResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Factory: run PM7+MOZYME on the conformation atoms.
    // net_charge: total charge of the system.
    // threads: OpenMP threads for MOZYME. 0 = auto (hardware_concurrency * 3/4).
    // Returns nullptr on failure (logged via OperationLog).
    static std::unique_ptr<MopacResult> Compute(
        ProteinConformation& conf,
        int net_charge = 0,
        int threads = 0);

    // --- Per-atom queries (O(1)) ---
    double ChargeAt(size_t atom_index) const;
    double SPopAt(size_t atom_index) const;
    double PPopAt(size_t atom_index) const;
    double ValencyAt(size_t atom_index) const;

    // --- Bond order queries ---

    // O(1) lookup by atom pair. Returns 0.0 if no MOPAC bond.
    // Symmetric: BondOrder(a,b) == BondOrder(b,a).
    double BondOrder(size_t atom_a, size_t atom_b) const;

    // O(1) lookup by topology bond index (parallel to protein.Bonds()).
    // Returns 0.0 if MOPAC did not report that bond pair.
    double TopologyBondOrder(size_t bond_index) const;
    const std::vector<double>& TopologyBondOrders() const { return topology_bond_orders_; }

    // Full bond order list (all pairs MOPAC reported).
    const std::vector<MopacBondOrder>& AllBondOrders() const { return bond_orders_; }

    // --- Molecule-level ---
    double HeatOfFormation() const { return heat_of_formation_; }
    Vec3 Dipole() const { return dipole_; }

    // --- Full-run capture accessors ---
    const MopacRunRecord& RunRecord() const { return run_record_; }
    const std::vector<MopacRawArtifact>& RawArtifacts() const { return run_record_.artifacts; }
    const std::vector<MopacTextSection>& OutputSections() const { return run_record_.sections; }
    const std::vector<MopacGenericTable>& ParsedTables() const { return run_record_.tables; }
    const std::vector<MopacAuxRecord>& AuxRecords() const { return run_record_.aux_records; }
    const std::vector<MopacAtomPopulation>& AtomPopulations() const { return run_record_.atom_populations; }
    const std::vector<MopacAOBasisFunction>& AOBasis() const { return run_record_.ao_basis; }
    const std::vector<MopacAtomicOrbitalPopulation>& AtomicOrbitalPopulations() const { return run_record_.atomic_orbital_populations; }
    const std::vector<MopacDipoleRow>& DipoleRows() const { return run_record_.dipole_rows; }
    const std::vector<MopacScalarTerm>& ScalarTerms() const { return run_record_.scalar_terms; }
    const std::vector<MopacPrintedBondEntry>& PrintedBondEntries() const { return run_record_.printed_bond_entries; }
    const std::vector<MopacUniqueBondOrder>& UniqueBondOrders() const { return run_record_.unique_bond_orders; }
    const std::vector<MopacTopologyBondOrderRecord>& TopologyBondRecords() const { return run_record_.topology_bond_records; }
    const std::vector<MopacMolecularOrbital>& MolecularOrbitals() const { return run_record_.molecular_orbitals; }
    const std::vector<MopacMatrixBlock>& MatrixBlocks() const { return run_record_.matrix_blocks; }

    const MopacAtomPopulation* AtomPopulationAt(size_t atom_index) const;
    const MopacAOBasisFunction* AOBasisAt(size_t ao_index) const;
    const MopacAuxRecord* AuxRecordByKey(const std::string& key) const;

    // --- Feature output ---
    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    // Per-atom (parallel to conf.AtomCount())
    std::vector<double> charges_;
    std::vector<double> s_pop_;
    std::vector<double> p_pop_;
    std::vector<double> valencies_;

    // Bond orders: full list + O(1) lookup
    std::vector<MopacBondOrder> bond_orders_;
    std::unordered_map<uint64_t, double> bond_order_map_;  // key: min(a,b)<<32 | max(a,b)

    // Parallel to protein.Bonds() for topology bridge
    std::vector<double> topology_bond_orders_;

    // Molecule-level
    double heat_of_formation_ = 0.0;
    Vec3 dipole_ = Vec3::Zero();
    MopacRunRecord run_record_;

    // Hash key for atom pair (symmetric)
    static uint64_t PairKey(size_t a, size_t b) {
        size_t lo = (a < b) ? a : b;
        size_t hi = (a < b) ? b : a;
        return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
    }
};

}  // namespace nmr
