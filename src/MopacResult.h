#pragma once
//
// MopacResult: diskless PM7/MOZYME single-point calculation.
//
// Compute calls the pinned libmopac mozyme_scf API in a crash-contained
// worker process, copies the complete direct electronic result, and stores
// the calculator-owned ConformationAtom fields. WriteFeatures only reads
// that stored result back to NPY.
//
// Dependencies: none. OperationRunner charge-gates this calculator in its
// fixed sequence after ChargeAssignmentResult has attached.
//


#include "ConformationResult.h"
#include "ProteinConformation.h"
#include "Types.h"

#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace nmr {

struct MopacBondOrder {
    size_t atom_a = 0;
    size_t atom_b = 0;
    double wiberg_order = 0.0;
};

class MopacResult : public ConformationResult {
public:
    std::string Name() const override { return "MopacResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }

    // Run PM7/MOZYME/1SCF at the supplied exact topology charge.
    // threads=0 selects the established 3/4-hardware policy. On failure,
    // returns null and, when supplied, fills error_out with the libmopac or
    // worker-process diagnostic.
    static std::unique_ptr<MopacResult> Compute(
        ProteinConformation& conf,
        int net_charge = 0,
        int threads = 0,
        std::string* error_out = nullptr);

    double ChargeAt(size_t atom_index) const;
    double SPopAt(size_t atom_index) const;
    double PPopAt(size_t atom_index) const;
    double ValencyAt(size_t atom_index) const;

    // Symmetric O(1) compatibility lookup. This reproduces the legacy
    // compact text-table parser (Fortran ordering, first six fields, F6.3,
    // parsed value >0.01). Zero means the pair was not visible there; the
    // complete direct sparse matrix remains in the direct NPY arrays.
    double BondOrder(size_t atom_a, size_t atom_b) const;
    double TopologyBondOrder(size_t bond_index) const;
    const std::vector<double>& TopologyBondOrders() const {
        return topology_bond_orders_;
    }
    const std::vector<MopacBondOrder>& AllBondOrders() const {
        return bond_orders_;
    }

    double HeatOfFormation() const {
        return compatibility_heat_of_formation_;
    }
    Vec3 Dipole() const { return compatibility_dipole_; }

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;

private:
    struct UniqueBond {
        std::int32_t atom_a = 0;
        std::int32_t atom_b = 0;
        double order = 0.0;       // a<b CSC entry: probe/archive convention
        double max_order = 0.0;   // exact maximum of directed API entries
        double mean_order = 0.0;
        double compatibility_max_order = 0.0;
        double compatibility_mean_order = 0.0;
        std::int32_t topology_bond_index = -1;
        std::int32_t directed_entry = -1;
    };

    // Per-atom direct quantities. These retain the exact libmopac doubles;
    // the separately stored compatibility projection below reproduces the
    // decimal fields read from the legacy MOPAC text tables.
    std::vector<double> charges_;
    std::vector<double> electron_populations_;
    std::vector<double> s_pop_;
    std::vector<double> p_pop_;
    std::vector<double> d_pop_;
    std::vector<double> valencies_;          // MOPAC CSC diagonal
    std::vector<double> compatibility_charges_;
    std::vector<double> compatibility_electron_populations_;
    std::vector<double> compatibility_s_pop_;
    std::vector<double> compatibility_p_pop_;
    std::vector<double> compatibility_d_pop_;
    std::vector<double> compatibility_valencies_;
    std::vector<double> compatibility_project_valencies_;
    std::vector<double> compatibility_atomic_orbital_populations_;
    std::vector<double> compatibility_atomic_orbital_population_totals_;

    // Sparse Wiberg matrix and its compatibility/query projections.
    std::vector<std::int32_t> bond_index_;
    std::vector<std::int32_t> bond_atom_;
    std::vector<double> bond_order_;
    std::vector<MopacBondOrder> bond_orders_;  // exact directed API entries
    std::vector<UniqueBond> unique_bonds_;
    // Legacy compact-table projection: unordered-pair maximum after MOPAC's
    // first-six/F6.3/>0.01 printed-row selection. This is intentionally
    // distinct from the complete API sparse pair set above.
    std::unordered_map<std::uint64_t, double> bond_order_map_;
    std::vector<double> topology_bond_orders_;
    std::vector<std::int32_t> topology_unique_pair_index_;

    // Molecular direct quantities.
    double heat_of_formation_ = 0.0;
    Vec3 dipole_ = Vec3::Zero();
    Vec3 dipole_point_charge_ = Vec3::Zero();
    Vec3 dipole_hybridization_ = Vec3::Zero();
    double compatibility_heat_of_formation_ = 0.0;
    Vec3 compatibility_dipole_ = Vec3::Zero();

    // AO density blocks. Stored in emitted C order: atom/entry, AO row,
    // AO column. The raw directed bond blocks remain parallel to the CSC
    // arrays; unique_bond_ao_density_ is the probe-compatible a<b view.
    std::int32_t ao_max_orbitals_ = 0;
    std::vector<std::int32_t> ao_orbitals_;
    std::vector<double> atom_ao_density_;
    std::vector<double> bond_ao_density_directed_;
    std::vector<double> unique_bond_ao_density_;
    std::vector<double> atomic_orbital_populations_;

    // MOZYME Lewis structure.
    std::vector<std::int32_t> lewis_bond_count_;
    std::vector<std::int32_t> lewis_bond_atoms_;

    // Probe-compatible packed live LMO views.
    std::vector<std::int32_t> lmo_occupied_atom_counts_;
    std::vector<std::int32_t> lmo_occupied_atoms_;
    std::vector<double> lmo_occupied_coefficients_;
    std::vector<std::int32_t> lmo_virtual_atom_counts_;
    std::vector<std::int32_t> lmo_virtual_atoms_;
    std::vector<double> lmo_virtual_coefficients_;
    std::vector<double> lmo_energy_levels_;

    // Exact noncompact native state storage and the libmopac offsets that
    // interpret its live slices.
    std::vector<std::int32_t> lmo_occupied_atom_offsets_native_;
    std::vector<std::int32_t> lmo_virtual_atom_offsets_native_;
    std::vector<std::int32_t> lmo_occupied_coefficient_offsets_native_;
    std::vector<std::int32_t> lmo_virtual_coefficient_offsets_native_;
    std::vector<std::int32_t> lmo_occupied_atom_storage_native_;
    std::vector<std::int32_t> lmo_virtual_atom_storage_native_;
    std::vector<double> lmo_occupied_coefficient_storage_native_;
    std::vector<double> lmo_virtual_coefficient_storage_native_;
    std::array<std::int32_t, 7> mozyme_state_dimensions_{};

    static std::uint64_t PairKey(size_t a, size_t b) {
        const size_t lo = (a < b) ? a : b;
        const size_t hi = (a < b) ? b : a;
        return (static_cast<std::uint64_t>(lo) << 32) |
               static_cast<std::uint64_t>(hi);
    }
};

}  // namespace nmr
