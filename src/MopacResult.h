#pragma once
//
// MopacResult: PM7+MOZYME semiempirical calculation via libmopac.
//
// Provides per-atom Mulliken charges, orbital populations (s, p),
// Wiberg bond orders (continuous, conformation-dependent), heat of
// formation, and dipole moment.
//
// Dependencies: none.
//
// Calls MOPAC via run_mopac_from_input() (linked libmopac, in-process).
// Writes a .mop input file to RuntimeEnvironment::TempFilePath,
// MOPAC writes the .out alongside it, we parse both structured API
// output and text output for orbital populations.
//
// Replaces XtbChargeResult. Same tier in OperationRunner (0.5).
//
// Future consumers:
//   MopacCoulombResult — ChargeAt(i) for QM-derived EFG
//   MopacBondAnisotropyResult — BondOrder(a,b) for delta-chi modulation
//   DipolarNearFieldFilter — mopac_p_pop for lobe offset calibration
//

#include "ConformationResult.h"
#include "ProteinConformation.h"
#include <vector>
#include <unordered_map>

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

    // Hash key for atom pair (symmetric)
    static uint64_t PairKey(size_t a, size_t b) {
        size_t lo = (a < b) ? a : b;
        size_t hi = (a < b) ? b : a;
        return (static_cast<uint64_t>(lo) << 32) | static_cast<uint64_t>(hi);
    }
};

}  // namespace nmr
