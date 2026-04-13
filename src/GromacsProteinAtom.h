#pragma once
//
// GromacsProteinAtom: per-atom accumulated data across trajectory frames.
//
// This is trajectory-level accumulation, NOT feature extraction.
// Feature extraction writes to ConformationAtom (per-frame, dies
// with the frame). This class accumulates statistics ACROSS frames
// by reading from ConformationAtom after each frame's calculators
// have run. The Welford accumulators live for the whole trajectory.
//
// Private constructor -- only GromacsProtein can create these.
// Atom identity (element, residue, name) goes through the Protein
// back-pointer, NOT duplicated here.
//
// Written to atom_catalog.csv by GromacsFinalResult at the end.
//
// Future: rename to TrajectoryProteinAtom.
//

#include <cmath>
#include <limits>

namespace nmr {

class GromacsProtein;

// Welford online accumulator: mean, variance, min/max with frame indices.
struct Welford {
    int    count = 0;
    double mean  = 0.0;
    double M2    = 0.0;
    double min_val =  std::numeric_limits<double>::max();
    double max_val = -std::numeric_limits<double>::max();
    int    min_frame = -1;
    int    max_frame = -1;

    void Update(double x, int frame) {
        ++count;
        double delta = x - mean;
        mean += delta / count;
        double delta2 = x - mean;
        M2 += delta * delta2;

        if (x < min_val) { min_val = x; min_frame = frame; }
        if (x > max_val) { max_val = x; max_frame = frame; }
    }

    double Variance() const {
        return (count > 1) ? M2 / (count - 1) : 0.0;
    }

    double Std() const {
        return std::sqrt(Variance());
    }

    double Range() const {
        return (count > 0) ? max_val - min_val : 0.0;
    }
};


// Frame-to-frame derivative tracker: Welford on (x_current - x_prev).
// Call UpdateDelta(x, frame) every frame; first frame is skipped.
struct DeltaTracker {
    Welford delta;
    double  prev = 0.0;
    bool    has_prev = false;

    void UpdateDelta(double x, int frame) {
        if (has_prev)
            delta.Update(x - prev, frame);
        prev = x;
        has_prev = true;
    }
};


// Transition counter for categorical or periodic quantities.
// Counts frames where the value changes by more than a threshold.
struct TransitionCounter {
    int    transitions = 0;
    int    frames = 0;
    double prev = 0.0;
    bool   has_prev = false;

    void Update(double x, double threshold) {
        ++frames;
        if (has_prev && std::abs(x - prev) > threshold)
            ++transitions;
        prev = x;
        has_prev = true;
    }

    double Rate() const {
        return (frames > 1) ? static_cast<double>(transitions) / (frames - 1) : 0.0;
    }
};


// Per-bond accumulated data across trajectory frames.
// One per covalent bond in the protein topology.
struct GromacsProteinBond {
    size_t atom_a = 0;
    size_t atom_b = 0;

    Welford length;              // bond length (A)
    DeltaTracker length_delta;   // bond length fluctuation rate
};


// Named pointer to a Welford accumulator — for data-driven iteration.
// WriteH5 and WriteCatalog iterate these instead of maintaining a
// parallel switch statement. Adding a Welford = add the field +
// one line in AllWelfords(). Nothing else changes.
struct NamedWelford {
    const char* name;
    const Welford* w;
};


class GromacsProteinAtom {
    friend class GromacsProtein;
public:
    // Index into the Protein's atom list (for identity lookup).
    size_t atom_index() const { return atom_index_; }

    // All Welford accumulators with their names — single source of truth
    // for CSV columns, H5 rollup, and SDK column names. Add new Welfords
    // here and they flow through to all outputs automatically.
    std::vector<NamedWelford> AllWelfords() const {
        return {
            // Ring current
            {"bs_T0", &bs_T0}, {"bs_T2mag", &bs_T2mag},
            {"hm_T0", &hm_T0}, {"hm_T2mag", &hm_T2mag},
            {"rs_T0", &rs_T0}, {"rs_T2mag", &rs_T2mag},
            // McConnell
            {"mc_T0", &mc_T0}, {"mc_T2mag", &mc_T2mag},
            {"mc_aromatic", &mc_aromatic}, {"mc_backbone", &mc_backbone},
            // PiQuad + Dispersion
            {"pq_T0", &pq_T0}, {"pq_T2mag", &pq_T2mag},
            {"disp_T0", &disp_T0}, {"disp_T2mag", &disp_T2mag},
            // H-bond
            {"hbond_inv_d3", &hbond_inv_d3}, {"hbond_count", &hbond_count},
            // APBS
            {"apbs_emag", &apbs_emag},
            {"apbs_efg_T0", &apbs_efg_T0}, {"apbs_efg_T2mag", &apbs_efg_T2mag},
            // AIMNet2
            {"aimnet2_charge", &aimnet2_charge},
            {"aimnet2_efg_T0", &aimnet2_efg_T0},
            {"aimnet2_efg_T2mag", &aimnet2_efg_T2mag},
            // SASA
            {"sasa", &sasa},
            {"sasa_normal_x", &sasa_normal_x},
            {"sasa_normal_y", &sasa_normal_y},
            {"sasa_normal_z", &sasa_normal_z},
            // Water
            {"water_n_first", &water_n_first}, {"water_n_second", &water_n_second},
            {"water_emag", &water_emag}, {"water_emag_first", &water_emag_first},
            {"water_efield_x", &water_efield_x},
            {"water_efield_y", &water_efield_y},
            {"water_efield_z", &water_efield_z},
            // Hydration
            {"half_shell", &half_shell}, {"dipole_cos", &dipole_cos},
            {"nearest_ion_dist", &nearest_ion_dist},
            // DSSP
            {"phi_cos", &phi_cos}, {"psi_cos", &psi_cos},
            {"dssp_hbond_energy", &dssp_hbond_energy},
            {"chi1_cos", &chi_cos[0]}, {"chi2_cos", &chi_cos[1]},
            {"chi3_cos", &chi_cos[2]}, {"chi4_cos", &chi_cos[3]},
            // Bond geometry
            {"mean_bond_angle_cos", &mean_bond_angle_cos},
            // Deltas (expose the inner Welford on the difference)
            {"bs_T0_delta", &bs_T0_delta.delta},
            {"aimnet2_charge_delta", &aimnet2_charge_delta.delta},
            {"sasa_delta", &sasa_delta.delta},
            {"water_n_first_delta", &water_n_first_delta.delta},
        };
    }

    // =================================================================
    // Positional dynamics
    // =================================================================
    Welford position_x;
    Welford position_y;
    Welford position_z;

    double RMSF() const {
        return std::sqrt(position_x.Variance()
                       + position_y.Variance()
                       + position_z.Variance());
    }

    // =================================================================
    // Ring current (BiotSavartResult, HaighMallionResult)
    // =================================================================
    Welford bs_T0;               // Biot-Savart isotropic total
    Welford bs_T2mag;            // Biot-Savart |T2| total
    Welford hm_T0;               // Haigh-Mallion isotropic total
    Welford hm_T2mag;            // Haigh-Mallion |T2| total
    Welford rs_T0;               // Ring susceptibility isotropic
    Welford rs_T2mag;            // Ring susceptibility |T2|

    DeltaTracker bs_T0_delta;    // frame-to-frame ring current change

    // =================================================================
    // Bond anisotropy (McConnellResult)
    // =================================================================
    Welford mc_T0;               // McConnell isotropic total
    Welford mc_T2mag;            // McConnell |T2| total
    Welford mc_aromatic;         // aromatic bond contribution (scalar)
    Welford mc_backbone;         // CO + CN contributions

    // =================================================================
    // Pi-quadrupole (PiQuadrupoleResult)
    // =================================================================
    Welford pq_T0;
    Welford pq_T2mag;

    // =================================================================
    // Dispersion (DispersionResult)
    // =================================================================
    Welford disp_T0;
    Welford disp_T2mag;

    // =================================================================
    // H-bond (HBondResult)
    // =================================================================
    Welford hbond_inv_d3;        // 1/d^3 to nearest H-bond partner
    Welford hbond_count;         // count within 3.5A

    // =================================================================
    // APBS solvated field (ApbsFieldResult)
    // =================================================================
    Welford apbs_emag;           // solvated E-field magnitude (V/A)
    Welford apbs_efg_T0;         // solvated EFG isotropic
    Welford apbs_efg_T2mag;      // solvated EFG |T2|

    // =================================================================
    // AIMNet2 charges (AIMNet2Result)
    // =================================================================
    Welford aimnet2_charge;      // Hirshfeld charge — variance IS polarisability
    Welford aimnet2_efg_T0;      // Coulomb EFG from AIMNet2 charges
    Welford aimnet2_efg_T2mag;

    DeltaTracker aimnet2_charge_delta;  // charge fluctuation rate

    // =================================================================
    // Solvent-accessible surface area (SasaResult)
    // =================================================================
    Welford sasa;
    Welford sasa_normal_x;       // surface normal x component
    Welford sasa_normal_y;       // surface normal y component
    Welford sasa_normal_z;       // surface normal z component

    DeltaTracker sasa_delta;     // breathing rate

    // =================================================================
    // Water environment (WaterFieldResult)
    // =================================================================
    Welford water_n_first;       // first-shell water count (within 3.5A)
    Welford water_n_second;      // second-shell water count (3.5-5.5A)
    Welford water_emag;          // water E-field magnitude (V/A)
    Welford water_emag_first;    // first-shell E-field magnitude
    Welford water_efield_x;      // E-field components (for anisotropy)
    Welford water_efield_y;
    Welford water_efield_z;

    DeltaTracker water_n_first_delta;  // water exchange rate

    // =================================================================
    // Hydration geometry (HydrationShellResult)
    // =================================================================
    Welford half_shell;          // half-shell asymmetry (0-1)
    Welford dipole_cos;          // mean water dipole orientation
    Welford nearest_ion_dist;    // distance to closest ion (A)

    // =================================================================
    // DSSP dynamics (DsspResult — per-residue, broadcast to atoms)
    // =================================================================
    // Chi1-4: circular quantities — track cos(chi) for Welford, transitions
    // for rotameric jumps. Threshold 0.5 on cos ≈ 60 degree change.
    // chi[k] is only populated for residues with >= k+1 chi angles.
    Welford chi_cos[4];                  // cos(chi1..4) — circular mean/variance
    TransitionCounter chi_transitions[4]; // rotamer flips per angle

    // Secondary structure: categorical, only transitions meaningful
    TransitionCounter ss8_transitions;   // SS state changes (threshold 0.5 on int cast)

    // H-bond energy: strongest acceptor energy per residue
    Welford dssp_hbond_energy;   // best (most negative) H-bond energy

    // Phi/psi backbone dihedrals
    Welford phi_cos;
    Welford psi_cos;

    // =================================================================
    // Bond geometry at this atom (for GNN edge messages)
    // =================================================================
    // Mean bond angle: average cos(angle) across all bond angles at this
    // atom (A-this-B for every pair of bonded neighbors A, B).
    // Variance tells you how rigid the local geometry is.
    Welford mean_bond_angle_cos;

    // Number of covalent neighbors (degree in bond graph).
    // Constant across frames but stored here for catalog output.
    int n_bonded_neighbors = 0;

    // =================================================================
    // Frame counting
    // =================================================================
    int n_frames_dry = 0;        // frames where water_n_first == 0 (buried)
    int n_frames_exposed = 0;    // frames where water_n_first >= 4 (exposed)

    // === GROMACS energy at this atom's extreme frames ===
    // (filled by cross-referencing .edr data at min/max frame indices)
    double energy_at_water_min = 0.0;
    double energy_at_water_max = 0.0;

private:
    explicit GromacsProteinAtom(size_t idx) : atom_index_(idx) {}
    const size_t atom_index_;
};

}  // namespace nmr
