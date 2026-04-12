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


class GromacsProteinAtom {
    friend class GromacsProtein;
public:
    // Index into the Protein's atom list (for identity lookup).
    size_t atom_index() const { return atom_index_; }

    // === Water environment accumulators ===
    Welford water_n_first;       // first-shell water count (within 3.5A)
    Welford water_n_second;      // second-shell water count (3.5-5.5A)
    Welford water_emag;          // water E-field magnitude (V/A)
    Welford water_emag_first;    // first-shell E-field magnitude

    // === Solvent exposure ===
    Welford sasa;                // Shrake-Rupley SASA (A^2)

    // === Hydration geometry ===
    Welford half_shell;          // half-shell asymmetry (0-1)
    Welford dipole_cos;          // mean water dipole orientation
    Welford nearest_ion_dist;    // distance to closest ion (A)

    // === Positional dynamics ===
    Welford position_x;          // for per-atom RMSF computation
    Welford position_y;
    Welford position_z;

    double RMSF() const {
        return std::sqrt(position_x.Variance()
                       + position_y.Variance()
                       + position_z.Variance());
    }

    // === Frame counting ===
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
