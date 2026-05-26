#pragma once
// libpq-backed loader for tripeptide DFT rows. QueryNearest rounds
// phi/psi/chi angles to the database grid and preserves frame_type so
// callers can apply method-specific calibration.

#include "LarsenResidue.h"
#include "Types.h"

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

struct pg_conn;
typedef struct pg_conn PGconn;

namespace nmr {


// Position and tensor share the DFT-output frame, so one Kabsch
// rotation aligns both onto the protein.
struct TripeptideDftAtom {
    int     atom_idx = 0;            // 1-based index in original DFT geometry
    Element element  = Element::Unknown;
    Vec3    position = Vec3::Zero();  // Angstroms, DFT-output frame

    Mat3   shielding_tensor = Mat3::Zero();   // full asymmetric 3x3, ppm
    double isotropic        = 0.0;            // T0 / 3, ppm
    double anisotropy       = 0.0;            // Haeberlen Δσ, ppm

    // sphericart layout matches SphericalTensor::T2.
    std::array<double, 5> t2_components = {};
};


// One row from raw_dft_calculations: a tripeptide at a (φ, ψ, χ*) pose.
struct TripeptideDftRecord {
    int  calc_id = 0;
    std::string tripeptide;          // "AAA", "AFA", "ASA", etc.

    int phi = 0, psi = 0;
    int chi1 = 0, chi2 = 0, chi3 = 0, chi4 = 0;

    std::string frame_type;

    int n_atoms = 0;
    std::vector<TripeptideDftAtom> atoms;

    // Perceived typed model of the 5-piece tripeptide. nullopt if
    // perception failed (logged via OperationLog at perception time);
    // calculators should decline records with no larsen.
    std::optional<LarsenTripeptide> larsen;

    bool IsHit() const { return calc_id != 0; }
};


class TripeptideDftTable {
public:
    explicit TripeptideDftTable(const std::string& conn_str);
    ~TripeptideDftTable();

    TripeptideDftTable(const TripeptideDftTable&) = delete;
    TripeptideDftTable& operator=(const TripeptideDftTable&) = delete;

    // Angles are degrees. n_chi_axes >= 0 limits the chi columns in
    // the WHERE clause; callers use that for explicit chi-depth
    // fallback. his_variant_hint applies only to HIS and uses
    // Residue::protonation_variant_index values 0=HID, 1=HIE, 2=HIP.
    TripeptideDftRecord QueryNearest(
        char residue_letter,
        double phi, double psi,
        double chi1 = 0.0, double chi2 = 0.0,
        double chi3 = 0.0, double chi4 = 0.0,
        int    n_chi_axes        = -1,
        int    his_variant_hint  = -1) const;

    // Transient connection breaks may not surface until QueryNearest.
    bool IsConnected() const;

private:
    PGconn* conn_ = nullptr;
};


}  // namespace nmr
