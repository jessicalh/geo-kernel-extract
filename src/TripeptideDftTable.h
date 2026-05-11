#pragma once
//
// TripeptideDftTable: libpq-backed loader for ProCS15 tripeptide DFT
// data from the local tensorcs15 Postgres replica.
//
// Owns a PGconn for the Session lifetime. Per-call QueryNearest returns
// a TripeptideDftRecord with the central-residue atom positions, full
// asymmetric Mat3 shielding tensors (T0+T1+T2 preserved, pre-decomposed
// at ingest), backbone atom indices into the record, and the
// frame_type method discriminator.
//
// frame_type per row encodes the DFT method (the SER discontinuity):
//   gaussian_standard_orientation : OPBE/6-31G(d,p) GIAO CPCM(water)
//                                    Larsen 2015 ProCS15, 19 residues
//   orca_input_orientation        : PBE/6-31G(d,p) GIAO CPCM(water)
//                                    project SER regen, 6,259 rows
// Per-(residue, atom_type, method) calibration absorbs the small
// PBE-vs-OPBE offset; downstream code reads frame_type per record and
// routes accordingly. See project_serine_pbe_discontinuity.
//
// Extension over the gotham reference (TripeptideAssembler/Database):
// chi3 + chi4 columns are queried for the K/R (chi3 + chi4 set) and
// E/M/Q (chi3 set) families. The gotham impl queries only chi1/chi2
// and leaves chi3/chi4 data on the floor. This port extends the
// lookup at the SQL level so the table picks up the richer geometry
// when the parent residue has chi3/chi4 defined and non-zero.
//
// Connection-string format (libpq URI / kv pair string):
//   "host=/var/run/postgresql dbname=tensorcs15 user=jessica"
// from ~/.nmr_tools.toml [databases].tensorcs15 (Session::LoadFromToml
// reads this and constructs the table at process start).
//

#include "LarsenResidue.h"
#include "Types.h"

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

// Forward declare libpq connection so callers don't need <libpq-fe.h>.
struct pg_conn;
typedef struct pg_conn PGconn;

namespace nmr {


// One atom from a tripeptide DFT calculation. Positions are in the
// DFT-output frame (Gaussian standard orientation for OPBE rows, ORCA
// input orientation for PBE rows). The tensor lives in the same frame
// as the positions, so a single rotation aligns both onto the protein.
struct TripeptideDftAtom {
    int     atom_idx = 0;            // 1-based index in original Gaussian/ORCA geometry
    Element element  = Element::Unknown;
    Vec3    position = Vec3::Zero();  // Angstroms, DFT-output frame

    // Pre-decomposed shielding tensor (irreps already split at ingest).
    Mat3   shielding_tensor = Mat3::Zero();   // full asymmetric 3x3, ppm
    double isotropic        = 0.0;            // T0 / 3, ppm
    double anisotropy       = 0.0;            // Haeberlen Δσ, ppm

    // 5-component T2 (sphericart layout matches our SphericalTensor::T2).
    std::array<double, 5> t2_components = {};
};


// One row from raw_dft_calculations: a tripeptide at a (φ, ψ, χ*) pose.
struct TripeptideDftRecord {
    int  calc_id = 0;
    std::string tripeptide;          // "AAA", "AFA", "ASA", etc.

    // Pose angles (degrees, integer grid points)
    int phi = 0, psi = 0;
    int chi1 = 0, chi2 = 0, chi3 = 0, chi4 = 0;

    // DFT method discriminator. Stashed per record so downstream
    // calibration code can route SER (orca_input_orientation, PBE) and
    // the rest (gaussian_standard_orientation, OPBE) separately.
    std::string frame_type;

    int n_atoms = 0;
    std::vector<TripeptideDftAtom> atoms;

    // Perceived typed model of the 5-piece tripeptide. nullopt if
    // perception failed (logged via OperationLog at perception time);
    // calculators should decline records with no larsen.
    //
    // Populated by `TripeptideDftTable::QueryNearest` after a successful
    // row fetch. Replaces the positional `central_*` heuristic above
    // with typed-identity-driven access to backbone slots and sidechain
    // atoms.
    std::optional<LarsenTripeptide> larsen;

    // True iff calc_id != 0 (i.e. the query found a row).
    bool IsHit() const { return calc_id != 0; }
};


class TripeptideDftTable {
public:
    explicit TripeptideDftTable(const std::string& conn_str);
    ~TripeptideDftTable();

    TripeptideDftTable(const TripeptideDftTable&) = delete;
    TripeptideDftTable& operator=(const TripeptideDftTable&) = delete;

    // Query the nearest tripeptide for a central residue at the given
    // angles. residue_letter is the 1-letter code (Markley convention,
    // e.g. 'F' for PHE, 'S' for SER). Angles in degrees.
    //
    // Grid spacing: 2° for ALA (the AAA reference baseline), 20° for
    // every other residue. Angles are rounded to grid before query.
    //
    // chi1/chi2 are used when the residue has them (everything except
    // ALA, GLY, PRO). chi3 is used when the residue has it (E, M, Q,
    // K, R). chi4 is used only for K and R.
    //
    // Returns an empty record (calc_id == 0) on miss; the call site
    // logs the miss with the rounded grid coordinates. Fallback on a
    // chi-specific miss is delegated to the call site (the calculator
    // re-tries with fewer chi axes), per the Bundle B fail-loud
    // discipline — no silent degradation in this layer.
    //
    // n_chi_axes (default -1) overrides the default chi-axis count
    // for the residue. When ≥ 0, exactly that many chi columns are
    // included in the WHERE clause (and the chi values beyond the
    // count are ignored). The calculator-side fallback uses this to
    // walk from the residue's natural chi depth down to phi/psi-only
    // when chi-specific lookups miss.
    //
    // his_variant_hint (default -1) is forwarded to
    // PerceiveLarsenTripeptide and applies only when residue_letter ==
    // 'H'. Canonical values: 0=HID, 1=HIE, 2=HIP. Calculator sites
    // should pass the protein residue's `protonation_variant_index`
    // so perception locks onto the matching variant and does not
    // silently mis-assign HD1/HE2 atoms.
    TripeptideDftRecord QueryNearest(
        char residue_letter,
        double phi, double psi,
        double chi1 = 0.0, double chi2 = 0.0,
        double chi3 = 0.0, double chi4 = 0.0,
        int    n_chi_axes        = -1,
        int    his_variant_hint  = -1) const;

    // Connection health. Reports the libpq connection state at the
    // last touch; transient breaks may not surface until the next
    // QueryNearest call.
    bool IsConnected() const;

private:
    PGconn* conn_ = nullptr;
};


}  // namespace nmr
