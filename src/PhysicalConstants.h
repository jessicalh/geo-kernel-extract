#pragma once
// Fixed reference constants. Tuneable calculator parameters live in TOML.

#include "Types.h"
#include <cmath>

namespace nmr {

constexpr double PI = 3.14159265358979323846;
constexpr double SQRT_2_OVER_PI = 0.79788456080286535588;  // √(2/π)

// Degree/radian conversion.  Used when reading GROMACS TPR parameters
// which store angles in degrees (equilibrium angles for harmonic bonds,
// proper/improper dihedral phase angles).
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / PI;

// Bohr radius.  CODATA 2018: a_0 = 0.529177210903(80) A.
// We use the 2014 value (0.52917721067 A) which is the one coded in
// dftd4, xTB, and DFTB+.  The difference (4e-11 A) is below any
// physical significance for EEQ charges.
constexpr double ANGSTROM_PER_BOHR = 0.52917721067;
constexpr double BOHR_PER_ANGSTROM = 1.0 / ANGSTROM_PER_BOHR;

// 2019 CODATA mu_0. The Biot-Savart wire model uses the binary-stable
// pre-2019 prefactor below instead of deriving from this reference value.
constexpr double VACUUM_PERMEABILITY = 1.25663706212e-6;   // T*m/A (mu_0, 2019 CODATA)

constexpr double ANGSTROMS_TO_METRES = 1.0e-10;
constexpr double NANOAMPERES_TO_AMPERES = 1.0e-9;
constexpr double PPM_FACTOR = 1.0e6;

// GROMACS works in nanometres; this library works in Angstroms.
// 1 nm = 10 A exactly (SI definitions).
constexpr double ANGSTROM_PER_NANOMETRE = 10.0;
constexpr double NANOMETRE_PER_ANGSTROM = 0.1;

// Electrostatics in {e, Angstrom, eV} units
// ke = e / (4 pi epsilon_0) = 14.3996 eV*A / e = 14.3996 V*A
// Converts E from e/A^2 (raw Coulomb sum) to V/A (physical E-field).
constexpr double COULOMB_KE = 14.3996;

// Thermal voltage at 298.15 K: kT/e = k_B * T / e = 0.025693 V
// Converts APBS potential/field from kT/e units to Volts.
constexpr double KT_OVER_E_298K = 0.025693;

// AMBER prmtop %FLAG CHARGE stores charges in internal units of
// e * sqrt(kcal/mol * A) = e * 18.2223; divide by this to recover
// elementary charge. (AMBER convention; mirrored in
// tools/amber/generate_ff14sb_pb_table.py.)
constexpr double AMBER_PRMTOP_CHARGE_FACTOR = 18.2223;

// Biot-Savart prefactor: mu_0/(4*pi) in SI units (T*m/A).
// Pre-2019 SI: exactly 1e-7. This is the value used in the JB wire model.
// The 2019 CODATA derivation (VACUUM_PERMEABILITY / 4*pi) gives ~1.00000000005e-7,
// which is NOT what the computation uses. Do not substitute.
constexpr double BIOT_SAVART_PREFACTOR = 1e-7;

// Defaults and legacy thresholds used by callers that do not read
// CalculatorConfig.
constexpr double MIN_DISTANCE = 0.1;            // Angstroms -- singularity cutoff
constexpr double NO_DATA_SENTINEL = 99.0;       // sentinel for missing data
constexpr double NEAR_ZERO_NORM = 1e-10;        // near-zero vector norm
constexpr double NEAR_ZERO_FIELD = 1e-15;       // near-zero field magnitude

constexpr double RING_COUNT_SHELL_1 = 3.0;      // Angstroms
constexpr double RING_COUNT_SHELL_2 = 5.0;
constexpr double RING_COUNT_SHELL_3 = 8.0;
constexpr double RING_COUNT_SHELL_4 = 12.0;

constexpr double RING_CALC_CUTOFF = 15.0;       // Angstroms -- ring current cutoff
constexpr double EXP_DECAY_LENGTH = 4.0;        // Angstroms
constexpr double PACKING_RADIUS = 8.0;          // Angstroms -- for heavy atom count

constexpr double HBOND_COUNT_RADIUS = 3.5;      // Angstroms
constexpr double HBOND_MAX_DIST = 50.0;         // Angstroms
constexpr double APBS_SANITY_LIMIT = 100.0;     // V/Angstrom

// Note: dispersion distance thresholds (R_SWITCH=4.3A, R_CUT=5.0A) are
// local to DispersionResult.cpp because they are specific to that
// CHARMM-form switching function.

constexpr int SEQUENTIAL_EXCLUSION_THRESHOLD = 2;

// Reference bond lengths (Angstroms)
//
// Allen, F. H., Kennard, O., Watson, D. G., Brammer, L., Orpen, A. G.,
// & Taylor, R. (1987). "Tables of bond lengths determined by X-ray and
// neutron diffraction. Part 1. Bond lengths in organic compounds."
// J. Chem. Soc., Perkin Trans. 2, S1-S19. DOI: 10.1039/P298700000S1.
//
// These are reference geometric scales used by synthetic test fixtures
// to place atoms at chemically realistic distances. Calculators don't
// consume them — bond detection is geometric (covalent-radius sum +
// tolerance) — but the test fixtures need defensible numbers rather
// than literal magic constants.

// Aromatic C-C bond length (benzene, X-ray + neutron average).
// Allen 1987 Table 1: aromatic C(sp2)-C(sp2) bond is 1.395 +- 0.003 A.
constexpr double BENZENE_CC_BOND_LENGTH = 1.40;

// Backbone peptide bond reference scales. Engh & Huber, Acta Cryst.
// A47, 392-400 (1991). DOI: 10.1107/S0108767391001071. The Engh-Huber
// values are the de facto standard for peptide ideal geometries.
constexpr double PEPTIDE_N_CA_BOND_LENGTH    = 1.458;  // A; Engh-Huber 1991 Table 1
constexpr double PEPTIDE_CA_C_BOND_LENGTH    = 1.525;  // A; Engh-Huber 1991 Table 1
constexpr double PEPTIDE_C_O_DOUBLE_BOND_LENGTH = 1.231;  // A; Engh-Huber 1991 Table 1

// Bondi van der Waals radii (Angstroms)
//
// Bondi, A. J. Phys. Chem. 68, 441-451 (1964).
// Used by SasaResult (Shrake-Rupley SASA).

inline double BondiVdwRadius(Element e) {
    switch (e) {
        case Element::H: return 1.20;
        case Element::C: return 1.70;
        case Element::N: return 1.55;
        case Element::O: return 1.52;
        case Element::S: return 1.80;
        default:         return 1.70;
    }
}

// D4 EEQ element parameters (atomic units: Hartree, Bohr)
//
// Caldeweyher, Ehlert, Hansen, Neugebauer, Spicher, Bannwarth & Grimme,
// J. Chem. Phys. 150, 154122 (2019).  DOI: 10.1063/1.5090222.
// Reference implementation: github.com/dftd4/dftd4 (Apache-2.0),
// src/dftd4/data/{en,hardness,rcov,rad}.f90.
//
// Fitted to reproduce Hirshfeld charges from DFT/def2-TZVP.
// Same parameters used in TURBOMOLE, ORCA D4, xTB, DFTB+.
//
// chi   — electronegativity [Hartree]
// gam   — chemical hardness (Hubbard U) [Hartree]
// kappa — CN-dependent electronegativity shift [Hartree]
// rcov  — covalent radius for CN counting [Bohr] (Pyykko 2009)
// rad   — Gaussian charge width for diagonal correction [Bohr]

struct D4EeqParams {
    double chi;
    double gam;
    double kappa;
    double rcov;
    double rad;
};

inline D4EeqParams D4EeqParamsFor(Element e) {
    switch (e) {
        //                      chi           gam          kappa         rcov      rad
        case Element::H:  return {-0.35015861, 0.47259288, -0.19793756, 0.80628, 1.61478};
        case Element::C:  return {-0.04726052, 0.25364654,  0.14216971, 1.51718, 2.49988};
        case Element::N:  return { 0.11527249, 0.28022740,  0.15169154, 1.42165, 2.23456};
        case Element::O:  return { 0.25136810, 0.36515829,  0.14510449, 1.24854, 1.89247};
        case Element::S:  return { 0.10789083, 0.25140725,  0.15916035, 2.00000, 3.29733};
        default:          return { 0.0,        0.30000000,  0.0,        1.50000, 2.50000};
    }
}

// Karplus 3J-coupling parameters
//
// Karplus form for the eight protein backbone / sidechain 3J coupling
// families emitted by JCouplingTimeSeriesTrajectoryResult (nine HDF5
// datasets because J(Halpha,Hbeta) is split into HB2 and HB3):
//
//     3J(alpha) = A * cos^2(alpha) + B * cos(alpha) + C   [Hz]
//
// where alpha is either (library_phi + theta_offset) for
// backbone phi channels, or the actual atomic 4-atom dihedral for
// chi1 channels.
//
// The published backbone Karplus papers express their fit as
//   J(phi_pub) = A * cos^2(phi_pub + theta_pub) +
//                B * cos(phi_pub + theta_pub) + C
// This library uses the opposite phi sign from the Wang-Bax/Vogeli
// plotting convention:
//
//   phi_library = -phi_pub
//   cos(phi_pub + theta_pub) = cos(phi_library - theta_pub)
//
// so stored backbone theta values are -theta_pub. Perez chi1-channel
// coefficients already include the relevant substituent offsets; those
// channels use the actual 4-atom dihedral with theta = 0.

// 3J(HN,Halpha), Vuister & Bax 1993.
// Wang-Bax theta_pub=-60 deg maps to library theta=+60 deg.
constexpr double KARPLUS_HN_HA_A =  6.51;
constexpr double KARPLUS_HN_HA_B = -1.76;
constexpr double KARPLUS_HN_HA_C =  1.60;
constexpr double KARPLUS_HN_HA_THETA =  PI / 3.0;  // library theta; Wang-Bax theta_pub=-60 deg

// 3J(HN,Halpha), Vogeli et al. 2007 rigid row; same H-N-CA-HA dihedral.
constexpr double KARPLUS_HN_HA_VOGELI_A =  7.97;
constexpr double KARPLUS_HN_HA_VOGELI_B = -1.26;
constexpr double KARPLUS_HN_HA_VOGELI_C =  0.63;
constexpr double KARPLUS_HN_HA_VOGELI_THETA = PI / 3.0;  // library theta;
        // Vogeli eq 5 eta_ik=-pi/3 in the published phi convention.

// 3J(HN,Cbeta), Wang & Bax 1996 Table 1 row 3.
// theta_pub=+60 deg maps to library theta=-60 deg.
constexpr double KARPLUS_HN_CB_A =  3.39;
constexpr double KARPLUS_HN_CB_B = -0.94;
constexpr double KARPLUS_HN_CB_C =  0.07;
constexpr double KARPLUS_HN_CB_THETA = -PI / 3.0;  // library theta;
        // Wang-Bax row 3 / Vogeli eta_ik=+pi/3 in the published convention.

// 3J(HN,C'), Wang & Bax 1996 Table 1 row 4. B is positive here.
// Vogeli eta=pi is the equivalent opposite-B parametrization.
constexpr double KARPLUS_HN_CP_A =  4.32;
constexpr double KARPLUS_HN_CP_B =  0.84;
constexpr double KARPLUS_HN_CP_C =  0.00;
constexpr double KARPLUS_HN_CP_THETA = 0.0;  // Wang-Bax row 4 (theta=0 deg).
        // Vogeli 2007 eq 5 gives eta_ik = pi for the same observable; the
        // two parametrizations are equivalent (cos(phi + pi) = -cos(phi) so
        // B flips sign) -- we ship the Wang-Bax theta=0 / B>0 form.

// 3J(Halpha,C'), Wang & Bax 1996 Table 1 row 2. The phi path is
// Halpha(i)-CA(i)-N(i)-C'(i-1); B is positive.
constexpr double KARPLUS_HA_CP_A =  3.75;
constexpr double KARPLUS_HA_CP_B =  2.19;
constexpr double KARPLUS_HA_CP_C =  1.28;
constexpr double KARPLUS_HA_CP_THETA = PI / 3.0;  // library theta; Wang-Bax row 2 theta_pub=-60 deg.

// 3J(N,Cgamma), Perez et al. 2001 Table 2 consensus row.
constexpr double KARPLUS_N_CG_A  =  1.29;
constexpr double KARPLUS_N_CG_B  = -0.49;
constexpr double KARPLUS_N_CG_C  =  0.37;
constexpr double KARPLUS_N_CG_THETA = 0.0;  // Perez 2001 uses chi1 = N-CA-CB-CG
        // directly; no phi-style offset.

// 3J(C',Cgamma), Perez et al. 2001 Table 2 consensus row.
// The coefficients include the C'-on-CA substituent offset.
constexpr double KARPLUS_CP_CG_A =  2.31;
constexpr double KARPLUS_CP_CG_B = -0.87;
constexpr double KARPLUS_CP_CG_C =  0.55;
constexpr double KARPLUS_CP_CG_THETA = 0.0;  // Perez 2001 internalizes the
        // C'-on-CA substituent offset (chi1+120 deg) in the per-coupling
        // (A, B, C); feeding the atomic C-CA-CB-CG dihedral matches the
        // Table 2 consensus row directly. See Table 2 footnote c.

// 3J(Halpha,Hbeta), Perez et al. 2001 Table 2 consensus row.
// Hbeta availability and HB2/HB3 splitting are handled by
// JCouplingTimeSeriesTrajectoryResult.
constexpr double KARPLUS_HA_HB_A =  7.23;
constexpr double KARPLUS_HA_HB_B = -1.37;
constexpr double KARPLUS_HA_HB_C =  2.22;
constexpr double KARPLUS_HA_HB_THETA = 0.0;  // Perez 2001 Table 2 footnote c.
        // The atomic dihedral HA-CA-CB-HB{2,3} differs from chi1 by the
        // Halpha and Hbeta substituent offsets (chi1 + Delta_chi1 ~ ±120°),
        // but the per-coupling (A, B, C) internalize Delta_chi1 -- feeding
        // the atomic dihedral matches the Table 2 consensus row directly.


}  // namespace nmr
