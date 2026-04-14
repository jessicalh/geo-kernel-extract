#pragma once
//
// PhysicalConstants.h — citable constants used by the calculation engine.
//
// Everything here has a literature source.  If a number comes from a
// paper, textbook, or standard and is compiled into a calculator, it
// lives here so the thesis can cite it and the examiner can audit it.
//
// Model-tuneable parameters (cutoffs, intensities) go in TOML via
// CalculatorConfig.  This file is for reference data that is fixed
// by the publication it came from.
//

#include "Types.h"
#include <cmath>

namespace nmr {

// ============================================================================
// Mathematical constants
// ============================================================================

constexpr double PI = 3.14159265358979323846;
constexpr double SQRT_2_OVER_PI = 0.79788456080286535588;  // √(2/π)

// Degree/radian conversion.  Used when reading GROMACS TPR parameters
// which store angles in degrees (equilibrium angles for harmonic bonds,
// proper/improper dihedral phase angles).
constexpr double DEG_TO_RAD = PI / 180.0;
constexpr double RAD_TO_DEG = 180.0 / PI;

// ============================================================================
// Unit conversions
// ============================================================================

// Bohr radius.  CODATA 2018: a_0 = 0.529177210903(80) A.
// We use the 2014 value (0.52917721067 A) which is the one coded in
// dftd4, xTB, and DFTB+.  The difference (4e-11 A) is below any
// physical significance for EEQ charges.
constexpr double ANGSTROM_PER_BOHR = 0.52917721067;
constexpr double BOHR_PER_ANGSTROM = 1.0 / ANGSTROM_PER_BOHR;


// ============================================================================
// Electromagnetic (SI)
// mu_0 = 4*pi*1e-7 T*m/A exactly in pre-2019 SI.
// 2019 CODATA redefined mu_0 = 1.25663706212e-6 T*m/A (measured, not exact).
// The Johnson-Bovey wire model in BiotSavartResult uses the pre-2019
// exact value. Changing this breaks binary reproduction of existing results.
constexpr double VACUUM_PERMEABILITY = 1.25663706212e-6;   // T*m/A (mu_0, 2019 CODATA)

// Unit conversions
constexpr double ANGSTROMS_TO_METRES = 1.0e-10;
constexpr double NANOAMPERES_TO_AMPERES = 1.0e-9;
constexpr double PPM_FACTOR = 1.0e6;

// Electrostatics in {e, Angstrom, eV} units
// ke = e / (4 pi epsilon_0) = 14.3996 eV*A / e = 14.3996 V*A
// Converts E from e/A^2 (raw Coulomb sum) to V/A (physical E-field).
constexpr double COULOMB_KE = 14.3996;

// Thermal voltage at 298.15 K: kT/e = k_B * T / e = 0.025693 V
// Converts APBS potential/field from kT/e units to Volts.
constexpr double KT_OVER_E_298K = 0.025693;

// Biot-Savart prefactor: mu_0/(4*pi) in SI units (T*m/A).
// Pre-2019 SI: exactly 1e-7. This is the value used in the JB wire model.
// The 2019 CODATA derivation (VACUUM_PERMEABILITY / 4*pi) gives ~1.00000000005e-7,
// which is NOT what the computation uses. Do not substitute.
constexpr double BIOT_SAVART_PREFACTOR = 1e-7;

// Constitution: numerical thresholds
// NOTE: Calculators read these from CalculatorConfig::Get() (TOML-configurable).
// These constexpr values remain for non-calculator consumers
// (ApbsFieldResult, MutationDeltaResult, tests).
constexpr double MIN_DISTANCE = 0.1;            // Angstroms -- singularity cutoff
constexpr double NO_DATA_SENTINEL = 99.0;       // sentinel for missing data
constexpr double NEAR_ZERO_NORM = 1e-10;        // near-zero vector norm
constexpr double NEAR_ZERO_FIELD = 1e-15;       // near-zero field magnitude

// Constitution: spatial shells for ring counting
constexpr double RING_COUNT_SHELL_1 = 3.0;      // Angstroms
constexpr double RING_COUNT_SHELL_2 = 5.0;
constexpr double RING_COUNT_SHELL_3 = 8.0;
constexpr double RING_COUNT_SHELL_4 = 12.0;

// Constitution: calculation cutoffs
constexpr double RING_CALC_CUTOFF = 15.0;       // Angstroms -- ring current cutoff
constexpr double EXP_DECAY_LENGTH = 4.0;        // Angstroms
constexpr double PACKING_RADIUS = 8.0;          // Angstroms -- for heavy atom count

// Constitution: H-bond thresholds
constexpr double HBOND_COUNT_RADIUS = 3.5;      // Angstroms
constexpr double HBOND_MAX_DIST = 50.0;         // Angstroms
constexpr double APBS_SANITY_LIMIT = 100.0;     // V/Angstrom

// Note: dispersion distance thresholds (R_SWITCH=4.3A, R_CUT=5.0A) are
// defined in DispersionResult.cpp with full physics documentation.
// They are NOT global constants because they are specific to the
// dispersion switching function (CHARMM form, Brooks et al. 1983).

// Constitution: sequence exclusion
constexpr int SEQUENTIAL_EXCLUSION_THRESHOLD = 2;


// ============================================================================
// Bondi van der Waals radii (Angstroms)
//
// Bondi, A. J. Phys. Chem. 68, 441-451 (1964).
// Used by SasaResult (Shrake-Rupley SASA).
// ============================================================================

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


// ============================================================================
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
// ============================================================================

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


}  // namespace nmr
