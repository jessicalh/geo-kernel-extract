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
// ============================================================================

// Aromatic C-C bond length (benzene, X-ray + neutron average).
// Allen 1987 Table 1: aromatic C(sp2)-C(sp2) bond is 1.395 +- 0.003 A.
constexpr double BENZENE_CC_BOND_LENGTH = 1.40;

// Backbone peptide bond reference scales. Engh & Huber, Acta Cryst.
// A47, 392-400 (1991). DOI: 10.1107/S0108767391001071. The Engh-Huber
// values are the de facto standard for peptide ideal geometries.
constexpr double PEPTIDE_N_CA_BOND_LENGTH    = 1.458;  // A; Engh-Huber 1991 Table 1
constexpr double PEPTIDE_CA_C_BOND_LENGTH    = 1.525;  // A; Engh-Huber 1991 Table 1
constexpr double PEPTIDE_C_O_DOUBLE_BOND_LENGTH = 1.231;  // A; Engh-Huber 1991 Table 1


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


// ============================================================================
// Karplus 3J-coupling parameters
//
// Karplus form for all four protein backbone / sidechain 3J couplings
// emitted by JCouplingTimeSeriesTrajectoryResult:
//
//     3J(theta) = A * cos^2(theta) + B * cos(theta) + C   [Hz]
//
// where theta is the relevant 3-bond H-X-Y-Z dihedral in RADIANS,
// measured via IUPAC signed atan2 directly from atom positions.
//
// Geometric identity used by the backbone phi observables:
// for L-amino acids in standard tetrahedral backbone geometry,
//   H-N-CA-HA  ~=  phi - 60 degrees   (basis for 3J(HN, Halpha))
//   H-N-CA-CB  ~=  phi + 60 degrees   (basis for 3J(HN, Cbeta))
// The Karplus parametrizations below were originally published as a
// function of (phi - 60) or (phi + 60); equivalently they take the
// actual atomic dihedral as input. We compute the dihedral directly
// from atom positions (never reconstructed via phi + offset), which
// is the modern usage explicitly documented by Li, Lee, Grishaev,
// Ying & Bax, ChemPhysChem 16, 572-578 (2015), DOI: 10.1002/cphc.
// 201402704.
// ============================================================================

// --- 3J(HN, Halpha) -- phi observable via H-N-CA-HA dihedral ---
//
// Vuister, G. W. & Bax, A. "Quantitative J correlation: a new approach
// for measuring homonuclear three-bond J(HNHalpha) coupling constants
// in 15N-enriched proteins." J. Am. Chem. Soc. 115, 7772-7777 (1993).
// DOI: 10.1021/ja00070a024.
//
// Most-cited Karplus parametrization for backbone phi from 3J(HN,Halpha).
// Note: Wang & Bax 1996 (DOI 10.1021/ja9535524) later reparametrized
// these to (6.98, -1.38, 1.72) on a refined ubiquitin NMR/X-ray phi
// set (Table 1 row 1); both are valid, the 1993 values are the
// thesis-default canonical reference. Reference PDF:
//   references/vuister-lecture-j-couplings.pdf (the Vuister teaching
//   lecture quotes A,B,C = 6.51, -1.76, 1.60 verbatim).
constexpr double KARPLUS_HN_HA_A =  6.51;
constexpr double KARPLUS_HN_HA_B = -1.76;
constexpr double KARPLUS_HN_HA_C =  1.60;

// --- 3J(HN, Cbeta) -- phi observable via H-N-CA-CB dihedral ---
//
// Wang, A. C. & Bax, A. "Determination of the backbone dihedral angles
// phi in human ubiquitin from reparametrized empirical Karplus
// equations." J. Am. Chem. Soc. 118, 2483-2494 (1996).
// DOI: 10.1021/ja9535524.
//
// Table 1 row 3, NMR/X-ray refined fit (page 2487):
//   theta = +60 degrees, A = 3.39 +- 0.07, B = -0.94 +- 0.07, C = 0.07
//   +- 0.02. Original published form is J = A * cos^2(phi + 60) + B *
//   cos(phi + 60) + C, where (phi + 60) ~= H-N-CA-CB dihedral in
//   standard L-amino-acid geometry; feeding the atomic dihedral
//   directly is equivalent and is the standard modern usage.
//
// Reference PDF: references/wang-bax-1996-karplus-phi-ubiquitin.pdf
// (Table 1, page 2487; byte-verified 2026-05-19 against the open
// Bax-group repository copy).
constexpr double KARPLUS_HN_CB_A =  3.39;
constexpr double KARPLUS_HN_CB_B = -0.94;
constexpr double KARPLUS_HN_CB_C =  0.07;

// --- 3J(N, Cgamma) -- chi1 observable via N-CA-CB-CG dihedral (= chi1) ---
//
// Perez, C., Lohr, F., Ruterjans, H. & Schmidt, J. M. "Self-consistent
// Karplus parametrization of 3J couplings depending on the polypeptide
// side-chain torsion chi1." J. Am. Chem. Soc. 123, 7081-7093 (2001).
// DOI: 10.1021/ja003724j.
//
// Widely-cited self-consistent Karplus fit for chi1 observables. The
// original paper is paywalled at ACS; the agent-level literature audit
// on 2026-05-19 confirmed the citation but could NOT byte-verify the
// coefficients against the published Table -- the values are circulated
// unchanged in downstream NMR software (TALOS-N, NMRViewJ) that cite
// Perez 2001. Byte-verification is pending institutional access; if
// the user obtains a copy and the values differ, update here and the
// JCouplingTimeSeriesTrajectoryResult attrs in the same commit.
constexpr double KARPLUS_N_CG_A  =  1.29;
constexpr double KARPLUS_N_CG_B  = -0.49;
constexpr double KARPLUS_N_CG_C  =  0.37;

// --- 3J(C', Cgamma) -- chi1 observable via C-CA-CB-CG dihedral ---
//
// Same paper: Perez, Lohr, Ruterjans & Schmidt, JACS 123:7081 (2001).
// DOI: 10.1021/ja003724j. Same paywall caveat as 3J(N, Cgamma).
//
// (Strictly: the C'-CA-CB-Cgamma dihedral differs from chi1 by ~120
// degrees around CA; the Perez self-consistent fit treats chi1 as
// the input variable with the substituent-position bookkeeping
// internalized in the per-coupling coefficients, so feeding the
// C-CA-CB-CG atomic dihedral directly is the correct modern usage.)
constexpr double KARPLUS_CP_CG_A =  1.74;
constexpr double KARPLUS_CP_CG_B = -0.57;
constexpr double KARPLUS_CP_CG_C =  0.25;


}  // namespace nmr
