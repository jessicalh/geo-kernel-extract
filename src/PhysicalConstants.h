#pragma once
//
// PhysicalConstants.h — fixed constants and legacy thresholds used by the
// calculation engine.
//
// Literature-sourced constants live here so the thesis can cite them and
// the examiner can audit them. A few compatibility thresholds also remain
// here for callers that have not moved to CalculatorConfig.
//
// Most model-tuneable parameters (cutoffs, intensities) go in TOML via
// CalculatorConfig.
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
// BiotSavartResult uses the pre-2019 exact prefactor below for the
// Johnson-Bovey wire model.
constexpr double VACUUM_PERMEABILITY = 1.25663706212e-6;   // T*m/A (mu_0, 2019 CODATA)

// Unit conversions
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

// Constitution: numerical thresholds
// Many calculator thresholds are now read from CalculatorConfig::Get()
// (TOML-configurable). These constexpr values remain for legacy or
// non-configured consumers.
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
constexpr double APBS_SANITY_LIMIT = 100.0;     // APBS-native kT/(e*Angstrom)

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
// probable source: Caldeweyher, Ehlert, Hansen, Neugebauer, Spicher,
// Bannwarth & Grimme, J. Chem. Phys. 150, 154122 (2019).
// DOI: 10.1063/1.5090222.
//
// Project-local EEQ parameter table consumed by EeqResult.
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
// Karplus form for the eight protein backbone / sidechain 3J coupling
// families emitted by JCouplingTimeSeriesTrajectoryResult (nine HDF5
// datasets because J(Halpha,Hbeta) is split into HB2 and HB3):
//
//     3J(alpha) = A * cos^2(alpha) + B * cos(alpha) + C   [Hz]
//
// where alpha is either (project_phi + project_theta_offset) for
// backbone phi channels, or the actual atomic 4-atom dihedral for
// chi1 channels.
//
// Convention -- per-channel, NOT a blanket claim:
//
// The published backbone Karplus papers express their fit as
//   J(phi_pub) = A * cos^2(phi_pub + theta_pub) +
//                B * cos(phi_pub + theta_pub) + C
// where theta_pub is a per-channel constant tabulated in the source:
//   - Wang & Bax 1996 JACS 118:2483 Table 1 (page 2487) tabulates
//     theta_pub for the four backbone channels: -60 deg for
//     3J(HN, Hα), +60 deg for 3J(HN, Cβ), -60 deg for 3J(Hα, C'),
//     0 deg for 3J(HN, C'). These correspond to the published
//     ideal-tetrahedral identity (alpha_atomic = phi + theta_offset).
//   - Vogeli et al. 2007 JACS 129:9377 eq 5 uses the same form,
//     restated with eta_ik in radians: -pi/3 for HN-Hα, +pi/3 for
//     HN-Cβ, pi for HN-C' (page 9383). Note Vogeli's eta = pi gives
//     B with the OPPOSITE sign to Wang-Bax's theta = 0 deg for the
//     same observable, because cos(phi + pi) = -cos(phi). Both
//     forms describe the same curve.
//   - Perez et al. 2001 JACS 123:7081 Table 2 (page 7086) tabulates
//     the chi1-related channels with the substituent rotation
//     book-kept inside the per-coupling (A, B, C); see footnote c
//     ("theta = chi1 + Delta_chi1, according to Figure 1, and
//     coefficients C0 already comprise incremental effects"). The
//     atomic dihedral on the relevant 4-atom path IS theta in
//     Perez's framework.
//   - Li, Lee, Grishaev, Ying & Bax 2015 ChemPhysChem 16:572 documents
//     the modern usage of evaluating the published curves on atomic
//     dihedrals computed directly from coordinates when coordinates
//     are available.
//
// The project DihedralTimeSeries uses the opposite sign from DSSP and
// from the Wang-Bax/Vogeli plotting convention (see
// DihedralTimeSeriesTrajectoryResult.cpp: phi_DSSP = -phi_IUPAC).
// Therefore the constants below store theta in the PROJECT convention:
//
//   phi_project = -phi_pub
//   cos(phi_pub + theta_pub) = cos(phi_project - theta_pub)
//
// so theta_project = -theta_pub. This is why HN-Halpha and
// Halpha-C' use +pi/3 here even though the Wang-Bax table prints
// theta=-60 deg, and HN-Cbeta uses -pi/3 even though the table prints
// theta=+60 deg. The 1UBQ probe confirms the sign convention:
// atomic H-N-CA-HA is phi_project + ~60.8 deg and atomic H-N-CA-CB is
// phi_project - ~57.5 deg.
//
// **For chi1 channels** (J(N,Cγ), J(C',Cγ), J(Hα,Hβ)) the parametriz-
// ation in Perez 2001 implicitly carries the rotation around C-alpha
// inside (A, B, C); feeding the atomic dihedral on the relevant
// 4-atom path is the Perez-intended usage. Channels referring to
// 3J(N,Cγ): chi1 terminal must be C; SER/CYS/THR carry O/S/O at
// chi[0].a[3] and the channel is NaN by element gate (Element::C
// at chi[0].a[3]).
//
// All channel families carry POSITIVE A (cos^2 amplitude). The B sign
// is channel-dependent. Karplus arithmetic bounds on
// J(alpha) = A * cos^2(alpha) + B * cos(alpha) + C over cos(alpha) in
// [-1, 1] (A > 0):
//   global extrema are at the endpoints f(+1) = A + B + C and f(-1)
//   = A - B + C, and at the vertex u* = -B/(2A) if u* lies in
//   (-1, 1), where f(u*) = C - B^2 / (4A). For B < 0 the vertex is
//   the MIN and f(-1) = A + |B| + C is the MAX. For B > 0 the vertex
//   is the MIN and f(+1) = A + B + C is the MAX.
// ============================================================================

// --- 3J(HN, Halpha) -- phi observable ---
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
// Wang-Bax fit form (published convention) has theta_pub=-60 deg.
// The project phi sign is opposite, so theta_project=+60 deg.
// Arithmetic range: cos(theta) in [-1, 1], A>0, B<0; max at f(-1) =
// A - B + C = 9.87 Hz; vertex u* = -B/(2A) = 0.135 in (-1, 1), f(u*)
// = C - B^2/(4A) = 1.48 Hz (MIN).
constexpr double KARPLUS_HN_HA_A =  6.51;
constexpr double KARPLUS_HN_HA_B = -1.76;
constexpr double KARPLUS_HN_HA_C =  1.60;
constexpr double KARPLUS_HN_HA_THETA =  PI / 3.0;  // project theta; Wang-Bax theta_pub=-60 deg

// --- 3J(HN, Halpha) Vogeli rigid -- phi observable via H-N-CA-HA dihedral ---
//
// Vogeli, B., Ying, J., Grishaev, A. & Bax, A. "Limits on variations
// in protein backbone dynamics from precise measurements of scalar
// couplings." J. Am. Chem. Soc. 129, 9377-9385 (2007).
// DOI: 10.1021/ja070324o.
//
// Table 1 (page 9383) "rigid" row: A=7.97, B=-1.26, C=0.63 with rmsd
// 0.42 Hz. Fit on GB3 RDC-refined structure (PDB 2OED) with hydrogens
// at idealized positions. Methods-accumulate alternate to Vuister-Bax
// 1993 (NOT a replacement -- see feedback_methods_accumulate). Same
// atomic dihedral as J_HN_Halpha (H-N-CA-HA).
// Reference PDF:
//   references/vogeli-2007-limits-backbone-dynamics-3j-couplings-gb3.pdf
//   Table 1 page 9383, byte-verified 2026-05-19.
// Arithmetic range: A>0, B<0; max at f(-1) = A - B + C = 9.86 Hz;
// vertex u* = -B/(2A) = 0.079 in (-1, 1), f(u*) = C - B^2/(4A) = 0.58
// Hz (MIN).
constexpr double KARPLUS_HN_HA_VOGELI_A =  7.97;
constexpr double KARPLUS_HN_HA_VOGELI_B = -1.26;
constexpr double KARPLUS_HN_HA_VOGELI_C =  0.63;
constexpr double KARPLUS_HN_HA_VOGELI_THETA = PI / 3.0;  // project theta;
        // Vogeli eq 5 eta_ik=-pi/3 in the published phi convention.

// --- 3J(HN, Cbeta) -- phi observable via H-N-CA-CB dihedral ---
//
// Wang, A. C. & Bax, A. "Determination of the backbone dihedral angles
// phi in human ubiquitin from reparametrized empirical Karplus
// equations." J. Am. Chem. Soc. 118, 2483-2494 (1996).
// DOI: 10.1021/ja9535524.
//
// Table 1 row 3, NMR/X-ray refined fit (page 2487):
//   theta = +60 degrees, A = 3.39 +- 0.07, B = -0.94 +- 0.07, C = 0.07
//   +- 0.02. Reference PDF:
//   references/wang-bax-1996-karplus-phi-ubiquitin.pdf (Table 1, page
//   2487; byte-verified 2026-05-19 against the open Bax-group PDF).
// Arithmetic range: A>0, B<0; max at f(-1) = A - B + C = 4.40 Hz;
// vertex u* = -B/(2A) = 0.139 in (-1, 1), f(u*) = C - B^2/(4A) =
// 0.005 Hz (MIN).
constexpr double KARPLUS_HN_CB_A =  3.39;
constexpr double KARPLUS_HN_CB_B = -0.94;
constexpr double KARPLUS_HN_CB_C =  0.07;
constexpr double KARPLUS_HN_CB_THETA = -PI / 3.0;  // project theta;
        // Wang-Bax row 3 / Vogeli eta_ik=+pi/3 in the published convention.

// --- 3J(HN, C') -- phi observable via H-N-CA-C dihedral ---
//
// Same paper: Wang & Bax, JACS 118:2483 (1996), DOI 10.1021/ja9535524.
// Table 1 row 4 (NMR/X-ray refined fit, page 2487):
//   theta = 0 deg, A = 4.32, B = +0.84, C = 0.00.
// IMPORTANT: B is POSITIVE for this coupling -- the bound derivation
// is different from A>0/B<0 channels. Reference PDF:
//   references/wang-bax-1996-karplus-phi-ubiquitin.pdf Table 1 page
//   2487, byte-verified 2026-05-19. The four-row mapping in Wang-Bax
//   Table 1 (page 2487 leftmost column, NMR/X-ray refined fit rows
//   in italics) is:
//     row 1 (theta=-60 deg): 3J(HN, Halpha)  -> (6.98, -1.38, 1.72)
//     row 2 (theta=-60 deg): 3J(Halpha, C')  -> (3.75, +2.19, 1.28)
//     row 3 (theta=+60 deg): 3J(HN, Cbeta)   -> (3.39, -0.94, 0.07)
//     row 4 (theta=  0 deg): 3J(HN, C')      -> (4.32, +0.84, 0.00)
//   The "Measurement of 3J_{HαC'}" text section heading directly
//   below the table (page 2487 right column) confirms row 2 is the
//   Halpha-C' channel. Vogeli 2007 page 9383 gives eta_ik = pi for
//   3J(HN, C') (B sign flips relative to Wang-Bax theta=0 because
//   cos(phi+pi) = -cos(phi); same curve, different parametrization).
// Arithmetic range: A>0, B>0; max at f(+1) = A + B + C = 5.16 Hz;
// vertex u* = -B/(2A) = -0.097 in (-1, 1), f(u*) = C - B^2/(4A) =
// -0.041 Hz (MIN, slightly negative -- physical, J can be negative).
constexpr double KARPLUS_HN_CP_A =  4.32;
constexpr double KARPLUS_HN_CP_B =  0.84;
constexpr double KARPLUS_HN_CP_C =  0.00;
constexpr double KARPLUS_HN_CP_THETA = 0.0;  // Wang-Bax row 4 (theta=0 deg).
        // Vogeli 2007 eq 5 gives eta_ik = pi for the same observable; the
        // two parametrizations are equivalent (cos(phi + pi) = -cos(phi) so
        // B flips sign) -- we ship the Wang-Bax theta=0 / B>0 form.

// --- 3J(Halpha, C') -- phi observable via Halpha-CA-N-C'(prev) dihedral ---
//
// Same paper: Wang & Bax, JACS 118:2483 (1996), DOI 10.1021/ja9535524.
// Table 1 row 2 (NMR/X-ray refined fit, page 2487):
//   theta = -60 deg (Wang-Bax sign convention), A = 3.75, B = +2.19,
//   C = 1.28.
// 4-atom path: there is only ONE 3-bond path from Halpha(i) to
// C'(i-1): Halpha(i) - CA(i) - N(i) - C'(i-1). The rotation is
// around the N-CA axis (phi axis), so this is a phi observable.
// The HA(i)-CA(i)-C(i)-N(i+1) path rotates around CA-C (psi axis),
// not phi. Vuister teaching lecture section 6.1 and Vogeli 2007
// page 9384 list 3J(C'(i-1), Halpha) among the six phi-related
// couplings.
// IMPORTANT: B is POSITIVE for this coupling -- bound derivation
// differs from A>0/B<0 channels. Reference PDF:
//   references/wang-bax-1996-karplus-phi-ubiquitin.pdf Table 1 page
//   2487, byte-verified 2026-05-19.
// Arithmetic range: A>0, B>0; max at f(+1) = A + B + C = 7.22 Hz;
// vertex u* = -B/(2A) = -0.292 in (-1, 1), f(u*) = C - B^2/(4A) =
// 0.96 Hz (MIN).
constexpr double KARPLUS_HA_CP_A =  3.75;
constexpr double KARPLUS_HA_CP_B =  2.19;
constexpr double KARPLUS_HA_CP_C =  1.28;
constexpr double KARPLUS_HA_CP_THETA = PI / 3.0;  // project theta; Wang-Bax row 2 theta_pub=-60 deg.

// --- 3J(N, Cgamma) -- chi1 observable via N-CA-CB-CG dihedral (= chi1) ---
//
// Perez, C., Lohr, F., Ruterjans, H. & Schmidt, J. M. "Self-consistent
// Karplus parametrization of 3J couplings depending on the polypeptide
// side-chain torsion chi1." J. Am. Chem. Soc. 123, 7081-7093 (2001).
// DOI: 10.1021/ja003724j.
//
// Table 2 (page 7086) "3J(N',Cgamma)" block, consensus row (cos power
// coefficients): A = 1.29, B = -0.49, C = 0.37. Byte-verified
// 2026-05-19 against the page-7086 PDF table. Reference PDF:
//   references/perez-2001-self-consistent-karplus-3j-chi1.pdf
// Arithmetic range: A>0, B<0; max at f(-1) = A - B + C = 2.15 Hz;
// vertex u* = -B/(2A) = 0.190 in (-1, 1), f(u*) = C - B^2/(4A) =
// 0.32 Hz (MIN).
constexpr double KARPLUS_N_CG_A  =  1.29;
constexpr double KARPLUS_N_CG_B  = -0.49;
constexpr double KARPLUS_N_CG_C  =  0.37;
constexpr double KARPLUS_N_CG_THETA = 0.0;  // Perez 2001 uses chi1 = N-CA-CB-CG
        // directly; no phi-style offset.

// --- 3J(C', Cgamma) -- chi1 observable via C-CA-CB-CG dihedral ---
//
// Same paper: Perez, Lohr, Ruterjans & Schmidt, JACS 123:7081 (2001).
// DOI: 10.1021/ja003724j.
//
// Table 2 (page 7086) "3J(C',Cgamma)" block, consensus row (cos power
// coefficients): A = 2.31, B = -0.87, C = 0.55. Byte-verified
// 2026-05-19 against the page-7086 PDF table. The C'-CA-CB-Cgamma
// dihedral differs from chi1 by ~120 deg around CA; the Perez
// self-consistent fit internalizes the substituent-position
// bookkeeping in the per-coupling coefficients, so feeding the
// C-CA-CB-CG atomic dihedral directly is the correct modern usage.
// Reference PDF: references/perez-2001-self-consistent-karplus-3j-
// chi1.pdf Table 2 page 7086.
// Arithmetic range: A>0, B<0; max at f(-1) = A - B + C = 3.73 Hz;
// vertex u* = -B/(2A) = 0.188 in (-1, 1), f(u*) = C - B^2/(4A) =
// 0.468 Hz (MIN).
constexpr double KARPLUS_CP_CG_A =  2.31;
constexpr double KARPLUS_CP_CG_B = -0.87;
constexpr double KARPLUS_CP_CG_C =  0.55;
constexpr double KARPLUS_CP_CG_THETA = 0.0;  // Perez 2001 internalizes the
        // C'-on-CA substituent offset (chi1+120 deg) in the per-coupling
        // (A, B, C); feeding the atomic C-CA-CB-CG dihedral matches the
        // Table 2 consensus row directly. See Table 2 footnote c.

// --- 3J(Halpha, Hbeta) -- chi1 observable via HA-CA-CB-HB dihedral ---
//
// Same paper: Perez, Lohr, Ruterjans & Schmidt, JACS 123:7081 (2001).
// DOI: 10.1021/ja003724j.
//
// Table 2 (page 7086) "3J(Halpha,Hbeta)" block, consensus row (cos
// power coefficients): A = 7.23, B = -1.37, C = 2.22. Byte-verified
// 2026-05-19 against the page-7086 PDF table. The Hbeta atoms are
// typically a prochiral methylene pair (HB2, HB3) on most residues;
// methine (single Hbeta) on Ile/Val/Thr; methyl (HB1/HB2/HB3) on Ala;
// absent on Gly. See JCouplingTimeSeriesTrajectoryResult Hbeta lookup
// for the per-residue policy.
// Reference PDF: references/perez-2001-self-consistent-karplus-3j-
// chi1.pdf Table 2 page 7086.
// Arithmetic range: A>0, B<0; max at f(-1) = A - B + C = 10.82 Hz;
// vertex u* = -B/(2A) = 0.095 in (-1, 1), f(u*) = C - B^2/(4A) =
// 2.155 Hz (MIN).
constexpr double KARPLUS_HA_HB_A =  7.23;
constexpr double KARPLUS_HA_HB_B = -1.37;
constexpr double KARPLUS_HA_HB_C =  2.22;
constexpr double KARPLUS_HA_HB_THETA = 0.0;  // Perez 2001 Table 2 footnote c.
        // The atomic dihedral HA-CA-CB-HB{2,3} differs from chi1 by the
        // Halpha and Hbeta substituent offsets (chi1 + Delta_chi1 ~ ±120°),
        // but the per-coupling (A, B, C) internalize Delta_chi1 -- feeding
        // the atomic dihedral matches the Table 2 consensus row directly.


// ============================================================================
// NMR relaxation constants (15N-1H dipolar + CSA)
//
// Consumed by ReorientationalDynamicsTrajectoryResult to turn Lipari-Szabo
// model-free parameters (S^2, tau_e, global tau_m) into 15N R1/R2/NOE via the
// standard dipolar + CSA rate equations:
//   Lipari & Szabo, J. Am. Chem. Soc. 104, 4546-4559 (1982),
//     DOI 10.1021/ja00381a009  (spectral-density / model-free form);
//   Kay, Torchia & Bax, Biochemistry 28, 8972-8979 (1989),
//     DOI 10.1021/bi00449a003  (the R1/R2/NOE rate expressions);
//   Palmer, Chem. Rev. 104, 3623-3640 (2004), DOI 10.1021/cr030413t (review).
//
// Gyromagnetic ratios and hbar are universal CODATA 2018 values. r_NH and the
// 15N CSA are cited fixed modelling choices (one effective bond length, one
// uniform CSA) -- not calibration-tuned, so they live here, not in the TOML.
// The reporting field B0 genuinely varies per experiment, so it is the one
// relaxation parameter read from CalculatorConfig (relaxation_field_tesla),
// NOT fixed here.
// ============================================================================

// Proton gyromagnetic ratio. CODATA 2018: gamma_p = 2.6752218744(11)e8.
constexpr double GAMMA_H = 2.6752218744e8;          // rad s^-1 T^-1

// 15N gyromagnetic ratio. NEGATIVE -- 15N has a negative magnetogyric ratio,
// which is why the steady-state 15N{1H} NOE is < 1 (often negative). CODATA
// 2018: gamma(15N) = -2.71261804(8)e7. The sign matters only in the NOE
// prefactor gamma_H/gamma_N; everywhere else gamma_N enters squared.
constexpr double GAMMA_N15 = -2.71261804e7;         // rad s^-1 T^-1

// Reduced Planck constant. CODATA 2018 (exact since the 2019 SI redefinition
// of h): hbar = 1.054571817e-34 J s.
constexpr double REDUCED_PLANCK = 1.054571817e-34;  // J s

// Effective N-H internuclear distance for dipolar relaxation. The
// vibrationally-averaged 1.02 A (longer than a ~0.99 A crystallographic N-H)
// is the standard relaxation value; Case, J. Biomol. NMR 15, 95-102 (1999).
// The dipolar coupling scales as r^-3, so this is load-bearing; it is a fixed
// modelling choice, not calibration-tuned.
constexpr double NH_DIPOLAR_BOND_LENGTH_A = 1.02;   // Angstroms

// 15N backbone-amide chemical shift anisotropy (axially-symmetric
// Dsigma = sigma_par - sigma_perp). -172 ppm; Fushman, Tjandra & Cowburn,
// J. Am. Chem. Soc. 120, 10947-10952 (1998). Stored in ppm; the rate
// equations multiply by 1e-6 for the dimensionless anisotropy. A single
// uniform CSA is a modelling choice (real CSA varies residue-to-residue); it
// enters R1/R2 squared, so its sign is immaterial without the dipole-CSA
// cross-correlation term (not computed in this v1).
constexpr double N15_CSA_PPM = -172.0;              // ppm


}  // namespace nmr
