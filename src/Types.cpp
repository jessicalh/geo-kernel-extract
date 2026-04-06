#include "Types.h"
#include "PhysicalConstants.h"

namespace nmr {

// ============================================================================
// SphericalTensor decomposition
// ============================================================================
//
// The decomposition of a 3x3 tensor into irreducible representations:
//
//   sigma = T0*I + epsilon*T1 + S
//
// where I is the identity, epsilon is the Levi-Civita tensor (mapping the
// antisymmetric pseudovector part back to a matrix), and S is the traceless
// symmetric part.
//
// T2 components use isometric normalization (real spherical harmonics).
// sum|T2_m|^2 = sum S_ij^2 (Frobenius norm preserved).
//

SphericalTensor SphericalTensor::Decompose(const Mat3& s) {
    SphericalTensor st;

    // T0: isotropic = trace / 3
    st.T0 = s.trace() / 3.0;

    // T1: antisymmetric pseudovector via Levi-Civita mapping.
    // For antisymmetric A_ij = (s_ij - s_ji)/2, the dual vector is:
    //   v_x = A_yz, v_y = A_zx, v_z = A_xy
    st.T1[0] = 0.5 * (s(1,2) - s(2,1));
    st.T1[1] = 0.5 * (s(2,0) - s(0,2));
    st.T1[2] = 0.5 * (s(0,1) - s(1,0));

    // Traceless symmetric part: S_ij = (s_ij + s_ji)/2 - (trace/3)*delta_ij
    double Sxx = s(0,0) - st.T0;
    double Syy = s(1,1) - st.T0;
    double Szz = s(2,2) - st.T0;
    double Sxy = 0.5 * (s(0,1) + s(1,0));
    double Sxz = 0.5 * (s(0,2) + s(2,0));
    double Syz = 0.5 * (s(1,2) + s(2,1));

    // T2: isometric real spherical harmonic basis.
    // Coefficients ensure sum|T2_m|^2 = sum S_ij^2 (Frobenius norm preserved).
    static const double SQRT2     = std::sqrt(2.0);
    static const double SQRT3_2   = std::sqrt(3.0 / 2.0);
    static const double INV_SQRT2 = 1.0 / SQRT2;

    st.T2[0] = SQRT2 * Sxy;                // m = -2:  sqrt(2) * Sxy
    st.T2[1] = SQRT2 * Syz;                // m = -1:  sqrt(2) * Syz
    st.T2[2] = SQRT3_2 * Szz;              // m =  0:  sqrt(3/2) * Szz
    st.T2[3] = SQRT2 * Sxz;                // m = +1:  sqrt(2) * Sxz
    st.T2[4] = INV_SQRT2 * (Sxx - Syy);    // m = +2:  (Sxx-Syy) / sqrt(2)

    return st;
}


Mat3 SphericalTensor::Reconstruct() const {
    // Invert the isometric normalization to recover Cartesian elements.
    static const double INV_SQRT2   = 1.0 / std::sqrt(2.0);
    static const double INV_SQRT3_2 = 1.0 / std::sqrt(3.0 / 2.0);
    static const double SQRT2       = std::sqrt(2.0);

    double Sxy = T2[0] * INV_SQRT2;
    double Syz = T2[1] * INV_SQRT2;
    double Szz = T2[2] * INV_SQRT3_2;
    double Sxz = T2[3] * INV_SQRT2;
    double Sxx_m_Syy = T2[4] * SQRT2;

    // Recover Sxx and Syy from Szz and (Sxx - Syy), using tracelessness:
    //   Sxx + Syy + Szz = 0
    //   Sxx - Syy = Sxx_m_Syy
    double Sxx = (-Szz + Sxx_m_Syy) / 2.0;
    double Syy = (-Szz - Sxx_m_Syy) / 2.0;

    // Reassemble: sigma_ij = T0*delta_ij + S_ij + A_ij
    // where A_ij is reconstructed from T1 via the Levi-Civita tensor.
    Mat3 result;
    result(0,0) = Sxx + T0;
    result(1,1) = Syy + T0;
    result(2,2) = Szz + T0;
    result(0,1) = Sxy + T1[2];     // S_xy + A_xy, where A_xy = +v_z
    result(1,0) = Sxy - T1[2];     // S_xy - A_xy
    result(0,2) = Sxz - T1[1];     // S_xz + A_xz, where A_xz = -v_y
    result(2,0) = Sxz + T1[1];     // S_xz - A_xz
    result(1,2) = Syz + T1[0];     // S_yz + A_yz, where A_yz = +v_x
    result(2,1) = Syz - T1[0];     // S_yz - A_yz
    return result;
}


double SphericalTensor::T2Magnitude() const {
    double sum = 0.0;
    for (double v : T2) sum += v * v;
    return std::sqrt(sum);
}


// ============================================================================
// AminoAcid convenience functions (delegate to AminoAcidType table)
// ============================================================================

// Forward-declared in Types.h, implemented here rather than requiring
// AminoAcidType.h in every translation unit.

static const char* AA_CODES[] = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY",
    "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER",
    "THR", "TRP", "TYR", "VAL", "UNK"
};

static const bool AA_AROMATIC[] = {
    false, false, false, false, false, false, false, false,
    true,  false, false, false, false, true,  false, false,
    false, true,  true,  false, false
};

AminoAcid AminoAcidFromThreeLetterCode(const std::string& code) {
    for (int i = 0; i < 20; ++i) {
        if (code == AA_CODES[i]) return static_cast<AminoAcid>(i);
    }
    // Handle common variants
    if (code == "HID" || code == "HIE" || code == "HIP" ||
        code == "HSD" || code == "HSE" || code == "HSP") return AminoAcid::HIS;
    if (code == "ASH") return AminoAcid::ASP;
    if (code == "GLH") return AminoAcid::GLU;
    if (code == "CYX" || code == "CYM") return AminoAcid::CYS;
    if (code == "LYN") return AminoAcid::LYS;
    if (code == "TYM") return AminoAcid::TYR;
    if (code == "MSE") return AminoAcid::MET;
    return AminoAcid::Unknown;
}

std::string ThreeLetterCodeForAminoAcid(AminoAcid aa) {
    int i = static_cast<int>(aa);
    if (i >= 0 && i <= 20) return AA_CODES[i];
    return "UNK";
}

bool IsAromaticAminoAcid(AminoAcid aa) {
    int i = static_cast<int>(aa);
    if (i >= 0 && i <= 20) return AA_AROMATIC[i];
    return false;
}

}  // namespace nmr
