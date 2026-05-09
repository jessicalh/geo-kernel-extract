#include "PlanarGeometryResult.h"

#include "AminoAcidType.h"
#include "Atom.h"
#include "ConformationAtom.h"
#include "EnrichmentResult.h"
#include "GeometryChoice.h"
#include "GeometryResult.h"
#include "LegacyAmberTopology.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "Ring.h"
#include "RingTopology.h"
#include "SemanticEnums.h"

#include <cmath>
#include <limits>

namespace nmr {

namespace {

// IEEE quiet NaN for "not applicable" rows in per-residue/per-ring vectors.
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();


// ─────────────────────────────────────────────────────────────────
// Dihedral angle from four positions (radians, range [-π, π])
// Standard formulation: project on the plane orthogonal to b2, take
// signed angle between b1 and b3.
// ─────────────────────────────────────────────────────────────────
double Dihedral(const Vec3& p1, const Vec3& p2,
                 const Vec3& p3, const Vec3& p4) {
    const Vec3 b1 = p2 - p1;
    const Vec3 b2 = p3 - p2;
    const Vec3 b3 = p4 - p3;
    const Vec3 n1 = b1.cross(b2);
    const Vec3 n2 = b2.cross(b3);
    const Vec3 m1 = n1.cross(b2.normalized());
    const double x = n1.dot(n2);
    const double y = m1.dot(n2);
    return std::atan2(y, x);
}


// Wrap an angle to (-π, π].
double WrapPi(double a) {
    while (a > M_PI)  a -= 2.0 * M_PI;
    while (a <= -M_PI) a += 2.0 * M_PI;
    return a;
}


// ─────────────────────────────────────────────────────────────────
// Pyramidalization at atom A whose three bonded neighbours are at
// positions B, C, D. Signed out-of-plane displacement of A from the
// plane through B, C, D, with sign by the right-hand rule on
// (B-centroid) × (C-centroid). Units: Å.
//
// The plane normal is constructed from (B-G) × (C-G) where G is the
// centroid of (B, C, D). For a perfectly planar sp2 site, A == G and
// the displacement is exactly zero.
// ─────────────────────────────────────────────────────────────────
double Pyramidalization(const Vec3& A, const Vec3& B,
                         const Vec3& C, const Vec3& D) {
    const Vec3 G = (B + C + D) / 3.0;
    const Vec3 n = (B - G).cross(C - G);
    const double n_norm = n.norm();
    if (n_norm < 1e-12) return 0.0;  // degenerate plane (collinear)
    const Vec3 n_hat = n / n_norm;
    return (A - G).dot(n_hat);
}


// ─────────────────────────────────────────────────────────────────
// Cremer-Pople pucker for a 5-ring.
//
// Cremer, D. & Pople, J. A. (1975) "A general definition of ring
// puckering coordinates." J. Am. Chem. Soc. 97, 1354–1358.
//
// Atoms in canonical cyclic order (length N = 5). Returns
// (Q, θ_degrees), with θ wrapped to [0, 360) and Q ≥ 0.
// θ mod 72° gives the envelope (E) / twist (T) configuration.
//
// For 5-rings the formulation has a single (Q₂, θ₂) pair (vs the
// (Q₂, θ₂, Q₃, …) family for larger rings). Equation 16 of
// Cremer & Pople 1975.
// ─────────────────────────────────────────────────────────────────
struct PuckerCP {
    double Q;        // amplitude (Å)
    double theta;    // phase (degrees), [0, 360)
};

PuckerCP CremerPople5Ring(const std::vector<Vec3>& positions) {
    if (positions.size() != 5) return {kNaN, kNaN};

    Vec3 G = Vec3::Zero();
    for (const auto& p : positions) G += p;
    G /= 5.0;

    // Mean-plane normal via SVD-style cross-product trick on edge
    // vectors. For a near-planar 5-ring, a robust normal is the
    // average of edge cross products.
    Vec3 n = Vec3::Zero();
    for (size_t j = 0; j < 5; ++j) {
        const Vec3 r_j   = positions[j] - G;
        const Vec3 r_jp1 = positions[(j + 1) % 5] - G;
        n += r_j.cross(r_jp1);
    }
    const double n_norm = n.norm();
    if (n_norm < 1e-12) return {kNaN, kNaN};
    const Vec3 n_hat = n / n_norm;

    // Per-atom z = displacement along normal. Cremer-Pople (Q₂ cos θ₂,
    // Q₂ sin θ₂) projection for the 2-fold (m=2) mode of a 5-ring:
    //
    //   Q₂ cos θ₂ = √(2/5) Σⱼ z_j cos(2 · 2π(j) / 5)
    //   Q₂ sin θ₂ = -√(2/5) Σⱼ z_j sin(2 · 2π(j) / 5)
    //
    // (j = 0..4, 0-indexed; equivalent to the 1..5 indexing in the
    // 1975 paper modulo phase choice).
    double cs = 0.0, sn = 0.0;
    for (size_t j = 0; j < 5; ++j) {
        const double z_j = (positions[j] - G).dot(n_hat);
        const double phi = 4.0 * M_PI * static_cast<double>(j) / 5.0;
        cs +=  z_j * std::cos(phi);
        sn += -z_j * std::sin(phi);
    }
    const double scale = std::sqrt(2.0 / 5.0);
    const double Qcos = scale * cs;
    const double Qsin = scale * sn;
    const double Q = std::sqrt(Qcos * Qcos + Qsin * Qsin);
    double theta = std::atan2(Qsin, Qcos) * 180.0 / M_PI;
    if (theta < 0.0) theta += 360.0;
    return {Q, theta};
}


// Three bonded neighbour atom indices for the pyramidalization plane
// at atom `ai`. Resolves the bond list (Atom::bond_indices is a list
// of BOND indices into the bonds_ table, not neighbour atom indices)
// and picks the other endpoint of each bond.
//
// Returns true and fills neighbours[0..2] when the atom has exactly
// three bonded neighbours; otherwise returns false and the calculator
// emits 0.0 for this atom (legitimate for ring-edge sp2 atoms whose
// bond graph degenerates in unusual structures).
bool ThreeBondedNeighbours(const Protein& protein, size_t ai,
                            std::array<size_t, 3>& neighbours) {
    const Atom& a = protein.AtomAt(ai);
    if (a.bond_indices.size() != 3) return false;
    const auto& bonds = protein.LegacyAmber().Bonds();
    for (size_t k = 0; k < 3; ++k) {
        const Bond& b = bonds.BondAt(a.bond_indices[k]);
        neighbours[k] = (b.atom_index_a == ai) ? b.atom_index_b
                                                : b.atom_index_a;
    }
    return true;
}


}  // anonymous namespace


std::vector<std::type_index> PlanarGeometryResult::Dependencies() const {
    return {
        std::type_index(typeid(GeometryResult)),
        std::type_index(typeid(EnrichmentResult))
    };
}


std::unique_ptr<PlanarGeometryResult> PlanarGeometryResult::Compute(
        ProteinConformation& conf) {

    OperationLog::Scope scope("PlanarGeometryResult::Compute",
        "atoms=" + std::to_string(conf.AtomCount()));

    const Protein& protein = conf.ProteinRef();
    const size_t N_atoms = conf.AtomCount();

    if (N_atoms == 0) {
        OperationLog::Error("PlanarGeometryResult::Compute",
            "Zero atoms — cannot compute planar geometry.");
        return nullptr;
    }
    if (!protein.LegacyAmber().HasAtomSemantic()) {
        OperationLog::Error("PlanarGeometryResult::Compute",
            "Substrate AtomSemanticTable not populated; "
            "PlanarGroupKind dispatch requires it.");
        return nullptr;
    }

    auto result_ptr = std::make_unique<PlanarGeometryResult>();
    result_ptr->conf_ = &conf;

    const LegacyAmberTopology& topo = protein.LegacyAmber();
    const size_t N_res = protein.ResidueCount();

    GeometryChoiceBuilder choices(conf);
    choices.Record(CalculatorId::PlanarGeometry, 0,
        "planar_geometry_per_atom_per_residue_per_ring",
        [&](GeometryChoice& gc) {
            AddNumber(gc, "atoms",        static_cast<double>(N_atoms),       "");
            AddNumber(gc, "residues",     static_cast<double>(N_res),         "");
            AddNumber(gc, "aromatic_rings", static_cast<double>(topo.AromaticRingCount()), "");
            AddNumber(gc, "saturated_rings", static_cast<double>(topo.SaturatedRingCount()), "");
        });

    // ──────────────────────────────────────────────────────────────
    // 1. Per-atom sp2 pyramidalization
    //
    // Signed out-of-plane displacement (Å) at every atom whose
    // typed planar_group != None. Atoms with planar_group == None
    // emit 0.0. Atoms whose bond graph does not have exactly three
    // neighbours (degenerate) also emit 0.0 — never a fail-loud case
    // because the substrate already carries the chemistry decision.
    // ──────────────────────────────────────────────────────────────
    int planar_atom_count = 0;
    double max_abs_pyr = 0.0;
    for (size_t ai = 0; ai < N_atoms; ++ai) {
        auto& ca = conf.MutableAtomAt(ai);
        ca.pyramidalization = 0.0;

        const AtomSemanticTable& sem = topo.SemanticAt(ai);
        if (sem.planar_group == PlanarGroupKind::None) continue;

        std::array<size_t, 3> nb;
        if (!ThreeBondedNeighbours(protein, ai, nb)) continue;

        const Vec3 A = conf.PositionAt(ai);
        const Vec3 B = conf.PositionAt(nb[0]);
        const Vec3 C = conf.PositionAt(nb[1]);
        const Vec3 D = conf.PositionAt(nb[2]);

        const double pyr = Pyramidalization(A, B, C, D);
        ca.pyramidalization = pyr;
        ++planar_atom_count;
        if (std::abs(pyr) > max_abs_pyr) max_abs_pyr = std::abs(pyr);
    }

    // ──────────────────────────────────────────────────────────────
    // 2. Per-residue ω (Cα(i)-C(i)-N(i+1)-Cα(i+1)) and Δω
    //
    // NaN at: C-terminus (no i+1), residues with bad atom indices,
    // and at the bond INTO Pro (i+1.type == Pro). Cα(i+1) cache holds
    // the next residue's CA index.
    // ──────────────────────────────────────────────────────────────
    result_ptr->omega_actual_.assign(N_res, kNaN);
    result_ptr->omega_deviation_.assign(N_res, kNaN);
    int omega_valid = 0;
    for (size_t ri = 0; ri + 1 < N_res; ++ri) {
        const Residue& res_i  = protein.ResidueAt(ri);
        const Residue& res_ip = protein.ResidueAt(ri + 1);

        if (res_ip.type == AminoAcid::PRO) continue;  // X-Pro: cis/trans isomerism, not deviation

        if (res_i.CA  == Residue::NONE) continue;
        if (res_i.C   == Residue::NONE) continue;
        if (res_ip.N  == Residue::NONE) continue;
        if (res_ip.CA == Residue::NONE) continue;

        const double omega = Dihedral(
            conf.PositionAt(res_i.CA),
            conf.PositionAt(res_i.C),
            conf.PositionAt(res_ip.N),
            conf.PositionAt(res_ip.CA));

        result_ptr->omega_actual_[ri] = omega;
        result_ptr->omega_deviation_[ri] = WrapPi(omega - M_PI);
        ++omega_valid;
    }

    // ──────────────────────────────────────────────────────────────
    // 3. Per-aromatic-ring χ₂ (parent residue's chi[1] dihedral)
    //
    // Per Akke & Weininger 2023 (M17): χ₂ IS the canonical ring-flip
    // observable. Per-frame value is *instantaneous*, NOT a flip rate.
    // ──────────────────────────────────────────────────────────────
    const size_t N_arom = topo.AromaticRingCount();
    result_ptr->aromatic_chi2_.assign(N_arom, kNaN);
    int chi2_valid = 0;
    for (size_t ri = 0; ri < N_arom; ++ri) {
        const Ring& ring = topo.AromaticRingAt(ri);
        const size_t parent = ring.parent_residue_index;
        if (parent >= N_res) continue;

        const Residue& res = protein.ResidueAt(parent);
        const Residue::ChiAtoms& chi = res.chi[1];
        if (!chi.Valid()) continue;

        const double chi2 = Dihedral(
            conf.PositionAt(chi.a[0]),
            conf.PositionAt(chi.a[1]),
            conf.PositionAt(chi.a[2]),
            conf.PositionAt(chi.a[3]));
        result_ptr->aromatic_chi2_[ri] = chi2;
        ++chi2_valid;
    }

    // ──────────────────────────────────────────────────────────────
    // 4. Per-saturated-ring Cremer-Pople pucker (Q, θ)
    //
    // 5-ring formulation per Cremer & Pople 1975. Atoms read from
    // the substrate-side canonical cyclic walk (already at
    // construction time per Bundle C / Slice B).
    // ──────────────────────────────────────────────────────────────
    const size_t N_sat = topo.SaturatedRingCount();
    result_ptr->pucker_Q_.assign(N_sat, kNaN);
    result_ptr->pucker_theta_.assign(N_sat, kNaN);
    int pucker_valid = 0;
    for (size_t ri = 0; ri < N_sat; ++ri) {
        const Ring& ring = topo.SaturatedRingAt(ri);
        if (ring.atom_indices.size() != 5) continue;  // 5-ring only

        std::vector<Vec3> ring_pos;
        ring_pos.reserve(5);
        for (size_t ai : ring.atom_indices) {
            ring_pos.push_back(conf.PositionAt(ai));
        }

        const PuckerCP cp = CremerPople5Ring(ring_pos);
        result_ptr->pucker_Q_[ri]     = cp.Q;
        result_ptr->pucker_theta_[ri] = cp.theta;
        if (!std::isnan(cp.Q)) ++pucker_valid;
    }

    OperationLog::Info(LogCalcOther, "PlanarGeometryResult::Compute",
        "pyramidalization: " + std::to_string(planar_atom_count) +
        " planar atoms, max |pyr|=" + std::to_string(max_abs_pyr) +
        " A; omega valid=" + std::to_string(omega_valid) +
        "/" + std::to_string(N_res) +
        "; aromatic_chi2 valid=" + std::to_string(chi2_valid) +
        "/" + std::to_string(N_arom) +
        "; pucker valid=" + std::to_string(pucker_valid) +
        "/" + std::to_string(N_sat));

    return result_ptr;
}


int PlanarGeometryResult::WriteFeatures(
        const ProteinConformation& conf,
        const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    int written = 0;

    // pyramidalization.npy (N,) float64 — per-atom signed out-of-plane Å
    {
        std::vector<double> data(N);
        for (size_t i = 0; i < N; ++i) {
            data[i] = conf.AtomAt(i).pyramidalization;
        }
        NpyWriter::WriteFloat64(output_dir + "/pyramidalization.npy",
                                data.data(), N);
        written++;
    }

    // omega_actual.npy (R,) float64 — per-residue ω (radians)
    {
        const auto& v = omega_actual_;
        NpyWriter::WriteFloat64(output_dir + "/omega_actual.npy",
                                v.data(), v.size());
        written++;
    }

    // omega_deviation.npy (R,) float64 — per-residue ω - π (radians, wrapped)
    {
        const auto& v = omega_deviation_;
        NpyWriter::WriteFloat64(output_dir + "/omega_deviation.npy",
                                v.data(), v.size());
        written++;
    }

    // aromatic_chi2.npy (R_aromatic,) float64 — per-aromatic-ring χ₂ (radians)
    {
        const auto& v = aromatic_chi2_;
        NpyWriter::WriteFloat64(output_dir + "/aromatic_chi2.npy",
                                v.data(), v.size());
        written++;
    }

    // pucker_Q.npy (R_saturated,) float64 — per-saturated-ring amplitude (Å)
    {
        const auto& v = pucker_Q_;
        NpyWriter::WriteFloat64(output_dir + "/pucker_Q.npy",
                                v.data(), v.size());
        written++;
    }

    // pucker_theta.npy (R_saturated,) float64 — per-saturated-ring phase (degrees)
    {
        const auto& v = pucker_theta_;
        NpyWriter::WriteFloat64(output_dir + "/pucker_theta.npy",
                                v.data(), v.size());
        written++;
    }

    return written;
}


}  // namespace nmr
