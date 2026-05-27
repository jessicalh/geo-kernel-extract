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

#include <algorithm>
#include <array>
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


// Wrap to [-π, π] via std::remainder (single-call, IEEE round-half-to-
// even at ±π). Matches DihedralTimeSeriesTrajectoryResult::WrapPi
// bit-identically so omega_deviation agrees across the two producers.
double WrapPi(double a) {
    if (!std::isfinite(a)) return a;
    return std::remainder(a, 2.0 * M_PI);
}


// ─────────────────────────────────────────────────────────────────
// Pyramidalization at atom A whose three bonded neighbours are at
// positions B, C, D. Signed out-of-plane displacement of A from the
// plane through B, C, D, with sign by the right-hand rule on
// (B-centroid) × (C-centroid). Units: Å.
//
// The plane normal is constructed from (B-centroid) × (C-centroid)
// where centroid is the centroid of (B, C, D). For a perfectly planar
// sp2 site A lies in the neighbour plane, so (A - centroid)·n̂ is zero.
// ─────────────────────────────────────────────────────────────────
double Pyramidalization(const Vec3& A, const Vec3& B,
                         const Vec3& C, const Vec3& D) {
    const Vec3 centroid = (B + C + D) / 3.0;
    const Vec3 normal = (B - centroid).cross(C - centroid);
    const double normal_norm = normal.norm();
    if (normal_norm < 1e-12) return 0.0;  // degenerate plane (collinear)
    const Vec3 normal_hat = normal / normal_norm;
    return (A - centroid).dot(normal_hat);
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

    // mean-plane normal
    //
    // Canonical Cremer-Pople 1975 construction, Eqs 11–12: orthogonal
    // sin/cos-weighted basis from the displacement vectors, with the
    // normal as the cross product.
    //
    //   R'₁ = Σⱼ (rⱼ - G) sin(2π j / N)
    //   R'₂ = Σⱼ (rⱼ - G) cos(2π j / N)
    //   n   = R'₁ × R'₂
    //
    // The simpler edge-cross-product accumulator
    //   n = Σⱼ rⱼ × rⱼ₊₁
    // gives a normal anti-parallel to the canonical direction for the
    // standard cyclic atom order of a pyrrolidine 5-ring, which inverts
    // the (Q, θ) phase by 180° and silently swaps envelope/twist
    // endo/exo labels.
    Vec3 sin_basis = Vec3::Zero();
    Vec3 cos_basis = Vec3::Zero();
    for (size_t j = 0; j < 5; ++j) {
        const Vec3 r_j = positions[j] - G;
        const double phi = 2.0 * M_PI * static_cast<double>(j) / 5.0;
        sin_basis += r_j * std::sin(phi);
        cos_basis += r_j * std::cos(phi);
    }
    const Vec3 mean_plane_normal = sin_basis.cross(cos_basis);
    const double mean_plane_normal_norm = mean_plane_normal.norm();
    if (mean_plane_normal_norm < 1e-12) return {kNaN, kNaN};
    const Vec3 mean_plane_normal_hat = mean_plane_normal / mean_plane_normal_norm;

    // (Q,θ) projection
    //
    // Per-atom z = displacement along normal. Cremer-Pople (Q₂ cos θ₂,
    // Q₂ sin θ₂) projection for the 2-fold (m=2) mode of a 5-ring:
    //
    //   Q₂ cos θ₂ = √(2/5) Σⱼ z_j cos(2 · 2π(j) / 5)
    //   Q₂ sin θ₂ = -√(2/5) Σⱼ z_j sin(2 · 2π(j) / 5)
    //
    // (j = 0..4, 0-indexed; equivalent to the 1..5 indexing in the
    // 1975 paper modulo phase choice).
    double q2_cos_sum = 0.0, q2_sin_sum = 0.0;
    for (size_t j = 0; j < 5; ++j) {
        const double z_j = (positions[j] - G).dot(mean_plane_normal_hat);
        // m=2 (2nd-harmonic) phase, hence 4π j / 5.
        const double phi = 4.0 * M_PI * static_cast<double>(j) / 5.0;
        q2_cos_sum +=  z_j * std::cos(phi);
        q2_sin_sum += -z_j * std::sin(phi);
    }
    const double scale = std::sqrt(2.0 / 5.0);
    const double q2_cos = scale * q2_cos_sum;
    const double q2_sin = scale * q2_sin_sum;
    const double Q = std::sqrt(q2_cos * q2_cos + q2_sin * q2_sin);
    // Sub-amplitude degeneracy guard: a perfectly planar pentagon makes
    // z_j ≈ 0 for every vertex, so q2_cos/q2_sin become floating-point
    // noise and theta = atan2(noise, noise) is meaningless. Typical
    // real proline pyrrolidine puckering is far above this cutoff
    // (probable source: Ho & Cornilescu 2000 JBNMR 18:155). Anything
    // below 1e-6 Å is structural noise at Å-scale coordinates. Return
    // NaN theta so downstream readers see
    // "degenerate" honestly (the H5 attr already says NaN means
    // degenerate).
    if (Q < 1e-6) return {Q, kNaN};
    double theta = std::atan2(q2_sin, q2_cos) * 180.0 / M_PI;
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
    // Sort by atom index so the pyramidalization sign is build-stable
    // and portable across structures (bond order in `bond_indices` is
    // not stable).
    std::sort(neighbours.begin(), neighbours.end());
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
    // 2. Per-residue ω (Cα(i)-C(i)-N(successor)-Cα(successor)) and Δω
    //
    // ω is emitted for every well-defined peptide bond, INCLUDING
    // X→Pro bonds. cis/trans isomerism at X-Pro is a real
    // conformational signal, not a "non-planar amide" deviation, but
    // the value itself is what the consumer wants. The
    // `omega_is_xpro` mask tags those rows so the BMRB-stratified
    // atlas can interpret deviation values appropriately at the
    // consumer side.
    //
    // Connectivity comes from `protein.BackboneSuccessor(ri)` (canonical
    // bond-graph walk). The query returns the residue index
    // whose N is covalently bonded to res(ri).C, or nullopt at chain
    // ends. This is the geometry-native substrate: the bond graph is
    // the authoritative answer, label-based heuristics (chain_id,
    // sequence_number, terminal_state, insertion_code) are banned per
    // the "Backbone connectivity discipline" in OBJECT_MODEL.md.
    //
    // The query is wrap-correct (cyclic peptide: Successor(N-1) = 0),
    // bonded across antibody insertion codes (100A -> 100B), correct
    // on residue numbering gaps with intact bonds (those ARE bonded;
    // the bond graph says so), and correct on ACE/NME caps (cap C/N
    // participate).
    //
    // NaN at: chain boundary (Successor returns nullopt), residues
    // with missing CA backbone-cache atoms (incomplete structure).
    // ──────────────────────────────────────────────────────────────
    result_ptr->omega_actual_.assign(N_res, kNaN);
    result_ptr->omega_deviation_.assign(N_res, kNaN);
    result_ptr->omega_is_xpro_.assign(N_res, 0);
    int omega_valid = 0;
    int omega_xpro = 0;
    for (size_t ri = 0; ri < N_res; ++ri) {
        auto next_idx = protein.BackboneSuccessor(ri);
        if (!next_idx) continue;
        // Successor guarantees res_i.C and res_next.N exist.
        const Residue& res_i    = protein.ResidueAt(ri);
        const Residue& res_next = protein.ResidueAt(*next_idx);
        if (res_i.CA    == Residue::NONE) continue;
        if (res_next.CA == Residue::NONE) continue;

        const double omega = Dihedral(
            conf.PositionAt(res_i.CA),
            conf.PositionAt(res_i.C),
            conf.PositionAt(res_next.N),
            conf.PositionAt(res_next.CA));

        result_ptr->omega_actual_[ri] = omega;
        result_ptr->omega_deviation_[ri] = WrapPi(omega - M_PI);
        if (res_next.type == AminoAcid::PRO) {
            result_ptr->omega_is_xpro_[ri] = 1;
            ++omega_xpro;
        }
        ++omega_valid;
    }

    // ──────────────────────────────────────────────────────────────
    // 3. Per-aromatic-ring χ₂ (parent residue's chi[1] dihedral)
    //
    // Per Akke & Weininger 2023 (M17), Phe/Tyr ring flips are 180°
    // rotamer transitions of χ₂. Per-frame value is *instantaneous*,
    // NOT a flip rate.
    // ──────────────────────────────────────────────────────────────
    const size_t N_arom = topo.AromaticRingCount();
    result_ptr->aromatic_chi2_.assign(N_arom, kNaN);
    int chi2_valid = 0;
    for (size_t arom_i = 0; arom_i < N_arom; ++arom_i) {
        const Ring& ring = topo.AromaticRingAt(arom_i);
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
        result_ptr->aromatic_chi2_[arom_i] = chi2;
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
    int pucker_finiteQ = 0;
    for (size_t sat_i = 0; sat_i < N_sat; ++sat_i) {
        const Ring& ring = topo.SaturatedRingAt(sat_i);
        if (ring.atom_indices.size() != 5) {
            // Only the 5-ring Cremer-Pople formulation is implemented
            // (Pro pyrrolidine today); a non-5 saturated ring is left
            // NaN. Note the skip so it is not silent.
            OperationLog::Info(LogCalcOther, "PlanarGeometryResult::Compute",
                "saturated ring " + std::to_string(sat_i) + " has " +
                std::to_string(ring.atom_indices.size()) +
                " atoms (not 5); pucker left NaN.");
            continue;
        }

        std::vector<Vec3> ring_pos;
        ring_pos.reserve(5);
        for (size_t ai : ring.atom_indices) {
            ring_pos.push_back(conf.PositionAt(ai));
        }

        const PuckerCP cp = CremerPople5Ring(ring_pos);
        result_ptr->pucker_Q_[sat_i]     = cp.Q;
        result_ptr->pucker_theta_[sat_i] = cp.theta;
        if (!std::isnan(cp.Q)) ++pucker_finiteQ;
    }

    OperationLog::Info(LogCalcOther, "PlanarGeometryResult::Compute",
        "pyramidalization: " + std::to_string(planar_atom_count) +
        " planar atoms, max |pyr|=" + std::to_string(max_abs_pyr) +
        " A; omega valid=" + std::to_string(omega_valid) +
        "/" + std::to_string(N_res) +
        " (xpro=" + std::to_string(omega_xpro) + ")" +
        "; aromatic_chi2 valid=" + std::to_string(chi2_valid) +
        "/" + std::to_string(N_arom) +
        "; pucker finite-Q=" + std::to_string(pucker_finiteQ) +
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

    // omega_is_xpro.npy (R,) int8 — per-residue mask: 1 where the
    // bond from residue i to its backbone successor is X→Pro
    // (cis/trans isomerism is real signal there, not a deviation),
    // 0 otherwise. Use this to interpret omega_deviation rows at the
    // consumer side.
    {
        const auto& v = omega_is_xpro_;
        NpyWriter::WriteInt8(output_dir + "/omega_is_xpro.npy",
                             reinterpret_cast<const int8_t*>(v.data()),
                             v.size());
        written++;
    }

    return written;
}


}  // namespace nmr
