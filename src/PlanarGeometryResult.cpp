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

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();


// Dihedral angle in radians, range [-pi, pi].
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


// Matches DihedralTimeSeriesTrajectoryResult::WrapPi bit-identically.
double WrapPi(double a) {
    if (!std::isfinite(a)) return a;
    return std::remainder(a, 2.0 * M_PI);
}


// Signed out-of-plane displacement of A from the B/C/D plane. Units: A.
double Pyramidalization(const Vec3& A, const Vec3& B,
                         const Vec3& C, const Vec3& D) {
    const Vec3 centroid = (B + C + D) / 3.0;
    const Vec3 normal = (B - centroid).cross(C - centroid);
    const double normal_norm = normal.norm();
    if (normal_norm < 1e-12) return 0.0;
    const Vec3 normal_hat = normal / normal_norm;
    return (A - centroid).dot(normal_hat);
}


// Cremer-Pople 5-ring pucker. Input atoms are in canonical cyclic order.
// Cremer & Pople, JACS 1975, 97, 1354-1358.
struct PuckerCP {
    double Q;        // amplitude (A)
    double theta;    // phase (degrees), [0, 360)
};

PuckerCP CremerPople5Ring(const std::vector<Vec3>& positions) {
    if (positions.size() != 5) return {kNaN, kNaN};

    Vec3 G = Vec3::Zero();
    for (const auto& p : positions) G += p;
    G /= 5.0;

    // Cremer-Pople's sin/cos-weighted normal fixes the phase convention;
    // the simpler edge-cross-product normal is antiparallel for the
    // standard pyrrolidine order and shifts theta by 180 degrees.
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

    double q2_cos_sum = 0.0, q2_sin_sum = 0.0;
    for (size_t j = 0; j < 5; ++j) {
        const double z_j = (positions[j] - G).dot(mean_plane_normal_hat);
        const double phi = 4.0 * M_PI * static_cast<double>(j) / 5.0;
        q2_cos_sum +=  z_j * std::cos(phi);
        q2_sin_sum += -z_j * std::sin(phi);
    }
    const double scale = std::sqrt(2.0 / 5.0);
    const double q2_cos = scale * q2_cos_sum;
    const double q2_sin = scale * q2_sin_sum;
    const double Q = std::sqrt(q2_cos * q2_cos + q2_sin * q2_sin);
    // Below this amplitude, theta is just atan2 of coordinate noise.
    if (Q < 1e-6) return {Q, kNaN};
    double theta = std::atan2(q2_sin, q2_cos) * 180.0 / M_PI;
    if (theta < 0.0) theta += 360.0;
    return {Q, theta};
}


// Three bonded neighbour atom indices for the pyramidalization plane
// at atom `ai`. Resolves the bond list (Atom::bond_indices is a list
// of BOND indices into the bonds_ table, not neighbour atom indices)
// and picks the other endpoint of each bond.
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
    // bond_indices order is not stable; the sign convention needs one.
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

    // BackboneSuccessor is bond-graph based, so cyclic peptides and
    // insertion-code or numbering gaps follow actual C-N connectivity.
    result_ptr->omega_actual_.assign(N_res, kNaN);
    result_ptr->omega_deviation_.assign(N_res, kNaN);
    result_ptr->omega_is_xpro_.assign(N_res, 0);
    int omega_valid = 0;
    int omega_xpro = 0;
    for (size_t ri = 0; ri < N_res; ++ri) {
        auto next_idx = protein.BackboneSuccessor(ri);
        if (!next_idx) continue;
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

    // Per-frame chi2 is an instantaneous dihedral, not a flip rate.
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

    const size_t N_sat = topo.SaturatedRingCount();
    result_ptr->pucker_Q_.assign(N_sat, kNaN);
    result_ptr->pucker_theta_.assign(N_sat, kNaN);
    int pucker_finiteQ = 0;
    for (size_t sat_i = 0; sat_i < N_sat; ++sat_i) {
        const Ring& ring = topo.SaturatedRingAt(sat_i);
        if (ring.atom_indices.size() != 5) {
            // Only the 5-ring Cremer-Pople formulation is implemented.
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
    // bond into i+1 is X→Pro (cis/trans isomerism is real signal
    // there, not a deviation), 0 otherwise. Use this to interpret
    // omega_deviation rows at the consumer side.
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
