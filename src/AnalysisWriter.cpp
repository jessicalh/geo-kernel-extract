#include "AnalysisWriter.h"
#include "ConformationAtom.h"
#include "Atom.h"
#include "Bond.h"
#include "Ring.h"
#include "Residue.h"
#include "AminoAcidType.h"
#include "Types.h"
#include "OperationLog.h"

#include <highfive/H5File.hpp>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5DataSpace.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <numeric>

namespace nmr {

// ── Construction ────────────────────────────────────────────────

AnalysisWriter::AnalysisWriter(const Protein& protein,
                               const std::string& protein_id)
    : protein_(protein), protein_id_(protein_id)
{}


// ── BeginTrajectory ─────────────────────────────────────────────

void AnalysisWriter::BeginTrajectory(size_t n_frames_estimate,
                                     size_t n_atoms, size_t stride) {
    n_atoms_ = n_atoms;
    n_frames_ = 0;
    stride_ = stride;

    const size_t T = n_frames_estimate;
    const size_t N = n_atoms;

    // Frame metadata
    frame_times_.reserve(T);
    frame_indices_.reserve(T);

    // Positions: T * N * 3
    positions_.reserve(T * N * 3);

    // Ring current — per-type: T * N * 8 (T0), T * N * 8 * 5 (T2)
    bs_T0_per_type_.reserve(T * N * 8);
    bs_T2_per_type_.reserve(T * N * 8 * 5);
    hm_T0_per_type_.reserve(T * N * 8);
    hm_T2_per_type_.reserve(T * N * 8 * 5);

    // Ring current — SphericalTensor totals: T * N * 9
    bs_shielding_.reserve(T * N * 9);
    hm_shielding_.reserve(T * N * 9);
    rs_shielding_.reserve(T * N * 9);

    // Ring proximity scalars: T * N
    n_rings_3A_.reserve(T * N);
    n_rings_5A_.reserve(T * N);
    n_rings_8A_.reserve(T * N);
    n_rings_12A_.reserve(T * N);
    mean_ring_dist_.reserve(T * N);
    nearest_ring_atom_.reserve(T * N);

    // Exponential-weighted sums
    G_iso_exp_sum_.reserve(T * N);
    G_T2_exp_sum_.reserve(T * N * 5);
    G_iso_var_8A_.reserve(T * N);

    // B-field vector: T * N * 3
    total_B_field_.reserve(T * N * 3);
}


// ── HarvestFrame ────────────────────────────────────────────────

void AnalysisWriter::HarvestFrame(const ProteinConformation& conf,
                                  size_t frame_idx, double time_ps) {
    assert(n_atoms_ > 0 && "BeginTrajectory not called");

    const size_t N = n_atoms_;
    assert(conf.AtomCount() == N);

    frame_times_.push_back(time_ps);
    frame_indices_.push_back(static_cast<int32_t>(frame_idx));

    for (size_t ai = 0; ai < N; ++ai) {
        const ConformationAtom& ca = conf.AtomAt(ai);

        // ── Positions ──────────────────────────────────────────
        Vec3 pos = ca.Position();
        positions_.push_back(pos.x());
        positions_.push_back(pos.y());
        positions_.push_back(pos.z());

        // ── Ring current: per-type BS ──────────────────────────
        for (size_t rt = 0; rt < 8; ++rt)
            bs_T0_per_type_.push_back(ca.per_type_G_T0_sum[rt]);
        for (size_t rt = 0; rt < 8; ++rt)
            for (size_t c = 0; c < 5; ++c)
                bs_T2_per_type_.push_back(ca.per_type_G_T2_sum[rt][c]);

        // ── Ring current: per-type HM ──────────────────────────
        for (size_t rt = 0; rt < 8; ++rt)
            hm_T0_per_type_.push_back(ca.per_type_hm_T0_sum[rt]);
        for (size_t rt = 0; rt < 8; ++rt)
            for (size_t c = 0; c < 5; ++c)
                hm_T2_per_type_.push_back(ca.per_type_hm_T2_sum[rt][c]);

        // ── Ring current: SphericalTensor totals ───────────────
        AppendSpherical(bs_shielding_, ca.bs_shielding_contribution);
        AppendSpherical(hm_shielding_, ca.hm_shielding_contribution);
        AppendSpherical(rs_shielding_, ca.ringchi_shielding_contribution);

        // ── Ring proximity scalars ─────────────────────────────
        n_rings_3A_.push_back(static_cast<int16_t>(ca.n_rings_within_3A));
        n_rings_5A_.push_back(static_cast<int16_t>(ca.n_rings_within_5A));
        n_rings_8A_.push_back(static_cast<int16_t>(ca.n_rings_within_8A));
        n_rings_12A_.push_back(static_cast<int16_t>(ca.n_rings_within_12A));
        mean_ring_dist_.push_back(ca.mean_ring_distance);
        nearest_ring_atom_.push_back(ca.nearest_ring_atom_distance);

        // ── Exponential-weighted sums ──────────────────────────
        G_iso_exp_sum_.push_back(ca.G_iso_exp_sum);
        for (size_t c = 0; c < 5; ++c)
            G_T2_exp_sum_.push_back(ca.G_T2_exp_sum[c]);
        G_iso_var_8A_.push_back(ca.G_iso_var_8A);

        // ── B-field vector ─────────────────────────────────────
        AppendVec3(total_B_field_, ca.total_B_field);

        // ── EFG group ─────────────────────────────────────────
        AppendSpherical(coulomb_efg_total_, ca.coulomb_EFG_total_spherical);
        AppendSpherical(coulomb_efg_backbone_, ca.coulomb_EFG_backbone_spherical);
        AppendSpherical(coulomb_efg_aromatic_, ca.coulomb_EFG_aromatic_spherical);
        AppendVec3(coulomb_E_total_, ca.coulomb_E_total);
        AppendVec3(coulomb_E_backbone_, ca.coulomb_E_backbone);
        AppendVec3(coulomb_E_sidechain_, ca.coulomb_E_sidechain);
        AppendVec3(coulomb_E_aromatic_, ca.coulomb_E_aromatic);
        AppendVec3(coulomb_E_solvent_, ca.coulomb_E_solvent);
        coulomb_E_magnitude_.push_back(ca.coulomb_E_magnitude);
        coulomb_E_bond_proj_.push_back(ca.coulomb_E_bond_proj);
        coulomb_E_bb_frac_.push_back(ca.coulomb_E_backbone_frac);
        AppendSpherical(apbs_efg_, ca.apbs_efg_spherical);
        AppendVec3(apbs_efield_, ca.apbs_efield);
        AppendSpherical(aimnet2_efg_total_, ca.aimnet2_EFG_total_spherical);
        AppendSpherical(aimnet2_efg_backbone_, ca.aimnet2_EFG_backbone_spherical);
        AppendSpherical(aimnet2_efg_aromatic_, ca.aimnet2_EFG_aromatic_spherical);
        AppendSpherical(aimnet2_shielding_, ca.aimnet2_shielding_contribution);

        // ── Bond anisotropy group ─────────────────────────────
        AppendSpherical(mc_shielding_, ca.mc_shielding_contribution);
        AppendSpherical(mc_T2_backbone_, ca.T2_backbone_total);
        AppendSpherical(mc_T2_sidechain_, ca.T2_sidechain_total);
        AppendSpherical(mc_T2_aromatic_, ca.T2_aromatic_total);
        AppendSpherical(mc_T2_CO_nearest_, ca.T2_CO_nearest);
        AppendSpherical(mc_T2_CN_nearest_, ca.T2_CN_nearest);
        mc_co_sum_.push_back(ca.mcconnell_co_sum);
        mc_cn_sum_.push_back(ca.mcconnell_cn_sum);
        mc_sidechain_sum_.push_back(ca.mcconnell_sidechain_sum);
        mc_aromatic_sum_.push_back(ca.mcconnell_aromatic_sum);
        mc_co_nearest_.push_back(ca.mcconnell_co_nearest);
        mc_nearest_CO_dist_.push_back(ca.nearest_CO_dist);
        mc_nearest_CN_dist_.push_back(ca.nearest_CN_dist);
        AppendVec3(mc_nearest_CO_mid_, ca.nearest_CO_midpoint);
        AppendVec3(mc_dir_nearest_CO_, ca.dir_nearest_CO);

        // ── Quadrupole group ──────────────────────────────────
        AppendSpherical(pq_shielding_, ca.piquad_shielding_contribution);
        for (size_t rt = 0; rt < 8; ++rt)
            pq_T0_per_type_.push_back(ca.per_type_pq_scalar_sum[rt]);
        for (size_t rt = 0; rt < 8; ++rt)
            for (size_t c = 0; c < 5; ++c)
                pq_T2_per_type_.push_back(ca.per_type_pq_T2_sum[rt][c]);

        // ── Dispersion group ──────────────────────────────────
        AppendSpherical(disp_shielding_, ca.disp_shielding_contribution);
        for (size_t rt = 0; rt < 8; ++rt)
            disp_T0_per_type_.push_back(ca.per_type_disp_scalar_sum[rt]);
        for (size_t rt = 0; rt < 8; ++rt)
            for (size_t c = 0; c < 5; ++c)
                disp_T2_per_type_.push_back(ca.per_type_disp_T2_sum[rt][c]);

        // ── H-bond group ──────────────────────────────────────
        AppendSpherical(hbond_shielding_, ca.hbond_shielding_contribution);
        AppendSpherical(hbond_nearest_sph_, ca.hbond_nearest_spherical);
        hbond_nearest_dist_.push_back(ca.hbond_nearest_dist);
        AppendVec3(hbond_nearest_dir_, ca.hbond_nearest_dir);
        hbond_inv_d3_.push_back(ca.hbond_inv_d3);
        hbond_count_.push_back(static_cast<int16_t>(ca.hbond_count_within_3_5A));
        hbond_is_donor_.push_back(ca.hbond_is_donor ? 1 : 0);
        hbond_is_acceptor_.push_back(ca.hbond_is_acceptor ? 1 : 0);
        hbond_is_backbone_.push_back(ca.hbond_is_backbone ? 1 : 0);

        // ── SASA group ────────────────────────────────────────
        sasa_.push_back(ca.atom_sasa);
        AppendVec3(sasa_normal_, ca.sasa_normal);

        // ── Water group ───────────────────────────────────────
        AppendVec3(water_efield_, ca.water_efield);
        AppendSpherical(water_efg_, ca.water_efg_spherical);
        AppendVec3(water_efield_first_, ca.water_efield_first);
        AppendSpherical(water_efg_first_, ca.water_efg_first_spherical);
        water_n_first_.push_back(static_cast<int16_t>(ca.water_n_first));
        water_n_second_.push_back(static_cast<int16_t>(ca.water_n_second));
        water_half_shell_.push_back(ca.half_shell_asymmetry);
        water_dipole_cos_.push_back(ca.mean_water_dipole_cos);
        water_nearest_ion_dist_.push_back(
            ca.nearest_ion_distance < 1e30 ? ca.nearest_ion_distance : -1.0);
        water_nearest_ion_charge_.push_back(ca.nearest_ion_charge);
        AppendVec3(water_dipole_vec_, ca.water_dipole_vector);
        AppendVec3(water_surface_normal_, ca.water_surface_normal);
        water_sasa_asym_.push_back(ca.sasa_half_shell_asymmetry);
        water_sasa_align_.push_back(ca.sasa_dipole_alignment);
        water_sasa_cohere_.push_back(ca.sasa_dipole_coherence);
        water_sasa_first_n_.push_back(static_cast<int16_t>(ca.sasa_first_shell_count));

        // ── Charges group ─────────────────────────────────────
        aimnet2_charge_.push_back(ca.aimnet2_charge);
        eeq_charge_.push_back(ca.eeq_charge);
        eeq_cn_.push_back(ca.eeq_cn);

        // ── AIMNet2 embedding ─────────────────────────────────
        for (size_t d = 0; d < AIMNET2_AIM_DIMS; ++d)
            aimnet2_aim_.push_back(ca.aimnet2_aim[d]);

        // ── Per-ring data (K=6 nearest) ──────────────────────
        {
            // Sort ring_neighbours by distance, take top K
            const auto& rn = ca.ring_neighbours;
            std::vector<size_t> order(rn.size());
            std::iota(order.begin(), order.end(), 0);
            std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
                return rn[a].distance_to_center < rn[b].distance_to_center;
            });

            for (size_t k = 0; k < K_RINGS; ++k) {
                if (k < order.size()) {
                    const auto& r = rn[order[k]];
                    // Geometry: distance, rho, z, theta, cos_phi, sin_phi
                    per_ring_geom_.push_back(r.distance_to_center);
                    per_ring_geom_.push_back(r.rho);
                    per_ring_geom_.push_back(r.z);
                    per_ring_geom_.push_back(r.theta);
                    per_ring_geom_.push_back(r.cos_phi);
                    per_ring_geom_.push_back(r.sin_phi);
                    // Ring type
                    per_ring_type_.push_back(static_cast<int8_t>(r.ring_type));
                    // BS T2
                    for (int c = 0; c < 5; ++c)
                        per_ring_bs_T2_.push_back(r.G_spherical.T2[c]);
                    // HM T2
                    for (int c = 0; c < 5; ++c)
                        per_ring_hm_T2_.push_back(r.hm_G_spherical.T2[c]);
                    // Chi T2
                    for (int c = 0; c < 5; ++c)
                        per_ring_chi_T2_.push_back(r.chi_spherical.T2[c]);
                    // PQ T2
                    for (int c = 0; c < 5; ++c)
                        per_ring_pq_T2_.push_back(r.quad_spherical.T2[c]);
                    // HM_H T2
                    for (int c = 0; c < 5; ++c)
                        per_ring_hm_H_T2_.push_back(r.hm_H_spherical.T2[c]);
                    // Dispersion scalar + T2
                    per_ring_disp_scalar_.push_back(r.disp_scalar);
                    for (int c = 0; c < 5; ++c)
                        per_ring_disp_T2_.push_back(r.disp_spherical.T2[c]);
                } else {
                    // Pad with zeros for atoms near fewer than K rings
                    for (int c = 0; c < 6; ++c) per_ring_geom_.push_back(0.0);
                    per_ring_type_.push_back(-1);
                    for (int c = 0; c < 5; ++c) per_ring_bs_T2_.push_back(0.0);
                    for (int c = 0; c < 5; ++c) per_ring_hm_T2_.push_back(0.0);
                    for (int c = 0; c < 5; ++c) per_ring_chi_T2_.push_back(0.0);
                    for (int c = 0; c < 5; ++c) per_ring_pq_T2_.push_back(0.0);
                    for (int c = 0; c < 5; ++c) per_ring_hm_H_T2_.push_back(0.0);
                    per_ring_disp_scalar_.push_back(0.0);
                    for (int c = 0; c < 5; ++c) per_ring_disp_T2_.push_back(0.0);
                }
            }
        }

        // ── Coulomb shielding SphericalTensor ────────────────
        AppendSpherical(coulomb_shielding_, ca.coulomb_shielding_contribution);
    }

    // ── Capture enrichment flags on first frame (frame-invariant) ──
    if (!enrichment_captured_) {
        enrichment_captured_ = true;
        is_amide_H_.resize(N);
        is_alpha_H_.resize(N);
        is_methyl_.resize(N);
        is_aromatic_H_.resize(N);
        is_on_aromatic_residue_.resize(N);
        is_hbond_donor_.resize(N);
        is_hbond_acceptor_.resize(N);
        parent_is_sp2_.resize(N);
        graph_dist_N_.resize(N);
        graph_dist_O_.resize(N);
        eneg_sum_1_.resize(N);
        eneg_sum_2_.resize(N);
        n_pi_bonds_3_.resize(N);
        bfs_to_nearest_ring_.resize(N);
        bfs_decay_.resize(N);
        partial_charge_.resize(N);
        vdw_radius_.resize(N);

        for (size_t i = 0; i < N; ++i) {
            const ConformationAtom& ca = conf.AtomAt(i);
            is_amide_H_[i] = ca.is_amide_H ? 1 : 0;
            is_alpha_H_[i] = ca.is_alpha_H ? 1 : 0;
            is_methyl_[i] = ca.is_methyl ? 1 : 0;
            is_aromatic_H_[i] = ca.is_aromatic_H ? 1 : 0;
            is_on_aromatic_residue_[i] = ca.is_on_aromatic_residue ? 1 : 0;
            is_hbond_donor_[i] = ca.is_hbond_donor ? 1 : 0;
            is_hbond_acceptor_[i] = ca.is_hbond_acceptor ? 1 : 0;
            parent_is_sp2_[i] = ca.parent_is_sp2 ? 1 : 0;
            graph_dist_N_[i] = ca.graph_dist_N;
            graph_dist_O_[i] = ca.graph_dist_O;
            eneg_sum_1_[i] = ca.eneg_sum_1;
            eneg_sum_2_[i] = ca.eneg_sum_2;
            n_pi_bonds_3_[i] = ca.n_pi_bonds_3;
            bfs_to_nearest_ring_[i] = ca.bfs_to_nearest_ring_atom;
            bfs_decay_[i] = ca.bfs_decay;
            partial_charge_[i] = ca.partial_charge;
            vdw_radius_[i] = ca.vdw_radius;
        }
    }

    // ── Per-frame ring geometry ──────────────────────────────────
    {
        const size_t R = protein_.RingCount();
        for (size_t ri = 0; ri < R; ++ri) {
            const Ring& ring = protein_.RingAt(ri);
            // Recompute ring geometry from current frame positions
            Vec3 center = Vec3::Zero();
            for (size_t ai : ring.atom_indices)
                center += conf.PositionAt(ai);
            center /= static_cast<double>(ring.atom_indices.size());

            // Normal from cross product of first two edge vectors
            // (matches GeometryResult's SVD for planar rings)
            Vec3 normal = Vec3::Zero();
            if (ring.atom_indices.size() >= 3) {
                Vec3 v1 = conf.PositionAt(ring.atom_indices[1]) -
                          conf.PositionAt(ring.atom_indices[0]);
                Vec3 v2 = conf.PositionAt(ring.atom_indices[2]) -
                          conf.PositionAt(ring.atom_indices[0]);
                normal = v1.cross(v2);
                double n = normal.norm();
                if (n > 1e-10) normal /= n;
            }

            double radius = 0.0;
            for (size_t ai : ring.atom_indices)
                radius += (conf.PositionAt(ai) - center).norm();
            radius /= static_cast<double>(ring.atom_indices.size());

            ring_geometry_.push_back(center.x());
            ring_geometry_.push_back(center.y());
            ring_geometry_.push_back(center.z());
            ring_geometry_.push_back(normal.x());
            ring_geometry_.push_back(normal.y());
            ring_geometry_.push_back(normal.z());
            ring_geometry_.push_back(radius);
        }
    }

    ++n_frames_;
}


// ── AppendSpherical ─────────────────────────────────────────────

void AnalysisWriter::AppendSpherical(std::vector<double>& buf,
                                     const SphericalTensor& st) {
    buf.push_back(st.T0);
    buf.push_back(st.T1[0]);
    buf.push_back(st.T1[1]);
    buf.push_back(st.T1[2]);
    buf.push_back(st.T2[0]);
    buf.push_back(st.T2[1]);
    buf.push_back(st.T2[2]);
    buf.push_back(st.T2[3]);
    buf.push_back(st.T2[4]);
}

void AnalysisWriter::AppendVec3(std::vector<double>& buf,
                                const Vec3& v) {
    buf.push_back(v.x());
    buf.push_back(v.y());
    buf.push_back(v.z());
}


// ── H5 write helpers (file-local, not class methods) ────────────
// These take HighFive::File& so the header never sees HighFive.

namespace {

void WriteAtomMetadata(HighFive::File& file,
                       const Protein& protein, size_t N) {
    std::vector<int32_t> element(N);
    std::vector<int32_t> residue_index(N);
    std::vector<std::string> atom_name(N);
    std::vector<int32_t> atom_role(N);
    std::vector<int32_t> hybridisation(N);
    std::vector<int32_t> n_bonded(N);
    std::vector<int32_t> graph_dist_ring(N);
    std::vector<int8_t> is_backbone(N);
    std::vector<int8_t> is_conjugated(N);

    // Read from the typed object model.
    // Atom: element (enum), residue_index, pdb_atom_name (display only),
    //       bond_indices.
    // ConformationAtom (from conformation 0): role (enum), hybridisation
    //       (enum), is_backbone, is_conjugated, graph_dist_ring.
    // These are all enrichment properties — constant across frames,
    // set once by EnrichmentResult on conformation 0.
    const auto& conf0 = protein.ConformationAt(0);

    for (size_t i = 0; i < N; ++i) {
        const Atom& atom = protein.AtomAt(i);
        const ConformationAtom& ca = conf0.AtomAt(i);

        element[i] = static_cast<int32_t>(AtomicNumberForElement(atom.element));
        residue_index[i] = static_cast<int32_t>(atom.residue_index);
        atom_name[i] = atom.pdb_atom_name;  // display/serialisation only
        atom_role[i] = static_cast<int32_t>(ca.role);
        hybridisation[i] = static_cast<int32_t>(ca.hybridisation);
        n_bonded[i] = static_cast<int32_t>(atom.bond_indices.size());
        graph_dist_ring[i] = ca.graph_dist_ring;
        is_backbone[i] = ca.is_backbone ? 1 : 0;
        is_conjugated[i] = ca.is_conjugated ? 1 : 0;
    }

    auto grp = file.createGroup("atoms");
    grp.createDataSet("element", element);
    grp.createDataSet("residue_index", residue_index);
    grp.createDataSet("atom_name", atom_name);
    grp.createDataSet("atom_role", atom_role);
    grp.createDataSet("hybridisation", hybridisation);
    grp.createDataSet("n_bonded", n_bonded);
    grp.createDataSet("graph_dist_ring", graph_dist_ring);
    grp.createDataSet("is_backbone", is_backbone);
    grp.createDataSet("is_conjugated", is_conjugated);
}

void WriteResidueMetadata(HighFive::File& file,
                          const Protein& protein) {
    const size_t R = protein.ResidueCount();

    std::vector<std::string> residue_name(R);
    std::vector<int32_t> residue_number(R);
    std::vector<std::string> chain_id(R);

    for (size_t i = 0; i < R; ++i) {
        const Residue& res = protein.ResidueAt(i);
        // Residue type is an AminoAcid enum; three-letter code comes
        // from the typed AminoAcidType, not from a string field.
        residue_name[i] = res.AminoAcidInfo().three_letter_code;
        residue_number[i] = res.sequence_number;
        chain_id[i] = res.chain_id;
    }

    auto grp = file.createGroup("residues");
    grp.createDataSet("residue_name", residue_name);
    grp.createDataSet("residue_number", residue_number);
    grp.createDataSet("chain_id", chain_id);
}

// Write topology group: bond graph + ring membership.
// This is what forces the Python consumer to build a protein model.
void WriteTopology(HighFive::File& file,
                   const Protein& protein) {
    auto grp = file.createGroup("topology");

    // Bonds: (B, 2) atom index pairs + (B,) category + (B,) order.
    // Sorted by (atom_a, atom_b) for inspectability.
    const size_t B = protein.BondCount();
    if (B > 0) {
        // Build sorted index
        std::vector<size_t> order(B);
        std::iota(order.begin(), order.end(), 0);
        std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
            const Bond& ba = protein.BondAt(a);
            const Bond& bb = protein.BondAt(b);
            if (ba.atom_index_a != bb.atom_index_a)
                return ba.atom_index_a < bb.atom_index_a;
            return ba.atom_index_b < bb.atom_index_b;
        });

        std::vector<int32_t> bond_atoms(B * 2);
        std::vector<int32_t> bond_category(B);
        std::vector<int32_t> bond_order(B);
        for (size_t i = 0; i < B; ++i) {
            const Bond& b = protein.BondAt(order[i]);
            bond_atoms[i * 2 + 0] = static_cast<int32_t>(b.atom_index_a);
            bond_atoms[i * 2 + 1] = static_cast<int32_t>(b.atom_index_b);
            bond_category[i] = static_cast<int32_t>(b.category);
            bond_order[i] = static_cast<int32_t>(b.order);
        }
        HighFive::DataSpace space({B, size_t{2}});
        auto ds = grp.createDataSet<int32_t>("bond_atoms", space);
        ds.write_raw(bond_atoms.data());
        grp.createDataSet("bond_category", bond_category);
        grp.createDataSet("bond_order", bond_order);
        grp.createAttribute("n_bonds", B);
    }

    // Rings: per-ring type index, parent residue, atom count.
    // Ring atom indices stored as a ragged array: flat list + offsets.
    const size_t R = protein.RingCount();
    if (R > 0) {
        std::vector<int32_t> ring_type(R);
        std::vector<int32_t> ring_residue(R);
        std::vector<int32_t> ring_fused(R);
        std::vector<int32_t> ring_offsets(R + 1);  // CSR-style
        std::vector<int32_t> ring_atom_indices;

        ring_offsets[0] = 0;
        for (size_t i = 0; i < R; ++i) {
            const Ring& ring = protein.RingAt(i);
            ring_type[i] = static_cast<int32_t>(ring.type_index);
            ring_residue[i] = static_cast<int32_t>(ring.parent_residue_index);
            ring_fused[i] = (ring.fused_partner_index == SIZE_MAX)
                ? -1 : static_cast<int32_t>(ring.fused_partner_index);
            for (size_t ai : ring.atom_indices)
                ring_atom_indices.push_back(static_cast<int32_t>(ai));
            ring_offsets[i + 1] = static_cast<int32_t>(ring_atom_indices.size());
        }
        grp.createDataSet("ring_type", ring_type);
        grp.createDataSet("ring_residue", ring_residue);
        grp.createDataSet("ring_fused_partner", ring_fused);
        grp.createDataSet("ring_offsets", ring_offsets);
        grp.createDataSet("ring_atom_indices", ring_atom_indices);
        grp.createAttribute("n_rings", R);
    }

    // Parent atom index: hydrogen→heavy atom mapping.
    // SIZE_MAX means not a hydrogen (or parent not assigned).
    const size_t N = protein.AtomCount();
    std::vector<int32_t> parent(N);
    for (size_t i = 0; i < N; ++i) {
        size_t p = protein.AtomAt(i).parent_atom_index;
        parent[i] = (p == SIZE_MAX) ? -1 : static_cast<int32_t>(p);
    }
    grp.createDataSet("parent_atom_index", parent);
}

// Write a flat double buffer as a shaped HDF5 dataset.
void WriteRawDouble(HighFive::Group& grp, const std::string& name,
                    const double* data, std::vector<size_t> shape) {
    HighFive::DataSpace space(shape);
    auto ds = grp.createDataSet<double>(name, space);
    ds.write_raw(data);
}

// Write a SphericalTensor dataset with layout documentation attribute.
void WriteSphericalTensor(HighFive::Group& grp, const std::string& name,
                          const double* data, size_t T, size_t N) {
    HighFive::DataSpace space({T, N, size_t{9}});
    auto ds = grp.createDataSet<double>(name, space);
    ds.write_raw(data);
    ds.createAttribute("layout", std::string("T0,T1[3],T2[5]"));
    ds.createAttribute("convention", std::string(
        "T0=isotropic, T1[0..2]=antisymmetric, T2[0..4]=traceless symmetric (sphericart)"));
}

void WriteRawInt16(HighFive::Group& grp, const std::string& name,
                   const int16_t* data, std::vector<size_t> shape) {
    HighFive::DataSpace space(shape);
    auto ds = grp.createDataSet<int16_t>(name, space);
    ds.write_raw(data);
}

void WriteRawInt8(HighFive::Group& grp, const std::string& name,
                  const int8_t* data, std::vector<size_t> shape) {
    HighFive::DataSpace space(shape);
    auto ds = grp.createDataSet<int8_t>(name, space);
    ds.write_raw(data);
}

}  // anonymous namespace


// ── WriteH5 ─────────────────────────────────────────────────────

void AnalysisWriter::WriteH5(const std::string& path) const {
    using namespace HighFive;

    const size_t T = n_frames_;
    const size_t N = n_atoms_;

    File file(path, File::Truncate);

    // ── /meta/ ──────────────────────────────────────────────────
    auto meta = file.createGroup("meta");
    meta.createDataSet("frame_times", frame_times_);
    meta.createDataSet("frame_indices", frame_indices_);
    meta.createAttribute("protein_id", protein_id_);
    meta.createAttribute("n_atoms", N);
    meta.createAttribute("n_frames", T);
    meta.createAttribute("n_residues", protein_.ResidueCount());
    meta.createAttribute("stride", stride_);

    // ── /atoms/ — typed object model, serialised for R/Python ───
    WriteAtomMetadata(file, protein_, N);

    // ── /residues/ ──────────────────────────────────────────────
    WriteResidueMetadata(file, protein_);

    // ── /topology/ — bond graph, ring membership, H→heavy ───────
    WriteTopology(file, protein_);

    // ── /positions/ ─────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto pos_grp = file.createGroup("positions");
        WriteRawDouble(pos_grp, "xyz", positions_.data(), {T, N, 3});
    }

    // ── /ring_current/ ──────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("ring_current");

        // Per-type arrays
        WriteRawDouble(grp, "bs_T0_per_type", bs_T0_per_type_.data(), {T, N, 8});
        WriteRawDouble(grp, "bs_T2_per_type", bs_T2_per_type_.data(), {T, N, 8, 5});
        WriteRawDouble(grp, "hm_T0_per_type", hm_T0_per_type_.data(), {T, N, 8});
        WriteRawDouble(grp, "hm_T2_per_type", hm_T2_per_type_.data(), {T, N, 8, 5});

        // SphericalTensor totals (9 = T0 + T1[3] + T2[5]) with layout attr
        WriteSphericalTensor(grp, "bs_shielding", bs_shielding_.data(), T, N);
        WriteSphericalTensor(grp, "hm_shielding", hm_shielding_.data(), T, N);
        WriteSphericalTensor(grp, "rs_shielding", rs_shielding_.data(), T, N);

        // Ring proximity scalars
        WriteRawInt16(grp, "n_rings_3A", n_rings_3A_.data(), {T, N});
        WriteRawInt16(grp, "n_rings_5A", n_rings_5A_.data(), {T, N});
        WriteRawInt16(grp, "n_rings_8A", n_rings_8A_.data(), {T, N});
        WriteRawInt16(grp, "n_rings_12A", n_rings_12A_.data(), {T, N});
        WriteRawDouble(grp, "mean_ring_dist", mean_ring_dist_.data(), {T, N});
        WriteRawDouble(grp, "nearest_ring_atom", nearest_ring_atom_.data(), {T, N});

        // Exponential-weighted sums
        WriteRawDouble(grp, "G_iso_exp_sum", G_iso_exp_sum_.data(), {T, N});
        WriteRawDouble(grp, "G_T2_exp_sum", G_T2_exp_sum_.data(), {T, N, 5});
        WriteRawDouble(grp, "G_iso_var_8A", G_iso_var_8A_.data(), {T, N});

        // B-field vector
        WriteRawDouble(grp, "total_B_field", total_B_field_.data(), {T, N, 3});
    }

    // ── /efg/ ────────────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("efg");
        WriteSphericalTensor(grp, "coulomb_total", coulomb_efg_total_.data(), T, N);
        WriteSphericalTensor(grp, "coulomb_backbone", coulomb_efg_backbone_.data(), T, N);
        WriteSphericalTensor(grp, "coulomb_aromatic", coulomb_efg_aromatic_.data(), T, N);
        WriteRawDouble(grp, "E_total", coulomb_E_total_.data(), {T, N, 3});
        WriteRawDouble(grp, "E_backbone", coulomb_E_backbone_.data(), {T, N, 3});
        WriteRawDouble(grp, "E_sidechain", coulomb_E_sidechain_.data(), {T, N, 3});
        WriteRawDouble(grp, "E_aromatic", coulomb_E_aromatic_.data(), {T, N, 3});
        WriteRawDouble(grp, "E_solvent", coulomb_E_solvent_.data(), {T, N, 3});
        WriteRawDouble(grp, "E_magnitude", coulomb_E_magnitude_.data(), {T, N});
        WriteRawDouble(grp, "E_bond_proj", coulomb_E_bond_proj_.data(), {T, N});
        WriteRawDouble(grp, "E_backbone_frac", coulomb_E_bb_frac_.data(), {T, N});
        WriteSphericalTensor(grp, "apbs_efg", apbs_efg_.data(), T, N);
        WriteRawDouble(grp, "apbs_efield", apbs_efield_.data(), {T, N, 3});
        WriteSphericalTensor(grp, "aimnet2_total", aimnet2_efg_total_.data(), T, N);
        WriteSphericalTensor(grp, "aimnet2_backbone", aimnet2_efg_backbone_.data(), T, N);
        WriteSphericalTensor(grp, "aimnet2_aromatic", aimnet2_efg_aromatic_.data(), T, N);
        WriteSphericalTensor(grp, "aimnet2_shielding", aimnet2_shielding_.data(), T, N);
    }

    // ── /bond_aniso/ ────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("bond_aniso");
        WriteSphericalTensor(grp, "mc_shielding", mc_shielding_.data(), T, N);
        WriteSphericalTensor(grp, "T2_backbone", mc_T2_backbone_.data(), T, N);
        WriteSphericalTensor(grp, "T2_sidechain", mc_T2_sidechain_.data(), T, N);
        WriteSphericalTensor(grp, "T2_aromatic", mc_T2_aromatic_.data(), T, N);
        WriteSphericalTensor(grp, "T2_CO_nearest", mc_T2_CO_nearest_.data(), T, N);
        WriteSphericalTensor(grp, "T2_CN_nearest", mc_T2_CN_nearest_.data(), T, N);
        WriteRawDouble(grp, "co_sum", mc_co_sum_.data(), {T, N});
        WriteRawDouble(grp, "cn_sum", mc_cn_sum_.data(), {T, N});
        WriteRawDouble(grp, "sidechain_sum", mc_sidechain_sum_.data(), {T, N});
        WriteRawDouble(grp, "aromatic_sum", mc_aromatic_sum_.data(), {T, N});
        WriteRawDouble(grp, "co_nearest", mc_co_nearest_.data(), {T, N});
        WriteRawDouble(grp, "nearest_CO_dist", mc_nearest_CO_dist_.data(), {T, N});
        WriteRawDouble(grp, "nearest_CN_dist", mc_nearest_CN_dist_.data(), {T, N});
        WriteRawDouble(grp, "nearest_CO_midpoint", mc_nearest_CO_mid_.data(), {T, N, 3});
        WriteRawDouble(grp, "dir_nearest_CO", mc_dir_nearest_CO_.data(), {T, N, 3});
    }

    // ── /quadrupole/ ────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("quadrupole");
        WriteSphericalTensor(grp, "pq_shielding", pq_shielding_.data(), T, N);
        WriteRawDouble(grp, "pq_T0_per_type", pq_T0_per_type_.data(), {T, N, 8});
        WriteRawDouble(grp, "pq_T2_per_type", pq_T2_per_type_.data(), {T, N, 8, 5});
    }

    // ── /dispersion/ ────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("dispersion");
        WriteSphericalTensor(grp, "disp_shielding", disp_shielding_.data(), T, N);
        WriteRawDouble(grp, "disp_T0_per_type", disp_T0_per_type_.data(), {T, N, 8});
        WriteRawDouble(grp, "disp_T2_per_type", disp_T2_per_type_.data(), {T, N, 8, 5});
    }

    // ── /hbond/ ─────────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("hbond");
        WriteSphericalTensor(grp, "hbond_shielding", hbond_shielding_.data(), T, N);
        WriteSphericalTensor(grp, "nearest_spherical", hbond_nearest_sph_.data(), T, N);
        WriteRawDouble(grp, "nearest_dist", hbond_nearest_dist_.data(), {T, N});
        WriteRawDouble(grp, "nearest_dir", hbond_nearest_dir_.data(), {T, N, 3});
        WriteRawDouble(grp, "inv_d3", hbond_inv_d3_.data(), {T, N});
        WriteRawInt16(grp, "count_3_5A", hbond_count_.data(), {T, N});
        WriteRawInt8(grp, "is_donor", hbond_is_donor_.data(), {T, N});
        WriteRawInt8(grp, "is_acceptor", hbond_is_acceptor_.data(), {T, N});
        WriteRawInt8(grp, "is_backbone", hbond_is_backbone_.data(), {T, N});
    }

    // ── /sasa/ ──────────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("sasa");
        WriteRawDouble(grp, "sasa", sasa_.data(), {T, N});
        WriteRawDouble(grp, "normal", sasa_normal_.data(), {T, N, 3});
    }

    // ── /water/ ─────────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("water");
        WriteRawDouble(grp, "efield", water_efield_.data(), {T, N, 3});
        WriteSphericalTensor(grp, "efg", water_efg_.data(), T, N);
        WriteRawDouble(grp, "efield_first", water_efield_first_.data(), {T, N, 3});
        WriteSphericalTensor(grp, "efg_first", water_efg_first_.data(), T, N);
        WriteRawInt16(grp, "n_first", water_n_first_.data(), {T, N});
        WriteRawInt16(grp, "n_second", water_n_second_.data(), {T, N});
        WriteRawDouble(grp, "half_shell_asymmetry", water_half_shell_.data(), {T, N});
        WriteRawDouble(grp, "dipole_cos", water_dipole_cos_.data(), {T, N});
        WriteRawDouble(grp, "nearest_ion_dist", water_nearest_ion_dist_.data(), {T, N});
        WriteRawDouble(grp, "nearest_ion_charge", water_nearest_ion_charge_.data(), {T, N});
        WriteRawDouble(grp, "dipole_vector", water_dipole_vec_.data(), {T, N, 3});
        WriteRawDouble(grp, "surface_normal", water_surface_normal_.data(), {T, N, 3});
        WriteRawDouble(grp, "sasa_asymmetry", water_sasa_asym_.data(), {T, N});
        WriteRawDouble(grp, "sasa_dipole_align", water_sasa_align_.data(), {T, N});
        WriteRawDouble(grp, "sasa_dipole_cohere", water_sasa_cohere_.data(), {T, N});
        WriteRawInt16(grp, "sasa_first_shell_n", water_sasa_first_n_.data(), {T, N});
    }

    // ── /charges/ ───────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("charges");
        WriteRawDouble(grp, "aimnet2_charge", aimnet2_charge_.data(), {T, N});
        WriteRawDouble(grp, "eeq_charge", eeq_charge_.data(), {T, N});
        WriteRawDouble(grp, "eeq_cn", eeq_cn_.data(), {T, N});
    }

    // ── /aimnet2_embedding/ ─────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("aimnet2_embedding");
        WriteRawDouble(grp, "aim", aimnet2_aim_.data(), {T, N, size_t{AIMNET2_AIM_DIMS}});
    }

    // ── /per_ring/ ────────────────────────────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.createGroup("per_ring");
        grp.createAttribute("K", K_RINGS);
        WriteRawDouble(grp, "geometry", per_ring_geom_.data(),
                       {T, N, K_RINGS, size_t{6}});
        grp.createDataSet("geometry_fields",
            std::vector<std::string>{"distance", "rho", "z", "theta", "cos_phi", "sin_phi"});
        WriteRawInt8(grp, "ring_type", per_ring_type_.data(), {T, N, K_RINGS});
        WriteRawDouble(grp, "bs_T2", per_ring_bs_T2_.data(), {T, N, K_RINGS, size_t{5}});
        WriteRawDouble(grp, "hm_T2", per_ring_hm_T2_.data(), {T, N, K_RINGS, size_t{5}});
        WriteRawDouble(grp, "chi_T2", per_ring_chi_T2_.data(), {T, N, K_RINGS, size_t{5}});
        WriteRawDouble(grp, "pq_T2", per_ring_pq_T2_.data(), {T, N, K_RINGS, size_t{5}});
        WriteRawDouble(grp, "hm_H_T2", per_ring_hm_H_T2_.data(), {T, N, K_RINGS, size_t{5}});
        WriteRawDouble(grp, "disp_scalar", per_ring_disp_scalar_.data(), {T, N, K_RINGS});
        WriteRawDouble(grp, "disp_T2", per_ring_disp_T2_.data(), {T, N, K_RINGS, size_t{5}});
    }

    // ── Additional /atoms/ datasets (enrichment + graph + charges) ──
    if (enrichment_captured_) {
        auto grp = file.getGroup("atoms");
        grp.createDataSet("is_amide_H", is_amide_H_);
        grp.createDataSet("is_alpha_H", is_alpha_H_);
        grp.createDataSet("is_methyl", is_methyl_);
        grp.createDataSet("is_aromatic_H", is_aromatic_H_);
        grp.createDataSet("is_on_aromatic_residue", is_on_aromatic_residue_);
        grp.createDataSet("is_hbond_donor", is_hbond_donor_);
        grp.createDataSet("is_hbond_acceptor", is_hbond_acceptor_);
        grp.createDataSet("parent_is_sp2", parent_is_sp2_);
        grp.createDataSet("graph_dist_N", graph_dist_N_);
        grp.createDataSet("graph_dist_O", graph_dist_O_);
        grp.createDataSet("eneg_sum_1", eneg_sum_1_);
        grp.createDataSet("eneg_sum_2", eneg_sum_2_);
        grp.createDataSet("n_pi_bonds_3", n_pi_bonds_3_);
        grp.createDataSet("bfs_to_nearest_ring", bfs_to_nearest_ring_);
        grp.createDataSet("bfs_decay", bfs_decay_);
        grp.createDataSet("partial_charge", partial_charge_);
        grp.createDataSet("vdw_radius", vdw_radius_);
    }

    // ── /ring_geometry/ (per-frame ring centers, normals, radii) ──
    {
        const size_t R = protein_.RingCount();
        if (T > 0 && R > 0) {
            auto grp = file.createGroup("ring_geometry");
            WriteRawDouble(grp, "data", ring_geometry_.data(), {T, R, size_t{7}});
            grp.createDataSet("fields",
                std::vector<std::string>{"center_x", "center_y", "center_z",
                                         "normal_x", "normal_y", "normal_z", "radius"});
        }
    }

    // ── Coulomb shielding in /efg/ group ────────────────────────
    if (T > 0 && N > 0) {
        auto grp = file.getGroup("efg");
        WriteSphericalTensor(grp, "coulomb_shielding",
                             coulomb_shielding_.data(), T, N);
    }

    OperationLog::Info(LogCalcOther, "AnalysisWriter::WriteH5",
        path + ": " + std::to_string(T) + " frames, " +
        std::to_string(N) + " atoms");
}

}  // namespace nmr
