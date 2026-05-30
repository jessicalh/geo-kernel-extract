// QtTrajectoryH5 — typed read boundary for trajectory.h5 (the
// per-TR-emits-its-own-group H5 format produced by
// `nmr_extract --trajectory` post-2026-05-13).
//
// Mirrors `ui/src/TrajectoryH5` in eager-load discipline + sparse-set
// tolerance + WarnShapeMismatch shape, but extended for the trajectory
// animator: every present TR group is eager-loaded in full (N, T, K)
// across all frames, not just frame-0. HighFive::File handle dies in
// the constructor; the buffers outlive the file.
//
// 56 TR groups in the 1P9J fixture (Agent 3's deep-dive). Each gets
// its own typed buffer accessor below. Sparse-tolerant via raw pointer
// (nullable) returns; absent groups return nullptr.
//
// Construction throws `std::runtime_error` on structural failure
// (missing /atoms, missing /trajectory/frames). Per-TR group failures
// (present but malformed) log Warn via OperationLog and leave the
// accessor nullptr. The "absent, not faked" discipline carries forward.
//
// Memory footprint on 1P9J fixture (846 atoms × 751 frames): ~1.6 GB
// once all TRs are loaded — close to the H5 file size, since the
// format overhead is the only compression. Fits 128 GB workstation
// comfortably; advisor laptops with <16 GB RAM should opt out of the
// embedding TR (256-dim × N × T × float32 = 800 MB by itself).

#pragma once

#include "../model/QtBondVectorBuffers.h"
#include "../model/QtPerAtomChannelBuffers.h"
#include "../model/QtPerResidueBuffers.h"
#include "../model/QtSelectionBag.h"
#include "../model/QtSpecialBuffers.h"
#include "../model/QtTimeSeriesBuffers.h"
#include "../model/QtWelfordBuffers.h"

#include <QString>
#include <QStringList>
#include <cstddef>
#include <memory>
#include <vector>

namespace h5reader::io {

class QtTrajectoryH5 {
public:
    // Throws std::runtime_error on structural failure. The h5_path
    // is the absolute path to trajectory.h5.
    explicit QtTrajectoryH5(const QString& h5_path);

    ~QtTrajectoryH5() = default;
    QtTrajectoryH5(const QtTrajectoryH5&) = delete;
    QtTrajectoryH5& operator=(const QtTrajectoryH5&) = delete;
    QtTrajectoryH5(QtTrajectoryH5&&) = default;
    QtTrajectoryH5& operator=(QtTrajectoryH5&&) = default;

    // ── Root metadata ─────────────────────────────────────────────
    std::size_t atomCount() const { return n_atoms_; }
    std::size_t frameCount() const { return n_frames_; }
    const QString& proteinId() const { return protein_id_; }

    // ── /atoms minimal identity (cross-check against sidecar) ─────
    const std::vector<int32_t>& atomElement() const { return atom_element_; }
    const std::vector<uint64_t>& atomResidueIndex() const { return atom_residue_index_; }
    const std::vector<QString>& atomPdbName() const { return atom_pdb_name_; }

    // ── /trajectory/frames ────────────────────────────────────────
    const std::vector<double>& frameTimes() const { return frame_times_; }
    const std::vector<uint64_t>& frameIndices() const { return frame_indices_; }

    // ── /trajectory attrs (source provenance) ─────────────────────
    const QString& xtcPath() const { return xtc_path_; }
    const QString& tprPath() const { return tpr_path_; }
    const QString& edrPath() const { return edr_path_; }
    const QString& configuration() const { return configuration_; }

    // ── Always-present TR (positions drives the animator) ─────────
    // Construction throws if /trajectory/positions is absent or does
    // not match the /atoms and /trajectory/frames dimensions.
    const h5reader::model::QtPositionsTimeSeries* positions() const { return positions_.get(); }

    // ── Shielding time-series family (13 TRs, all (N, T, 9)) ──────
    const h5reader::model::QtShieldingTimeSeries* bsShielding() const { return bs_shielding_.get(); }
    const h5reader::model::QtShieldingTimeSeries* hmShielding() const { return hm_shielding_.get(); }
    const h5reader::model::QtShieldingTimeSeries* mcShielding() const { return mc_shielding_.get(); }
    const h5reader::model::QtShieldingTimeSeries* piQuadShielding() const { return piquad_shielding_.get(); }
    const h5reader::model::QtShieldingTimeSeries* ringChiShielding() const { return ringchi_shielding_.get(); }
    const h5reader::model::QtShieldingTimeSeries* dispShielding() const { return disp_shielding_.get(); }
    const h5reader::model::QtShieldingTimeSeries* hbondShielding() const { return hbond_shielding_.get(); }
    const h5reader::model::QtT2TimeSeries* mopacCoulombShielding() const { return mopac_coulomb_shielding_.get(); }
    const h5reader::model::QtShieldingTimeSeries* mopacMcShielding() const { return mopac_mc_shielding_.get(); }
    const h5reader::model::QtScalarTimeSeries* mopacVsFf14sbReconciliation() const {
        return mopac_vs_ff14sb_reconciliation_.get();
    }
    const h5reader::model::QtShieldingTimeSeries* tripeptideBbShielding() const { return tripeptide_bb_shielding_.get(); }
    const h5reader::model::QtShieldingTimeSeries* tripeptideNeighborShielding() const {
        return tripeptide_neighbor_shielding_.get();
    }
    const h5reader::model::QtShieldingTimeSeries* larsenHBond1pHBShielding() const { return larsen_1pHB_.get(); }
    const h5reader::model::QtShieldingTimeSeries* larsenHBond1pHaBShielding() const { return larsen_1pHaB_.get(); }
    const h5reader::model::QtShieldingTimeSeries* larsenHBond2pHBShielding() const { return larsen_2pHB_.get(); }
    const h5reader::model::QtShieldingTimeSeries* larsenHBond2pHaBShielding() const { return larsen_2pHaB_.get(); }

    // ── Scalar time-series (N, T) ─────────────────────────────────
    const h5reader::model::QtScalarTimeSeries* sasa() const { return sasa_.get(); }
    const h5reader::model::QtScalarTimeSeries* aimnet2Charge() const { return aimnet2_charge_.get(); }
    const h5reader::model::QtScalarTimeSeries* jCoupling() const { return j_coupling_.get(); }
    const h5reader::model::QtScalarTimeSeries* larsenHBondCount() const { return larsen_hbond_count_.get(); }
    const h5reader::model::QtScalarTimeSeries* larsenHBondWaterTerm() const { return larsen_hbond_water_term_.get(); }

    // ── Vec3 time-series (N, T, 3) ────────────────────────────────
    const h5reader::model::QtVec3TimeSeries* apbsEfield() const { return apbs_efield_.get(); }
    const h5reader::model::QtVec3TimeSeries* tripeptideBbResidualVec() const { return tripeptide_bb_residual_vec_.get(); }
    const h5reader::model::QtVec3TimeSeries* tripeptideNeighborResidualVecPrev() const {
        return tripeptide_neighbor_residual_vec_prev_.get();
    }
    const h5reader::model::QtVec3TimeSeries* tripeptideNeighborResidualVecNext() const {
        return tripeptide_neighbor_residual_vec_next_.get();
    }

    // ── Special shapes ────────────────────────────────────────────
    const h5reader::model::QtT2TimeSeries* apbsEfg() const { return apbs_efg_.get(); }
    const h5reader::model::QtTagTimeSeries* tripeptideBbMethodTag() const { return tripeptide_bb_method_tag_.get(); }
    const h5reader::model::QtEmbeddingTimeSeries* aimnet2Embedding() const { return aimnet2_embedding_.get(); }
    const h5reader::model::QtAimnet2ChargeResponseGradientTimeSeries* aimnet2ChargeResponseGradient() const {
        return aimnet2_crg_.get();
    }

    // ── Bonded energy (per-atom split) ────────────────────────────
    const h5reader::model::QtScalarTimeSeries* bondedEnergyTotal() const { return bonded_energy_total_.get(); }

    // ── Per-residue TRs ───────────────────────────────────────────
    const h5reader::model::QtDihedralTimeSeries* dihedrals() const { return dihedrals_.get(); }
    const h5reader::model::QtDssp8TimeSeries* dssp8() const { return dssp8_.get(); }
    const h5reader::model::QtRingPuckerTimeSeries* ringPucker() const { return ring_pucker_.get(); }
    const h5reader::model::QtJCouplingTimeSeries* jCouplingTimeSeries() const { return j_coupling_full_.get(); }

    // ── Special-shape TRs ─────────────────────────────────────────
    const h5reader::model::QtRingNeighbourhoodTimeSeries* ringNeighbourhood() const { return ring_neighbourhood_.get(); }
    const h5reader::model::QtBondLengthStats* bondLengthStats() const { return bond_length_stats_.get(); }
    const h5reader::model::QtSystemEnergyTimeSeries* gromacsEnergy() const { return gromacs_energy_.get(); }
    const h5reader::model::QtRmsdTracking* rmsdTracking() const { return rmsd_tracking_.get(); }

    // ── Welford rollups (atom-axis, no T) ─────────────────────────
    const h5reader::model::QtShieldingWelford* bsWelford() const { return bs_welford_.get(); }
    const h5reader::model::QtShieldingWelford* hmWelford() const { return hm_welford_.get(); }
    const h5reader::model::QtShieldingWelford* mcWelford() const { return mc_welford_.get(); }
    const h5reader::model::QtScalarWelford* sasaWelford() const { return sasa_welford_.get(); }
    const h5reader::model::QtScalarWelford* eeqWelford() const { return eeq_welford_.get(); }
    const h5reader::model::QtScalarWelford* hbondCountWelford() const { return hbond_count_welford_.get(); }
    const h5reader::model::QtScalarWelford* mopacChargeWelford() const { return mopac_charge_welford_.get(); }
    const h5reader::model::QtBondOrderWelford* mopacBondOrderWelford() const { return mopac_bond_order_welford_.get(); }
    const h5reader::model::QtVec3Welford* waterFieldWelford() const { return water_field_welford_.get(); }
    const h5reader::model::QtVec3Welford* aimnet2ChargeResponseGradientWelford() const { return aimnet2_crg_welford_.get(); }
    const h5reader::model::QtHydrationWelford* hydrationShellWelford() const { return hydration_shell_welford_.get(); }
    const h5reader::model::QtHydrationWelford* hydrationGeometryWelford() const { return hydration_geometry_welford_.get(); }
    const h5reader::model::QtAutocorrelation* bsT0Autocorrelation() const { return bs_t0_autocorrelation_.get(); }

    // ── Bond-vector axis (TR-specific identity tables) ────────────
    const h5reader::model::QtIRedOrderParameters* iredOrderParameters() const {
        return ired_order_parameters_.get();
    }

    // ── Per-atom × per-channel composites ─────────────────────────
    const h5reader::model::QtKernelDynamics* kernelDynamics() const {
        return kernel_dynamics_.get();
    }
    const h5reader::model::QtReorientationalDynamics* reorientationalDynamics() const {
        return reorientational_dynamics_.get();
    }
    const h5reader::model::QtDihedralAutocorrelation* dihedralAutocorrelation() const {
        return dihedral_autocorrelation_.get();
    }
    const h5reader::model::QtKernelCoherence* kernelCoherence() const {
        return kernel_coherence_.get();
    }

    // ── Transition trackers ───────────────────────────────────────
    const h5reader::model::QtDssp8Transitions* dssp8Transitions() const { return dssp8_transitions_.get(); }
    const h5reader::model::QtDihedralBinTransitions* dihedralBinTransitions() const { return dihedral_bin_transitions_.get(); }

    // ── Hydration + water field per-frame (composite buffers) ─────
    const h5reader::model::QtHydrationShellTimeSeries* hydrationShellTimeSeries() const { return hydration_shell_ts_.get(); }
    const h5reader::model::QtHydrationGeometryTimeSeries* hydrationGeometryTimeSeries() const {
        return hydration_geometry_ts_.get();
    }
    const h5reader::model::QtWaterFieldTimeSeries* waterFieldTimeSeries() const { return water_field_ts_.get(); }

    // ── Selections (typed event bag) ──────────────────────────────
    const h5reader::model::QtSelectionBag& selections() const { return selections_; }

    // ── Sparse-set introspection ──────────────────────────────────
    const QStringList& groupsPresent() const { return groups_present_; }
    bool hasGroup(const QString& name) const { return groups_present_.contains(name); }

private:
    std::size_t n_atoms_ = 0;
    std::size_t n_frames_ = 0;
    QString protein_id_;

    std::vector<int32_t> atom_element_;
    std::vector<uint64_t> atom_residue_index_;
    std::vector<QString> atom_pdb_name_;

    std::vector<double> frame_times_;
    std::vector<uint64_t> frame_indices_;

    QString xtc_path_, tpr_path_, edr_path_, configuration_;

    // Inventory of /trajectory child group names actually present.
    QStringList groups_present_;

    // ── Per-TR buffer slots ───────────────────────────────────────
    std::unique_ptr<h5reader::model::QtPositionsTimeSeries> positions_;

    // Shielding family
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> bs_shielding_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> hm_shielding_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> mc_shielding_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> piquad_shielding_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> ringchi_shielding_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> disp_shielding_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> hbond_shielding_;
    std::unique_ptr<h5reader::model::QtT2TimeSeries> mopac_coulomb_shielding_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> mopac_mc_shielding_;
    std::unique_ptr<h5reader::model::QtScalarTimeSeries> mopac_vs_ff14sb_reconciliation_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> tripeptide_bb_shielding_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> tripeptide_neighbor_shielding_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> larsen_1pHB_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> larsen_1pHaB_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> larsen_2pHB_;
    std::unique_ptr<h5reader::model::QtShieldingTimeSeries> larsen_2pHaB_;
    std::unique_ptr<h5reader::model::QtWaterFieldTimeSeries> water_field_ts_;

    // Scalar TS family
    std::unique_ptr<h5reader::model::QtScalarTimeSeries> sasa_;
    std::unique_ptr<h5reader::model::QtScalarTimeSeries> aimnet2_charge_;
    std::unique_ptr<h5reader::model::QtScalarTimeSeries> j_coupling_;
    std::unique_ptr<h5reader::model::QtScalarTimeSeries> larsen_hbond_count_;
    std::unique_ptr<h5reader::model::QtScalarTimeSeries> larsen_hbond_water_term_;
    std::unique_ptr<h5reader::model::QtScalarTimeSeries> bonded_energy_total_;
    std::unique_ptr<h5reader::model::QtHydrationShellTimeSeries> hydration_shell_ts_;
    std::unique_ptr<h5reader::model::QtHydrationGeometryTimeSeries> hydration_geometry_ts_;

    // Vec3 TS family
    std::unique_ptr<h5reader::model::QtVec3TimeSeries> apbs_efield_;
    std::unique_ptr<h5reader::model::QtVec3TimeSeries> tripeptide_bb_residual_vec_;
    std::unique_ptr<h5reader::model::QtVec3TimeSeries> tripeptide_neighbor_residual_vec_prev_;
    std::unique_ptr<h5reader::model::QtVec3TimeSeries> tripeptide_neighbor_residual_vec_next_;

    // Special TS shapes
    std::unique_ptr<h5reader::model::QtT2TimeSeries> apbs_efg_;
    std::unique_ptr<h5reader::model::QtTagTimeSeries> tripeptide_bb_method_tag_;
    std::unique_ptr<h5reader::model::QtEmbeddingTimeSeries> aimnet2_embedding_;
    std::unique_ptr<h5reader::model::QtAimnet2ChargeResponseGradientTimeSeries> aimnet2_crg_;

    // Per-residue
    std::unique_ptr<h5reader::model::QtDihedralTimeSeries> dihedrals_;
    std::unique_ptr<h5reader::model::QtDssp8TimeSeries> dssp8_;
    std::unique_ptr<h5reader::model::QtRingPuckerTimeSeries> ring_pucker_;
    std::unique_ptr<h5reader::model::QtJCouplingTimeSeries> j_coupling_full_;

    // Special-shape TRs
    std::unique_ptr<h5reader::model::QtRingNeighbourhoodTimeSeries> ring_neighbourhood_;
    std::unique_ptr<h5reader::model::QtBondLengthStats> bond_length_stats_;
    std::unique_ptr<h5reader::model::QtSystemEnergyTimeSeries> gromacs_energy_;
    std::unique_ptr<h5reader::model::QtRmsdTracking> rmsd_tracking_;

    // Welford
    std::unique_ptr<h5reader::model::QtShieldingWelford> bs_welford_;
    std::unique_ptr<h5reader::model::QtShieldingWelford> hm_welford_;
    std::unique_ptr<h5reader::model::QtShieldingWelford> mc_welford_;
    std::unique_ptr<h5reader::model::QtScalarWelford> sasa_welford_;
    std::unique_ptr<h5reader::model::QtScalarWelford> eeq_welford_;
    std::unique_ptr<h5reader::model::QtScalarWelford> hbond_count_welford_;
    std::unique_ptr<h5reader::model::QtScalarWelford> mopac_charge_welford_;
    std::unique_ptr<h5reader::model::QtBondOrderWelford> mopac_bond_order_welford_;
    std::unique_ptr<h5reader::model::QtVec3Welford> water_field_welford_;
    std::unique_ptr<h5reader::model::QtVec3Welford> aimnet2_crg_welford_;
    std::unique_ptr<h5reader::model::QtHydrationWelford> hydration_shell_welford_;
    std::unique_ptr<h5reader::model::QtHydrationWelford> hydration_geometry_welford_;
    std::unique_ptr<h5reader::model::QtAutocorrelation> bs_t0_autocorrelation_;

    // Bond-vector axis (each TR owns its identity table; see QtBondVectorBuffers.h).
    std::unique_ptr<h5reader::model::QtIRedOrderParameters> ired_order_parameters_;

    // Per-atom × per-channel composites.
    std::unique_ptr<h5reader::model::QtKernelDynamics> kernel_dynamics_;
    std::unique_ptr<h5reader::model::QtReorientationalDynamics> reorientational_dynamics_;
    std::unique_ptr<h5reader::model::QtDihedralAutocorrelation> dihedral_autocorrelation_;
    std::unique_ptr<h5reader::model::QtKernelCoherence> kernel_coherence_;

    // Transitions
    std::unique_ptr<h5reader::model::QtDssp8Transitions> dssp8_transitions_;
    std::unique_ptr<h5reader::model::QtDihedralBinTransitions> dihedral_bin_transitions_;

    // Selections
    h5reader::model::QtSelectionBag selections_;
};

}  // namespace h5reader::io
