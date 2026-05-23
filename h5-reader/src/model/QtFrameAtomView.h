// QtFrameAtomView — materialised typed view of one (atom, frame) pair.
//
// Where QtFrame's accessors return default-zero on absent TRs (the
// "absent-is-zero" backward-compat policy), QtFrameAtomView returns
// std::optional<T> for every per-TR field so the inspector dock can
// distinguish "no source attached" from "source attached, value is
// zero".
//
// Used by the atom-inspector dock to display one atom at the current
// frame. Materialised on-demand via `QtFrame::atomView(atomIdx)` (not
// pre-built — N×T views would explode memory).
//
// Composes typed substrate (via QtAtom from QtProtein) + per-frame
// values (via QtFrame's underlying QtTrajectoryH5 buffers) into one
// snapshot for the inspector tree.

#pragma once

#include "Types.h"

#include <QString>
#include <cstddef>
#include <optional>

namespace h5reader::model {

class QtAtom;
class QtAtomNames;

struct QtFrameAtomView {
    // Identity (from QtProtein — never changes per frame)
    const QtAtom* atom = nullptr;
    const QtAtomNames* names = nullptr;

    // Frame identity
    std::size_t frame_index = 0;
    double time_ps = 0.0;

    // Per-frame position (always present from positions TR)
    Vec3 position = Vec3::Zero();

    // Per-frame computed values — std::nullopt when source TR absent
    // OR conditional-attach mask says source wasn't attached at this frame.
    std::optional<SphericalTensor> bs_shielding;
    std::optional<SphericalTensor> hm_shielding;
    std::optional<SphericalTensor> mc_shielding;
    std::optional<SphericalTensor> piquad_shielding;
    std::optional<SphericalTensor> ringchi_shielding;
    std::optional<SphericalTensor> disp_shielding;
    std::optional<SphericalTensor> hbond_shielding;
    std::optional<SphericalTensor> mopac_coulomb_shielding;
    std::optional<SphericalTensor> mopac_mc_shielding;
    std::optional<SphericalTensor> mopac_vs_ff14sb_reconciliation;
    std::optional<SphericalTensor> tripeptide_bb_shielding;
    std::optional<SphericalTensor> tripeptide_neighbor_shielding;
    std::optional<SphericalTensor> larsen_1pHB_shielding;
    std::optional<SphericalTensor> larsen_1pHaB_shielding;
    std::optional<SphericalTensor> larsen_2pHB_shielding;
    std::optional<SphericalTensor> larsen_2pHaB_shielding;
    std::optional<SphericalTensor> water_field_shielding;

    std::optional<Vec3> apbs_efield;
    std::optional<double> apbs_efg_t2_magnitude;
    std::optional<double> sasa;
    std::optional<double> aimnet2_charge;
    std::optional<Vec3> aimnet2_charge_response_gradient_vec;
    std::optional<double> aimnet2_charge_response_gradient_scalar;
    std::optional<double> j_coupling;
    std::optional<double> larsen_hbond_count;
    std::optional<double> larsen_hbond_water_term;
    std::optional<double> bonded_energy_total;

    std::optional<Vec3> tripeptide_bb_residual_vec;
    std::optional<Vec3> tripeptide_neighbor_residual_vec_prev;
    std::optional<Vec3> tripeptide_neighbor_residual_vec_next;

    std::optional<int> tripeptide_bb_method_tag;

    // Source-attached masks: true when the TR was present AND
    // source_attached_per_frame[t] != 0 at this frame.
    std::optional<bool> sasa_source_attached;
    std::optional<bool> tripeptide_source_attached;
};

}  // namespace h5reader::model
