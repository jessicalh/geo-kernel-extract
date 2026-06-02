// AllAtomEquivariant -- all-atom, all-category lab-frame geometry emit.
//
// This is the corrected e3nn substrate: raw source geometry plus the DFT
// shielding target in one consistent molecular/lab frame. It intentionally does
// not build or apply per-atom local frames.
//
// ORIENTATION PARITY CONTRACT (#51 C6 resolution — read before using
// `orientation_a` as an e3nn irrep). The source `orientation_a` (disp_*, r, etc.
// are genuine polar 1o vectors; orientation_a is NOT, for two source kinds):
//   - BOND sources (mechanism="bond"): orientation_a is the unit bond axis
//     (posB − posA) in the producer's index-oriented endpoint order (min/max
//     atom index), so its per-bond SIGN IS ARBITRARY. Do NOT consume it as a
//     polar 1o vector. It is an UNORIENTED axis: use its even form (the outer
//     product n⊗n, a 0e+2e irrep) or tag it index-oriented. The McConnell tensor
//     (3cos²θ−1)/r³ and `dipolar` are already even in the axis sign and are the
//     correct way to use the bond direction; orientation_a is carried only for
//     provenance / an even-form construction downstream.
//   - RING sources (mechanism="ring"): orientation_a is the ring NORMAL, an
//     AXIAL pseudovector (the producer's canonical normal flip fixes a
//     convention but the l=2 ring-current form is even under n → −n). Treat it
//     as a pseudovector / use its even outer-product form, NOT a polar 1o.
// disp_* (target→source) is the one genuine polar 1o orientation on every
// source row; cos_theta/dipolar are even (axis-sign-invariant) by construction.

#pragma once

#include "AllAtomEquivariantSink.h"
#include "AnalysisBody.h"

#include <array>
#include <cstddef>

namespace h5reader::rediscover {

struct AllAtomEquivariantConfig {
    double ring_cutoff_A = 8.0;
    double bond_cutoff_A = 10.0;
    double charge_cutoff_A = 6.0;
    double mc_near_field_ratio = 0.5;
};

struct AllAtomEquivariantStats {
    std::size_t atom_count = 0;
    std::size_t dft_rows = 0;
    std::size_t target_rows = 0;
    std::size_t dft_present = 0;
    std::size_t source_rows = 0;
    std::size_t ring_rows = 0;
    std::size_t bond_rows = 0;
    std::size_t charge_ff14sb_rows = 0;
    std::size_t charge_aimnet2_rows = 0;
    std::size_t charge_mopac_rows = 0;
    std::size_t apbs_efield_rows = 0;
    std::size_t apbs_efg_rows = 0;
    std::size_t aimnet2_atom_rows = 0;
    std::size_t aimnet2_embedding_present = 0;
    std::size_t mopac_coulomb_shielding_rows = 0;
    std::size_t mopac_mc_shielding_rows = 0;

    std::array<std::size_t, 9> ring_type_rows = {};
    std::array<std::size_t, 8> bond_category_rows = {};
};

AllAtomEquivariantStats RunAllAtomEquivariantEmit(const Body& body,
                                                  AllAtomEquivariantSink& sink,
                                                  const AllAtomEquivariantConfig& config);

}  // namespace h5reader::rediscover
