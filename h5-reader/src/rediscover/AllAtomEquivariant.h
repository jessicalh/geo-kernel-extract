// AllAtomEquivariant -- all-atom, all-category lab-frame geometry emit.
//
// This is the corrected e3nn substrate: raw source geometry plus the DFT
// shielding target in one consistent molecular/lab frame. It intentionally does
// not build or apply per-atom local frames.

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
    std::size_t apbs_efield_rows = 0;
    std::size_t apbs_efg_rows = 0;
    std::size_t aimnet2_atom_rows = 0;
    std::size_t aimnet2_embedding_present = 0;

    std::array<std::size_t, 9> ring_type_rows = {};
    std::array<std::size_t, 8> bond_category_rows = {};
};

AllAtomEquivariantStats RunAllAtomEquivariantEmit(const Body& body,
                                                  AllAtomEquivariantSink& sink,
                                                  const AllAtomEquivariantConfig& config);

}  // namespace h5reader::rediscover
