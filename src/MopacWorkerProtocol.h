#pragma once
//
// Private wire protocol between MopacResult and a fresh libmopac worker.
// The worker is a separate executable so a libmopac abort/segfault cannot
// terminate the extractor and so no post-fork allocator/OpenMP state is used.
//

#include <array>
#include <cstdint>
#include <string>
#include <vector>

namespace nmr::mopac_worker {

struct Input {
    std::int32_t net_charge = 0;
    std::vector<std::int32_t> atoms;
    std::vector<double> coordinates;
};

struct Output {
    std::int32_t natom = 0;
    std::int32_t ao_max_orbitals = 0;
    // [numat,noccupied,nvirtual,icocc_dim,icvir_dim,cocc_dim,cvir_dim]
    std::array<std::int32_t, 7> state_dimensions{};

    double heat = 0.0;
    std::array<double, 3> dipole{};
    std::array<double, 3> dipole_point_charge{};
    std::array<double, 3> dipole_hybridization{};

    std::vector<double> charge;
    std::vector<std::int32_t> bond_index;
    std::vector<std::int32_t> bond_atom;
    std::vector<double> bond_order;
    std::vector<std::int32_t> ao_orbitals;
    std::vector<double> atom_ao_density_fortran;
    std::vector<double> bond_ao_density_fortran;

    std::vector<std::int32_t> nbonds;
    std::vector<std::int32_t> ibonds_fortran;
    std::vector<std::int32_t> iorbs;
    std::vector<std::int32_t> ncf;
    std::vector<std::int32_t> nce;
    std::vector<std::int32_t> icocc;
    std::vector<std::int32_t> icvir;
    std::vector<double> cocc;
    std::vector<double> cvir;

    std::vector<double> lmo_energy;
    std::vector<std::int32_t> occupied_atom_offsets;
    std::vector<std::int32_t> virtual_atom_offsets;
    std::vector<std::int32_t> occupied_coefficient_offsets;
    std::vector<std::int32_t> virtual_coefficient_offsets;
};

bool WriteInput(int fd, const Input& input);
bool ReadInput(int fd, Input& input, std::string& error);

bool WriteSuccess(int fd, const Output& output);
bool WriteFailure(int fd, const std::string& error);

// A true return means a complete success or logical-failure record was read.
// On logical failure, child_error is nonempty. A false return means the wire
// record itself was truncated, malformed, or dimensionally impossible.
bool ReadOutput(int fd, size_t expected_natom, Output& output,
                std::string& child_error, std::string& protocol_error);

}  // namespace nmr::mopac_worker
