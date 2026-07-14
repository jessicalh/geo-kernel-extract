#include "MopacWorkerProtocol.h"

#include <mopac.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <limits>
#include <string>
#include <vector>

#include <unistd.h>
#include <limits.h>
#include <dlfcn.h>
#include <sys/resource.h>

namespace {

constexpr double kMopacTolerance = 1.0;
constexpr int kMopacMaxTimeSeconds = 7200;
constexpr const char* kPinnedMopacLibrary =
    "/shared/2026Thesis/mopac2/local/lib/libmopac.so.2";

bool VerifyPinnedMopac(std::string& error) {
    Dl_info loaded{};
    if (::dladdr(reinterpret_cast<const void*>(&mozyme_scf), &loaded) == 0 ||
        !loaded.dli_fname) {
        error = "cannot identify the loaded libmopac object";
        return false;
    }

    std::array<char, PATH_MAX> loaded_real{};
    std::array<char, PATH_MAX> pinned_real{};
    if (!::realpath(loaded.dli_fname, loaded_real.data()) ||
        !::realpath(kPinnedMopacLibrary, pinned_real.data())) {
        error = "cannot resolve the loaded and pinned libmopac paths";
        return false;
    }
    if (std::string(loaded_real.data()) != std::string(pinned_real.data())) {
        error = "refusing non-pinned libmopac object: " +
                std::string(loaded_real.data());
        return false;
    }
    return true;
}

bool CheckedProduct(std::initializer_list<size_t> factors, size_t& value) {
    value = 1;
    for (size_t factor : factors) {
        if (factor != 0 &&
            value > std::numeric_limits<size_t>::max() / factor) {
            return false;
        }
        value *= factor;
    }
    return true;
}

template <typename Out, typename In>
bool CopyApiArray(const In* pointer, size_t count, const char* label,
                  std::vector<Out>& output, std::string& error) {
    if (count == 0) {
        output.clear();
        return true;
    }
    if (!pointer) {
        error = std::string("libmopac omitted requested ") + label;
        return false;
    }
    output.resize(count);
    for (size_t i = 0; i < count; ++i) {
        output[i] = static_cast<Out>(pointer[i]);
    }
    return true;
}

std::string MopacErrorText(const mopac_properties& properties) {
    if (properties.nerror <= 0 || !properties.error_msg) {
        return "libmopac API error count " +
               std::to_string(properties.nerror);
    }
    std::string result;
    for (int i = 0; i < properties.nerror; ++i) {
        if (!result.empty()) result += "; ";
        result += properties.error_msg[i] ? properties.error_msg[i]
                                          : "<null libmopac error>";
    }
    return result;
}

struct MopacAllocationGuard {
    mopac_properties* properties = nullptr;
    mozyme_state* state = nullptr;
    ~MopacAllocationGuard() {
        if (properties) destroy_mopac_properties(properties);
        if (state) destroy_mozyme_state(state);
    }
};

bool CallMozyme(const nmr::mopac_worker::Input& input,
                nmr::mopac_worker::Output& data,
                std::string& error) {
    std::vector<int> atoms(input.atoms.begin(), input.atoms.end());

    mopac_system system{};
    system.natom = static_cast<int>(atoms.size());
    system.natom_move = 0;
    system.charge = input.net_charge;
    system.spin = 0;
    system.model = 0;  // PM7
    system.epsilon = 1.0;
    system.atom = atoms.data();
    system.coord = const_cast<double*>(input.coordinates.data());
    system.nlattice = 0;
    system.nlattice_move = 0;
    system.pressure = 0.0;
    system.lattice = nullptr;
    system.tolerance = kMopacTolerance;
    system.max_time = kMopacMaxTimeSeconds;

    mozyme_state state{};       // numat=0: Lewis-structure initial guess
    mopac_properties properties{};
    mozyme_scf(&system, &state, &properties);
    MopacAllocationGuard guard{&properties, &state};

    if (properties.nerror != 0) {
        error = MopacErrorText(properties);
        return false;
    }

    const size_t natom = atoms.size();
    if (state.numat != static_cast<int>(natom)) {
        error = "libmopac returned MOZYME numat=" +
                std::to_string(state.numat) + " for " +
                std::to_string(natom) + " input atoms";
        return false;
    }
    if (properties.ao_max_orbitals <= 0 ||
        properties.ao_max_orbitals > 64) {
        error = "libmopac returned invalid ao_max_orbitals=" +
                std::to_string(properties.ao_max_orbitals);
        return false;
    }
    if (state.noccupied < 0 || state.nvirtual < 0 ||
        state.icocc_dim < 0 || state.icvir_dim < 0 ||
        state.cocc_dim < 0 || state.cvir_dim < 0) {
        error = "libmopac returned a negative MOZYME state dimension";
        return false;
    }

    data.natom = static_cast<std::int32_t>(natom);
    data.ao_max_orbitals = properties.ao_max_orbitals;
    data.state_dimensions = {
        state.numat, state.noccupied, state.nvirtual,
        state.icocc_dim, state.icvir_dim, state.cocc_dim, state.cvir_dim
    };
    data.heat = properties.heat;
    std::copy_n(properties.dipole, 3, data.dipole.begin());
    std::copy_n(properties.dipole_point_charge, 3,
                data.dipole_point_charge.begin());
    std::copy_n(properties.dipole_hybridization, 3,
                data.dipole_hybridization.begin());

    if (!CopyApiArray(properties.charge, natom, "atomic charges",
                      data.charge, error) ||
        !CopyApiArray(properties.bond_index, natom + 1,
                      "bond_index", data.bond_index, error)) {
        return false;
    }
    if (data.bond_index.empty() || data.bond_index.front() != 0) {
        error = "libmopac returned an invalid CSC bond_index origin";
        return false;
    }
    for (size_t i = 1; i < data.bond_index.size(); ++i) {
        if (data.bond_index[i] < data.bond_index[i - 1]) {
            error = "libmopac returned non-monotonic CSC bond_index";
            return false;
        }
    }
    const std::int32_t entry_count_i32 = data.bond_index.back();
    if (entry_count_i32 < 0 ||
        static_cast<std::uint64_t>(entry_count_i32) >
            static_cast<std::uint64_t>(natom) * natom) {
        error = "libmopac returned an invalid CSC entry count";
        return false;
    }
    const size_t entry_count = static_cast<size_t>(entry_count_i32);
    const size_t width = static_cast<size_t>(data.ao_max_orbitals);
    size_t atom_density_count = 0;
    size_t bond_density_count = 0;
    if (!CheckedProduct({natom, width, width}, atom_density_count) ||
        !CheckedProduct({entry_count, width, width}, bond_density_count)) {
        error = "libmopac AO-density dimensions overflow size_t";
        return false;
    }

    const size_t occupied = static_cast<size_t>(state.noccupied);
    const size_t virtual_count = static_cast<size_t>(state.nvirtual);
    if (occupied > std::numeric_limits<size_t>::max() - virtual_count) {
        error = "libmopac LMO count overflows size_t";
        return false;
    }
    const size_t lmo_count = occupied + virtual_count;

    return
        CopyApiArray(properties.bond_atom, entry_count, "bond_atom",
                     data.bond_atom, error) &&
        CopyApiArray(properties.bond_order, entry_count, "bond_order",
                     data.bond_order, error) &&
        CopyApiArray(properties.ao_orbitals, natom, "ao_orbitals",
                     data.ao_orbitals, error) &&
        CopyApiArray(properties.atom_ao_density, atom_density_count,
                     "atom_ao_density", data.atom_ao_density_fortran,
                     error) &&
        CopyApiArray(properties.bond_ao_density, bond_density_count,
                     "bond_ao_density", data.bond_ao_density_fortran,
                     error) &&
        CopyApiArray(state.nbonds, natom, "MOZYME nbonds",
                     data.nbonds, error) &&
        CopyApiArray(state.ibonds, 9 * natom, "MOZYME ibonds",
                     data.ibonds_fortran, error) &&
        CopyApiArray(state.iorbs, natom, "MOZYME iorbs",
                     data.iorbs, error) &&
        CopyApiArray(state.ncf, occupied, "MOZYME ncf",
                     data.ncf, error) &&
        CopyApiArray(state.nce, virtual_count, "MOZYME nce",
                     data.nce, error) &&
        CopyApiArray(state.icocc, static_cast<size_t>(state.icocc_dim),
                     "MOZYME icocc", data.icocc, error) &&
        CopyApiArray(state.icvir, static_cast<size_t>(state.icvir_dim),
                     "MOZYME icvir", data.icvir, error) &&
        CopyApiArray(state.cocc, static_cast<size_t>(state.cocc_dim),
                     "MOZYME cocc", data.cocc, error) &&
        CopyApiArray(state.cvir, static_cast<size_t>(state.cvir_dim),
                     "MOZYME cvir", data.cvir, error) &&
        CopyApiArray(properties.lmo_energy, lmo_count, "LMO energies",
                     data.lmo_energy, error) &&
        CopyApiArray(properties.lmo_occupied_atom_offset, occupied,
                     "occupied-LMO atom offsets",
                     data.occupied_atom_offsets, error) &&
        CopyApiArray(properties.lmo_virtual_atom_offset, virtual_count,
                     "virtual-LMO atom offsets",
                     data.virtual_atom_offsets, error) &&
        CopyApiArray(properties.lmo_occupied_coefficient_offset, occupied,
                     "occupied-LMO coefficient offsets",
                     data.occupied_coefficient_offsets, error) &&
        CopyApiArray(properties.lmo_virtual_coefficient_offset,
                     virtual_count, "virtual-LMO coefficient offsets",
                     data.virtual_coefficient_offsets, error);
}

}  // namespace

int main() {
    const rlimit no_core{0, 0};
    (void)::setrlimit(RLIMIT_CORE, &no_core);

    std::string error;
    if (!VerifyPinnedMopac(error)) {
        const bool wrote = nmr::mopac_worker::WriteFailure(
            STDOUT_FILENO, error);
        return wrote ? 2 : 3;
    }

    nmr::mopac_worker::Input input;
    if (!nmr::mopac_worker::ReadInput(STDIN_FILENO, input, error)) {
        const bool wrote = nmr::mopac_worker::WriteFailure(
            STDOUT_FILENO, error);
        return wrote ? 2 : 3;
    }

    nmr::mopac_worker::Output output;
    bool success = false;
    try {
        success = CallMozyme(input, output, error);
    } catch (const std::exception& ex) {
        error = std::string("MOPAC worker exception: ") + ex.what();
    } catch (...) {
        error = "MOPAC worker raised an unknown exception";
    }

    const bool wrote = success
        ? nmr::mopac_worker::WriteSuccess(STDOUT_FILENO, output)
        : nmr::mopac_worker::WriteFailure(STDOUT_FILENO, error);
    if (!wrote) return 3;
    return success ? 0 : 2;
}
