#include "MopacResult.h"

#include "NpyWriter.h"
#include "OperationLog.h"
#include "Protein.h"

#include <mopac.h>
#include <omp.h>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <map>
#include <numeric>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <dlfcn.h>
#include <limits.h>

namespace nmr {
namespace {

constexpr double kMopacTolerance = 1.0;
constexpr int kMopacMaxTimeSeconds = 7200;
constexpr const char* kPinnedMopacLibrary =
    "/shared/2026Thesis/mopac2/local/lib/libmopac.so.2";

struct MozymeOutput {
    std::int32_t natom = 0;
    std::int32_t ao_max_orbitals = 0;
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

double QuietNaN() {
    return std::numeric_limits<double>::quiet_NaN();
}

// Reproduce the numeric value recovered by the legacy parser after a fixed
// precision Fortran F field was written. Field width only adds whitespace for
// the finite values accepted here, so the decimal precision determines the
// parsed double. Keeping the conversion textual also preserves signed zero and
// the runtime library's round-to-nearest behaviour at decimal boundaries.
bool QuantizeLegacyPrinted(double value, int precision, double& quantized) {
    if (!std::isfinite(value) || precision < 0) return false;
    char text[128];
    const int count =
        std::snprintf(text, sizeof(text), "%.*f", precision, value);
    if (count <= 0 || static_cast<size_t>(count) >= sizeof(text)) return false;
    char* end = nullptr;
    errno = 0;
    quantized = std::strtod(text, &end);
    return errno == 0 && end != text && *end == '\0' &&
           std::isfinite(quantized);
}

// The legacy subprocess writer used %14.8f coordinates. cif++ exposes its
// coordinate storage through float-expanded doubles, so passing the raw
// doubles changes the frozen benchmark at about 1e-5. Reproducing the same
// decimal quantization in memory keeps the physics/input identical without
// writing a MOPAC input file.
bool QuantizeLegacyCoordinate(double value, double& quantized) {
    if (!std::isfinite(value)) return false;
    char text[64];
    const int count = std::snprintf(text, sizeof(text), "%.8f", value);
    if (count <= 0 || static_cast<size_t>(count) >= sizeof(text)) return false;
    char* end = nullptr;
    errno = 0;
    quantized = std::strtod(text, &end);
    return errno == 0 && end != text && *end == '\0' &&
           std::isfinite(quantized);
}

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

struct OpenMpThreadGuard {
    int previous_threads = 1;

    explicit OpenMpThreadGuard(int threads)
        : previous_threads(omp_get_max_threads()) {
        omp_set_num_threads(threads);
    }
    ~OpenMpThreadGuard() { omp_set_num_threads(previous_threads); }
};

}  // namespace

std::unique_ptr<MopacResult> MopacResult::Compute(
        ProteinConformation& conf,
        int net_charge,
        int threads,
        std::string* error_out) {
    if (error_out) error_out->clear();
    const Protein& protein = conf.ProteinRef();
    const size_t natoms = conf.AtomCount();

    auto fail = [&](const std::string& detail) -> std::unique_ptr<MopacResult> {
        const std::string message = "PM7/MOZYME/1SCF failed: " + detail;
        if (error_out) *error_out = message;
        OperationLog::Error("MopacResult::Compute", message);
        return nullptr;
    };

    if (natoms == 0) return fail("the conformation has no atoms");
    if (natoms > static_cast<size_t>(std::numeric_limits<std::int32_t>::max()) ||
        natoms > static_cast<size_t>(std::numeric_limits<std::uint32_t>::max())) {
        return fail("atom count exceeds the libmopac/index representation");
    }

    std::vector<std::int32_t> atoms(natoms);
    std::vector<double> coordinates(3 * natoms);
    for (size_t i = 0; i < natoms; ++i) {
        const int atomic_number =
            AtomicNumberForElement(protein.AtomAt(i).element);
        if (atomic_number <= 0) {
            return fail("unsupported typed element at atom " +
                        std::to_string(i));
        }
        atoms[i] = atomic_number;
        const Vec3& position = conf.PositionAt(i);
        for (size_t axis = 0; axis < 3; ++axis) {
            if (!QuantizeLegacyCoordinate(position[axis],
                                          coordinates[3 * i + axis])) {
                return fail("non-finite or unrepresentable coordinate at atom " +
                            std::to_string(i));
            }
        }
    }

    if (threads <= 0) {
        const int hardware =
            static_cast<int>(std::thread::hardware_concurrency());
        threads = std::max(4, (hardware * 3) / 4);
    }
    if (threads <= 0) return fail("invalid OpenMP thread count");
    OperationLog::Scope scope(
        "MopacResult::Compute",
        "PM7/MOZYME/1SCF libmopac.so.2 atoms=" +
        std::to_string(natoms) + " charge=" + std::to_string(net_charge) +
        " threads=" + std::to_string(threads));

    std::string api_error;
    if (!VerifyPinnedMopac(api_error)) return fail(api_error);

    MozymeOutput api;
    {
        std::vector<int> mopac_atoms(atoms.begin(), atoms.end());

        mopac_system system{};
        system.natom = static_cast<int>(mopac_atoms.size());
        system.natom_move = 0;
        system.charge = net_charge;
        system.spin = 0;
        system.model = 0;  // PM7
        system.epsilon = 1.0;
        system.atom = mopac_atoms.data();
        system.coord = coordinates.data();
        system.nlattice = 0;
        system.nlattice_move = 0;
        system.pressure = 0.0;
        system.lattice = nullptr;
        system.tolerance = kMopacTolerance;
        system.max_time = kMopacMaxTimeSeconds;

        mozyme_state state{};       // numat=0: Lewis-structure initial guess
        mopac_properties properties{};
        OpenMpThreadGuard thread_guard(threads);
        mozyme_scf(&system, &state, &properties);
        MopacAllocationGuard allocation_guard{&properties, &state};

        if (properties.nerror != 0) {
            return fail(MopacErrorText(properties));
        }

        const size_t natom = mopac_atoms.size();
        if (state.numat != static_cast<int>(natom)) {
            return fail("libmopac returned MOZYME numat=" +
                        std::to_string(state.numat) + " for " +
                        std::to_string(natom) + " input atoms");
        }
        if (properties.ao_max_orbitals <= 0 ||
            properties.ao_max_orbitals > 64) {
            return fail("libmopac returned invalid ao_max_orbitals=" +
                        std::to_string(properties.ao_max_orbitals));
        }
        if (state.noccupied < 0 || state.nvirtual < 0 ||
            state.icocc_dim < 0 || state.icvir_dim < 0 ||
            state.cocc_dim < 0 || state.cvir_dim < 0) {
            return fail("libmopac returned a negative MOZYME state dimension");
        }

        api.natom = static_cast<std::int32_t>(natom);
        api.ao_max_orbitals = properties.ao_max_orbitals;
        api.state_dimensions = {
            state.numat, state.noccupied, state.nvirtual,
            state.icocc_dim, state.icvir_dim, state.cocc_dim, state.cvir_dim
        };
        api.heat = properties.heat;
        std::copy_n(properties.dipole, 3, api.dipole.begin());
        std::copy_n(properties.dipole_point_charge, 3,
                    api.dipole_point_charge.begin());
        std::copy_n(properties.dipole_hybridization, 3,
                    api.dipole_hybridization.begin());

        if (!CopyApiArray(properties.charge, natom, "atomic charges",
                          api.charge, api_error) ||
            !CopyApiArray(properties.bond_index, natom + 1,
                          "bond_index", api.bond_index, api_error)) {
            return fail(api_error);
        }
        if (api.bond_index.empty() || api.bond_index.front() != 0) {
            return fail("libmopac returned an invalid CSC bond_index origin");
        }
        for (size_t i = 1; i < api.bond_index.size(); ++i) {
            if (api.bond_index[i] < api.bond_index[i - 1]) {
                return fail(
                    "libmopac returned non-monotonic CSC bond_index");
            }
        }
        const std::int32_t entry_count_i32 = api.bond_index.back();
        if (entry_count_i32 < 0 ||
            static_cast<std::uint64_t>(entry_count_i32) >
                static_cast<std::uint64_t>(natom) * natom) {
            return fail("libmopac returned an invalid CSC entry count");
        }
        const size_t entry_count = static_cast<size_t>(entry_count_i32);
        const size_t width = static_cast<size_t>(api.ao_max_orbitals);
        size_t atom_density_count = 0;
        size_t bond_density_count = 0;
        if (!CheckedProduct({natom, width, width}, atom_density_count) ||
            !CheckedProduct({entry_count, width, width},
                            bond_density_count)) {
            return fail("libmopac AO-density dimensions overflow size_t");
        }

        const size_t occupied = static_cast<size_t>(state.noccupied);
        const size_t virtual_count = static_cast<size_t>(state.nvirtual);
        if (occupied > std::numeric_limits<size_t>::max() - virtual_count) {
            return fail("libmopac LMO count overflows size_t");
        }
        const size_t lmo_count = occupied + virtual_count;
        size_t orbital_capacity = 0;
        if (!CheckedProduct({natom, width}, orbital_capacity) ||
            lmo_count > orbital_capacity) {
            return fail("libmopac returned impossible LMO dimensions");
        }

        if (!CopyApiArray(properties.bond_atom, entry_count, "bond_atom",
                          api.bond_atom, api_error) ||
            !CopyApiArray(properties.bond_order, entry_count, "bond_order",
                          api.bond_order, api_error) ||
            !CopyApiArray(properties.ao_orbitals, natom, "ao_orbitals",
                          api.ao_orbitals, api_error) ||
            !CopyApiArray(properties.atom_ao_density, atom_density_count,
                          "atom_ao_density", api.atom_ao_density_fortran,
                          api_error) ||
            !CopyApiArray(properties.bond_ao_density, bond_density_count,
                          "bond_ao_density", api.bond_ao_density_fortran,
                          api_error) ||
            !CopyApiArray(state.nbonds, natom, "MOZYME nbonds",
                          api.nbonds, api_error) ||
            !CopyApiArray(state.ibonds, 9 * natom, "MOZYME ibonds",
                          api.ibonds_fortran, api_error) ||
            !CopyApiArray(state.iorbs, natom, "MOZYME iorbs",
                          api.iorbs, api_error) ||
            !CopyApiArray(state.ncf, occupied, "MOZYME ncf",
                          api.ncf, api_error) ||
            !CopyApiArray(state.nce, virtual_count, "MOZYME nce",
                          api.nce, api_error) ||
            !CopyApiArray(state.icocc, static_cast<size_t>(state.icocc_dim),
                          "MOZYME icocc", api.icocc, api_error) ||
            !CopyApiArray(state.icvir, static_cast<size_t>(state.icvir_dim),
                          "MOZYME icvir", api.icvir, api_error) ||
            !CopyApiArray(state.cocc, static_cast<size_t>(state.cocc_dim),
                          "MOZYME cocc", api.cocc, api_error) ||
            !CopyApiArray(state.cvir, static_cast<size_t>(state.cvir_dim),
                          "MOZYME cvir", api.cvir, api_error) ||
            !CopyApiArray(properties.lmo_energy, lmo_count, "LMO energies",
                          api.lmo_energy, api_error) ||
            !CopyApiArray(properties.lmo_occupied_atom_offset, occupied,
                          "occupied-LMO atom offsets",
                          api.occupied_atom_offsets, api_error) ||
            !CopyApiArray(properties.lmo_virtual_atom_offset, virtual_count,
                          "virtual-LMO atom offsets",
                          api.virtual_atom_offsets, api_error) ||
            !CopyApiArray(properties.lmo_occupied_coefficient_offset,
                          occupied, "occupied-LMO coefficient offsets",
                          api.occupied_coefficient_offsets, api_error) ||
            !CopyApiArray(properties.lmo_virtual_coefficient_offset,
                          virtual_count,
                          "virtual-LMO coefficient offsets",
                          api.virtual_coefficient_offsets, api_error)) {
            return fail(api_error);
        }
    }

    auto result = std::make_unique<MopacResult>();
    result->charges_ = std::move(api.charge);
    result->bond_index_ = std::move(api.bond_index);
    result->bond_atom_ = std::move(api.bond_atom);
    result->bond_order_ = std::move(api.bond_order);
    result->ao_max_orbitals_ = api.ao_max_orbitals;
    // The benchmark names mozyme_state.iorbs as the live LMO coefficient
    // width. The appended properties.ao_orbitals field must be identical;
    // retain the state field and fail if the two API views ever diverge.
    result->ao_orbitals_ = std::move(api.iorbs);
    result->heat_of_formation_ = api.heat;
    result->dipole_ = Vec3(api.dipole[0], api.dipole[1], api.dipole[2]);
    result->dipole_point_charge_ = Vec3(
        api.dipole_point_charge[0], api.dipole_point_charge[1],
        api.dipole_point_charge[2]);
    result->dipole_hybridization_ = Vec3(
        api.dipole_hybridization[0], api.dipole_hybridization[1],
        api.dipole_hybridization[2]);

    if (!std::isfinite(result->heat_of_formation_) ||
        !result->dipole_.allFinite() ||
        !result->dipole_point_charge_.allFinite() ||
        !result->dipole_hybridization_.allFinite()) {
        return fail("libmopac returned a non-finite molecular property");
    }
    for (double charge : result->charges_) {
        if (!std::isfinite(charge)) {
            return fail("libmopac returned a non-finite atomic charge");
        }
    }

    const size_t width = static_cast<size_t>(result->ao_max_orbitals_);
    const size_t block_size = width * width;
    const size_t entries = result->bond_atom_.size();
    for (size_t i = 0; i < natoms; ++i) {
        if (result->ao_orbitals_[i] <= 0 ||
            result->ao_orbitals_[i] > result->ao_max_orbitals_ ||
            result->ao_orbitals_[i] != api.ao_orbitals[i]) {
            return fail("libmopac AO-width fields disagree at atom " +
                        std::to_string(i));
        }
    }

    // Convert the Fortran [W,W,N/E] buffers to emitted C [N/E,W,W].
    result->atom_ao_density_.resize(natoms * block_size);
    for (size_t atom = 0; atom < natoms; ++atom) {
        for (size_t row = 0; row < width; ++row) {
            for (size_t col = 0; col < width; ++col) {
                const double value = api.atom_ao_density_fortran[
                    row + width * (col + width * atom)];
                if (!std::isfinite(value)) {
                    return fail("non-finite atom AO density at atom " +
                                std::to_string(atom));
                }
                result->atom_ao_density_[
                    (atom * width + row) * width + col] = value;
            }
        }
    }
    result->bond_ao_density_directed_.resize(entries * block_size);
    for (size_t entry = 0; entry < entries; ++entry) {
        for (size_t row = 0; row < width; ++row) {
            for (size_t col = 0; col < width; ++col) {
                const double value = api.bond_ao_density_fortran[
                    row + width * (col + width * entry)];
                if (!std::isfinite(value)) {
                    return fail("non-finite bond AO density at CSC entry " +
                                std::to_string(entry));
                }
                result->bond_ao_density_directed_[
                    (entry * width + row) * width + col] = value;
            }
        }
    }

    result->atomic_orbital_populations_.resize(natoms * width);
    result->electron_populations_.assign(natoms, 0.0);
    result->s_pop_.assign(natoms, 0.0);
    result->p_pop_.assign(natoms, 0.0);
    result->d_pop_.assign(natoms, 0.0);
    for (size_t atom = 0; atom < natoms; ++atom) {
        for (size_t ao = 0; ao < width; ++ao) {
            const double population = result->atom_ao_density_[
                (atom * width + ao) * width + ao];
            result->atomic_orbital_populations_[atom * width + ao] =
                population;
            result->electron_populations_[atom] += population;
        }
        result->s_pop_[atom] =
            result->atomic_orbital_populations_[atom * width];
        for (size_t ao = 1; ao < std::min<size_t>(4, width); ++ao) {
            result->p_pop_[atom] +=
                result->atomic_orbital_populations_[atom * width + ao];
        }
        for (size_t ao = 4; ao < std::min<size_t>(9, width); ++ao) {
            result->d_pop_[atom] +=
                result->atomic_orbital_populations_[atom * width + ao];
        }
    }

    // Materialize the legacy printed-table projection during Compute. These
    // values are calculator state, not a WriteFeatures-time transformation:
    // F15.6 charge, F14.4 total population, 3F12.5 shell populations, and
    // 9F10.5 AO populations. Only each atom's live AO width was printed.
    result->compatibility_charges_.resize(natoms);
    result->compatibility_electron_populations_.resize(natoms);
    result->compatibility_s_pop_.resize(natoms);
    result->compatibility_p_pop_.resize(natoms);
    result->compatibility_d_pop_.resize(natoms);
    result->compatibility_atomic_orbital_populations_.assign(
        natoms * 9, QuietNaN());
    result->compatibility_atomic_orbital_population_totals_.assign(
        natoms * 3, QuietNaN());
    for (size_t atom = 0; atom < natoms; ++atom) {
        if (!QuantizeLegacyPrinted(result->charges_[atom], 6,
                                   result->compatibility_charges_[atom]) ||
            !QuantizeLegacyPrinted(
                result->electron_populations_[atom], 4,
                result->compatibility_electron_populations_[atom]) ||
            !QuantizeLegacyPrinted(result->s_pop_[atom], 5,
                                   result->compatibility_s_pop_[atom]) ||
            !QuantizeLegacyPrinted(result->p_pop_[atom], 5,
                                   result->compatibility_p_pop_[atom]) ||
            !QuantizeLegacyPrinted(result->d_pop_[atom], 5,
                                   result->compatibility_d_pop_[atom])) {
            return fail("cannot reproduce a legacy atomic print field");
        }
        const size_t live_width = std::min<size_t>(
            9, static_cast<size_t>(result->ao_orbitals_[atom]));
        for (size_t ao = 0; ao < live_width; ++ao) {
            if (!QuantizeLegacyPrinted(
                    result->atomic_orbital_populations_[atom * width + ao],
                    5,
                    result->compatibility_atomic_orbital_populations_[
                        atom * 9 + ao])) {
                return fail("cannot reproduce a legacy AO print field");
            }
        }
        result->compatibility_atomic_orbital_population_totals_[
            atom * 3] =
            result->compatibility_atomic_orbital_populations_[atom * 9];
        if (live_width >= 4) {
            result->compatibility_atomic_orbital_population_totals_[
                atom * 3 + 1] =
                result->compatibility_atomic_orbital_populations_[
                    atom * 9 + 1] +
                result->compatibility_atomic_orbital_populations_[
                    atom * 9 + 2] +
                result->compatibility_atomic_orbital_populations_[
                    atom * 9 + 3];
        }
        if (live_width >= 9) {
            double d_total = 0.0;
            for (size_t ao = 4; ao < 9; ++ao) {
                d_total +=
                    result->compatibility_atomic_orbital_populations_[
                        atom * 9 + ao];
            }
            result->compatibility_atomic_orbital_population_totals_[
                atom * 3 + 2] = d_total;
        }
    }

    if (!QuantizeLegacyPrinted(result->heat_of_formation_, 5,
                               result->compatibility_heat_of_formation_)) {
        return fail("cannot reproduce the legacy heat print field");
    }
    for (size_t axis = 0; axis < 3; ++axis) {
        if (!QuantizeLegacyPrinted(result->dipole_[axis], 3,
                                   result->compatibility_dipole_[axis])) {
            return fail("cannot reproduce a legacy dipole print field");
        }
    }

    struct PairAccum {
        double sum = 0.0;
        double max = -std::numeric_limits<double>::infinity();
        size_t count = 0;
        bool has_low_entry = false;
        double low_order = 0.0;
        std::int32_t low_entry = -1;
    };
    std::map<std::uint64_t, PairAccum> pair_accum;
    std::vector<std::uint64_t> low_pair_order;
    result->valencies_.assign(natoms, QuietNaN());
    result->compatibility_valencies_.assign(natoms, QuietNaN());
    result->compatibility_project_valencies_.assign(natoms, 0.0);

    for (size_t atom_a = 0; atom_a < natoms; ++atom_a) {
        const size_t start = static_cast<size_t>(result->bond_index_[atom_a]);
        const size_t stop = static_cast<size_t>(result->bond_index_[atom_a + 1]);
        if (stop > entries) return fail("CSC offset exceeds bond arrays");
        size_t diagonal_count = 0;
        std::vector<MopacBondOrder> compact_candidates;
        compact_candidates.reserve(stop - start);
        for (size_t entry = start; entry < stop; ++entry) {
            const std::int32_t atom_b_i32 = result->bond_atom_[entry];
            const double order = result->bond_order_[entry];
            if (atom_b_i32 < 0 || static_cast<size_t>(atom_b_i32) >= natoms ||
                !std::isfinite(order)) {
                return fail("invalid CSC bond entry " + std::to_string(entry));
            }
            const size_t atom_b = static_cast<size_t>(atom_b_i32);
            if (atom_a == atom_b) {
                result->valencies_[atom_a] = order;
                ++diagonal_count;
                continue;
            }

            result->bond_orders_.push_back({atom_a, atom_b, order});
            compact_candidates.push_back({atom_a, atom_b, order});
            const std::uint64_t key = PairKey(atom_a, atom_b);
            auto& accum = pair_accum[key];
            accum.sum += order;
            accum.max = std::max(accum.max, order);
            ++accum.count;
            if (atom_a < atom_b) {
                if (accum.has_low_entry) {
                    return fail("duplicate a<b CSC bond entry");
                }
                accum.has_low_entry = true;
                accum.low_order = order;
                accum.low_entry = static_cast<std::int32_t>(entry);
                low_pair_order.push_back(key);
            }
        }
        if (diagonal_count != 1 ||
            !std::isfinite(result->valencies_[atom_a])) {
            return fail("libmopac omitted the CSC valency diagonal for atom " +
                        std::to_string(atom_a));
        }
        if (!QuantizeLegacyPrinted(
                result->valencies_[atom_a], 3,
                result->compatibility_valencies_[atom_a])) {
            return fail("cannot reproduce a legacy valency print field");
        }

        // bonds_for_MOZYME visits candidate atoms in ascending index order,
        // then print_bonds_compact applies this selection sort. Its 1e-4
        // inflation deliberately makes near ties stable. The old parser read
        // only the physical row-start line (six entries), parsed F6.3, and
        // retained values strictly greater than 0.01.
        std::sort(compact_candidates.begin(), compact_candidates.end(),
                  [](const MopacBondOrder& lhs,
                     const MopacBondOrder& rhs) {
                      return lhs.atom_b < rhs.atom_b;
                  });
        for (size_t j = 0; j < compact_candidates.size(); ++j) {
            size_t selected = j;
            double running = 0.0;
            bool found = false;
            for (size_t k = j; k < compact_candidates.size(); ++k) {
                if (compact_candidates[k].wiberg_order > running) {
                    selected = k;
                    running = compact_candidates[k].wiberg_order *
                              (1.0 + 1.0e-4);
                    found = true;
                }
            }
            if (!found) {
                return fail("invalid non-positive compact bond candidate");
            }
            std::swap(compact_candidates[j], compact_candidates[selected]);
        }
        const size_t printed_count =
            std::min<size_t>(6, compact_candidates.size());
        for (size_t j = 0; j < printed_count; ++j) {
            double printed_order = 0.0;
            if (!QuantizeLegacyPrinted(
                    compact_candidates[j].wiberg_order, 3,
                    printed_order)) {
                return fail("cannot reproduce a legacy bond-order field");
            }
            if (!(printed_order > 0.01)) continue;
            result->compatibility_project_valencies_[atom_a] +=
                printed_order;
            const std::uint64_t key = PairKey(
                atom_a, compact_candidates[j].atom_b);
            const auto inserted =
                result->bond_order_map_.emplace(key, printed_order);
            if (!inserted.second && printed_order > inserted.first->second) {
                inserted.first->second = printed_order;
            }
        }
    }

    std::unordered_map<std::uint64_t, size_t> topology_by_pair;
    topology_by_pair.reserve(protein.BondCount());
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const Bond& bond = protein.BondAt(bi);
        topology_by_pair[PairKey(bond.atom_index_a, bond.atom_index_b)] = bi;
    }

    result->unique_bonds_.reserve(low_pair_order.size());
    result->unique_bond_ao_density_.reserve(
        low_pair_order.size() * block_size);
    if (low_pair_order.size() != pair_accum.size()) {
        return fail("libmopac returned an asymmetric Wiberg pair set");
    }
    for (std::uint64_t key : low_pair_order) {
        const auto it = pair_accum.find(key);
        if (it == pair_accum.end() || !it->second.has_low_entry ||
            it->second.count != 2) {
            return fail("incomplete symmetric Wiberg projection");
        }
        const double other_order = it->second.sum - it->second.low_order;
        const double symmetry_scale = std::max(
            {1.0, std::abs(it->second.low_order), std::abs(other_order)});
        if (std::abs(other_order - it->second.low_order) >
            1.0e-12 * symmetry_scale) {
            return fail("directed Wiberg values disagree for one atom pair");
        }
        const size_t atom_a = static_cast<size_t>(key >> 32);
        const size_t atom_b = static_cast<size_t>(key & 0xffffffffULL);
        UniqueBond unique;
        unique.atom_a = static_cast<std::int32_t>(atom_a);
        unique.atom_b = static_cast<std::int32_t>(atom_b);
        unique.order = it->second.low_order;
        unique.max_order = it->second.max;
        unique.mean_order = it->second.sum /
                            static_cast<double>(it->second.count);
        double printed_low_order = 0.0;
        double printed_other_order = 0.0;
        if (!QuantizeLegacyPrinted(it->second.low_order, 3,
                                   printed_low_order) ||
            !QuantizeLegacyPrinted(other_order, 3,
                                   printed_other_order)) {
            return fail("cannot reproduce a legacy unique-bond field");
        }
        unique.compatibility_max_order =
            std::max(printed_low_order, printed_other_order);
        unique.compatibility_mean_order =
            (printed_low_order + printed_other_order) / 2.0;
        unique.directed_entry = it->second.low_entry;
        const auto topology_it = topology_by_pair.find(key);
        if (topology_it != topology_by_pair.end()) {
            unique.topology_bond_index =
                static_cast<std::int32_t>(topology_it->second);
        }

        const size_t entry = static_cast<size_t>(unique.directed_entry);
        const auto block_begin =
            result->bond_ao_density_directed_.begin() + entry * block_size;
        result->unique_bond_ao_density_.insert(
            result->unique_bond_ao_density_.end(), block_begin,
            block_begin + block_size);
        double reconstructed = 0.0;
        for (size_t k = 0; k < block_size; ++k) {
            const double value = *(block_begin + k);
            reconstructed += value * value;
        }
        if (std::abs(reconstructed - unique.order) > 1.0e-8) {
            return fail("AO-density/Wiberg mismatch for atom pair " +
                        std::to_string(atom_a) + "," +
                        std::to_string(atom_b));
        }
        result->unique_bonds_.push_back(unique);
    }

    result->topology_bond_orders_.assign(protein.BondCount(), 0.0);
    result->topology_unique_pair_index_.assign(protein.BondCount(), -1);
    for (size_t i = 0; i < result->unique_bonds_.size(); ++i) {
        const UniqueBond& unique = result->unique_bonds_[i];
        if (unique.topology_bond_index >= 0) {
            const size_t bi =
                static_cast<size_t>(unique.topology_bond_index);
            result->topology_unique_pair_index_[bi] =
                static_cast<std::int32_t>(i);
        }
    }
    // Query/trajectory consumers historically saw the compact-table parser,
    // not the complete ALLBONDS table. Populate them from that source-exact
    // compatibility map while retaining the direct unique indices for the
    // reduced mopac_bond_orders_unique interpretation.
    for (const auto& item : result->bond_order_map_) {
        const auto topology_it = topology_by_pair.find(item.first);
        if (topology_it != topology_by_pair.end()) {
            result->topology_bond_orders_[topology_it->second] = item.second;
        }
    }

    // Lewis structure: expose zero-based live entries and explicit -1 in
    // unused slots. This is lossless for the API's 1-based/zero-unused array.
    result->lewis_bond_count_ = std::move(api.nbonds);
    result->lewis_bond_atoms_.assign(9 * natoms, -1);
    for (size_t atom = 0; atom < natoms; ++atom) {
        const std::int32_t count = result->lewis_bond_count_[atom];
        if (count < 0 || count > 9) {
            return fail("MOZYME Lewis-bond count is outside 0..9");
        }
        for (std::int32_t slot = 0; slot < count; ++slot) {
            const std::int32_t native =
                api.ibonds_fortran[static_cast<size_t>(slot) + 9 * atom];
            if (native <= 0 || static_cast<size_t>(native) > natoms) {
                return fail("MOZYME Lewis-bond atom is out of range");
            }
            result->lewis_bond_atoms_[9 * atom +
                                      static_cast<size_t>(slot)] = native - 1;
        }
    }

    result->mozyme_state_dimensions_ = api.state_dimensions;
    result->lmo_occupied_atom_counts_ = std::move(api.ncf);
    result->lmo_virtual_atom_counts_ = std::move(api.nce);
    result->lmo_energy_levels_ = std::move(api.lmo_energy);
    result->lmo_occupied_atom_offsets_native_ =
        std::move(api.occupied_atom_offsets);
    result->lmo_virtual_atom_offsets_native_ =
        std::move(api.virtual_atom_offsets);
    result->lmo_occupied_coefficient_offsets_native_ =
        std::move(api.occupied_coefficient_offsets);
    result->lmo_virtual_coefficient_offsets_native_ =
        std::move(api.virtual_coefficient_offsets);
    result->lmo_occupied_atom_storage_native_ = std::move(api.icocc);
    result->lmo_virtual_atom_storage_native_ = std::move(api.icvir);
    result->lmo_occupied_coefficient_storage_native_ = std::move(api.cocc);
    result->lmo_virtual_coefficient_storage_native_ = std::move(api.cvir);

    const auto all_finite = [](const std::vector<double>& values) {
        return std::all_of(values.begin(), values.end(),
                           [](double value) { return std::isfinite(value); });
    };
    if (!all_finite(result->lmo_energy_levels_) ||
        !all_finite(result->lmo_occupied_coefficient_storage_native_) ||
        !all_finite(result->lmo_virtual_coefficient_storage_native_)) {
        return fail("libmopac returned non-finite LMO data");
    }

    auto pack_lmos = [&](const char* label,
                         const std::vector<std::int32_t>& counts,
                         const std::vector<std::int32_t>& atom_offsets,
                         const std::vector<std::int32_t>& coefficient_offsets,
                         const std::vector<std::int32_t>& native_atoms,
                         const std::vector<double>& native_coefficients,
                         std::vector<std::int32_t>& packed_atoms,
                         std::vector<double>& packed_coefficients) -> bool {
        for (size_t lmo = 0; lmo < counts.size(); ++lmo) {
            const std::int32_t count_i32 = counts[lmo];
            const std::int32_t atom_offset_i32 = atom_offsets[lmo];
            const std::int32_t coefficient_offset_i32 =
                coefficient_offsets[lmo];
            if (count_i32 < 0 || atom_offset_i32 < 0 ||
                coefficient_offset_i32 < 0) {
                api_error = std::string(label) +
                    " LMO has a negative count/offset";
                return false;
            }
            const size_t count = static_cast<size_t>(count_i32);
            const size_t atom_offset = static_cast<size_t>(atom_offset_i32);
            const size_t coefficient_offset =
                static_cast<size_t>(coefficient_offset_i32);
            if (atom_offset > native_atoms.size() ||
                count > native_atoms.size() - atom_offset) {
                api_error = std::string(label) +
                    " LMO atom offset is outside native storage";
                return false;
            }
            size_t coefficient_count = 0;
            for (size_t j = 0; j < count; ++j) {
                const std::int32_t native_atom =
                    native_atoms[atom_offset + j];
                if (native_atom <= 0 ||
                    static_cast<size_t>(native_atom) > natoms) {
                    api_error = std::string(label) +
                        " LMO atom index is outside the molecule";
                    return false;
                }
                const size_t atom = static_cast<size_t>(native_atom - 1);
                const size_t ao_count =
                    static_cast<size_t>(result->ao_orbitals_[atom]);
                if (coefficient_count >
                    std::numeric_limits<size_t>::max() - ao_count) {
                    api_error = std::string(label) +
                        " LMO coefficient count overflow";
                    return false;
                }
                coefficient_count += ao_count;
                packed_atoms.push_back(native_atom - 1);
            }
            if (coefficient_offset > native_coefficients.size() ||
                coefficient_count >
                    native_coefficients.size() - coefficient_offset) {
                api_error = std::string(label) +
                    " LMO coefficient offset is outside native storage";
                return false;
            }
            packed_coefficients.insert(
                packed_coefficients.end(),
                native_coefficients.begin() + coefficient_offset,
                native_coefficients.begin() + coefficient_offset +
                    coefficient_count);
        }
        return true;
    };

    if (!pack_lmos("occupied", result->lmo_occupied_atom_counts_,
                   result->lmo_occupied_atom_offsets_native_,
                   result->lmo_occupied_coefficient_offsets_native_,
                   result->lmo_occupied_atom_storage_native_,
                   result->lmo_occupied_coefficient_storage_native_,
                   result->lmo_occupied_atoms_,
                   result->lmo_occupied_coefficients_) ||
        !pack_lmos("virtual", result->lmo_virtual_atom_counts_,
                   result->lmo_virtual_atom_offsets_native_,
                   result->lmo_virtual_coefficient_offsets_native_,
                   result->lmo_virtual_atom_storage_native_,
                   result->lmo_virtual_coefficient_storage_native_,
                   result->lmo_virtual_atoms_,
                   result->lmo_virtual_coefficients_)) {
        return fail(api_error);
    }

    std::vector<std::vector<MopacBondNeighbour>> neighbours(natoms);
    for (const auto& item : result->bond_order_map_) {
        const size_t atom_a = static_cast<size_t>(item.first >> 32);
        const size_t atom_b = static_cast<size_t>(
            item.first & 0xffffffffULL);
        const auto topology_it = topology_by_pair.find(item.first);
        const size_t topology_index = topology_it == topology_by_pair.end()
            ? SIZE_MAX
            : topology_it->second;
        neighbours[atom_a].push_back(
            {atom_b, item.second, topology_index});
        neighbours[atom_b].push_back(
            {atom_a, item.second, topology_index});
    }
    for (auto& atom_neighbours : neighbours) {
        std::sort(atom_neighbours.begin(), atom_neighbours.end(),
            [](const MopacBondNeighbour& lhs,
               const MopacBondNeighbour& rhs) {
                return lhs.wiberg_order > rhs.wiberg_order;
            });
    }

    // Commit calculator-owned ConformationAtom fields only after the entire
    // API payload and the per-atom neighbour vectors are ready.
    for (size_t atom = 0; atom < natoms; ++atom) {
        auto& conformation_atom = conf.MutableAtomAt(atom);
        conformation_atom.mopac_charge =
            result->compatibility_charges_[atom];
        conformation_atom.mopac_s_pop =
            result->compatibility_s_pop_[atom];
        conformation_atom.mopac_p_pop =
            result->compatibility_p_pop_[atom];
        // Preserve the consumed compatibility surface: legacy MopacResult
        // summed the quantized bonds visible on each compact-table row. The
        // exact API CSC diagonal and sparse entries remain stored separately.
        conformation_atom.mopac_valency =
            result->compatibility_project_valencies_[atom];
        conformation_atom.mopac_bond_neighbours =
            std::move(neighbours[atom]);
    }

    OperationLog::Log(
        OperationLog::Level::Info, LogMopac, "MopacResult::Compute",
        "heat=" + std::to_string(result->heat_of_formation_) +
        " kcal/mol, sparse_entries=" + std::to_string(entries) +
        ", direct_unique_bonds=" +
        std::to_string(result->unique_bonds_.size()) +
        ", compatibility_bonds=" +
        std::to_string(result->bond_order_map_.size()) +
        ", ao_width=" + std::to_string(result->ao_max_orbitals_));
    return result;
}

double MopacResult::ChargeAt(size_t atom_index) const {
    return atom_index < compatibility_charges_.size()
        ? compatibility_charges_[atom_index]
        : 0.0;
}

double MopacResult::SPopAt(size_t atom_index) const {
    return atom_index < compatibility_s_pop_.size()
        ? compatibility_s_pop_[atom_index]
        : 0.0;
}

double MopacResult::PPopAt(size_t atom_index) const {
    return atom_index < compatibility_p_pop_.size()
        ? compatibility_p_pop_[atom_index]
        : 0.0;
}

double MopacResult::ValencyAt(size_t atom_index) const {
    return atom_index < compatibility_project_valencies_.size()
        ? compatibility_project_valencies_[atom_index]
        : 0.0;
}

double MopacResult::BondOrder(size_t atom_a, size_t atom_b) const {
    const auto it = bond_order_map_.find(PairKey(atom_a, atom_b));
    return it == bond_order_map_.end() ? 0.0 : it->second;
}

double MopacResult::TopologyBondOrder(size_t bond_index) const {
    return bond_index < topology_bond_orders_.size()
        ? topology_bond_orders_[bond_index]
        : 0.0;
}

int MopacResult::WriteFeatures(const ProteinConformation& conf,
                               const std::string& output_dir) const {
    const size_t atom_count = conf.AtomCount();
    const size_t unique_count = unique_bonds_.size();
    const size_t compatibility_bond_count = bond_order_map_.size();
    const size_t entry_count = bond_atom_.size();
    const size_t width = static_cast<size_t>(ao_max_orbitals_);
    int written = 0;
    auto record = [&](bool success) {
        if (success) ++written;
    };

    // Existing compatibility surface.
    {
        std::vector<double> data(atom_count);
        for (size_t i = 0; i < atom_count; ++i) {
            data[i] = conf.AtomAt(i).mopac_charge;
        }
        record(NpyWriter::WriteFloat64(output_dir + "/mopac_charges.npy",
                                       data.data(), atom_count));
    }
    {
        std::vector<double> data(atom_count * 4);
        for (size_t i = 0; i < atom_count; ++i) {
            data[4 * i + 0] = conf.AtomAt(i).mopac_charge;
            data[4 * i + 1] = conf.AtomAt(i).mopac_s_pop;
            data[4 * i + 2] = conf.AtomAt(i).mopac_p_pop;
            data[4 * i + 3] = conf.AtomAt(i).mopac_valency;
        }
        record(NpyWriter::WriteFloat64(output_dir + "/mopac_scalars.npy",
                                       data.data(), atom_count, 4));
    }
    {
        std::vector<double> data(compatibility_bond_count * 3);
        size_t row = 0;
        // Preserve the historical unordered-map iteration as well as its
        // compact-table membership and values; trained consumers saw this
        // exact row ordering rather than atom-sorted direct API pairs.
        for (const auto& item : bond_order_map_) {
            data[3 * row + 0] =
                static_cast<double>(item.first >> 32);
            data[3 * row + 1] = static_cast<double>(
                item.first & 0xffffffffULL);
            data[3 * row + 2] = item.second;
            ++row;
        }
        record(NpyWriter::WriteFloat64(
            output_dir + "/mopac_bond_orders.npy",
            data.empty() ? nullptr : data.data(),
            compatibility_bond_count, 3));
    }
    {
        size_t rows = 0;
        for (size_t i = 0; i < atom_count; ++i) {
            rows += conf.AtomAt(i).mopac_bond_neighbours.size();
        }
        std::vector<double> data(rows * 4);
        size_t row = 0;
        for (size_t atom = 0; atom < atom_count; ++atom) {
            for (const auto& neighbour :
                 conf.AtomAt(atom).mopac_bond_neighbours) {
                data[4 * row + 0] = static_cast<double>(atom);
                data[4 * row + 1] =
                    static_cast<double>(neighbour.other_atom);
                data[4 * row + 2] = neighbour.wiberg_order;
                data[4 * row + 3] = neighbour.topology_bond_index == SIZE_MAX
                    ? -1.0
                    : static_cast<double>(neighbour.topology_bond_index);
                ++row;
            }
        }
        record(NpyWriter::WriteFloat64(
            output_dir + "/mopac_bond_neighbors.npy",
            data.empty() ? nullptr : data.data(), rows, 4));
    }
    {
        const double data[4] = {
            compatibility_heat_of_formation_,
            compatibility_dipole_.x(), compatibility_dipole_.y(),
            compatibility_dipole_.z()
        };
        record(NpyWriter::WriteFloat64(output_dir + "/mopac_global.npy",
                                       data, 4));
    }
    {
        std::vector<double> data(atom_count * 12, QuietNaN());
        for (size_t i = 0; i < atom_count; ++i) {
            data[12 * i + 0] = compatibility_charges_[i];
            data[12 * i + 1] = compatibility_electron_populations_[i];
            data[12 * i + 2] = compatibility_s_pop_[i];
            const size_t live_width =
                static_cast<size_t>(ao_orbitals_[i]);
            if (live_width >= 4) {
                data[12 * i + 3] = compatibility_p_pop_[i];
            }
            if (live_width >= 9) {
                data[12 * i + 4] = compatibility_d_pop_[i];
            }
            // 5 (printed f population) and 6:10 (printed per-atom dipole)
            // have no API source and remain explicit NaN/N/A.
            data[12 * i + 10] = compatibility_valencies_[i];
            data[12 * i + 11] =
                compatibility_project_valencies_[i];
        }
        record(NpyWriter::WriteFloat64(
            output_dir + "/mopac_atom_populations.npy",
            data.data(), atom_count, 12));
    }
    {
        record(NpyWriter::WriteFloat64(
            output_dir + "/mopac_atomic_orbital_populations.npy",
            compatibility_atomic_orbital_populations_.data(),
            atom_count, 9));
    }
    {
        record(NpyWriter::WriteFloat64(
            output_dir + "/mopac_atomic_orbital_population_totals.npy",
            compatibility_atomic_orbital_population_totals_.data(),
            atom_count, 3));
    }
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_bond_valencies.npy",
        compatibility_valencies_.data(), atom_count));
    {
        std::vector<double> data(unique_count * 8, QuietNaN());
        for (size_t i = 0; i < unique_count; ++i) {
            data[8 * i + 0] = unique_bonds_[i].atom_a;
            data[8 * i + 1] = unique_bonds_[i].atom_b;
            data[8 * i + 2] =
                unique_bonds_[i].compatibility_max_order;
            data[8 * i + 3] =
                unique_bonds_[i].compatibility_mean_order;
            // 4:7 were printed-entry count/indices and are retired as NaN.
            if (unique_bonds_[i].topology_bond_index >= 0) {
                data[8 * i + 7] =
                    unique_bonds_[i].topology_bond_index;
            }
        }
        record(NpyWriter::WriteFloat64(
            output_dir + "/mopac_bond_orders_unique.npy",
            data.empty() ? nullptr : data.data(), unique_count, 8));
    }
    {
        const Protein& protein = conf.ProteinRef();
        const size_t topology_count = protein.BondCount();
        std::vector<double> data(topology_count * 8, QuietNaN());
        for (size_t i = 0; i < topology_count; ++i) {
            const Bond& bond = protein.BondAt(i);
            const bool present = topology_unique_pair_index_[i] >= 0;
            data[8 * i + 0] = static_cast<double>(i);
            data[8 * i + 1] = static_cast<double>(bond.atom_index_a);
            data[8 * i + 2] = static_cast<double>(bond.atom_index_b);
            if (present) {
                const size_t unique_index = static_cast<size_t>(
                    topology_unique_pair_index_[i]);
                data[8 * i + 3] =
                    unique_bonds_[unique_index].compatibility_max_order;
            }
            data[8 * i + 4] = present ? 1.0 : 0.0;
            if (present) {
                data[8 * i + 5] = topology_unique_pair_index_[i];
            }
            data[8 * i + 6] = present ? 0.0 : 1.0;
            // 7 was the text printed-entry count and is retired as NaN.
        }
        record(NpyWriter::WriteFloat64(
            output_dir + "/mopac_topology_bond_orders_full.npy",
            data.empty() ? nullptr : data.data(), topology_count, 8));
    }

    // Direct libmopac quantities and the index arrays needed to interpret
    // their padded/sparse/native storage.
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_charges_full_precision.npy",
        charges_.data(), atom_count));
    {
        std::vector<double> data(unique_count * 3);
        for (size_t i = 0; i < unique_count; ++i) {
            data[3 * i + 0] = unique_bonds_[i].atom_a;
            data[3 * i + 1] = unique_bonds_[i].atom_b;
            // Exact a<b CSC entry, matching the direct probe/archive view.
            data[3 * i + 2] = unique_bonds_[i].order;
        }
        record(NpyWriter::WriteFloat64(
            output_dir + "/mopac_bond_orders_full_precision.npy",
            data.empty() ? nullptr : data.data(), unique_count, 3));
    }
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_bond_valencies_full_precision.npy",
        valencies_.data(), atom_count));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_heat_kcal_mol.npy",
        &heat_of_formation_, 1));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_dipole_debye.npy", dipole_.data(), 3));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_dipole_point_charge_debye.npy",
        dipole_point_charge_.data(), 3));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_dipole_hybridization_debye.npy",
        dipole_hybridization_.data(), 3));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_bond_index.npy", bond_index_.data(),
        bond_index_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_bond_atom.npy",
        bond_atom_.empty() ? nullptr : bond_atom_.data(), entry_count));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_bond_order.npy",
        bond_order_.empty() ? nullptr : bond_order_.data(), entry_count));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_ao_max_orbitals.npy", &ao_max_orbitals_, 1));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_ao_orbitals_per_atom.npy",
        ao_orbitals_.data(), atom_count));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_atom_ao_density.npy",
        atom_ao_density_.data(), {atom_count, width, width}));
    record(NpyWriter::WriteFloat64(
        output_dir +
            "/mopac_atomic_orbital_populations_full_precision.npy",
        atomic_orbital_populations_.data(), atom_count, width));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_bond_ao_density_directed.npy",
        bond_ao_density_directed_.empty()
            ? nullptr : bond_ao_density_directed_.data(),
        {entry_count, width, width}));
    {
        std::vector<std::int32_t> pairs(unique_count * 2);
        for (size_t i = 0; i < unique_count; ++i) {
            pairs[2 * i + 0] = unique_bonds_[i].atom_a;
            pairs[2 * i + 1] = unique_bonds_[i].atom_b;
        }
        record(NpyWriter::WriteInt32(
            output_dir + "/mopac_bond_density_pairs.npy",
            pairs.empty() ? nullptr : pairs.data(), unique_count, 2));
    }
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_bond_ao_density.npy",
        unique_bond_ao_density_.empty()
            ? nullptr : unique_bond_ao_density_.data(),
        {unique_count, width, width}));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_atom_electron_population.npy",
        electron_populations_.data(), atom_count));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_atom_s_population.npy",
        s_pop_.data(), atom_count));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_atom_p_population.npy",
        p_pop_.data(), atom_count));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_atom_d_population.npy",
        d_pop_.data(), atom_count));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lewis_bond_count.npy",
        lewis_bond_count_.data(), atom_count));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lewis_bond_atoms.npy",
        lewis_bond_atoms_.data(), atom_count, 9));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_lmo_energy_levels.npy",
        lmo_energy_levels_.empty() ? nullptr : lmo_energy_levels_.data(),
        lmo_energy_levels_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_occupied_atom_counts.npy",
        lmo_occupied_atom_counts_.empty()
            ? nullptr : lmo_occupied_atom_counts_.data(),
        lmo_occupied_atom_counts_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_occupied_atoms.npy",
        lmo_occupied_atoms_.empty() ? nullptr : lmo_occupied_atoms_.data(),
        lmo_occupied_atoms_.size()));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_lmo_occupied_coefficients.npy",
        lmo_occupied_coefficients_.empty()
            ? nullptr : lmo_occupied_coefficients_.data(),
        lmo_occupied_coefficients_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_virtual_atom_counts.npy",
        lmo_virtual_atom_counts_.empty()
            ? nullptr : lmo_virtual_atom_counts_.data(),
        lmo_virtual_atom_counts_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_virtual_atoms.npy",
        lmo_virtual_atoms_.empty() ? nullptr : lmo_virtual_atoms_.data(),
        lmo_virtual_atoms_.size()));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_lmo_virtual_coefficients.npy",
        lmo_virtual_coefficients_.empty()
            ? nullptr : lmo_virtual_coefficients_.data(),
        lmo_virtual_coefficients_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_occupied_atom_offsets_native.npy",
        lmo_occupied_atom_offsets_native_.empty()
            ? nullptr : lmo_occupied_atom_offsets_native_.data(),
        lmo_occupied_atom_offsets_native_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_virtual_atom_offsets_native.npy",
        lmo_virtual_atom_offsets_native_.empty()
            ? nullptr : lmo_virtual_atom_offsets_native_.data(),
        lmo_virtual_atom_offsets_native_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_occupied_coefficient_offsets_native.npy",
        lmo_occupied_coefficient_offsets_native_.empty()
            ? nullptr : lmo_occupied_coefficient_offsets_native_.data(),
        lmo_occupied_coefficient_offsets_native_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_virtual_coefficient_offsets_native.npy",
        lmo_virtual_coefficient_offsets_native_.empty()
            ? nullptr : lmo_virtual_coefficient_offsets_native_.data(),
        lmo_virtual_coefficient_offsets_native_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_occupied_atom_storage_native.npy",
        lmo_occupied_atom_storage_native_.empty()
            ? nullptr : lmo_occupied_atom_storage_native_.data(),
        lmo_occupied_atom_storage_native_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_lmo_virtual_atom_storage_native.npy",
        lmo_virtual_atom_storage_native_.empty()
            ? nullptr : lmo_virtual_atom_storage_native_.data(),
        lmo_virtual_atom_storage_native_.size()));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_lmo_occupied_coefficient_storage_native.npy",
        lmo_occupied_coefficient_storage_native_.empty()
            ? nullptr : lmo_occupied_coefficient_storage_native_.data(),
        lmo_occupied_coefficient_storage_native_.size()));
    record(NpyWriter::WriteFloat64(
        output_dir + "/mopac_lmo_virtual_coefficient_storage_native.npy",
        lmo_virtual_coefficient_storage_native_.empty()
            ? nullptr : lmo_virtual_coefficient_storage_native_.data(),
        lmo_virtual_coefficient_storage_native_.size()));
    record(NpyWriter::WriteInt32(
        output_dir + "/mopac_mozyme_state_dimensions.npy",
        mozyme_state_dimensions_.data(), mozyme_state_dimensions_.size()));

    return written;
}

}  // namespace nmr
