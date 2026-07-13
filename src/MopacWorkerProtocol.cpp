#include "MopacWorkerProtocol.h"

#include <cerrno>
#include <cstring>
#include <limits>
#include <type_traits>

#include <unistd.h>
#include <sys/socket.h>

namespace nmr::mopac_worker {
namespace {

constexpr std::uint64_t kInputMagic = 0x4d4f50414332494eULL;   // MOPAC2IN
constexpr std::uint64_t kOutputMagic = 0x4d4f504143324f55ULL;  // MOPAC2OU
constexpr std::uint64_t kWireTrailer = 0x454e444d4f504143ULL;  // ENDMOPAC
constexpr std::uint32_t kWireVersion = 1;
constexpr std::int32_t kWireSuccess = 0;
constexpr std::int32_t kWireFailure = 1;
constexpr std::uint64_t kMaxWireStringBytes = 1024 * 1024;
// This is a corruption guard with headroom over the source-derived allocation
// bound for the validated 11,405-atom system, not a curation threshold.
constexpr std::uint64_t kMaxWirePayloadBytes = 2ULL * 1024 * 1024 * 1024;
constexpr std::uint64_t kMaxWireAtoms =
    static_cast<std::uint64_t>(std::numeric_limits<std::int32_t>::max());

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

bool AddPayloadBytes(size_t count, size_t element_size,
                     std::uint64_t& total) {
    if (count != 0 && element_size >
            std::numeric_limits<std::uint64_t>::max() / count) {
        return false;
    }
    const std::uint64_t bytes =
        static_cast<std::uint64_t>(count) * element_size;
    if (bytes > kMaxWirePayloadBytes ||
        total > kMaxWirePayloadBytes - bytes) {
        return false;
    }
    total += bytes;
    return true;
}

bool WriteAll(int fd, const void* bytes, size_t count) {
    const auto* cursor = static_cast<const unsigned char*>(bytes);
    while (count > 0) {
        const ssize_t written = ::send(fd, cursor, count, MSG_NOSIGNAL);
        if (written < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        if (written == 0) return false;
        cursor += written;
        count -= static_cast<size_t>(written);
    }
    return true;
}

bool ReadAll(int fd, void* bytes, size_t count) {
    auto* cursor = static_cast<unsigned char*>(bytes);
    while (count > 0) {
        const ssize_t received = ::read(fd, cursor, count);
        if (received < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        if (received == 0) return false;
        cursor += received;
        count -= static_cast<size_t>(received);
    }
    return true;
}

template <typename T>
bool WriteScalar(int fd, const T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    return WriteAll(fd, &value, sizeof(value));
}

template <typename T>
bool ReadScalar(int fd, T& value) {
    static_assert(std::is_trivially_copyable_v<T>);
    return ReadAll(fd, &value, sizeof(value));
}

template <typename T>
bool WriteVector(int fd, const std::vector<T>& values) {
    const std::uint64_t count = static_cast<std::uint64_t>(values.size());
    return WriteScalar(fd, count) &&
           (values.empty() ||
            WriteAll(fd, values.data(), values.size() * sizeof(T)));
}

template <typename T>
bool ReadVectorExact(int fd, size_t expected, std::vector<T>& values) {
    std::uint64_t count = 0;
    if (!ReadScalar(fd, count) || count != expected) return false;
    values.resize(expected);
    return expected == 0 ||
           ReadAll(fd, values.data(), expected * sizeof(T));
}

template <typename T>
bool ReadVectorBounded(int fd, std::uint64_t maximum,
                       std::vector<T>& values) {
    std::uint64_t count = 0;
    if (!ReadScalar(fd, count) || count > maximum ||
        count > std::numeric_limits<size_t>::max() / sizeof(T)) {
        return false;
    }
    values.resize(static_cast<size_t>(count));
    return values.empty() ||
           ReadAll(fd, values.data(), values.size() * sizeof(T));
}

bool WriteString(int fd, const std::string& value) {
    const std::uint64_t size = static_cast<std::uint64_t>(value.size());
    return WriteScalar(fd, size) &&
           (value.empty() || WriteAll(fd, value.data(), value.size()));
}

bool ReadString(int fd, std::string& value) {
    std::uint64_t size = 0;
    if (!ReadScalar(fd, size) || size > kMaxWireStringBytes) return false;
    value.resize(static_cast<size_t>(size));
    return value.empty() || ReadAll(fd, value.data(), value.size());
}

}  // namespace

bool WriteInput(int fd, const Input& input) {
    return WriteScalar(fd, kInputMagic) &&
           WriteScalar(fd, kWireVersion) &&
           WriteScalar(fd, input.net_charge) &&
           WriteVector(fd, input.atoms) &&
           WriteVector(fd, input.coordinates) &&
           WriteScalar(fd, kWireTrailer);
}

bool ReadInput(int fd, Input& input, std::string& error) {
    std::uint64_t magic = 0;
    std::uint32_t version = 0;
    if (!ReadScalar(fd, magic) || !ReadScalar(fd, version) ||
        magic != kInputMagic || version != kWireVersion) {
        error = "invalid or truncated MOPAC worker input header";
        return false;
    }
    if (!ReadScalar(fd, input.net_charge) ||
        !ReadVectorBounded(fd, kMaxWireAtoms, input.atoms)) {
        error = "invalid or truncated MOPAC worker atom input";
        return false;
    }
    const std::uint64_t coordinate_count =
        static_cast<std::uint64_t>(input.atoms.size()) * 3;
    if (!ReadVectorExact(fd, static_cast<size_t>(coordinate_count),
                         input.coordinates)) {
        error = "invalid or truncated MOPAC worker coordinates";
        return false;
    }
    std::uint64_t trailer = 0;
    if (!ReadScalar(fd, trailer) || trailer != kWireTrailer) {
        error = "MOPAC worker input trailer is missing";
        return false;
    }
    return true;
}

bool WriteSuccess(int fd, const Output& output) {
    if (!WriteScalar(fd, kOutputMagic) ||
        !WriteScalar(fd, kWireVersion) ||
        !WriteScalar(fd, kWireSuccess) ||
        !WriteScalar(fd, output.natom) ||
        !WriteScalar(fd, output.ao_max_orbitals) ||
        !WriteAll(fd, output.state_dimensions.data(),
                  output.state_dimensions.size() * sizeof(std::int32_t)) ||
        !WriteScalar(fd, output.heat) ||
        !WriteAll(fd, output.dipole.data(), 3 * sizeof(double)) ||
        !WriteAll(fd, output.dipole_point_charge.data(), 3 * sizeof(double)) ||
        !WriteAll(fd, output.dipole_hybridization.data(),
                  3 * sizeof(double))) {
        return false;
    }

    const bool vectors_ok =
        WriteVector(fd, output.charge) &&
        WriteVector(fd, output.bond_index) &&
        WriteVector(fd, output.bond_atom) &&
        WriteVector(fd, output.bond_order) &&
        WriteVector(fd, output.ao_orbitals) &&
        WriteVector(fd, output.atom_ao_density_fortran) &&
        WriteVector(fd, output.bond_ao_density_fortran) &&
        WriteVector(fd, output.nbonds) &&
        WriteVector(fd, output.ibonds_fortran) &&
        WriteVector(fd, output.iorbs) &&
        WriteVector(fd, output.ncf) &&
        WriteVector(fd, output.nce) &&
        WriteVector(fd, output.icocc) &&
        WriteVector(fd, output.icvir) &&
        WriteVector(fd, output.cocc) &&
        WriteVector(fd, output.cvir) &&
        WriteVector(fd, output.lmo_energy) &&
        WriteVector(fd, output.occupied_atom_offsets) &&
        WriteVector(fd, output.virtual_atom_offsets) &&
        WriteVector(fd, output.occupied_coefficient_offsets) &&
        WriteVector(fd, output.virtual_coefficient_offsets);
    return vectors_ok && WriteScalar(fd, kWireTrailer);
}

bool WriteFailure(int fd, const std::string& error) {
    return WriteScalar(fd, kOutputMagic) &&
           WriteScalar(fd, kWireVersion) &&
           WriteScalar(fd, kWireFailure) &&
           WriteString(fd, error) &&
           WriteScalar(fd, kWireTrailer);
}

bool ReadOutput(int fd, size_t expected_natom, Output& output,
                std::string& child_error, std::string& protocol_error) {
    std::uint64_t magic = 0;
    std::uint32_t version = 0;
    std::int32_t status = -1;
    if (!ReadScalar(fd, magic) || !ReadScalar(fd, version) ||
        !ReadScalar(fd, status)) {
        protocol_error = "MOPAC worker produced a truncated protocol header";
        return false;
    }
    if (magic != kOutputMagic || version != kWireVersion) {
        protocol_error = "MOPAC worker produced an incompatible protocol header";
        return false;
    }
    if (status == kWireFailure) {
        std::uint64_t trailer = 0;
        if (!ReadString(fd, child_error) || !ReadScalar(fd, trailer) ||
            trailer != kWireTrailer) {
            protocol_error = "MOPAC worker failure diagnostic was truncated";
            return false;
        }
        if (child_error.empty()) child_error = "MOPAC worker reported failure";
        return true;
    }
    if (status != kWireSuccess) {
        protocol_error = "MOPAC worker returned an unknown protocol status";
        return false;
    }

    if (!ReadScalar(fd, output.natom) ||
        !ReadScalar(fd, output.ao_max_orbitals) ||
        !ReadAll(fd, output.state_dimensions.data(),
                 output.state_dimensions.size() * sizeof(std::int32_t)) ||
        !ReadScalar(fd, output.heat) ||
        !ReadAll(fd, output.dipole.data(), 3 * sizeof(double)) ||
        !ReadAll(fd, output.dipole_point_charge.data(), 3 * sizeof(double)) ||
        !ReadAll(fd, output.dipole_hybridization.data(), 3 * sizeof(double))) {
        protocol_error = "MOPAC worker produced truncated scalar metadata";
        return false;
    }
    if (output.natom < 0 ||
        static_cast<size_t>(output.natom) != expected_natom ||
        output.ao_max_orbitals <= 0 || output.ao_max_orbitals > 64) {
        protocol_error = "MOPAC worker returned invalid scalar dimensions";
        return false;
    }
    for (std::int32_t extent : output.state_dimensions) {
        if (extent < 0) {
            protocol_error = "MOPAC worker returned a negative state dimension";
            return false;
        }
    }
    if (output.state_dimensions[0] != output.natom) {
        protocol_error = "MOPAC worker atom/state dimensions disagree";
        return false;
    }

    const size_t natom = expected_natom;
    const size_t width = static_cast<size_t>(output.ao_max_orbitals);
    if (!ReadVectorExact(fd, natom, output.charge) ||
        !ReadVectorExact(fd, natom + 1, output.bond_index)) {
        protocol_error = "MOPAC worker truncated atomic/CSC metadata";
        return false;
    }
    if (output.bond_index.empty() || output.bond_index.front() != 0) {
        protocol_error = "MOPAC worker returned invalid CSC offsets";
        return false;
    }
    for (size_t i = 1; i < output.bond_index.size(); ++i) {
        if (output.bond_index[i] < output.bond_index[i - 1]) {
            protocol_error = "MOPAC worker returned non-monotonic CSC offsets";
            return false;
        }
    }
    if (output.bond_index.back() < 0) {
        protocol_error = "MOPAC worker returned a negative CSC size";
        return false;
    }
    const size_t entries = static_cast<size_t>(output.bond_index.back());
    if (static_cast<std::uint64_t>(entries) >
        static_cast<std::uint64_t>(natom) * natom) {
        protocol_error = "MOPAC worker returned an impossible CSC size";
        return false;
    }
    size_t atom_density_count = 0;
    size_t bond_density_count = 0;
    if (!CheckedProduct({natom, width, width}, atom_density_count) ||
        !CheckedProduct({entries, width, width}, bond_density_count)) {
        protocol_error = "MOPAC worker density dimensions overflow size_t";
        return false;
    }

    const size_t occupied = static_cast<size_t>(output.state_dimensions[1]);
    const size_t virtual_count =
        static_cast<size_t>(output.state_dimensions[2]);
    const size_t icocc_dim = static_cast<size_t>(output.state_dimensions[3]);
    const size_t icvir_dim = static_cast<size_t>(output.state_dimensions[4]);
    const size_t cocc_dim = static_cast<size_t>(output.state_dimensions[5]);
    const size_t cvir_dim = static_cast<size_t>(output.state_dimensions[6]);
    if (occupied > std::numeric_limits<size_t>::max() - virtual_count) {
        protocol_error = "MOPAC worker LMO dimensions overflow size_t";
        return false;
    }
    size_t orbital_capacity = 0;
    if (!CheckedProduct({natom, width}, orbital_capacity) ||
        occupied + virtual_count > orbital_capacity) {
        protocol_error = "MOPAC worker returned impossible LMO dimensions";
        return false;
    }

    std::uint64_t payload_bytes = 0;
    const auto add_i32 = [&](size_t count) {
        return AddPayloadBytes(count, sizeof(std::int32_t), payload_bytes);
    };
    const auto add_f64 = [&](size_t count) {
        return AddPayloadBytes(count, sizeof(double), payload_bytes);
    };
    if (!add_f64(natom) || !add_i32(natom + 1) ||
        !add_i32(entries) || !add_f64(entries) || !add_i32(natom) ||
        !add_f64(atom_density_count) || !add_f64(bond_density_count) ||
        !add_i32(natom) || !add_i32(9 * natom) || !add_i32(natom) ||
        !add_i32(occupied) || !add_i32(virtual_count) ||
        !add_i32(icocc_dim) || !add_i32(icvir_dim) ||
        !add_f64(cocc_dim) || !add_f64(cvir_dim) ||
        !add_f64(occupied + virtual_count) ||
        !add_i32(2 * occupied) || !add_i32(2 * virtual_count)) {
        protocol_error = "MOPAC worker payload exceeds the protocol safety ceiling";
        return false;
    }

    const bool vectors_ok =
        ReadVectorExact(fd, entries, output.bond_atom) &&
        ReadVectorExact(fd, entries, output.bond_order) &&
        ReadVectorExact(fd, natom, output.ao_orbitals) &&
        ReadVectorExact(fd, atom_density_count,
                        output.atom_ao_density_fortran) &&
        ReadVectorExact(fd, bond_density_count,
                        output.bond_ao_density_fortran) &&
        ReadVectorExact(fd, natom, output.nbonds) &&
        ReadVectorExact(fd, 9 * natom, output.ibonds_fortran) &&
        ReadVectorExact(fd, natom, output.iorbs) &&
        ReadVectorExact(fd, occupied, output.ncf) &&
        ReadVectorExact(fd, virtual_count, output.nce) &&
        ReadVectorExact(fd, icocc_dim, output.icocc) &&
        ReadVectorExact(fd, icvir_dim, output.icvir) &&
        ReadVectorExact(fd, cocc_dim, output.cocc) &&
        ReadVectorExact(fd, cvir_dim, output.cvir) &&
        ReadVectorExact(fd, occupied + virtual_count, output.lmo_energy) &&
        ReadVectorExact(fd, occupied, output.occupied_atom_offsets) &&
        ReadVectorExact(fd, virtual_count, output.virtual_atom_offsets) &&
        ReadVectorExact(fd, occupied,
                        output.occupied_coefficient_offsets) &&
        ReadVectorExact(fd, virtual_count,
                        output.virtual_coefficient_offsets);
    if (!vectors_ok) {
        protocol_error = "MOPAC worker produced a truncated data payload";
        return false;
    }
    std::uint64_t trailer = 0;
    if (!ReadScalar(fd, trailer) || trailer != kWireTrailer) {
        protocol_error = "MOPAC worker payload trailer is missing";
        return false;
    }
    return true;
}

}  // namespace nmr::mopac_worker
