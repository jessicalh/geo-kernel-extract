#include "TestEnvironment.h"
#include "DirectionalTestHelpers.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "ApbsFieldResult.h"
#include "ApbsEfieldTimeSeriesTrajectoryResult.h"
#include "ApbsEfgTimeSeriesTrajectoryResult.h"
#include "CoulombResult.h"
#include "SpatialIndexResult.h"
#include "CalculatorConfig.h"
#include "PhysicalConstants.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include <highfive/H5File.hpp>
#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <cstring>
#include <utility>
#include <unistd.h>

#ifndef NMR_TEST_PYTHON_EXECUTABLE
#error "NMR_TEST_PYTHON_EXECUTABLE must be defined"
#endif
#ifndef NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT
#error "NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT must be defined"
#endif

using namespace nmr;

namespace {

template <typename T>
std::vector<T> ReadApbsH5Flat(
        const std::filesystem::path& path,
        const std::string& dataset,
        std::vector<std::size_t>* dimensions = nullptr) {
    HighFive::File file(path.string(), HighFive::File::ReadOnly);
    auto data_set = file.getDataSet(dataset);
    const auto dims = data_set.getSpace().getDimensions();
    if (dimensions) *dimensions = dims;
    std::size_t count = 1;
    for (const std::size_t dim : dims) count *= dim;
    std::vector<T> values(count);
    if (!values.empty()) data_set.read(values.data());
    return values;
}

std::string ReadApbsH5StringAttribute(
        const std::filesystem::path& path,
        const std::string& object,
        const std::string& attribute) {
    HighFive::File file(path.string(), HighFive::File::ReadOnly);
    std::string value;
    if (file.exist(object) && file.getDataSet(object).hasAttribute(attribute))
        file.getDataSet(object).getAttribute(attribute).read(value);
    return value;
}

bool CheckApbsClampDirectionalMetadata(
        const std::filesystem::path& path) {
    const std::string root = "/trajectory/apbs_efield_time_series/";
    const std::string mask = root + "clamp_mask";
    const std::string scale = root + "clamp_scale";
    const std::string xyz = root + "xyz";
    const std::string t2 =
        "/trajectory/apbs_efg_time_series/t2";
    const std::string frame = "lab_fixed_apbs_finite_difference_grid";
    const std::string mask_law =
        "continuum rotation/translation/reflection-invariant scalar "
        "threshold diagnostic derived from |E|; the live axis-aligned "
        "finite-difference APBS solve has no exact O(3) law";
    const std::string scale_law =
        "continuum rotation/translation/reflection-invariant scalar "
        "derived from |E| and the configured clamp threshold; the live "
        "axis-aligned finite-difference APBS solve has no exact O(3) law";
    const std::string xyz_law =
        "continuum polar_vector: v'=R v; translation invariant. The live "
        "axis-aligned finite-difference APBS solve has no exact O(3) law; "
        "transformed production reruns use the recorded 1.8e-2 V/Angstrom "
        "absolute + 5e-2 relative finite-grid envelope";
    const std::string t2_law =
        "continuum even_rank2: T'=R T R^T; emitted values are project-native "
        "T2 coefficients. The live axis-aligned finite-difference APBS solve "
        "has no exact O(3) law; transformed production reruns use the recorded "
        "4e-2 V/Angstrom^2 absolute + 5e-2 relative finite-grid envelope";
    return ReadApbsH5StringAttribute(path, xyz, "transformation") == xyz_law
        && ReadApbsH5StringAttribute(path, t2, "transformation") == t2_law
        && ReadApbsH5StringAttribute(path, mask, "units") == "dimensionless"
        && ReadApbsH5StringAttribute(path, scale, "units") == "dimensionless"
        && ReadApbsH5StringAttribute(path, mask, "coordinate_frame") == frame
        && ReadApbsH5StringAttribute(path, scale, "coordinate_frame") == frame
        && ReadApbsH5StringAttribute(path, mask, "parity") == "mixed"
        && ReadApbsH5StringAttribute(path, scale, "parity") == "mixed"
        && ReadApbsH5StringAttribute(path, mask, "transformation") == mask_law
        && ReadApbsH5StringAttribute(path, scale, "transformation") == scale_law
        && ReadApbsH5StringAttribute(path, mask, "validity").find(
               "1 iff the canonical reaction E-field row was") !=
               std::string::npos
        && ReadApbsH5StringAttribute(path, scale, "validity").find(
               "absent-source frames are NaN") != std::string::npos;
}

void WriteApbsTestConfig(const std::filesystem::path& path,
                         double pdie,
                         double sdie,
                         double ionic_strength,
                         double clamp,
                         double manual_padding_A = 70.0,
                         double manual_min_dim_A = 70.0) {
    std::ofstream out(path);
    // 161 is the production grid and avoids APBS's coarse-grid direct-solver
    // fallback, whose legacy teardown is not safe at 65^3 for this protein.
    out << "apbs_grid_dim = 161\n"
        << "apbs_manual_grid_padding_A = " << manual_padding_A << "\n"
        << "apbs_manual_grid_min_dim_A = " << manual_min_dim_A << "\n"
        << "apbs_protein_dielectric = " << pdie << "\n"
        << "apbs_solvent_dielectric = " << sdie << "\n"
        << "apbs_temperature_K = 298.15\n"
        << "apbs_ionic_strength_M = " << ionic_strength << "\n"
        << "efield_magnitude_sanity_clamp = " << clamp << "\n";
}

std::unique_ptr<Protein> BuildSingleChargeProbeProtein() {
    auto protein = std::make_unique<Protein>();
    Residue residue;
    residue.type = AminoAcid::ALA;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const size_t ri = protein->AddResidue(std::move(residue));

    auto source = Atom::Create(Element::C);
    source->pdb_atom_name = "CA";
    source->residue_index = ri;
    const size_t source_i = protein->AddAtom(std::move(source));
    auto probe = Atom::Create(Element::C);
    probe->pdb_atom_name = "CB";
    probe->residue_index = ri;
    const size_t probe_i = protein->AddAtom(std::move(probe));
    protein->MutableResidueAt(ri).atom_indices = {source_i, probe_i};

    // 70 A / (161 - 1) = 0.4375 A, so both sites lie exactly on the
    // configured grid.  Only atom 0 is charged; atom 1 is a neutral field
    // probe at r=3.5 A along +x from the source.
    const std::vector<Vec3> positions = {
        Vec3(-1.75, 0.0, 0.0), Vec3(1.75, 0.0, 0.0)};
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "single-charge APBS sign probe");
    return protein;
}

bool WriteApbsTrajectoryFrame(
        const ProteinConformation& source,
        const std::filesystem::path& path,
        std::size_t frame_index,
        double time_ps) {
    auto storage = TrajectoryProtein::CreateForTesting(
        BuildSingleChargeProbeProtein());
    if (!storage || storage->AtomCount() != source.AtomCount()) return false;

    auto efield = ApbsEfieldTimeSeriesTrajectoryResult::Create(*storage);
    auto efg = ApbsEfgTimeSeriesTrajectoryResult::Create(*storage);
    if (!efield || !efg ||
        !storage->AttachResult(std::move(efield)) ||
        !storage->AttachResult(std::move(efg))) {
        return false;
    }

    Trajectory dummy("", "", "");
    storage->DispatchCompute(source, dummy, frame_index, time_ps);
    storage->FinalizeAllResults(dummy);
    HighFive::File file(path.string(), HighFive::File::Truncate);
    storage->WriteH5(file);
    return file.exist("/trajectory/apbs_efield_time_series/xyz") &&
           file.exist("/trajectory/apbs_efg_time_series/t2");
}

bool CheckApbsGridProvenance(
        const ProteinConformation& source,
        const std::filesystem::path& path) {
    if (!source.HasResult<ApbsFieldResult>()) return false;
    const auto& apbs = source.Result<ApbsFieldResult>();
    const std::vector<std::uint64_t> expected_dims(
        apbs.GridDims().begin(), apbs.GridDims().end());
    for (const std::string& root : {
            std::string("/trajectory/apbs_efield_time_series/"),
            std::string("/trajectory/apbs_efg_time_series/")}) {
        std::vector<std::size_t> dimensions;
        const auto emitted = ReadApbsH5Flat<std::uint64_t>(
            path, root + "apbs_grid_dims_per_frame", &dimensions);
        if (dimensions != std::vector<std::size_t>{1u, 3u} ||
            emitted != expected_dims) {
            return false;
        }
    }
    const std::array<std::pair<const char*, std::array<double, 3>>, 3>
        channels{{
            {"apbs_grid_lengths_A_per_frame", apbs.GridLengthsA()},
            {"apbs_grid_origin_A_per_frame", apbs.GridOriginA()},
            {"apbs_grid_spacing_A_per_frame", apbs.GridSpacingA()},
        }};
    for (const auto& [name, expected] : channels) {
        std::vector<std::size_t> efield_dims;
        std::vector<std::size_t> efg_dims;
        const auto efield = ReadApbsH5Flat<double>(
            path, "/trajectory/apbs_efield_time_series/" +
                      std::string(name),
            &efield_dims);
        const auto efg = ReadApbsH5Flat<double>(
            path, "/trajectory/apbs_efg_time_series/" +
                      std::string(name),
            &efg_dims);
        if (efield_dims != std::vector<std::size_t>{1u, 3u} ||
            efg_dims != efield_dims || efield.size() != 3u ||
            efg.size() != 3u) {
            return false;
        }
        for (std::size_t component = 0; component < 3; ++component) {
            if (!std::isfinite(efield[component]) ||
                !std::isfinite(efg[component]) ||
                efield[component] != expected[component] ||
                efg[component] != expected[component]) {
                return false;
            }
        }
    }
    return true;
}

bool CheckApbsFrameProvenance(
        const std::filesystem::path& path,
        std::size_t expected_frame_index,
        double expected_time_ps) {
    for (const std::string& root : {
            std::string("/trajectory/apbs_efield_time_series/"),
            std::string("/trajectory/apbs_efg_time_series/")}) {
        std::vector<std::size_t> frame_index_dims;
        std::vector<std::size_t> frame_time_dims;
        const auto frame_indices = ReadApbsH5Flat<std::uint64_t>(
            path, root + "frame_indices", &frame_index_dims);
        const auto frame_times = ReadApbsH5Flat<double>(
            path, root + "frame_times", &frame_time_dims);
        if (frame_index_dims != std::vector<std::size_t>{1u} ||
            frame_time_dims != frame_index_dims ||
            frame_indices != std::vector<std::uint64_t>{
                static_cast<std::uint64_t>(expected_frame_index)} ||
            frame_times != std::vector<double>{expected_time_ps}) {
            return false;
        }
    }
    return true;
}

bool AttachSingleCharge(ProteinConformation& conf, double charge_e) {
    std::vector<AtomChargeRadius> values = {
        {charge_e, 1.5, ChargeAssignmentStatus::Matched},
        {0.0, 1.5, ChargeAssignmentStatus::Matched},
    };
    PreloadedChargeSource source(std::move(values), ForceField::Amber_ff14SB);
    auto charges = ChargeAssignmentResult::Compute(conf, source);
    return charges != nullptr && conf.AttachResult(std::move(charges));
}

std::vector<double> ReadFloat64NpyPayload(
        const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return {};
    const std::vector<unsigned char> bytes(
        (std::istreambuf_iterator<char>(in)),
        std::istreambuf_iterator<char>());
    if (bytes.size() < 10 || std::memcmp(bytes.data(), "\x93NUMPY", 6) != 0)
        return {};

    size_t header_len = 0;
    size_t payload_offset = 0;
    if (bytes[6] == 1) {
        header_len = static_cast<size_t>(bytes[8]) |
                     (static_cast<size_t>(bytes[9]) << 8);
        payload_offset = 10 + header_len;
    } else if (bytes[6] == 2 || bytes[6] == 3) {
        if (bytes.size() < 12) return {};
        header_len = static_cast<size_t>(bytes[8]) |
                     (static_cast<size_t>(bytes[9]) << 8) |
                     (static_cast<size_t>(bytes[10]) << 16) |
                     (static_cast<size_t>(bytes[11]) << 24);
        payload_offset = 12 + header_len;
    } else {
        return {};
    }
    if (payload_offset > bytes.size() ||
        (bytes.size() - payload_offset) % sizeof(double) != 0)
        return {};

    std::vector<double> data(
        (bytes.size() - payload_offset) / sizeof(double));
    std::memcpy(data.data(), bytes.data() + payload_offset,
                data.size() * sizeof(double));
    return data;
}

std::vector<std::uint8_t> ReadUInt8NpyPayload(
        const std::filesystem::path& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) return {};
    const std::vector<unsigned char> bytes(
        (std::istreambuf_iterator<char>(in)),
        std::istreambuf_iterator<char>());
    if (bytes.size() < 10 || std::memcmp(bytes.data(), "\x93NUMPY", 6) != 0)
        return {};

    size_t header_len = 0;
    size_t payload_offset = 0;
    if (bytes[6] == 1) {
        header_len = static_cast<size_t>(bytes[8]) |
                     (static_cast<size_t>(bytes[9]) << 8);
        payload_offset = 10 + header_len;
    } else if (bytes[6] == 2 || bytes[6] == 3) {
        if (bytes.size() < 12) return {};
        header_len = static_cast<size_t>(bytes[8]) |
                     (static_cast<size_t>(bytes[9]) << 8) |
                     (static_cast<size_t>(bytes[10]) << 16) |
                     (static_cast<size_t>(bytes[11]) << 24);
        payload_offset = 12 + header_len;
    } else {
        return {};
    }
    if (payload_offset > bytes.size()) return {};
    return std::vector<std::uint8_t>(bytes.begin() + payload_offset,
                                     bytes.end());
}

void RemoveApbsNpys(const std::filesystem::path& output_dir) {
    for (const char* filename : {
            "apbs_E.npy", "apbs_efg.npy", "apbs_phi.npy",
            "apbs_E_clamp_mask.npy", "apbs_E_clamp_scale.npy",
            "apbs_nonfinite_sanitizer_mask.npy",
            "apbs_E_total_diagnostic.npy",
            "apbs_efg_total_diagnostic.npy"}) {
        std::filesystem::remove(output_dir / filename);
    }
}

void RemoveCoulombNpys(const std::filesystem::path& output_dir) {
    for (const char* filename : {
            "coulomb_efg.npy", "coulomb_efg_t2.npy", "coulomb_E.npy",
            "coulomb_E_backbone.npy", "coulomb_E_sidechain.npy",
            "coulomb_E_aromatic.npy", "coulomb_efg_backbone.npy",
            "coulomb_efg_sidechain.npy", "coulomb_efg_aromatic.npy",
            "coulomb_scalars.npy", "coulomb_aromatic_E_proj.npy",
            "coulomb_aromatic_n_src.npy", "coulomb_E_solvent.npy",
            "coulomb_efg_solvent.npy"}) {
        std::filesystem::remove(output_dir / filename);
    }
}

void RunHomogeneousReferenceForcingFunctionInChild() {
    const auto nonce = std::chrono::steady_clock::now()
                           .time_since_epoch().count();
    const auto config_path = std::filesystem::temp_directory_path() /
        ("apbs_homogeneous_reference_" + std::to_string(::getpid()) +
         "_" + std::to_string(nonce) + ".toml");
    const auto output_dir = std::filesystem::temp_directory_path() /
        ("apbs_homogeneous_reference_" + std::to_string(::getpid()) +
         "_" + std::to_string(nonce));
    WriteApbsTestConfig(config_path, 1.0, 1.0, 0.0, 100.0);
    CalculatorConfig::Load(config_path.string());

    auto protein = BuildSingleChargeProbeProtein();
    bool ok = protein != nullptr;
    if (ok) {
        auto& conf = protein->Conformation();
        ok = AttachSingleCharge(conf, 1.0);
        auto apbs = ok ? ApbsFieldResult::Compute(conf) : nullptr;
        ok = ok && apbs != nullptr;

        double max_reaction = 0.0;
        double max_total_diagnostic = 0.0;
        if (ok) {
            for (size_t i = 0; i < conf.AtomCount(); ++i) {
                const auto& atom = conf.AtomAt(i);
                max_reaction = std::max(max_reaction,
                    atom.apbs_efield.norm() + atom.apbs_efg.norm() +
                    std::abs(atom.apbs_phi));
                max_total_diagnostic = std::max(max_total_diagnostic,
                    atom.apbs_efield_total_diagnostic.norm() +
                    atom.apbs_efg_total_diagnostic.norm());
                ok = ok && atom.apbs_efield_clamp_mask == 0u
                        && atom.apbs_efield_clamp_scale == 1.0;
            }
            // With production physics set equal to the independently defined
            // homogeneous reference, the two deterministic grids must cancel.
            ok = ok && max_reaction < 1e-8
                    && max_total_diagnostic > 1e-8;

            // Independently forced unclamped serialization: exact
            // total-reference cancellation cannot cross the 100 V/A clamp.
            // Read the actual NPY payloads so swapped/default mask or scale
            // arrays cannot hide behind the in-memory checks above.
            ok = ok && std::filesystem::create_directory(output_dir);
            ok = ok && apbs->WriteFeatures(conf, output_dir.string()) == 8;
            const auto emitted_mask = ReadUInt8NpyPayload(
                output_dir / "apbs_E_clamp_mask.npy");
            const auto emitted_scale = ReadFloat64NpyPayload(
                output_dir / "apbs_E_clamp_scale.npy");
            ok = ok && emitted_mask.size() == conf.AtomCount()
                    && emitted_scale.size() == conf.AtomCount();
            if (ok) {
                for (size_t i = 0; i < conf.AtomCount(); ++i) {
                    ok = ok && emitted_mask[i] == 0u
                            && emitted_scale[i] == 1.0;
                }
            }
        }
    }

    RemoveApbsNpys(output_dir);
    std::filesystem::remove(output_dir);
    std::filesystem::remove(config_path);
    _exit(ok ? 0 : 1);
}

void RunDirectionalCovarianceForcingFunctionInChild(
        const std::filesystem::path& metrics_path) {
    using nmr::test::directional::EvenRank2;
    using nmr::test::directional::Polar;
    using nmr::test::directional::Positions;
    using nmr::test::directional::RotateNativeT2;
    using nmr::test::directional::SeededTransform;

    const auto nonce = std::chrono::steady_clock::now()
                           .time_since_epoch().count();
    const auto config_path = std::filesystem::temp_directory_path() /
        ("apbs_directional_" + std::to_string(::getpid()) + "_" +
         std::to_string(nonce) + ".toml");
    const auto output_root = std::filesystem::temp_directory_path() /
        ("apbs_directional_npy_" + std::to_string(::getpid()) + "_" +
         std::to_string(nonce));
    // A cubic domain isolates the finite Cartesian stencil from the separate
    // orientation-dependence of axis-wise bounding-box grid lengths.
    WriteApbsTestConfig(config_path, 1.0, 78.54, 0.0, 100.0,
                        /*manual_padding_A=*/0.0,
                        /*manual_min_dim_A=*/70.0);
    CalculatorConfig::Load(config_path.string());

    auto protein = BuildSingleChargeProbeProtein();
    bool ok = protein != nullptr;
    if (!ok) _exit(1);
    // Do not privilege the reference solve by placing its charge pair exactly
    // on one Cartesian grid axis.  Both the reference and transformed solves
    // start from generic, deterministically seeded sub-grid orientations.
    const auto initial_orientation =
        SeededTransform(0xA9B5F13CULL, /*improper=*/false);
    auto& original = protein->AddConformation(
        Positions(initial_orientation,
                  protein->Conformation().Positions()),
        "APBS generic reference orientation");
    ok = ok && AttachSingleCharge(original, 1.0);
    auto original_result = ok ? ApbsFieldResult::Compute(original) : nullptr;
    ok = ok && original_result != nullptr;
    if (ok) ok = original.AttachResult(std::move(original_result));
    if (ok) ok = original.AttachResult(GeometryResult::Compute(original));
    if (ok) ok = original.AttachResult(SpatialIndexResult::Compute(original));
    auto original_coulomb = ok ? CoulombResult::Compute(original) : nullptr;
    ok = ok && original_coulomb != nullptr;
    if (ok) ok = original.AttachResult(std::move(original_coulomb));

    const auto original_dir = output_root / "original";
    ok = ok && std::filesystem::create_directories(original_dir);
    if (ok) {
        ok = original.Result<ApbsFieldResult>().WriteFeatures(
                 original, original_dir.string()) == 8
            && original.Result<CoulombResult>().WriteFeatures(
                 original, original_dir.string()) == 14;
    }
    if (ok) {
        ok = nmr::test::directional::RunNumpyAllowPickleFalse(
                 NMR_TEST_PYTHON_EXECUTABLE,
                 NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                 {original_dir / "apbs_E.npy",
                  original_dir / "apbs_efg.npy",
                  original_dir / "apbs_phi.npy",
                  original_dir / "apbs_E_clamp_mask.npy",
                  original_dir / "apbs_E_clamp_scale.npy",
                  original_dir / "apbs_nonfinite_sanitizer_mask.npy",
                  original_dir / "apbs_E_total_diagnostic.npy",
                  original_dir / "apbs_efg_total_diagnostic.npy",
                  original_dir / "coulomb_E_solvent.npy",
                  original_dir / "coulomb_efg_solvent.npy"}) == 0;
    }
    const auto original_E = ReadFloat64NpyPayload(
        original_dir / "apbs_E.npy");
    const auto original_E_total = ReadFloat64NpyPayload(
        original_dir / "apbs_E_total_diagnostic.npy");
    const auto original_efg = ReadFloat64NpyPayload(
        original_dir / "apbs_efg.npy");
    const auto original_efg_total = ReadFloat64NpyPayload(
        original_dir / "apbs_efg_total_diagnostic.npy");
    const auto original_phi = ReadFloat64NpyPayload(
        original_dir / "apbs_phi.npy");
    const auto original_sanitizer = ReadUInt8NpyPayload(
        original_dir / "apbs_nonfinite_sanitizer_mask.npy");
    const auto original_clamp = ReadUInt8NpyPayload(
        original_dir / "apbs_E_clamp_mask.npy");
    const auto original_clamp_scale = ReadFloat64NpyPayload(
        original_dir / "apbs_E_clamp_scale.npy");
    const auto original_solvent_E = ReadFloat64NpyPayload(
        original_dir / "coulomb_E_solvent.npy");
    const auto original_solvent_efg = ReadFloat64NpyPayload(
        original_dir / "coulomb_efg_solvent.npy");
    ok = ok && original_E.size() == original.AtomCount() * 3
            && original_E_total.size() == original.AtomCount() * 3
            && original_efg.size() == original.AtomCount() * 5
            && original_efg_total.size() == original.AtomCount() * 5
            && original_phi.size() == original.AtomCount()
            && original_sanitizer.size() == original.AtomCount()
            && original_clamp.size() == original.AtomCount()
            && original_clamp_scale.size() == original.AtomCount()
            && original_solvent_E == original_E
            && original_solvent_efg == original_efg;

    // Feed the same attached production owner into the real trajectory
    // accumulators and cross the exact H5 boundary.  APBS grid diagnostics
    // are lab-axis provenance, not vectors: CheckApbsGridProvenance pins
    // shape/finiteness/exact owner serialization without imposing Q laws.
    const auto original_h5 = original_dir / "apbs_directional.h5";
    if (ok) {
        ok = WriteApbsTrajectoryFrame(
            original, original_h5, /*frame_index=*/73u, /*time_ps=*/4.25);
    }
    std::vector<double> original_h5_E;
    std::vector<double> original_h5_efg;
    std::vector<std::uint8_t> original_h5_clamp;
    std::vector<double> original_h5_clamp_scale;
    if (ok) {
        std::vector<std::size_t> e_dims;
        std::vector<std::size_t> q_dims;
        std::vector<std::size_t> clamp_dims;
        std::vector<std::size_t> clamp_scale_dims;
        std::vector<std::size_t> e_mask_dims;
        std::vector<std::size_t> q_mask_dims;
        original_h5_E = ReadApbsH5Flat<double>(
            original_h5,
            "/trajectory/apbs_efield_time_series/xyz", &e_dims);
        original_h5_efg = ReadApbsH5Flat<double>(
            original_h5,
            "/trajectory/apbs_efg_time_series/t2", &q_dims);
        original_h5_clamp = ReadApbsH5Flat<std::uint8_t>(
            original_h5,
            "/trajectory/apbs_efield_time_series/clamp_mask",
            &clamp_dims);
        original_h5_clamp_scale = ReadApbsH5Flat<double>(
            original_h5,
            "/trajectory/apbs_efield_time_series/clamp_scale",
            &clamp_scale_dims);
        const auto e_mask = ReadApbsH5Flat<std::uint8_t>(
            original_h5,
            "/trajectory/apbs_efield_time_series/"
            "source_attached_per_frame",
            &e_mask_dims);
        const auto q_mask = ReadApbsH5Flat<std::uint8_t>(
            original_h5,
            "/trajectory/apbs_efg_time_series/"
            "source_attached_per_frame",
            &q_mask_dims);
        ok = e_dims == std::vector<std::size_t>{
                 original.AtomCount(), 1u, 3u}
            && q_dims == std::vector<std::size_t>{
                 original.AtomCount(), 1u, 5u}
            && clamp_dims == std::vector<std::size_t>{
                 original.AtomCount(), 1u}
            && clamp_scale_dims == clamp_dims
            && original_h5_clamp == original_clamp
            && original_h5_clamp_scale == original_clamp_scale
            && e_mask_dims == std::vector<std::size_t>{1u}
            && q_mask_dims == e_mask_dims
            && e_mask == std::vector<std::uint8_t>{1u}
            && q_mask == e_mask
            && CheckApbsClampDirectionalMetadata(original_h5)
            && CheckApbsGridProvenance(original, original_h5)
            && CheckApbsFrameProvenance(
                original_h5, /*expected_frame_index=*/73u,
                /*expected_time_ps=*/4.25);
    }
    const auto invariant_grid_dims =
        original.Result<ApbsFieldResult>().GridDims();

    // APBS is rerun on a fresh, axis-aligned finite-difference grid after
    // each rigid transform.  Unlike analytic Coulomb kernels, the sampled
    // PDE solution is only approximately O(3)-closed at 161^3.  The neutral
    // probe's recorded proper+improper maxima are 0.0146388 V/A (27.87%)
    // for E, 0.033567 V/A^2 (35.66%) for EFG, and 9.24672e-4 V for phi.
    // These bounds retain a small deterministic finite-grid margin while a
    // sign/parity error at the nonzero probe remains far outside tolerance.
    constexpr double kVectorAbs = 1.8e-2;
    constexpr double kVectorRel = 5.0e-2;
    constexpr double kTensorAbs = 4.0e-2;
    constexpr double kTensorRel = 5.0e-2;
    constexpr double kPhiAbs = 1.2e-3;
    double max_vector_error = 0.0;
    double max_vector_relative = 0.0;
    double max_tensor_error = 0.0;
    double max_tensor_relative = 0.0;
    double max_phi_error = 0.0;

    for (const bool improper : {false, true}) {
        const auto transform = SeededTransform(0xA9B5F13DULL, improper);
        auto& moved = protein->AddConformation(
            Positions(transform, original.Positions()),
            improper ? "APBS improper rerun" : "APBS proper rerun");
        ok = ok && AttachSingleCharge(moved, 1.0);
        auto moved_result = ok ? ApbsFieldResult::Compute(moved) : nullptr;
        ok = ok && moved_result != nullptr;
        if (ok) ok = moved.AttachResult(std::move(moved_result));
        if (ok) ok = moved.AttachResult(GeometryResult::Compute(moved));
        if (ok) ok = moved.AttachResult(SpatialIndexResult::Compute(moved));
        auto moved_coulomb = ok ? CoulombResult::Compute(moved) : nullptr;
        ok = ok && moved_coulomb != nullptr;
        if (ok) ok = moved.AttachResult(std::move(moved_coulomb));
        if (!ok) break;

        const auto moved_dir = output_root /
            (improper ? "improper" : "proper");
        ok = ok && std::filesystem::create_directories(moved_dir);
        if (ok) {
            ok = moved.Result<ApbsFieldResult>().WriteFeatures(
                     moved, moved_dir.string()) == 8
                && moved.Result<CoulombResult>().WriteFeatures(
                     moved, moved_dir.string()) == 14;
        }
        if (ok) {
            ok = nmr::test::directional::RunNumpyAllowPickleFalse(
                     NMR_TEST_PYTHON_EXECUTABLE,
                     NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                     {moved_dir / "apbs_E.npy",
                      moved_dir / "apbs_efg.npy",
                      moved_dir / "apbs_phi.npy",
                      moved_dir / "apbs_E_clamp_mask.npy",
                      moved_dir / "apbs_E_clamp_scale.npy",
                      moved_dir / "apbs_nonfinite_sanitizer_mask.npy",
                      moved_dir / "apbs_E_total_diagnostic.npy",
                      moved_dir / "apbs_efg_total_diagnostic.npy",
                      moved_dir / "coulomb_E_solvent.npy",
                      moved_dir / "coulomb_efg_solvent.npy"}) == 0;
        }
        const auto moved_E = ReadFloat64NpyPayload(
            moved_dir / "apbs_E.npy");
        const auto moved_E_total = ReadFloat64NpyPayload(
            moved_dir / "apbs_E_total_diagnostic.npy");
        const auto moved_efg = ReadFloat64NpyPayload(
            moved_dir / "apbs_efg.npy");
        const auto moved_efg_total = ReadFloat64NpyPayload(
            moved_dir / "apbs_efg_total_diagnostic.npy");
        const auto moved_phi = ReadFloat64NpyPayload(
            moved_dir / "apbs_phi.npy");
        const auto moved_sanitizer = ReadUInt8NpyPayload(
            moved_dir / "apbs_nonfinite_sanitizer_mask.npy");
        const auto moved_clamp = ReadUInt8NpyPayload(
            moved_dir / "apbs_E_clamp_mask.npy");
        const auto moved_clamp_scale = ReadFloat64NpyPayload(
            moved_dir / "apbs_E_clamp_scale.npy");
        const auto moved_solvent_E = ReadFloat64NpyPayload(
            moved_dir / "coulomb_E_solvent.npy");
        const auto moved_solvent_efg = ReadFloat64NpyPayload(
            moved_dir / "coulomb_efg_solvent.npy");
        ok = ok && moved_E.size() == original_E.size()
                && moved_E_total.size() == original_E_total.size()
                && moved_efg.size() == original_efg.size()
                && moved_efg_total.size() == original_efg_total.size()
                && moved_phi.size() == original_phi.size()
                && moved_sanitizer.size() == original_sanitizer.size()
                && moved_clamp.size() == original_clamp.size()
                && moved_clamp_scale.size() == original_clamp_scale.size()
                && moved_solvent_E == moved_E
                && moved_solvent_efg == moved_efg;
        if (!ok) break;

        const auto moved_h5 = moved_dir / "apbs_directional.h5";
        ok = WriteApbsTrajectoryFrame(
            moved, moved_h5, /*frame_index=*/73u, /*time_ps=*/4.25);
        if (!ok) break;
        std::vector<std::size_t> moved_e_dims;
        std::vector<std::size_t> moved_q_dims;
        std::vector<std::size_t> moved_clamp_dims;
        std::vector<std::size_t> moved_clamp_scale_dims;
        std::vector<std::size_t> moved_e_mask_dims;
        std::vector<std::size_t> moved_q_mask_dims;
        const auto moved_h5_E = ReadApbsH5Flat<double>(
            moved_h5, "/trajectory/apbs_efield_time_series/xyz",
            &moved_e_dims);
        const auto moved_h5_efg = ReadApbsH5Flat<double>(
            moved_h5, "/trajectory/apbs_efg_time_series/t2",
            &moved_q_dims);
        const auto moved_h5_clamp = ReadApbsH5Flat<std::uint8_t>(
            moved_h5,
            "/trajectory/apbs_efield_time_series/clamp_mask",
            &moved_clamp_dims);
        const auto moved_h5_clamp_scale = ReadApbsH5Flat<double>(
            moved_h5,
            "/trajectory/apbs_efield_time_series/clamp_scale",
            &moved_clamp_scale_dims);
        const auto moved_e_mask = ReadApbsH5Flat<std::uint8_t>(
            moved_h5,
            "/trajectory/apbs_efield_time_series/"
            "source_attached_per_frame",
            &moved_e_mask_dims);
        const auto moved_q_mask = ReadApbsH5Flat<std::uint8_t>(
            moved_h5,
            "/trajectory/apbs_efg_time_series/"
            "source_attached_per_frame",
            &moved_q_mask_dims);
        ok = moved_e_dims == std::vector<std::size_t>{
                 moved.AtomCount(), 1u, 3u}
            && moved_q_dims == std::vector<std::size_t>{
                 moved.AtomCount(), 1u, 5u}
            && moved_clamp_dims == std::vector<std::size_t>{
                 moved.AtomCount(), 1u}
            && moved_clamp_scale_dims == moved_clamp_dims
            && moved_h5_clamp == moved_clamp
            && moved_h5_clamp_scale == moved_clamp_scale
            && moved_e_mask_dims == std::vector<std::size_t>{1u}
            && moved_q_mask_dims == moved_e_mask_dims
            && moved_e_mask == std::vector<std::uint8_t>{1u}
            && moved_q_mask == moved_e_mask
            && moved_h5_E.size() == original_h5_E.size()
            && moved_h5_efg.size() == original_h5_efg.size()
            && moved.Result<ApbsFieldResult>().GridDims() ==
               invariant_grid_dims
            && CheckApbsClampDirectionalMetadata(moved_h5)
            && CheckApbsGridProvenance(moved, moved_h5)
            && CheckApbsFrameProvenance(
                moved_h5, /*expected_frame_index=*/73u,
                /*expected_time_ps=*/4.25);
        if (!ok) break;

        // Atom 0 is the deposited point charge.  A Cartesian charge-
        // assignment stencil does not provide a converged directional field
        // at its own singular support: sub-grid rotations change that
        // self-stencil by O(1).  Atom 1 is the neutral, finite-distance probe
        // and is the scientifically meaningful covariance observation.  The
        // writer payload and masks still include both rows and are size-
        // checked above; no source row is dropped from production output.
        for (std::size_t atom = 1; atom < original.AtomCount(); ++atom) {
            const auto& a = original.AtomAt(atom);
            const auto& b = moved.AtomAt(atom);
            for (const auto pair : {
                    std::pair<const Vec3*, const Vec3*>{
                        &a.apbs_efield, &b.apbs_efield},
                    std::pair<const Vec3*, const Vec3*>{
                        &a.apbs_efield_total_diagnostic,
                        &b.apbs_efield_total_diagnostic}}) {
                const Vec3 expected = Polar(transform, *pair.first);
                const double error = (*pair.second - expected).norm();
                max_vector_error = std::max(max_vector_error, error);
                max_vector_relative = std::max(
                    max_vector_relative,
                    error / std::max(1.0e-12, expected.norm()));
                ok = ok && error <= kVectorAbs + kVectorRel * expected.norm();
            }
            for (const auto pair : {
                    std::pair<const Mat3*, const Mat3*>{
                        &a.apbs_efg, &b.apbs_efg},
                    std::pair<const Mat3*, const Mat3*>{
                        &a.apbs_efg_total_diagnostic,
                        &b.apbs_efg_total_diagnostic}}) {
                const Mat3 expected = EvenRank2(transform, *pair.first);
                const double error = (*pair.second - expected).norm();
                max_tensor_error = std::max(max_tensor_error, error);
                max_tensor_relative = std::max(
                    max_tensor_relative,
                    error / std::max(1.0e-12, expected.norm()));
                ok = ok && error <= kTensorAbs + kTensorRel * expected.norm();
            }
            max_phi_error = std::max(
                max_phi_error, std::abs(b.apbs_phi - a.apbs_phi));
            ok = ok && std::abs(b.apbs_phi - a.apbs_phi) <= kPhiAbs;
            ok = ok && b.apbs_nonfinite_sanitizer_mask ==
                           a.apbs_nonfinite_sanitizer_mask;
            ok = ok && b.apbs_efield_clamp_mask == a.apbs_efield_clamp_mask;
            ok = ok && b.apbs_efield_clamp_scale ==
                           a.apbs_efield_clamp_scale;
            ok = ok && moved_clamp_scale[atom] ==
                           original_clamp_scale[atom];
            ok = ok && moved_h5_clamp[atom] == original_h5_clamp[atom]
                    && moved_h5_clamp_scale[atom] ==
                           original_h5_clamp_scale[atom];
            ok = ok && std::abs(b.apbs_efg_spherical.T0) <= 2.0e-10;
            ok = ok && std::abs(b.apbs_efg_total_diagnostic_spherical.T0)
                           <= 2.0e-10;
            for (int component = 0; component < 3; ++component) {
                ok = ok && std::abs(b.apbs_efg_spherical.T1[component])
                               <= 2.0e-10;
                ok = ok && std::abs(
                    b.apbs_efg_total_diagnostic_spherical.T1[component])
                               <= 2.0e-10;
            }

            // The same laws must hold for the bytes emitted by the owning
            // writer, including the project-native five-component T2 basis.
            const Vec3 serialized_E0(
                original_E[atom * 3], original_E[atom * 3 + 1],
                original_E[atom * 3 + 2]);
            const Vec3 serialized_E1(
                moved_E[atom * 3], moved_E[atom * 3 + 1],
                moved_E[atom * 3 + 2]);
            const Vec3 serialized_total_E0(
                original_E_total[atom * 3],
                original_E_total[atom * 3 + 1],
                original_E_total[atom * 3 + 2]);
            const Vec3 serialized_total_E1(
                moved_E_total[atom * 3], moved_E_total[atom * 3 + 1],
                moved_E_total[atom * 3 + 2]);
            ok = ok && (serialized_E1 - Polar(transform, serialized_E0)).norm()
                           <= kVectorAbs + kVectorRel * serialized_E0.norm();
            const Vec3 serialized_solvent_E0(
                original_solvent_E[atom * 3],
                original_solvent_E[atom * 3 + 1],
                original_solvent_E[atom * 3 + 2]);
            const Vec3 serialized_solvent_E1(
                moved_solvent_E[atom * 3], moved_solvent_E[atom * 3 + 1],
                moved_solvent_E[atom * 3 + 2]);
            ok = ok && (serialized_solvent_E1 -
                        Polar(transform, serialized_solvent_E0)).norm()
                           <= kVectorAbs +
                              kVectorRel * serialized_solvent_E0.norm();
            ok = ok && (serialized_total_E1 -
                        Polar(transform, serialized_total_E0)).norm()
                           <= kVectorAbs +
                              kVectorRel * serialized_total_E0.norm();

            const Vec3 h5_source_E(
                original_h5_E[atom * 3],
                original_h5_E[atom * 3 + 1],
                original_h5_E[atom * 3 + 2]);
            const Vec3 h5_moved_E(
                moved_h5_E[atom * 3],
                moved_h5_E[atom * 3 + 1],
                moved_h5_E[atom * 3 + 2]);
            ok = ok &&
                (h5_moved_E - Polar(transform, h5_source_E)).norm()
                    <= kVectorAbs + kVectorRel * h5_source_E.norm();

            for (const auto pair : {
                    std::pair<const std::vector<double>*,
                              const std::vector<double>*>{
                        &original_efg, &moved_efg},
                    std::pair<const std::vector<double>*,
                              const std::vector<double>*>{
                        &original_efg_total, &moved_efg_total}}) {
                SphericalTensor source;
                for (std::size_t component = 0; component < 5; ++component) {
                    source.T2[component] =
                        (*pair.first)[atom * 5 + component];
                }
                const SphericalTensor expected =
                    RotateNativeT2(transform, source);
                for (std::size_t component = 0; component < 5; ++component) {
                    ok = ok && std::abs(
                        (*pair.second)[atom * 5 + component] -
                        expected.T2[component]) <=
                        kTensorAbs + kTensorRel *
                            std::abs(expected.T2[component]);
                }
            }
            {
                SphericalTensor source;
                for (std::size_t component = 0; component < 5; ++component)
                    source.T2[component] =
                        original_solvent_efg[atom * 5 + component];
                const SphericalTensor expected =
                    RotateNativeT2(transform, source);
                for (std::size_t component = 0; component < 5; ++component) {
                    ok = ok && std::abs(
                        moved_solvent_efg[atom * 5 + component] -
                        expected.T2[component]) <=
                        kTensorAbs + kTensorRel *
                            std::abs(expected.T2[component]);
                }
            }
            {
                SphericalTensor source;
                for (std::size_t component = 0; component < 5; ++component) {
                    source.T2[component] =
                        original_h5_efg[atom * 5 + component];
                }
                const SphericalTensor expected =
                    RotateNativeT2(transform, source);
                for (std::size_t component = 0; component < 5; ++component) {
                    ok = ok && std::abs(
                        moved_h5_efg[atom * 5 + component] -
                        expected.T2[component]) <=
                        kTensorAbs + kTensorRel *
                            std::abs(expected.T2[component]);
                }
            }
            ok = ok && std::abs(moved_phi[atom] - original_phi[atom])
                           <= kPhiAbs;
            ok = ok && original_sanitizer[atom] == 0u
                    && moved_sanitizer[atom] == original_sanitizer[atom]
                    && moved_clamp[atom] == original_clamp[atom];
        }
    }

    std::fprintf(stderr,
        "APBS directional rerun maxima: vector_abs=%.12g vector_rel=%.12g "
        "tensor_abs=%.12g tensor_rel=%.12g phi_abs=%.12g\n",
        max_vector_error, max_vector_relative,
        max_tensor_error, max_tensor_relative, max_phi_error);
    {
        std::ofstream metrics(metrics_path);
        metrics << "APBS directional rerun maxima: vector_abs="
                << max_vector_error << " vector_rel=" << max_vector_relative
                << " tensor_abs=" << max_tensor_error
                << " tensor_rel=" << max_tensor_relative
                << " phi_abs=" << max_phi_error;
    }
    for (const char* subdir : {"original", "proper", "improper"}) {
        const auto path = output_root / subdir;
        RemoveApbsNpys(path);
        RemoveCoulombNpys(path);
        std::filesystem::remove(path / "apbs_directional.h5");
        std::filesystem::remove(path);
    }
    std::filesystem::remove(output_root);
    std::filesystem::remove(config_path);
    _exit(ok ? 0 : 1);
}

void RunAnalyticSignAndClampForcingFunctionInChild() {
    constexpr double kClampVPerA = 100.0;
    constexpr double kChargeE = 500.0;
    constexpr double kDistanceA = 3.5;
    constexpr double kDielectric = 4.0;
    const auto nonce = std::chrono::steady_clock::now()
                           .time_since_epoch().count();
    const auto config_path = std::filesystem::temp_directory_path() /
        ("apbs_analytic_sign_clamp_" + std::to_string(::getpid()) +
         "_" + std::to_string(nonce) + ".toml");
    const auto output_dir = std::filesystem::temp_directory_path() /
        ("apbs_analytic_sign_clamp_" + std::to_string(::getpid()) +
         "_" + std::to_string(nonce));
    WriteApbsTestConfig(config_path, kDielectric, kDielectric, 0.0,
                        kClampVPerA);
    CalculatorConfig::Load(config_path.string());

    auto protein = BuildSingleChargeProbeProtein();
    bool ok = protein != nullptr;
    if (ok) {
        auto& conf = protein->Conformation();
        ok = AttachSingleCharge(conf, kChargeE);
        auto apbs = ok ? ApbsFieldResult::Compute(conf) : nullptr;
        ok = ok && apbs != nullptr;
        if (ok) {
            constexpr size_t kProbe = 1;
            const auto& atom = conf.AtomAt(kProbe);

            // Independent continuum truth for one +q source in a homogeneous
            // dielectric: phi=q*k/(eps*r), E points away from the charge,
            // and Hessian(phi)=q*k*(3 rr^T-r^2 I)/(eps*r^5).
            // Canonical APBS is total(eps=4) - reference(eps=1), so every
            // canonical sign is opposite the raw total diagnostic here.
            const double expected_total_E =
                COULOMB_KE * kChargeE /
                (kDielectric * kDistanceA * kDistanceA);
            const double expected_reaction_E =
                expected_total_E - COULOMB_KE * kChargeE /
                (kDistanceA * kDistanceA);
            const double expected_reaction_phi =
                COULOMB_KE * kChargeE / kDistanceA *
                (1.0 / kDielectric - 1.0);
            const double expected_total_xx =
                2.0 * COULOMB_KE * kChargeE /
                (kDielectric * kDistanceA * kDistanceA * kDistanceA);
            const double expected_reaction_xx =
                expected_total_xx - 2.0 * COULOMB_KE * kChargeE /
                (kDistanceA * kDistanceA * kDistanceA);

            const auto relative_close = [](double actual, double expected) {
                return std::abs(actual - expected) <=
                    0.25 * std::abs(expected);
            };
            std::fprintf(stderr,
                "APBS analytic probe: phi=%g E=(%g,%g,%g) |E|=%g "
                "scale=%g mask=%u Etotal=(%g,%g,%g) "
                "EFGdiag=(%g,%g,%g) EFGtotaldiag=(%g,%g,%g)\n",
                atom.apbs_phi,
                atom.apbs_efield.x(), atom.apbs_efield.y(),
                atom.apbs_efield.z(), atom.apbs_efield.norm(),
                atom.apbs_efield_clamp_scale,
                static_cast<unsigned>(atom.apbs_efield_clamp_mask),
                atom.apbs_efield_total_diagnostic.x(),
                atom.apbs_efield_total_diagnostic.y(),
                atom.apbs_efield_total_diagnostic.z(),
                atom.apbs_efg(0, 0), atom.apbs_efg(1, 1),
                atom.apbs_efg(2, 2),
                atom.apbs_efg_total_diagnostic(0, 0),
                atom.apbs_efg_total_diagnostic(1, 1),
                atom.apbs_efg_total_diagnostic(2, 2));
            ok = ok
                && atom.apbs_phi < 0.0
                && relative_close(atom.apbs_phi, expected_reaction_phi)
                && atom.apbs_efield.x() < 0.0
                && std::abs(atom.apbs_efield.norm() - kClampVPerA) < 1e-10
                && std::abs(atom.apbs_efield.y()) < 1e-5
                && std::abs(atom.apbs_efield.z()) < 1e-5
                && atom.apbs_efield_clamp_mask == 1u
                && relative_close(atom.apbs_efield_clamp_scale,
                                  kClampVPerA /
                                      std::abs(expected_reaction_E))
                && atom.apbs_efield_total_diagnostic.x() > 0.0
                && relative_close(atom.apbs_efield_total_diagnostic.x(),
                                  expected_total_E)
                && atom.apbs_efg(0, 0) < 0.0
                && atom.apbs_efg(1, 1) > 0.0
                && atom.apbs_efg(2, 2) > 0.0
                && relative_close(atom.apbs_efg(0, 0),
                                  expected_reaction_xx)
                && atom.apbs_efg_total_diagnostic(0, 0) > 0.0
                && atom.apbs_efg_total_diagnostic(1, 1) < 0.0
                && atom.apbs_efg_total_diagnostic(2, 2) < 0.0
                && relative_close(
                    atom.apbs_efg_total_diagnostic(0, 0),
                    expected_total_xx);

            ok = ok && std::filesystem::create_directory(output_dir);
            ok = ok && apbs->WriteFeatures(conf, output_dir.string()) == 8;
            const auto emitted_E =
                ReadFloat64NpyPayload(output_dir / "apbs_E.npy");
            const auto emitted_efg =
                ReadFloat64NpyPayload(output_dir / "apbs_efg.npy");
            const auto emitted_phi =
                ReadFloat64NpyPayload(output_dir / "apbs_phi.npy");
            const auto emitted_total_E = ReadFloat64NpyPayload(
                output_dir / "apbs_E_total_diagnostic.npy");
            const auto emitted_total_efg = ReadFloat64NpyPayload(
                output_dir / "apbs_efg_total_diagnostic.npy");
            const auto emitted_clamp_mask = ReadUInt8NpyPayload(
                output_dir / "apbs_E_clamp_mask.npy");
            const auto emitted_clamp_scale = ReadFloat64NpyPayload(
                output_dir / "apbs_E_clamp_scale.npy");
            ok = ok && emitted_E.size() == 6
                    && emitted_efg.size() == 10
                    && emitted_phi.size() == 2
                    && emitted_total_E.size() == 6
                    && emitted_total_efg.size() == 10
                    && emitted_clamp_mask.size() == 2
                    && emitted_clamp_scale.size() == 2;
            if (ok) {
                const size_t e = kProbe * 3;
                const size_t t2 = kProbe * 5;
                // First pin the serialization mapping exactly against every
                // independently computed in-memory audit value.  This catches
                // row swaps as well as an accidentally constant mask/scale.
                for (size_t i = 0; i < conf.AtomCount(); ++i) {
                    ok = ok
                        && emitted_clamp_mask[i] ==
                            conf.AtomAt(i).apbs_efield_clamp_mask
                        && emitted_clamp_scale[i] ==
                            conf.AtomAt(i).apbs_efield_clamp_scale;
                }
                // These are direct pins on the frozen payloads, derived from
                // the continuum configuration above rather than from the
                // in-memory writer source fields.
                ok = ok
                    && emitted_E[e] < 0.0
                    && std::abs(emitted_E[e] + kClampVPerA) < 1e-10
                    && emitted_phi[kProbe] < 0.0
                    && emitted_efg[t2 + 2] > 0.0
                    && emitted_efg[t2 + 4] < 0.0
                    && emitted_total_E[e] > 0.0
                    && emitted_total_efg[t2 + 2] < 0.0
                    && emitted_total_efg[t2 + 4] > 0.0
                    && emitted_clamp_mask[kProbe] == 1u
                    && emitted_clamp_scale[kProbe] > 0.0
                    && emitted_clamp_scale[kProbe] < 1.0
                    && relative_close(
                        emitted_clamp_scale[kProbe],
                        kClampVPerA / std::abs(expected_reaction_E));
            }
        }
    }

    RemoveApbsNpys(output_dir);
    std::filesystem::remove(output_dir);
    std::filesystem::remove(config_path);
    _exit(ok ? 0 : 1);
}

}  // namespace


TEST(ApbsFieldResultStandalone, ThermalVoltageUsesConfiguredTemperature) {
    EXPECT_DOUBLE_EQ(ApbsThermalVoltage(298.15), KT_OVER_E_298K);
    EXPECT_DOUBLE_EQ(ApbsThermalVoltage(310.0),
                     KT_OVER_E_298K * 310.0 / 298.15);
}

TEST(ApbsFieldResultStandalone, HomogeneousReferenceCancelsCanonicalField) {
    EXPECT_EXIT({ RunHomogeneousReferenceForcingFunctionInChild(); },
                ::testing::ExitedWithCode(0), "");
}

TEST(ApbsFieldResultStandalone,
     AnalyticSingleChargePinsReactionSignAndLiteral100VPerAngstromClamp) {
    EXPECT_EXIT({ RunAnalyticSignAndClampForcingFunctionInChild(); },
                ::testing::ExitedWithCode(0), "");
}

TEST(ApbsFieldResultStandalone,
     DirectionalFieldsRerunProductionSolveOnO3TransformedInput) {
    const auto metrics_path = std::filesystem::temp_directory_path() /
        ("apbs_directional_metrics_" + std::to_string(::getpid()) + ".txt");
    std::filesystem::remove(metrics_path);
    ASSERT_EXIT(
        { RunDirectionalCovarianceForcingFunctionInChild(metrics_path); },
        ::testing::ExitedWithCode(0), "");
    std::ifstream metrics(metrics_path);
    ASSERT_TRUE(metrics.is_open());
    const std::string line((std::istreambuf_iterator<char>(metrics)),
                           std::istreambuf_iterator<char>());
    ASSERT_FALSE(line.empty());
    std::cerr << line << "\n";
    std::filesystem::remove(metrics_path);
}


class ApbsFieldResultTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!std::filesystem::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ PDB not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << r.error;
        protein = std::move(r.protein);
    }
    std::unique_ptr<Protein> protein;
};

TEST_F(ApbsFieldResultTest, ComputeSucceeds) {
    auto& conf = protein->Conformation();

    // Attach ChargeAssignmentResult (dependency)
    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    ASSERT_NE(charges, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(charges)));

    // Compute and attach ApbsFieldResult
    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(apbs)));

    const auto& result = conf.Result<ApbsFieldResult>();
    const auto& dims = result.GridDims();
    const auto& lengths = result.GridLengthsA();
    const auto& origin = result.GridOriginA();
    const auto& spacing = result.GridSpacingA();
    Vec3 bbox_min = conf.PositionAt(0);
    Vec3 bbox_max = bbox_min;
    for (size_t i = 1; i < conf.AtomCount(); ++i) {
        bbox_min = bbox_min.cwiseMin(conf.PositionAt(i));
        bbox_max = bbox_max.cwiseMax(conf.PositionAt(i));
    }
    const Vec3 extent = bbox_max - bbox_min;
    const auto configured_dim = static_cast<std::uint64_t>(
        CalculatorConfig::Get("apbs_grid_dim"));
    for (size_t d = 0; d < 3; ++d) {
        EXPECT_EQ(dims[d], configured_dim);
        const double expected_length = std::max(
            extent(static_cast<Eigen::Index>(d)) +
                CalculatorConfig::Get("apbs_manual_grid_padding_A"),
            CalculatorConfig::Get("apbs_manual_grid_min_dim_A"));
        EXPECT_DOUBLE_EQ(lengths[d], expected_length);
        EXPECT_TRUE(std::isfinite(origin[d]));
        EXPECT_NEAR(spacing[d] * static_cast<double>(dims[d] - 1),
                    lengths[d], 1e-12);
    }
    EXPECT_DOUBLE_EQ(result.ManualGridPaddingA(),
                     CalculatorConfig::Get("apbs_manual_grid_padding_A"));
    EXPECT_DOUBLE_EQ(result.ManualGridMinDimA(),
                     CalculatorConfig::Get("apbs_manual_grid_min_dim_A"));
    EXPECT_DOUBLE_EQ(result.TemperatureK(),
                     CalculatorConfig::Get("apbs_temperature_K"));
    EXPECT_DOUBLE_EQ(result.ThermalVoltageV(),
                     ApbsThermalVoltage(result.TemperatureK()));
    EXPECT_DOUBLE_EQ(result.EfieldClampThreshold(),
                     CalculatorConfig::Get("efield_magnitude_sanity_clamp"));

    // Distinct bit-pattern sentinels make the new audit emission a real
    // field-mapping test: writing the clamp mask or an all-zero substitute
    // cannot satisfy the uint8 payload assertions below.
    ASSERT_GE(conf.AtomCount(), 2u);
    conf.MutableAtomAt(0).apbs_nonfinite_sanitizer_mask = 0x05u;
    conf.MutableAtomAt(1).apbs_nonfinite_sanitizer_mask = 0x0au;

    const auto output_dir = std::filesystem::temp_directory_path() /
        ("apbs_static_schema_" + std::to_string(::getpid()));
    std::filesystem::create_directories(output_dir);
    EXPECT_EQ(result.WriteFeatures(conf, output_dir.string()), 8);
    for (const char* filename : {
            "apbs_E.npy", "apbs_efg.npy", "apbs_phi.npy",
            "apbs_E_clamp_mask.npy", "apbs_E_clamp_scale.npy",
            "apbs_nonfinite_sanitizer_mask.npy",
            "apbs_E_total_diagnostic.npy",
            "apbs_efg_total_diagnostic.npy"}) {
        EXPECT_TRUE(std::filesystem::exists(output_dir / filename))
            << filename;
    }
    std::ifstream mask_file(output_dir / "apbs_E_clamp_mask.npy",
                            std::ios::binary);
    const std::string mask_bytes(
        (std::istreambuf_iterator<char>(mask_file)),
        std::istreambuf_iterator<char>());
    EXPECT_NE(mask_bytes.find("|u1"), std::string::npos);
    mask_file.close();
    std::ifstream sanitizer_file(
        output_dir / "apbs_nonfinite_sanitizer_mask.npy",
        std::ios::binary);
    const std::string sanitizer_bytes(
        (std::istreambuf_iterator<char>(sanitizer_file)),
        std::istreambuf_iterator<char>());
    EXPECT_NE(sanitizer_bytes.find("|u1"), std::string::npos);
    ASSERT_GE(sanitizer_bytes.size(), conf.AtomCount());
    const size_t sanitizer_payload =
        sanitizer_bytes.size() - conf.AtomCount();
    EXPECT_EQ(static_cast<unsigned char>(
                  sanitizer_bytes[sanitizer_payload + 0]), 0x05u);
    EXPECT_EQ(static_cast<unsigned char>(
                  sanitizer_bytes[sanitizer_payload + 1]), 0x0au);
    sanitizer_file.close();
    for (const char* filename : {
            "apbs_E.npy", "apbs_efg.npy", "apbs_phi.npy",
            "apbs_E_clamp_mask.npy", "apbs_E_clamp_scale.npy",
            "apbs_nonfinite_sanitizer_mask.npy",
            "apbs_E_total_diagnostic.npy",
            "apbs_efg_total_diagnostic.npy"}) {
        std::filesystem::remove(output_dir / filename);
    }
    std::filesystem::remove(output_dir);
}

TEST_F(ApbsFieldResultTest, EFieldNonZeroForSomeAtoms) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    const auto& result = conf.Result<ApbsFieldResult>();

    int nonzero = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Vec3 E = result.ElectricFieldAt(ai);
        if (E.norm() > 1e-10) nonzero++;
    }
    // A normal heterogeneous PB solve must retain a nonzero reaction field
    // for a substantial fraction of atoms.
    EXPECT_GT(nonzero, static_cast<int>(conf.AtomCount()) / 2);
}

TEST_F(ApbsFieldResultTest, NoNanOrInf) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    const auto& result = conf.Result<ApbsFieldResult>();

    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        Vec3 E = result.ElectricFieldAt(ai);
        Mat3 EFG = result.FieldGradientAt(ai);
        const auto& atom = conf.AtomAt(ai);

        for (int d = 0; d < 3; ++d) {
            EXPECT_FALSE(std::isnan(E(d))) << "NaN in E-field, atom " << ai;
            EXPECT_FALSE(std::isinf(E(d))) << "Inf in E-field, atom " << ai;
        }

        for (int a = 0; a < 3; ++a) {
            for (int b = 0; b < 3; ++b) {
                EXPECT_FALSE(std::isnan(EFG(a,b)))
                    << "NaN in EFG, atom " << ai << " (" << a << "," << b << ")";
                EXPECT_FALSE(std::isinf(EFG(a,b)))
                    << "Inf in EFG, atom " << ai << " (" << a << "," << b << ")";
            }
        }
        EXPECT_TRUE(std::isfinite(atom.apbs_phi));
        EXPECT_TRUE(atom.apbs_efield_total_diagnostic.allFinite());
        EXPECT_TRUE(atom.apbs_efg_total_diagnostic.allFinite());
        EXPECT_TRUE(atom.apbs_efield_clamp_mask == 0u ||
                    atom.apbs_efield_clamp_mask == 1u);
        if (atom.apbs_efield_clamp_mask == 0u)
            EXPECT_DOUBLE_EQ(atom.apbs_efield_clamp_scale, 1.0);
        else {
            EXPECT_GT(atom.apbs_efield_clamp_scale, 0.0);
            EXPECT_LT(atom.apbs_efield_clamp_scale, 1.0);
        }
    }
}

TEST_F(ApbsFieldResultTest, EFGDecomposesCorrectly) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    const auto& result = conf.Result<ApbsFieldResult>();

    // SphericalTensor roundtrip: decompose and reconstruct should match
    // Test on a few atoms
    for (size_t ai = 0; ai < std::min(conf.AtomCount(), size_t(20)); ++ai) {
        Mat3 EFG = result.FieldGradientAt(ai);
        SphericalTensor st = result.FieldGradientSphericalAt(ai);

        // Reconstruct from spherical
        Mat3 reconstructed = st.Reconstruct();

        // Should match within numerical tolerance
        double diff = (EFG - reconstructed).norm();
        EXPECT_LT(diff, 1e-10)
            << "SphericalTensor roundtrip failed at atom " << ai
            << " diff=" << diff;
    }
}

TEST_F(ApbsFieldResultTest, EFGIsTraceless) {
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    const auto& result = conf.Result<ApbsFieldResult>();

    // EFG tensor should be traceless (trace projected out)
    for (size_t ai = 0; ai < std::min(conf.AtomCount(), size_t(50)); ++ai) {
        Mat3 EFG = result.FieldGradientAt(ai);
        double trace = EFG.trace();
        EXPECT_NEAR(trace, 0.0, 1e-10)
            << "Non-zero trace at atom " << ai << " trace=" << trace;
    }
}

TEST_F(ApbsFieldResultTest, FullTensorStored) {
    // Both Mat3 AND SphericalTensor must be stored on ConformationAtom
    auto& conf = protein->Conformation();

    auto charges = ChargeAssignmentResult::Compute(conf, nmr::test::TestEnvironment::Ff14sbParams().c_str());
    conf.AttachResult(std::move(charges));
    auto apbs = ApbsFieldResult::Compute(conf);
    conf.AttachResult(std::move(apbs));

    // Check that ConformationAtom fields are populated
    for (size_t ai = 0; ai < std::min(conf.AtomCount(), size_t(10)); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        EXPECT_EQ(ca.apbs_nonfinite_sanitizer_mask, 0u);

        // At least one of these should be nonzero in the normal PB solve.
        bool efield_ok = ca.apbs_efield.norm() > 1e-15;

        if (efield_ok) {
            // SphericalTensor should also be set
            SphericalTensor st = ca.apbs_efg_spherical;
            Mat3 roundtrip = st.Reconstruct();
            double diff = (ca.apbs_efg - roundtrip).norm();
            EXPECT_LT(diff, 1e-10);

            const Mat3 total_roundtrip =
                ca.apbs_efg_total_diagnostic_spherical.Reconstruct();
            EXPECT_LT((ca.apbs_efg_total_diagnostic - total_roundtrip).norm(),
                      1e-10);
        }
    }
}

TEST_F(ApbsFieldResultTest, DependencyEnforced) {
    auto& conf = protein->Conformation();

    // Try to attach ApbsFieldResult WITHOUT ChargeAssignmentResult
    auto apbs = ApbsFieldResult::Compute(conf);
    ASSERT_NE(apbs, nullptr);

    // AttachResult should fail because ChargeAssignmentResult is missing
    EXPECT_FALSE(conf.AttachResult(std::move(apbs)));
}
