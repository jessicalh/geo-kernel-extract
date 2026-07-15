//
// test_gromacs_energy_time_series: discipline + integration for
// GromacsEnergyTimeSeriesTrajectoryResult. System-scalar per-frame TR;
// breaks the per-atom DenseBuffer pattern that the BS/HM/Sasa TS TRs
// follow. Storage is a flat std::vector<GromacsEnergy> indexed by frame.
//

#include "CalculatorConfig.h"
#include "BondedEnergyResult.h"
#include "BondedEnergyTimeSeriesTrajectoryResult.h"
#include "FullSystemReader.h"
#include "GeometryResult.h"
#include "GromacsEnergyResult.h"
#include "GromacsEnergyTimeSeriesTrajectoryResult.h"
#include "OperationLog.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "Types.h"
#include "DirectionalTestHelpers.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif
#ifndef NMR_TEST_PYTHON_EXECUTABLE
#error "NMR_TEST_PYTHON_EXECUTABLE must be defined"
#endif
#ifndef NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT
#error "NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT must be defined"
#endif

namespace {

constexpr const char* kFixtureProtein = "1P9J_5801";

std::string TrrPathFor(const std::string& p) {
    return fs::path(p).replace_extension(".trr").string();
}
std::string ProductionDirFor(const std::string& p) {
    return fs::path(p).parent_path().string();
}
bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() && fs::exists(fix.tpr_path)
        && fs::exists(TrrPathFor(fix.tpr_path)) && fs::exists(fix.edr_path);
}

std::unique_ptr<nmr::Protein> MakeDirectionalProtein(
        const std::vector<nmr::Vec3>& positions) {
    auto protein = std::make_unique<nmr::Protein>();
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::ALA;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const size_t residue_index = protein->AddResidue(std::move(residue));
    auto atom = nmr::Atom::Create(nmr::Element::C);
    atom->residue_index = residue_index;
    const size_t atom_index = protein->AddAtom(std::move(atom));
    protein->MutableResidueAt(residue_index).atom_indices.push_back(atom_index);
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "GROMACS external tensor source");
    return protein;
}

std::unique_ptr<nmr::Protein> MakeMixedImproperProtein(
        const std::vector<nmr::Vec3>& positions) {
    auto protein = std::make_unique<nmr::Protein>();
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::ALA;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const size_t residue_index = protein->AddResidue(std::move(residue));
    for (std::size_t i = 0; i < positions.size(); ++i) {
        auto atom = nmr::Atom::Create(nmr::Element::C);
        atom->residue_index = residue_index;
        const size_t atom_index = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(residue_index).atom_indices.push_back(
            atom_index);
    }
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "mixed GROMACS improper fixture");
    return protein;
}

double Sum(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0);
}

// GromacsEnergy::vir/pres and both serialized forms are explicitly row-major:
// XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ.  Do not use Eigen::Map here: Mat3's default
// column-major storage would transpose a genuinely nonsymmetric source and
// let a test accidentally bless the wrong external-source layout.
nmr::Mat3 LoadRowMajor(const double* values) {
    nmr::Mat3 matrix;
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            matrix(row, col) = values[3 * row + col];
        }
    }
    return matrix;
}

void StoreRowMajor(const nmr::Mat3& matrix, double* values) {
    for (int row = 0; row < 3; ++row) {
        for (int col = 0; col < 3; ++col) {
            values[3 * row + col] = matrix(row, col);
        }
    }
}

struct DoubleNpy {
    std::vector<std::size_t> shape;
    std::vector<double> values;
};

std::string Trim(std::string value) {
    const auto is_space = [](unsigned char c) {
        return std::isspace(c) != 0;
    };
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(),
                             [&](char c) { return !is_space(c); }));
    value.erase(std::find_if(value.rbegin(), value.rend(),
                             [&](char c) { return !is_space(c); }).base(),
                value.end());
    return value;
}

DoubleNpy ReadFloat64Npy(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    EXPECT_TRUE(input.is_open()) << path;
    DoubleNpy array;
    if (!input.is_open()) return array;

    char magic[6] = {};
    input.read(magic, sizeof(magic));
    EXPECT_EQ(std::string(magic, sizeof(magic)),
              std::string("\x93NUMPY", sizeof(magic))) << path;
    std::uint8_t major = 0;
    std::uint8_t minor = 0;
    input.read(reinterpret_cast<char*>(&major), sizeof(major));
    input.read(reinterpret_cast<char*>(&minor), sizeof(minor));
    EXPECT_EQ(major, 1u) << path;
    EXPECT_EQ(minor, 0u) << path;
    std::uint16_t header_length = 0;
    input.read(reinterpret_cast<char*>(&header_length), sizeof(header_length));
    std::string header(header_length, '\0');
    input.read(header.data(), static_cast<std::streamsize>(header.size()));
    EXPECT_NE(header.find("'descr': '<f8'"), std::string::npos) << header;
    EXPECT_NE(header.find("'fortran_order': False"), std::string::npos)
        << header;

    const auto shape_key = header.find("'shape': (");
    if (shape_key == std::string::npos) {
        ADD_FAILURE() << header;
        return array;
    }
    const auto shape_begin = header.find('(', shape_key);
    const auto shape_end = header.find(')', shape_begin);
    if (shape_begin == std::string::npos || shape_end == std::string::npos) {
        ADD_FAILURE() << header;
        return array;
    }
    std::stringstream shape_stream(
        header.substr(shape_begin + 1, shape_end - shape_begin - 1));
    std::string token;
    std::size_t count = 1;
    while (std::getline(shape_stream, token, ',')) {
        token = Trim(token);
        if (!token.empty()) {
            const auto dimension =
                static_cast<std::size_t>(std::stoull(token));
            array.shape.push_back(dimension);
            count *= dimension;
        }
    }
    array.values.resize(count);
    input.read(reinterpret_cast<char*>(array.values.data()),
               static_cast<std::streamsize>(count * sizeof(double)));
    EXPECT_EQ(input.gcount(),
              static_cast<std::streamsize>(count * sizeof(double))) << path;
    return array;
}

fs::path FreshDirectionalDirectory(const std::string& suffix) {
    const fs::path path = fs::temp_directory_path() /
        ("gromacs_directional_" + std::to_string(::getpid()) + "_" + suffix);
    std::error_code ec;
    fs::remove(path / "gromacs_energy.npy", ec);
    ec.clear();
    fs::remove(path / "gromacs_energy.h5", ec);
    ec.clear();
    fs::remove(path, ec);
    ec.clear();
    EXPECT_TRUE(fs::create_directories(path, ec)) << ec.message();
    return path;
}

void RemoveDirectionalDirectory(const fs::path& path) {
    std::error_code ec;
    EXPECT_TRUE(fs::remove(path / "gromacs_energy.npy", ec)) << ec.message();
    ec.clear();
    EXPECT_TRUE(fs::remove(path / "gromacs_energy.h5", ec)) << ec.message();
    ec.clear();
    EXPECT_TRUE(fs::remove(path, ec)) << ec.message();
}

std::vector<double> ExpectedNpyRow(const nmr::GromacsEnergy& energy) {
    return {
        energy.coulomb_sr, energy.coulomb_recip, energy.coulomb_14,
        energy.bond, energy.angle, energy.urey_bradley,
        energy.proper_dih, energy.improper_dih,
        energy.periodic_improper_dih, energy.cmap_dih,
        energy.lj_sr, energy.lj_14, energy.disper_corr,
        energy.potential, energy.kinetic, energy.total_energy,
        energy.enthalpy, energy.temperature, energy.pressure,
        energy.volume, energy.density,
        energy.box_x, energy.box_y, energy.box_z,
        energy.vir[0], energy.vir[1], energy.vir[2],
        energy.vir[3], energy.vir[4], energy.vir[5],
        energy.vir[6], energy.vir[7], energy.vir[8],
        energy.pres[0], energy.pres[1], energy.pres[2],
        energy.pres[3], energy.pres[4], energy.pres[5],
        energy.pres[6], energy.pres[7], energy.pres[8],
        energy.T_protein, energy.T_non_protein,
    };
}
}  // namespace


TEST(GromacsImproperSources,
     MixedFixtureKeepsExternalAndPerAtomIdentitiesDistinct) {
    const fs::path fixture_dir = fs::path(NMR_TEST_DATA_DIR)
        .parent_path().parent_path() /
        "data/test_fixtures/gromacs_mixed_improper";
    const fs::path tpr_path = fixture_dir / "mixed.tpr";
    const fs::path edr_path = fixture_dir / "mixed.edr";
    ASSERT_TRUE(fs::exists(tpr_path)) << tpr_path;
    ASSERT_TRUE(fs::exists(edr_path)) << edr_path;

    // Raw `gmx energy` values from this checked-in EDR. These constants
    // are the external-tool oracle; the production evaluator is not used
    // to derive them.
    constexpr double kGromacsHarmonicImproper = 85.660431;
    constexpr double kGromacsPeriodicImproper = 13.136967;
    constexpr double kComparisonTolerance = 1.0e-3;

    nmr::Trajectory trajectory("", tpr_path, edr_path);
    ASSERT_TRUE(trajectory.HasEdr());
    const nmr::GromacsEnergy* external = trajectory.EnergyAtTime(0.0);
    ASSERT_NE(external, nullptr);
    ASSERT_TRUE(std::isfinite(external->improper_dih));
    ASSERT_TRUE(std::isfinite(external->periodic_improper_dih));
    EXPECT_NEAR(external->improper_dih,
                kGromacsHarmonicImproper, 5.0e-6);
    EXPECT_NEAR(external->periodic_improper_dih,
                kGromacsPeriodicImproper, 5.0e-6);
    EXPECT_NE(external->improper_dih, external->periodic_improper_dih);

    nmr::FullSystemReader reader;
    ASSERT_TRUE(reader.ReadTopology(tpr_path.string())) << reader.error();
    const nmr::BondedParameters& imported = reader.BondedParams();
    nmr::BondedParameters harmonic_only;
    nmr::BondedParameters periodic_only;
    nmr::BondedParameters both;
    for (const auto& interaction : imported.interactions) {
        if (interaction.type == nmr::BondedInteraction::ImproperDih) {
            harmonic_only.interactions.push_back(interaction);
            both.interactions.push_back(interaction);
        } else if (interaction.type ==
                   nmr::BondedInteraction::PeriodicImproperDih) {
            periodic_only.interactions.push_back(interaction);
            both.interactions.push_back(interaction);
        }
    }
    ASSERT_EQ(imported.interactions.size(), 2u);
    ASSERT_EQ(harmonic_only.interactions.size(), 1u);
    ASSERT_EQ(periodic_only.interactions.size(), 1u);
    ASSERT_EQ(both.interactions.size(), 2u);

    const std::vector<nmr::Vec3> positions = {
        nmr::Vec3(1.0, 2.0, 3.0),
        nmr::Vec3(2.2, 1.5, 2.6),
        nmr::Vec3(3.1, 2.7, 1.9),
        nmr::Vec3(4.3, 2.1, 3.6),
    };
    auto protein = MakeMixedImproperProtein(positions);
    nmr::ProteinConformation& conf = protein->Conformation();

    auto harmonic_result =
        nmr::BondedEnergyResult::Compute(conf, harmonic_only);
    auto periodic_result =
        nmr::BondedEnergyResult::Compute(conf, periodic_only);
    ASSERT_NE(harmonic_result, nullptr);
    ASSERT_NE(periodic_result, nullptr);

    EXPECT_NEAR(Sum(harmonic_result->ImproperDihEnergy()),
                kGromacsHarmonicImproper, kComparisonTolerance);
    EXPECT_DOUBLE_EQ(Sum(harmonic_result->PeriodicImproperDihEnergy()), 0.0);
    EXPECT_DOUBLE_EQ(Sum(periodic_result->ImproperDihEnergy()), 0.0);
    EXPECT_NEAR(Sum(periodic_result->PeriodicImproperDihEnergy()),
                kGromacsPeriodicImproper, kComparisonTolerance);

    auto bonded_result = nmr::BondedEnergyResult::Compute(conf, both);
    ASSERT_NE(bonded_result, nullptr);
    EXPECT_NEAR(Sum(bonded_result->ImproperDihEnergy()),
                kGromacsHarmonicImproper, kComparisonTolerance);
    EXPECT_NEAR(Sum(bonded_result->PeriodicImproperDihEnergy()),
                kGromacsPeriodicImproper, kComparisonTolerance);
    EXPECT_NEAR(Sum(bonded_result->TotalBonded()),
                kGromacsHarmonicImproper + kGromacsPeriodicImproper,
                2.0 * kComparisonTolerance);

    const fs::path output_dir =
        FreshDirectionalDirectory("mixed_improper_sources");
    auto gromacs_result =
        nmr::GromacsEnergyResult::Compute(conf, *external);
    ASSERT_NE(gromacs_result, nullptr);
    ASSERT_EQ(gromacs_result->WriteFeatures(conf, output_dir.string()), 1);
    ASSERT_EQ(bonded_result->WriteFeatures(conf, output_dir.string()), 1);

    const DoubleNpy gromacs_npy =
        ReadFloat64Npy(output_dir / "gromacs_energy.npy");
    ASSERT_EQ(gromacs_npy.shape,
              (std::vector<std::size_t>{1u, 44u}));
    EXPECT_NEAR(gromacs_npy.values[7],
                kGromacsHarmonicImproper, 5.0e-6);
    EXPECT_NEAR(gromacs_npy.values[8],
                kGromacsPeriodicImproper, 5.0e-6);

    const DoubleNpy bonded_npy =
        ReadFloat64Npy(output_dir / "bonded_energy.npy");
    ASSERT_EQ(bonded_npy.shape,
              (std::vector<std::size_t>{positions.size(), 8u}));
    double npy_harmonic = 0.0;
    double npy_periodic = 0.0;
    double npy_total = 0.0;
    for (std::size_t atom = 0; atom < positions.size(); ++atom) {
        npy_harmonic += bonded_npy.values[atom * 8 + 4];
        npy_periodic += bonded_npy.values[atom * 8 + 5];
        npy_total += bonded_npy.values[atom * 8 + 7];
    }
    EXPECT_NEAR(npy_harmonic,
                kGromacsHarmonicImproper, kComparisonTolerance);
    EXPECT_NEAR(npy_periodic,
                kGromacsPeriodicImproper, kComparisonTolerance);
    EXPECT_NEAR(npy_total,
                kGromacsHarmonicImproper + kGromacsPeriodicImproper,
                2.0 * kComparisonTolerance);

    ASSERT_TRUE(conf.AttachResult(std::move(gromacs_result)));
    ASSERT_TRUE(conf.AttachResult(std::move(bonded_result)));
    auto tp = nmr::TrajectoryProtein::CreateForTesting(std::move(protein));
    ASSERT_NE(tp, nullptr);
    auto gromacs_time_series =
        nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(*tp);
    auto bonded_time_series =
        nmr::BondedEnergyTimeSeriesTrajectoryResult::Create(*tp);
    ASSERT_NE(gromacs_time_series, nullptr);
    ASSERT_NE(bonded_time_series, nullptr);
    gromacs_time_series->Compute(conf, *tp, trajectory, 0, 0.0);
    bonded_time_series->Compute(conf, *tp, trajectory, 0, 0.0);
    gromacs_time_series->Finalize(*tp, trajectory);
    bonded_time_series->Finalize(*tp, trajectory);

    const fs::path h5_path = output_dir / "gromacs_energy.h5";
    {
        HighFive::File file(h5_path.string(), HighFive::File::Truncate);
        gromacs_time_series->WriteH5Group(*tp, file);
        bonded_time_series->WriteH5Group(*tp, file);
    }
    {
        HighFive::File file(h5_path.string(), HighFive::File::ReadOnly);
        auto gromacs_group = file.getGroup(
            "/trajectory/gromacs_energy_time_series");
        std::vector<double> external_harmonic;
        std::vector<double> external_periodic;
        gromacs_group.getDataSet("improper_dih").read(external_harmonic);
        gromacs_group.getDataSet("periodic_improper_dih")
            .read(external_periodic);
        ASSERT_EQ(external_harmonic.size(), 1u);
        ASSERT_EQ(external_periodic.size(), 1u);
        EXPECT_NEAR(external_harmonic[0],
                    kGromacsHarmonicImproper, 5.0e-6);
        EXPECT_NEAR(external_periodic[0],
                    kGromacsPeriodicImproper, 5.0e-6);

        auto bonded_group = file.getGroup(
            "/trajectory/bonded_energy_time_series");
        std::vector<std::vector<double>> per_atom_harmonic;
        std::vector<std::vector<double>> per_atom_periodic;
        std::vector<std::vector<double>> per_atom_total;
        bonded_group.getDataSet("improper_dih").read(per_atom_harmonic);
        bonded_group.getDataSet("periodic_improper_dih")
            .read(per_atom_periodic);
        bonded_group.getDataSet("total").read(per_atom_total);
        ASSERT_EQ(per_atom_harmonic.size(), positions.size());
        ASSERT_EQ(per_atom_periodic.size(), positions.size());
        ASSERT_EQ(per_atom_total.size(), positions.size());
        double h5_harmonic = 0.0;
        double h5_periodic = 0.0;
        double h5_total = 0.0;
        for (std::size_t atom = 0; atom < positions.size(); ++atom) {
            ASSERT_EQ(per_atom_harmonic[atom].size(), 1u);
            ASSERT_EQ(per_atom_periodic[atom].size(), 1u);
            ASSERT_EQ(per_atom_total[atom].size(), 1u);
            h5_harmonic += per_atom_harmonic[atom][0];
            h5_periodic += per_atom_periodic[atom][0];
            h5_total += per_atom_total[atom][0];
        }
        EXPECT_NEAR(h5_harmonic,
                    kGromacsHarmonicImproper, kComparisonTolerance);
        EXPECT_NEAR(h5_periodic,
                    kGromacsPeriodicImproper, kComparisonTolerance);
        EXPECT_NEAR(h5_total,
                    kGromacsHarmonicImproper + kGromacsPeriodicImproper,
                    2.0 * kComparisonTolerance);
    }

    std::error_code ec;
    EXPECT_TRUE(fs::remove(output_dir / "bonded_energy.npy", ec))
        << ec.message();
    RemoveDirectionalDirectory(output_dir);
}


TEST(GromacsEnergyDirectionalCovariance,
     ProductionRerunAndExactNpyH5UnderO3) {
    const std::vector<nmr::Vec3> positions{nmr::Vec3(0.31, -1.7, 2.4)};
    auto original_protein = MakeDirectionalProtein(positions);
    nmr::ProteinConformation& original = original_protein->Conformation();

    nmr::GromacsEnergy source;
    source.time_ps = 12.5;
    source.coulomb_sr = -113.0;
    source.coulomb_recip = -29.0;
    source.coulomb_14 = 4.0;
    source.bond = 5.0;
    source.angle = 6.0;
    source.urey_bradley = 7.0;
    source.proper_dih = 8.0;
    source.improper_dih = 9.0;
    source.periodic_improper_dih = 9.5;
    source.cmap_dih = 10.0;
    source.lj_sr = -41.0;
    source.lj_14 = -3.0;
    source.disper_corr = -2.0;
    source.potential = -139.0;
    source.kinetic = 51.0;
    source.total_energy = -88.0;
    source.enthalpy = -87.5;
    source.temperature = 301.25;
    source.pressure = 1.07;
    source.volume = 125.0;
    source.density = 998.0;
    // Only diagonal box lengths are retained by this source schema; an
    // isotropic cell is the O(3)-closed case.  A general triclinic cell
    // cannot be reconstructed from these three scalars.
    source.box_x = source.box_y = source.box_z = 5.0;
    source.T_protein = 300.5;
    source.T_non_protein = 301.7;

    // Deliberately nonsymmetric.  A transposed row/column-major load cannot
    // hide behind the usual symmetric-stress special case.
    nmr::Mat3 virial;
    virial << 8.0, -1.5,  0.7,
              2.4,  5.0,  2.1,
             -0.8,  3.2,  3.0;
    nmr::Mat3 pressure;
    pressure << 1.2,  0.11, -0.07,
               -0.31, 0.9,   0.04,
                0.26, 0.18,  1.1;
    ASSERT_GT((virial - virial.transpose()).norm(), 1.0);
    ASSERT_GT((pressure - pressure.transpose()).norm(), 0.1);
    StoreRowMajor(virial, source.vir);
    StoreRowMajor(pressure, source.pres);

    auto original_result =
        nmr::GromacsEnergyResult::Compute(original, source);
    ASSERT_NE(original_result, nullptr);

    constexpr std::uint64_t kTransformSeed = 0x47524f4d41435345ULL;
    constexpr double kTensorAbsTolerance = 3.0e-13;
    constexpr double kTensorRelTolerance = 3.0e-14;
    nmr::Trajectory dummy("", "", "");
    ASSERT_TRUE(original.AttachResult(std::move(original_result)));
    auto original_tp = nmr::TrajectoryProtein::CreateForTesting(
        std::move(original_protein));
    ASSERT_NE(original_tp, nullptr);
    auto original_trajectory_result =
        nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(*original_tp);
    ASSERT_NE(original_trajectory_result, nullptr);
    original_trajectory_result->Compute(
        original, *original_tp, dummy, 19, 4.5);
    original_trajectory_result->Finalize(*original_tp, dummy);
    const fs::path source_h5_path = fs::temp_directory_path() /
        ("gromacs_directional_" + std::to_string(::getpid()) +
         "_source.h5");
    fs::remove(source_h5_path);
    {
        HighFive::File file(
            source_h5_path.string(), HighFive::File::Truncate);
        original_trajectory_result->WriteH5Group(*original_tp, file);
    }
    std::vector<std::size_t> source_frame_indices;
    std::vector<double> source_frame_times;
    {
        HighFive::File file(
            source_h5_path.string(), HighFive::File::ReadOnly);
        const auto group = file.getGroup(
            "/trajectory/gromacs_energy_time_series");
        group.getDataSet("frame_indices").read(source_frame_indices);
        group.getDataSet("frame_times").read(source_frame_times);
    }
    ASSERT_EQ(source_frame_indices, (std::vector<std::size_t>{19u}));
    ASSERT_EQ(source_frame_times, (std::vector<double>{4.5}));
    for (const bool improper : {false, true}) {
        const auto transform = nmr::test::directional::SeededTransform(
            kTransformSeed, improper);
        const auto moved_positions =
            nmr::test::directional::Positions(transform, positions);
        auto moved_tp = nmr::TrajectoryProtein::CreateForTesting(
            MakeDirectionalProtein(moved_positions));
        ASSERT_NE(moved_tp, nullptr);
        auto moved = moved_tp->TickConformation(moved_positions);
        ASSERT_NE(moved, nullptr);

        nmr::GromacsEnergy moved_source = source;
        const nmr::Mat3 expected_virial =
            nmr::test::directional::EvenRank2(transform, virial);
        const nmr::Mat3 expected_pressure =
            nmr::test::directional::EvenRank2(transform, pressure);
        StoreRowMajor(expected_virial, moved_source.vir);
        StoreRowMajor(expected_pressure, moved_source.pres);

        auto moved_result =
            nmr::GromacsEnergyResult::Compute(*moved, moved_source);
        ASSERT_NE(moved_result, nullptr);
        const auto& emitted = moved_result->Energy();
        const nmr::Mat3 emitted_virial = LoadRowMajor(emitted.vir);
        const nmr::Mat3 emitted_pressure = LoadRowMajor(emitted.pres);
        EXPECT_TRUE(nmr::test::directional::NearMatrix(
            emitted_virial, expected_virial,
            kTensorAbsTolerance, kTensorRelTolerance));
        EXPECT_TRUE(nmr::test::directional::NearMatrix(
            emitted_pressure, expected_pressure,
            kTensorAbsTolerance, kTensorRelTolerance));

        // Representative scalar channels and the isotropic box remain
        // rotation/translation invariant when the source is transformed
        // consistently before the owning Compute() call.
        EXPECT_DOUBLE_EQ(emitted.coulomb_sr, source.coulomb_sr);
        EXPECT_DOUBLE_EQ(emitted.total_energy, source.total_energy);
        EXPECT_DOUBLE_EQ(emitted.temperature, source.temperature);
        EXPECT_DOUBLE_EQ(emitted.volume, source.volume);
        EXPECT_DOUBLE_EQ(emitted.box_x, source.box_x);
        EXPECT_DOUBLE_EQ(emitted.box_y, source.box_y);
        EXPECT_DOUBLE_EQ(emitted.box_z, source.box_z);
        EXPECT_DOUBLE_EQ(emitted.T_protein, source.T_protein);
        EXPECT_DOUBLE_EQ(emitted.T_non_protein, source.T_non_protein);

        // Static serialization boundary: production writer, exact filename,
        // exact (1,44) shape and explicit row-major tensor reconstruction.
        const fs::path output_dir = FreshDirectionalDirectory(
            improper ? "improper" : "proper");
        ASSERT_EQ(moved_result->WriteFeatures(*moved, output_dir.string()), 1);
        ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                      NMR_TEST_PYTHON_EXECUTABLE,
                      NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                      {output_dir / "gromacs_energy.npy"}),
                  0);
        const DoubleNpy npy =
            ReadFloat64Npy(output_dir / "gromacs_energy.npy");
        ASSERT_EQ(npy.shape,
                  (std::vector<std::size_t>{1u, 44u}));
        const auto expected_row = ExpectedNpyRow(moved_source);
        ASSERT_EQ(npy.values.size(), expected_row.size());
        for (std::size_t column = 0; column < expected_row.size(); ++column) {
            SCOPED_TRACE("gromacs_energy.npy column=" +
                         std::to_string(column));
            EXPECT_DOUBLE_EQ(npy.values[column], expected_row[column]);
        }
        EXPECT_TRUE(nmr::test::directional::NearMatrix(
            LoadRowMajor(npy.values.data() + 24), expected_virial,
            kTensorAbsTolerance, kTensorRelTolerance));
        EXPECT_TRUE(nmr::test::directional::NearMatrix(
            LoadRowMajor(npy.values.data() + 33), expected_pressure,
            kTensorAbsTolerance, kTensorRelTolerance));
        // All non-tensor columns, including isotropic Box-X/Y/Z, remain
        // exact invariants under both proper and improper transforms.
        const auto source_row = ExpectedNpyRow(source);
        for (std::size_t column = 0; column < 24; ++column) {
            EXPECT_DOUBLE_EQ(npy.values[column], source_row[column])
                << "gromacs_energy.npy invariant column=" << column;
        }
        EXPECT_DOUBLE_EQ(npy.values[42], source_row[42]);
        EXPECT_DOUBLE_EQ(npy.values[43], source_row[43]);

        // Attach the freshly rerun production result and exercise the real
        // trajectory rollup/H5 writer on the transformed input.
        ASSERT_TRUE(moved->AttachResult(std::move(moved_result)));
        auto trajectory_result =
            nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(*moved_tp);
        ASSERT_NE(trajectory_result, nullptr);
        trajectory_result->Compute(*moved, *moved_tp, dummy, 19, 4.5);
        trajectory_result->Finalize(*moved_tp, dummy);
        const fs::path h5_path = output_dir / "gromacs_energy.h5";
        {
            HighFive::File file(h5_path.string(), HighFive::File::Truncate);
            trajectory_result->WriteH5Group(*moved_tp, file);
        }
        {
            HighFive::File h5(h5_path.string(), HighFive::File::ReadOnly);
            auto group = h5.getGroup(
                "/trajectory/gromacs_energy_time_series");
            std::vector<std::vector<double>> h5_virial;
            std::vector<std::vector<double>> h5_pressure;
            group.getDataSet("virial").read(h5_virial);
            group.getDataSet("pressure_tensor").read(h5_pressure);
            ASSERT_EQ(h5_virial.size(), 1u);
            ASSERT_EQ(h5_pressure.size(), 1u);
            ASSERT_EQ(h5_virial[0].size(), 9u);
            ASSERT_EQ(h5_pressure[0].size(), 9u);
            EXPECT_TRUE(nmr::test::directional::NearMatrix(
                LoadRowMajor(h5_virial[0].data()), expected_virial,
                kTensorAbsTolerance, kTensorRelTolerance));
            EXPECT_TRUE(nmr::test::directional::NearMatrix(
                LoadRowMajor(h5_pressure[0].data()), expected_pressure,
                kTensorAbsTolerance, kTensorRelTolerance));
            std::vector<double> box_x;
            std::vector<double> box_y;
            std::vector<double> box_z;
            group.getDataSet("box_x").read(box_x);
            group.getDataSet("box_y").read(box_y);
            group.getDataSet("box_z").read(box_z);
            ASSERT_EQ(box_x.size(), 1u);
            ASSERT_EQ(box_y.size(), 1u);
            ASSERT_EQ(box_z.size(), 1u);
            EXPECT_DOUBLE_EQ(box_x[0], source.box_x);
            EXPECT_DOUBLE_EQ(box_y[0], source.box_y);
            EXPECT_DOUBLE_EQ(box_z[0], source.box_z);
            for (const char* name : {"box_x", "box_y", "box_z"}) {
                auto dataset = group.getDataSet(name);
                std::string frame;
                std::string parity;
                std::string law;
                dataset.getAttribute("coordinate_frame").read(frame);
                dataset.getAttribute("parity").read(parity);
                dataset.getAttribute("transformation").read(law);
                EXPECT_EQ(frame, "gromacs_simulation_box_axes") << name;
                EXPECT_EQ(parity, "mixed") << name;
                EXPECT_NE(law.find("no closed O(3) law"),
                          std::string::npos) << name;
            }
            for (const char* name : {"virial", "pressure_tensor"}) {
                auto dataset = group.getDataSet(name);
                std::string basis;
                std::string order;
                std::string frame;
                std::string parity;
                std::string law;
                std::string alias_frame;
                std::string alias_parity;
                std::string alias_law;
                std::string structural_zeros;
                std::string e3nn_export;
                dataset.getAttribute("tensor_basis").read(basis);
                dataset.getAttribute("tensor_component_order").read(order);
                dataset.getAttribute("tensor_frame").read(frame);
                dataset.getAttribute("tensor_parity").read(parity);
                dataset.getAttribute("tensor_transformation").read(law);
                dataset.getAttribute("coordinate_frame").read(alias_frame);
                dataset.getAttribute("parity").read(alias_parity);
                dataset.getAttribute("transformation").read(alias_law);
                dataset.getAttribute("structural_zero_components").read(
                    structural_zeros);
                dataset.getAttribute("e3nn_export").read(e3nn_export);
                EXPECT_EQ(basis, "cartesian_matrix_row_major") << name;
                EXPECT_EQ(order, "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ") << name;
                EXPECT_EQ(frame, "gromacs_simulation_cartesian_xyz") << name;
                EXPECT_EQ(parity, "even") << name;
                EXPECT_EQ(law, "even_rank2: T'=R T R^T") << name;
                EXPECT_EQ(alias_frame, frame) << name;
                EXPECT_EQ(alias_parity, parity) << name;
                EXPECT_EQ(alias_law, law) << name;
                EXPECT_EQ(structural_zeros, "none") << name;
                EXPECT_EQ(
                    e3nn_export,
                    "decompose the row-major Cartesian matrix into "
                    "project-native T0/T1/T2 before explicit conversion "
                    "to e3nn") << name;
            }
            const std::vector<std::pair<const char*, double>>
                scalar_invariants = {
                    {"coulomb_sr", source.coulomb_sr},
                    {"coulomb_recip", source.coulomb_recip},
                    {"coulomb_14", source.coulomb_14},
                    {"bond", source.bond},
                    {"angle", source.angle},
                    {"urey_bradley", source.urey_bradley},
                    {"proper_dih", source.proper_dih},
                    {"improper_dih", source.improper_dih},
                    {"periodic_improper_dih",
                     source.periodic_improper_dih},
                    {"cmap_dih", source.cmap_dih},
                    {"lj_sr", source.lj_sr},
                    {"lj_14", source.lj_14},
                    {"disper_corr", source.disper_corr},
                    {"potential", source.potential},
                    {"kinetic", source.kinetic},
                    {"total_energy", source.total_energy},
                    {"enthalpy", source.enthalpy},
                    {"temperature", source.temperature},
                    {"pressure", source.pressure},
                    {"volume", source.volume},
                    {"density", source.density},
                    {"box_x", source.box_x},
                    {"box_y", source.box_y},
                    {"box_z", source.box_z},
                    {"T_protein", source.T_protein},
                    {"T_non_protein", source.T_non_protein},
                    {"energy_frame_times_ps", source.time_ps},
                };
            for (const auto& [name, expected] : scalar_invariants) {
                std::vector<double> values;
                group.getDataSet(name).read(values);
                ASSERT_EQ(values.size(), 1u) << name;
                EXPECT_DOUBLE_EQ(values[0], expected) << name;
            }
            std::vector<std::uint8_t> attached;
            group.getDataSet("source_attached_per_frame").read(attached);
            EXPECT_EQ(attached, (std::vector<std::uint8_t>{1u}));
            std::vector<std::size_t> frame_indices;
            std::vector<double> frame_times;
            group.getDataSet("frame_indices").read(frame_indices);
            group.getDataSet("frame_times").read(frame_times);
            EXPECT_EQ(frame_indices, source_frame_indices);
            EXPECT_EQ(frame_times, source_frame_times);
        }
        RemoveDirectionalDirectory(output_dir);
    }
    EXPECT_TRUE(fs::remove(source_h5_path));
}


// ============================================================================
// SYNTHETIC: attach a fake GromacsEnergyResult on each of three frames and
// confirm the TR rolls up the snapshots in order with the correct values.
// System-scalar TR cannot use the per-atom-field shortcut other TRs do;
// the source is GromacsEnergyResult itself, so we attach it.
// ============================================================================

TEST(GromacsEnergyTimeSeries, SyntheticThreeFrames) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    constexpr std::size_t kFrames = 3;
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");

        // Synthesize a GromacsEnergy snapshot with a recognizable pattern.
        nmr::GromacsEnergy ge;
        ge.time_ps      = static_cast<double>(t);
        ge.coulomb_sr   = -1000.0 - t;
        ge.coulomb_recip = -500.0 - t;
        ge.potential    = -5000.0 - 10.0 * t;
        ge.kinetic      = 2000.0 + t;
        ge.total_energy = ge.potential + ge.kinetic;
        ge.temperature  = 300.0 + 0.1 * t;
        ge.pressure     = 1.0 + 0.01 * t;
        ge.volume       = 100.0 + t;
        ge.density      = 1000.0 - t;
        ge.box_x        = ge.box_y = ge.box_z = std::cbrt(ge.volume);
        for (int k = 0; k < 9; ++k) ge.vir[k]  = static_cast<double>(t * 100 + k);
        for (int k = 0; k < 9; ++k) ge.pres[k] = static_cast<double>(t * 10  + k);

        conf->AttachResult(nmr::GromacsEnergyResult::Compute(*conf, ge));
        tr->Compute(*conf, tp, traj, t, ge.time_ps);
    }
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), kFrames);

    // H5 round-trip the snapshots and read back the channels.
    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_unit_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/gromacs_energy_time_series"));
    auto grp = reopen.getGroup("/trajectory/gromacs_energy_time_series");

    std::vector<double> potential;
    grp.getDataSet("potential").read(potential);
    ASSERT_EQ(potential.size(), kFrames);
    for (std::size_t t = 0; t < kFrames; ++t)
        EXPECT_DOUBLE_EQ(potential[t], -5000.0 - 10.0 * t);

    std::vector<std::vector<double>> vir;
    grp.getDataSet("virial").read(vir);
    ASSERT_EQ(vir.size(), kFrames);
    ASSERT_EQ(vir[0].size(), 9u);
    EXPECT_DOUBLE_EQ(vir[2][3], 2 * 100 + 3);  // (frame=2, k=3)

    std::string units, tensor_layout;
    grp.getAttribute("units").read(units);
    grp.getAttribute("tensor_layout").read(tensor_layout);
    EXPECT_EQ(units, "kJ/mol");
    EXPECT_EQ(tensor_layout, "XX,XY,XZ,YX,YY,YZ,ZX,ZY,ZZ");

    fs::remove(h5_path);
}


// ============================================================================
// SYNTHETIC: source-absent path. Mix 2 attached + 2 absent frames, verify
// the source-attached gate behavior — H5 group lists source_attached_count=2
// and absent frames carry NaN on the energy datasets. R2 review 2026-05-18:
// the gate logic was uncovered until this test landed.
// ============================================================================

TEST(GromacsEnergyTimeSeries, SyntheticSourceAbsentFrames) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    constexpr std::size_t kFrames = 4;
    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < kFrames; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        // Attach GromacsEnergyResult only on EVEN frames (0, 2). Odd
        // frames (1, 3) are source-absent.
        if (t % 2 == 0) {
            nmr::GromacsEnergy ge;
            ge.time_ps     = static_cast<double>(t);
            ge.potential   = -1000.0 - 10.0 * t;
            ge.temperature = 300.0;
            conf->AttachResult(nmr::GromacsEnergyResult::Compute(*conf, ge));
        }
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);

    // H5 round-trip + verify NaN on absent frames.
    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_absent_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/gromacs_energy_time_series");

    std::size_t source_attached_count = 0;
    grp.getAttribute("source_attached_count").read(source_attached_count);
    EXPECT_EQ(source_attached_count, 2u)
        << "Expected 2/4 frames source-attached (even frames only)";

    std::vector<std::uint8_t> mask;
    grp.getDataSet("source_attached_per_frame").read(mask);
    ASSERT_EQ(mask.size(), kFrames);
    EXPECT_EQ(mask[0], 1u); EXPECT_EQ(mask[1], 0u);
    EXPECT_EQ(mask[2], 1u); EXPECT_EQ(mask[3], 0u);

    std::vector<double> potential;
    grp.getDataSet("potential").read(potential);
    ASSERT_EQ(potential.size(), kFrames);
    // Frame 0, 2: finite values (attached). Frame 1, 3: NaN.
    EXPECT_TRUE(std::isfinite(potential[0]));
    EXPECT_TRUE(std::isnan(potential[1]));
    EXPECT_TRUE(std::isfinite(potential[2]));
    EXPECT_TRUE(std::isnan(potential[3]));
    EXPECT_DOUBLE_EQ(potential[0], -1000.0);
    EXPECT_DOUBLE_EQ(potential[2], -1020.0);

    // energy_frame_times_ps NaN-filled too.
    std::vector<double> et;
    grp.getDataSet("energy_frame_times_ps").read(et);
    EXPECT_TRUE(std::isfinite(et[0]));
    EXPECT_TRUE(std::isnan(et[1]));

    fs::remove(h5_path);
}


// ============================================================================
// SYNTHETIC: ALL frames source-absent → WriteH5Group skips emission entirely.
// ============================================================================

TEST(GromacsEnergyTimeSeries, SyntheticAllAbsentSkipsGroup) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    const std::size_t N = tp.AtomCount();
    auto tr = nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp);
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);

    std::vector<nmr::Vec3> positions(N, nmr::Vec3::Zero());
    for (std::size_t t = 0; t < 3; ++t) {
        auto conf = std::make_unique<nmr::ProteinConformation>(
            &tp.ProteinRef(), positions, "synthetic frame");
        // No AttachResult — source is absent for every frame.
        tr->Compute(*conf, tp, traj, t, static_cast<double>(t));
    }
    tr->Finalize(tp, traj);
    EXPECT_EQ(tr->NumFrames(), 3u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_allabsent_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr->WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    EXPECT_FALSE(reopen.exist("/trajectory/gromacs_energy_time_series"))
        << "All-absent run should skip group emission entirely.";

    fs::remove(h5_path);
}


// ============================================================================
// DISCIPLINE: Frame-0 semantics — single-frame fixture exercises Compute once.
// ============================================================================

TEST(GromacsEnergyTimeSeries, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::GromacsEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);
    EXPECT_EQ(traj.FrameCount(), 1u);

    const auto& tr = tp.Result<nmr::GromacsEnergyTimeSeriesTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


// ============================================================================
// DISCIPLINE: Finalize idempotency — second Finalize must not corrupt state.
// ============================================================================

TEST(GromacsEnergyTimeSeries, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::GromacsEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::GromacsEnergyTimeSeriesTrajectoryResult>();
    const std::size_t T_first = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T_first);
}


// ============================================================================
// DISCIPLINE: H5 round-trip — group + all expected datasets land.
// ============================================================================

TEST(GromacsEnergyTimeSeries, H5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::GromacsEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(99999);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::GromacsEnergyTimeSeriesTrajectoryResult>();
    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_h5_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/gromacs_energy_time_series"));
    auto grp = reopen.getGroup("/trajectory/gromacs_energy_time_series");

    std::string units;
    grp.getAttribute("units").read(units);
    EXPECT_EQ(units, "kJ/mol");

    // Expected datasets — sample 4 representative ones across categories.
    EXPECT_TRUE(grp.exist("potential"));
    EXPECT_TRUE(grp.exist("temperature"));
    EXPECT_TRUE(grp.exist("virial"));
    EXPECT_TRUE(grp.exist("frame_times"));

    // Virial / pressure_tensor shape: (T, 9).
    const auto vir_dims = grp.getDataSet("virial").getSpace().getDimensions();
    ASSERT_EQ(vir_dims.size(), 2u);
    EXPECT_EQ(vir_dims[1], 9u);

    fs::remove(h5_path);
}


// ============================================================================
// INTEGRATION: real .edr energies through Trajectory::Run, multi-frame.
// ============================================================================

TEST(GromacsEnergyTimeSeries, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true; opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::GromacsEnergyResult));
    config.AddTrajectoryResultFactory([](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::GromacsEnergyTimeSeriesTrajectoryResult::Create(tp_in);
    });
    config.SetStride(300);

    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path))) << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::GromacsEnergyTimeSeriesTrajectoryResult>();
    EXPECT_GE(tr.NumFrames(), 2u);

    // H5 round-trip to inspect populated channels — confirm finite values
    // on a sample of the physically-meaningful channels.
    const std::string h5_path = (fs::temp_directory_path() /
        ("gromacs_energy_ts_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate); tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/gromacs_energy_time_series");

    std::vector<double> total, temperature, pressure;
    grp.getDataSet("total_energy").read(total);
    grp.getDataSet("temperature").read(temperature);
    grp.getDataSet("pressure").read(pressure);
    ASSERT_EQ(total.size(), tr.NumFrames());

    for (double v : total)       EXPECT_TRUE(std::isfinite(v));
    for (double T : temperature) EXPECT_TRUE(std::isfinite(T));
    for (double P : pressure)    EXPECT_TRUE(std::isfinite(P));

    // Loose physical-sanity: thermostatted MD around 300 K within ±50 K.
    // Per feedback_log_overages_dont_assert this is logged, not asserted.
    for (double T : temperature) {
        if (T < 250.0 || T > 400.0) {
            std::cerr << "WARN: temperature out of [250,400] K: " << T << "\n";
        }
    }
    std::cout << "GromacsEnergyTimeSeries: " << tr.NumFrames() << " frames; "
              << "total_energy[0]=" << total[0]
              << " T[0]=" << temperature[0]
              << " P[0]=" << pressure[0] << "\n";

    fs::remove(h5_path);
}
