#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <filesystem>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <unistd.h>

#include "MutationDeltaResult.h"
#include "OrcaRunLoader.h"
#include "OrcaShieldingResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "GeometryResult.h"
#include "EnrichmentResult.h"
#include "SpatialIndexResult.h"
#include "DsspResult.h"
#include "ApbsFieldResult.h"
#include "MopacResult.h"
#include "MolecularGraphResult.h"
#include "ProtonationDetectionResult.h"
#include "DirectionalTestHelpers.h"

namespace fs = std::filesystem;
using namespace nmr;

#ifndef NMR_TEST_PYTHON_EXECUTABLE
#error "NMR_TEST_PYTHON_EXECUTABLE must be defined"
#endif
#ifndef NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT
#error "NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT must be defined"
#endif

namespace {

std::vector<double> ReadFloat64Npy(const fs::path& path,
                                   size_t expected_values) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    if (!in) return {};
    char magic[6] = {};
    char version[2] = {};
    std::uint16_t header_len = 0;
    in.read(magic, 6);
    in.read(version, 2);
    in.read(reinterpret_cast<char*>(&header_len), sizeof(header_len));
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 0);
    std::string header(header_len, '\0');
    in.read(header.data(), static_cast<std::streamsize>(header_len));
    EXPECT_NE(header.find("'descr': '<f8'"), std::string::npos);
    std::vector<char> bytes{std::istreambuf_iterator<char>(in),
                            std::istreambuf_iterator<char>()};
    EXPECT_EQ(bytes.size(), expected_values * sizeof(double));
    if (bytes.size() != expected_values * sizeof(double)) return {};
    std::vector<double> values(expected_values);
    std::memcpy(values.data(), bytes.data(), bytes.size());
    return values;
}

void WriteOrcaTensorMatrix(std::ofstream& out, const Mat3& tensor) {
    for (int row = 0; row < 3; ++row) {
        out << ' ' << tensor(row, 0) << ' ' << tensor(row, 1) << ' '
            << tensor(row, 2) << '\n';
    }
}

fs::path WriteTransformedOrcaSource(
        const ProteinConformation& source,
        const nmr::test::directional::OrthogonalTransform& transform,
        const std::string& tag) {
    const fs::path path = fs::temp_directory_path() /
        ("directional_orca_" + tag + "_" + std::to_string(::getpid()) +
         ".out");
    std::ofstream out(path);
    EXPECT_TRUE(out.is_open()) << path;
    out << std::setprecision(17);
    out << "CHEMICAL SHIELDINGS (ppm)\n";
    for (size_t i = 0; i < source.AtomCount(); ++i) {
        const auto& atom = source.ProteinRef().AtomAt(i);
        const auto& ca = source.AtomAt(i);
        out << " Nucleus " << i << SymbolForElement(atom.element) << " :\n";
        out << "Diamagnetic contribution to the shielding tensor (ppm) :\n";
        WriteOrcaTensorMatrix(out, nmr::test::directional::EvenRank2(
            transform, ca.orca_shielding_diamagnetic));
        out << "Paramagnetic contribution to the shielding tensor (ppm):\n";
        WriteOrcaTensorMatrix(out, nmr::test::directional::EvenRank2(
            transform, ca.orca_shielding_paramagnetic));
        out << "Total shielding tensor (ppm):\n";
        WriteOrcaTensorMatrix(out, nmr::test::directional::EvenRank2(
            transform, ca.orca_shielding_total));
    }
    return path;
}

void CopyTransformedMutationInputs(
        ProteinConformation& destination,
        const ProteinConformation& source,
        const nmr::test::directional::OrthogonalTransform& transform,
        bool include_apbs) {
    ASSERT_EQ(destination.AtomCount(), source.AtomCount());
    for (size_t i = 0; i < source.AtomCount(); ++i) {
        auto& dst = destination.MutableAtomAt(i);
        const auto& src = source.AtomAt(i);
        dst.partial_charge = src.partial_charge;
        dst.mopac_charge = src.mopac_charge;
        dst.role = src.role;
        dst.is_backbone = src.is_backbone;
        if (include_apbs) {
            dst.apbs_efield =
                nmr::test::directional::Polar(transform, src.apbs_efield);
            dst.apbs_efg = nmr::test::directional::EvenRank2(
                transform, src.apbs_efg);
            dst.apbs_efg_spherical =
                SphericalTensor::Decompose(dst.apbs_efg);
        }
    }
    if (include_apbs) {
        destination.ForceAttachResultForTesting(
            std::make_unique<ApbsFieldResult>());
    }
}

}  // namespace



// ============================================================================
// Helper: load an ORCA protein and run the full Layer 0 pipeline
// ============================================================================

struct LoadedProtein {
    std::unique_ptr<Protein> protein;
    bool ok = false;
};

static LoadedProtein LoadAndPrepare(const std::string& prefix) {
    OrcaRunFiles files;
    files.pdb_path    = std::string(nmr::test::TestEnvironment::OrcaDir()) + prefix + ".pdb";
    files.xyz_path    = std::string(nmr::test::TestEnvironment::OrcaDir()) + prefix + ".xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + prefix + ".prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        return {nullptr, false};

    auto load = BuildFromOrca(files);
    if (!load.Ok()) return {nullptr, false};

    auto& conf = load.protein->Conformation();

    // Full Layer 0 pipeline
    conf.AttachResult(ProtonationDetectionResult::Compute(conf));
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(EnrichmentResult::Compute(conf));

    PrmtopChargeSource charge_source(files.prmtop_path);
    conf.AttachResult(ChargeAssignmentResult::Compute(conf, charge_source));

    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto dssp = DsspResult::Compute(conf);
    if (dssp) conf.AttachResult(std::move(dssp));

    conf.AttachResult(ApbsFieldResult::Compute(conf));
    conf.AttachResult(MolecularGraphResult::Compute(conf));

    // ORCA shielding tensors
    std::string nmr_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + prefix + "_nmr.out";
    if (fs::exists(nmr_path)) {
        auto orca = OrcaShieldingResult::Compute(conf, nmr_path);
        if (orca) conf.AttachResult(std::move(orca));
    }

    return {std::move(load.protein), true};
}


// ============================================================================
// Test fixture
// ============================================================================

class MutationDeltaTest : public ::testing::Test {
protected:
    void SetUp() override {
        wt_ = LoadAndPrepare("A0A7C5FAR6_WT");
        ala_ = LoadAndPrepare("A0A7C5FAR6_ALA");

        if (!wt_.ok || !ala_.ok)
            GTEST_SKIP() << "ORCA test data not available";

        auto& wt_conf = wt_.protein->Conformation();
        auto& ala_conf = ala_.protein->Conformation();

        if (!wt_conf.HasResult<OrcaShieldingResult>() ||
            !ala_conf.HasResult<OrcaShieldingResult>())
            GTEST_SKIP() << "ORCA shielding tensors not loaded";
    }

    LoadedProtein wt_;
    LoadedProtein ala_;
};


// ============================================================================
// Tests
// ============================================================================

TEST_F(MutationDeltaTest, ComputeSucceeds) {
    auto& wt_conf = wt_.protein->Conformation();
    const auto& ala_conf = ala_.protein->Conformation();

    auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
    ASSERT_NE(delta, nullptr);
    ASSERT_TRUE(wt_conf.AttachResult(std::move(delta)));
}


TEST_F(MutationDeltaTest,
       OrcaAndMutationDirectionalSourcesRerunUnderProperAndImproperO3) {
    auto& wt_original = wt_.protein->Conformation();
    auto& mut_original = ala_.protein->Conformation();
    auto baseline = MutationDeltaResult::Compute(wt_original, mut_original);
    ASSERT_NE(baseline, nullptr);
    ASSERT_GT(baseline->MatchedAtomCount(), 0u);

    const bool include_apbs =
        wt_original.HasResult<ApbsFieldResult>() &&
        mut_original.HasResult<ApbsFieldResult>();
    ASSERT_TRUE(include_apbs)
        << "delta_apbs.npy covariance proof requires both canonical APBS "
           "sources; do not silently downgrade this exact-name check";

    constexpr double kTensorAbsTolerance = 2.0e-9;
    constexpr double kTensorRelTolerance = 3.0e-12;
    constexpr double kVectorAbsTolerance = 2.0e-10;
    constexpr double kVectorRelTolerance = 3.0e-12;
    constexpr double kGeometryAbsTolerance = 2.0e-8;
    constexpr double kStructuralZeroTolerance = 2.0e-9;

    using MatrixMember = Mat3 MatchedAtomData::*;
    struct ShieldingArray {
        const char* name;
        MatrixMember matrix;
    };
    const std::array<ShieldingArray, 7> mutation_shielding_arrays{{
        {"delta_shielding.npy", &MatchedAtomData::delta_shielding},
        {"wt_shielding_diamagnetic.npy",
         &MatchedAtomData::wt_shielding_diamagnetic},
        {"wt_shielding_paramagnetic.npy",
         &MatchedAtomData::wt_shielding_paramagnetic},
        {"mut_shielding_diamagnetic.npy",
         &MatchedAtomData::mut_shielding_diamagnetic},
        {"mut_shielding_paramagnetic.npy",
         &MatchedAtomData::mut_shielding_paramagnetic},
        {"delta_shielding_diamagnetic.npy",
         &MatchedAtomData::delta_shielding_diamagnetic},
        {"delta_shielding_paramagnetic.npy",
         &MatchedAtomData::delta_shielding_paramagnetic},
    }};

    auto expect_full_tensor = [&](const Mat3& actual,
                                  const Mat3& source,
                                  const SphericalTensor& actual_spherical,
                                  const SphericalTensor& source_spherical,
                                  const auto& transform,
                                  size_t atom_index,
                                  const char* label) {
        const Mat3 expected =
            nmr::test::directional::EvenRank2(transform, source);
        EXPECT_TRUE(nmr::test::directional::NearMatrix(
            actual, expected, kTensorAbsTolerance, kTensorRelTolerance))
            << label << " atom=" << atom_index;
        EXPECT_TRUE(nmr::test::directional::Near(
            actual_spherical.T0, source_spherical.T0,
            kTensorAbsTolerance, kTensorRelTolerance))
            << label << " T0 atom=" << atom_index;
        EXPECT_TRUE(nmr::test::directional::NearVector(
            nmr::test::directional::T1Vector(actual_spherical),
            nmr::test::directional::Axial(
                transform,
                nmr::test::directional::T1Vector(source_spherical)),
            kTensorAbsTolerance, kTensorRelTolerance))
            << label << " T1 atom=" << atom_index;
        const SphericalTensor expected_t2 =
            nmr::test::directional::RotateNativeT2(
                transform, source_spherical);
        for (size_t component = 0; component < 5; ++component) {
            EXPECT_TRUE(nmr::test::directional::Near(
                actual_spherical.T2[component],
                expected_t2.T2[component],
                kTensorAbsTolerance, kTensorRelTolerance))
                << label << " T2=" << component << " atom=" << atom_index;
        }
    };

    for (const bool improper : {false, true}) {
        const auto transform = nmr::test::directional::SeededTransform(
            0x4f5243414d555441ULL, improper);
        ProteinConformation& wt_moved = wt_.protein->AddConformation(
            nmr::test::directional::Positions(transform,
                                              wt_original.Positions()),
            improper ? "mutation WT improper" : "mutation WT proper");
        ProteinConformation& mut_moved = ala_.protein->AddConformation(
            nmr::test::directional::Positions(transform,
                                              mut_original.Positions()),
            improper ? "mutation mutant improper" :
                       "mutation mutant proper");

        CopyTransformedMutationInputs(
            wt_moved, wt_original, transform, include_apbs);
        CopyTransformedMutationInputs(
            mut_moved, mut_original, transform, include_apbs);
        ASSERT_TRUE(wt_moved.AttachResult(GeometryResult::Compute(wt_moved)));
        ASSERT_TRUE(mut_moved.AttachResult(GeometryResult::Compute(mut_moved)));

        const std::string suffix = improper ? "improper" : "proper";
        const fs::path wt_orca_source = WriteTransformedOrcaSource(
            wt_original, transform, "wt_" + suffix);
        const fs::path mut_orca_source = WriteTransformedOrcaSource(
            mut_original, transform, "mut_" + suffix);
        auto wt_orca = OrcaShieldingResult::Compute(
            wt_moved, wt_orca_source.string());
        auto mut_orca = OrcaShieldingResult::Compute(
            mut_moved, mut_orca_source.string());
        ASSERT_NE(wt_orca, nullptr);
        ASSERT_NE(mut_orca, nullptr);
        ASSERT_TRUE(wt_moved.AttachResult(std::move(wt_orca)));
        ASSERT_TRUE(mut_moved.AttachResult(std::move(mut_orca)));

        for (size_t i = 0; i < wt_original.AtomCount(); ++i) {
            const auto& a = wt_original.AtomAt(i);
            const auto& b = wt_moved.AtomAt(i);
            expect_full_tensor(
                b.orca_shielding_total, a.orca_shielding_total,
                b.orca_shielding_total_spherical,
                a.orca_shielding_total_spherical,
                transform, i, "orca_total.npy");
            expect_full_tensor(
                b.orca_shielding_diamagnetic,
                a.orca_shielding_diamagnetic,
                b.orca_shielding_diamagnetic_spherical,
                a.orca_shielding_diamagnetic_spherical,
                transform, i, "orca_diamagnetic.npy");
            expect_full_tensor(
                b.orca_shielding_paramagnetic,
                a.orca_shielding_paramagnetic,
                b.orca_shielding_paramagnetic_spherical,
                a.orca_shielding_paramagnetic_spherical,
                transform, i, "orca_paramagnetic.npy");
        }

        auto moved = MutationDeltaResult::Compute(wt_moved, mut_moved);
        ASSERT_NE(moved, nullptr);
        ASSERT_EQ(moved->MatchedAtomCount(), baseline->MatchedAtomCount());
        for (size_t i = 0; i < wt_original.AtomCount(); ++i) {
            ASSERT_EQ(moved->HasMatch(i), baseline->HasMatch(i));
            if (!baseline->HasMatch(i)) continue;
            const auto& a = baseline->MatchedDataAt(i);
            const auto& b = moved->MatchedDataAt(i);
            ASSERT_EQ(b.mut_index, a.mut_index);
            EXPECT_NEAR(b.match_distance, a.match_distance,
                        kGeometryAbsTolerance);

            for (const auto& spec : mutation_shielding_arrays) {
                const Mat3& source = a.*(spec.matrix);
                const Mat3& actual = b.*(spec.matrix);
                EXPECT_TRUE(nmr::test::directional::NearMatrix(
                    actual,
                    nmr::test::directional::EvenRank2(transform, source),
                    kTensorAbsTolerance, kTensorRelTolerance))
                    << spec.name << " atom=" << i;
            }

            if (include_apbs) {
                EXPECT_TRUE(nmr::test::directional::NearVector(
                    b.delta_efield,
                    nmr::test::directional::Polar(
                        transform, a.delta_efield),
                    kVectorAbsTolerance, kVectorRelTolerance))
                    << "delta_apbs.npy E atom=" << i;
                EXPECT_TRUE(nmr::test::directional::NearMatrix(
                    b.delta_efg,
                    nmr::test::directional::EvenRank2(
                        transform, a.delta_efg),
                    kTensorAbsTolerance, kTensorRelTolerance))
                    << "delta_apbs.npy EFG atom=" << i;
                EXPECT_NEAR(b.delta_efg_spherical.T0, 0.0,
                            kStructuralZeroTolerance);
                for (double component : b.delta_efg_spherical.T1) {
                    EXPECT_NEAR(component, 0.0,
                                kStructuralZeroTolerance);
                }
                const auto expected_t2 =
                    nmr::test::directional::RotateNativeT2(
                        transform, a.delta_efg_spherical);
                for (size_t component = 0; component < 5; ++component) {
                    EXPECT_TRUE(nmr::test::directional::Near(
                        b.delta_efg_spherical.T2[component],
                        expected_t2.T2[component],
                        kTensorAbsTolerance, kTensorRelTolerance));
                }
            }

            ASSERT_EQ(b.removed_ring_proximity.size(),
                      a.removed_ring_proximity.size());
            for (size_t ring = 0;
                 ring < a.removed_ring_proximity.size(); ++ring) {
                const auto& ar = a.removed_ring_proximity[ring];
                const auto& br = b.removed_ring_proximity[ring];
                EXPECT_NEAR(br.distance, ar.distance,
                            kGeometryAbsTolerance);
                EXPECT_NEAR(br.z,
                            transform.Determinant() * ar.z,
                            kGeometryAbsTolerance);
                EXPECT_NEAR(br.rho, ar.rho, kGeometryAbsTolerance);
                EXPECT_NEAR(br.theta,
                            std::atan2(ar.rho,
                                       transform.Determinant() * ar.z),
                            kGeometryAbsTolerance);
                EXPECT_NEAR(br.mcconnell_factor, ar.mcconnell_factor,
                            kGeometryAbsTolerance);
                EXPECT_NEAR(br.exp_decay, ar.exp_decay,
                            kGeometryAbsTolerance);
            }
        }

        const fs::path out_dir = fs::temp_directory_path() /
            ("directional_mutation_" + suffix + "_" +
             std::to_string(::getpid()));
        fs::create_directories(out_dir);
        ASSERT_EQ(wt_moved.Result<OrcaShieldingResult>().WriteFeatures(
                      wt_moved, out_dir.string()), 3);
        for (const auto& spec : {
                 std::pair{"orca_total.npy",
                           &ConformationAtom::orca_shielding_total},
                 std::pair{"orca_diamagnetic.npy",
                           &ConformationAtom::orca_shielding_diamagnetic},
                 std::pair{"orca_paramagnetic.npy",
                           &ConformationAtom::orca_shielding_paramagnetic}}) {
            const auto values = ReadFloat64Npy(
                out_dir / spec.first, wt_moved.AtomCount() * 9);
            ASSERT_EQ(values.size(), wt_moved.AtomCount() * 9);
            for (size_t i = 0; i < wt_moved.AtomCount(); ++i) {
                double expected[9] = {};
                SphericalTensor::Decompose(
                    wt_moved.AtomAt(i).*(spec.second)).PackFull9(expected);
                for (size_t c = 0; c < 9; ++c) {
                    EXPECT_NEAR(values[i * 9 + c], expected[c],
                                kTensorAbsTolerance)
                        << spec.first << " atom=" << i << " column=" << c;
                }
            }
        }

        ASSERT_GT(moved->WriteFeatures(wt_moved, out_dir.string()), 0);
        ASSERT_EQ(nmr::test::directional::RunNumpyAllowPickleFalse(
                      NMR_TEST_PYTHON_EXECUTABLE,
                      NMR_NPY_ALLOW_PICKLE_FALSE_SCRIPT,
                      {out_dir / "orca_total.npy",
                       out_dir / "orca_diamagnetic.npy",
                       out_dir / "orca_paramagnetic.npy",
                       out_dir / "delta_shielding.npy",
                       out_dir / "wt_shielding_diamagnetic.npy",
                       out_dir / "wt_shielding_paramagnetic.npy",
                       out_dir / "mut_shielding_diamagnetic.npy",
                       out_dir / "mut_shielding_paramagnetic.npy",
                       out_dir / "delta_shielding_diamagnetic.npy",
                       out_dir / "delta_shielding_paramagnetic.npy",
                       out_dir / "delta_apbs.npy",
                       out_dir / "delta_ring_proximity.npy"}),
                  0);
        for (const auto& spec : mutation_shielding_arrays) {
            const auto values = ReadFloat64Npy(
                out_dir / spec.name, wt_moved.AtomCount() * 9);
            ASSERT_EQ(values.size(), wt_moved.AtomCount() * 9);
            for (size_t i = 0; i < wt_moved.AtomCount(); ++i) {
                double expected[9] = {};
                if (baseline->HasMatch(i)) {
                    const Mat3 transformed =
                        nmr::test::directional::EvenRank2(
                            transform,
                            baseline->MatchedDataAt(i).*(spec.matrix));
                    SphericalTensor::Decompose(transformed)
                        .PackFull9(expected);
                }
                for (size_t c = 0; c < 9; ++c) {
                    EXPECT_NEAR(values[i * 9 + c], expected[c],
                                kTensorAbsTolerance)
                        << spec.name << " atom=" << i << " column=" << c;
                }
            }
        }

        if (include_apbs) {
            const auto values = ReadFloat64Npy(
                out_dir / "delta_apbs.npy", wt_moved.AtomCount() * 12);
            ASSERT_EQ(values.size(), wt_moved.AtomCount() * 12);
            for (size_t i = 0; i < wt_moved.AtomCount(); ++i) {
                double expected[12] = {};
                if (baseline->HasMatch(i)) {
                    const auto& base = baseline->MatchedDataAt(i);
                    const Vec3 e = nmr::test::directional::Polar(
                        transform, base.delta_efield);
                    expected[0] = e.x();
                    expected[1] = e.y();
                    expected[2] = e.z();
                    SphericalTensor::Decompose(
                        nmr::test::directional::EvenRank2(
                            transform, base.delta_efg))
                        .PackFull9(expected + 3);
                }
                for (size_t c = 0; c < 12; ++c) {
                    EXPECT_NEAR(values[i * 12 + c], expected[c],
                                kTensorAbsTolerance)
                        << "delta_apbs.npy atom=" << i
                        << " column=" << c;
                }
            }
        }

        size_t removed_rings = 0;
        for (size_t i = 0; i < wt_original.AtomCount(); ++i) {
            if (baseline->HasMatch(i)) {
                removed_rings = baseline->MatchedDataAt(i)
                                    .removed_ring_proximity.size();
                break;
            }
        }
        ASSERT_GT(removed_rings, 0u);
        const size_t ring_columns = removed_rings * 6;
        const auto ring_values = ReadFloat64Npy(
            out_dir / "delta_ring_proximity.npy",
            wt_moved.AtomCount() * ring_columns);
        ASSERT_EQ(ring_values.size(),
                  wt_moved.AtomCount() * ring_columns);
        for (size_t i = 0; i < wt_moved.AtomCount(); ++i) {
            if (!baseline->HasMatch(i)) {
                for (size_t c = 0; c < ring_columns; ++c) {
                    EXPECT_DOUBLE_EQ(ring_values[i * ring_columns + c], 0.0);
                }
                continue;
            }
            const auto& base = baseline->MatchedDataAt(i);
            ASSERT_EQ(base.removed_ring_proximity.size(), removed_rings);
            for (size_t ring = 0; ring < removed_rings; ++ring) {
                const auto& rp = base.removed_ring_proximity[ring];
                const double expected[6] = {
                    rp.distance,
                    transform.Determinant() * rp.z,
                    rp.rho,
                    std::atan2(rp.rho,
                               transform.Determinant() * rp.z),
                    rp.mcconnell_factor,
                    rp.exp_decay,
                };
                for (size_t c = 0; c < 6; ++c) {
                    EXPECT_NEAR(
                        ring_values[i * ring_columns + ring * 6 + c],
                        expected[c], kGeometryAbsTolerance)
                        << "delta_ring_proximity.npy atom=" << i
                        << " ring=" << ring << " column=" << c;
                }
            }
        }

        for (const auto& entry : fs::directory_iterator(out_dir)) {
            EXPECT_EQ(std::remove(entry.path().string().c_str()), 0)
                << entry.path();
        }
        EXPECT_EQ(std::remove(out_dir.string().c_str()), 0) << out_dir;
        fs::remove(wt_orca_source);
        fs::remove(mut_orca_source);
    }
}


TEST_F(MutationDeltaTest, FourMutationSitesDetected) {
    auto& wt_conf = wt_.protein->Conformation();
    const auto& ala_conf = ala_.protein->Conformation();

    auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
    ASSERT_NE(delta, nullptr);

    const auto& sites = delta->MutationSites();
    EXPECT_EQ(sites.size(), 4u);

    for (const auto& site : sites) {
        EXPECT_EQ(site.mut_type, AminoAcid::ALA);
        std::cout << "  Mutation site " << site.residue_index << ": "
                  << ThreeLetterCodeForAminoAcid(site.wt_type) << " -> "
                  << ThreeLetterCodeForAminoAcid(site.mut_type)
                  << " (" << site.wt_ring_indices.size() << " rings)\n";
    }

    // Check all aromatic types present
    std::set<AminoAcid> wt_types;
    for (const auto& s : sites) wt_types.insert(s.wt_type);
    EXPECT_TRUE(wt_types.count(AminoAcid::TRP) > 0);
    EXPECT_TRUE(wt_types.count(AminoAcid::TYR) > 0);
    EXPECT_TRUE(wt_types.count(AminoAcid::HIS) > 0);
    EXPECT_TRUE(wt_types.count(AminoAcid::PHE) > 0);

    // TRP should have 3 rings (benzene + pyrrole + perimeter), others 1
    for (const auto& s : sites) {
        if (s.wt_type == AminoAcid::TRP)
            EXPECT_EQ(s.wt_ring_indices.size(), 3u);
        else
            EXPECT_EQ(s.wt_ring_indices.size(), 1u);
    }

    wt_conf.AttachResult(std::move(delta));
}


TEST_F(MutationDeltaTest, AtomMatchingReasonable) {
    auto& wt_conf = wt_.protein->Conformation();
    const auto& ala_conf = ala_.protein->Conformation();

    auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
    ASSERT_NE(delta, nullptr);

    std::cout << "  WT=" << wt_conf.AtomCount()
              << " ALA=" << ala_conf.AtomCount()
              << " matched=" << delta->MatchedAtomCount()
              << " unmatched=" << delta->UnmatchedWtAtomCount() << "\n";

    EXPECT_GT(delta->MatchedAtomCount(), 400u);
    EXPECT_GT(delta->UnmatchedWtAtomCount(), 0u);
    EXPECT_EQ(delta->MatchedAtomCount() + delta->UnmatchedWtAtomCount(),
              wt_conf.AtomCount());
}


TEST_F(MutationDeltaTest, BackboneAtomsMatch) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    const Protein& p = wt_conf.ProteinRef();
    size_t res0_N = p.ResidueAt(0).N;
    if (res0_N != Residue::NONE) {
        EXPECT_TRUE(delta->HasMatch(res0_N));
        EXPECT_LT(delta->MatchDistanceAt(res0_N), 0.1);
    }
}


TEST_F(MutationDeltaTest, DeltaShieldingFinite) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    int checked = 0;
    for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
        if (!delta->HasMatch(ai)) continue;
        const auto& st = delta->DeltaShieldingSphericalAt(ai);
        EXPECT_FALSE(std::isnan(st.T0)) << "NaN T0 at atom " << ai;
        for (int i = 0; i < 5; ++i)
            EXPECT_FALSE(std::isnan(st.T2[i])) << "NaN T2[" << i << "] at " << ai;
        checked++;
    }
    EXPECT_GT(checked, 400);
}

TEST_F(MutationDeltaTest,
       ApbsDeltaPayloadUsesCanonicalReactionFieldsAndStructuralZeros) {
    auto& wt_conf = wt_.protein->Conformation();
    auto& ala_conf = ala_.protein->Conformation();

    auto mapping = MutationDeltaResult::Compute(wt_conf, ala_conf);
    ASSERT_NE(mapping, nullptr);
    size_t pinned_wt = SIZE_MAX;
    for (size_t wi = 0; wi < wt_conf.AtomCount(); ++wi) {
        if (mapping->HasMatch(wi)) {
            pinned_wt = wi;
            break;
        }
    }
    ASSERT_NE(pinned_wt, SIZE_MAX);
    const size_t pinned_mut = mapping->MutantAtomFor(pinned_wt);
    ASSERT_LT(pinned_mut, ala_conf.AtomCount());

    // Independent sentinels pin WT-minus-mutant sign and, critically, the
    // canonical reaction fields rather than the deliberately different raw
    // total-PB diagnostics.
    auto& wt_atom = wt_conf.MutableAtomAt(pinned_wt);
    auto& mut_atom = ala_conf.MutableAtomAt(pinned_mut);
    wt_atom.apbs_efield = Vec3(1.0, -2.0, 3.0);
    mut_atom.apbs_efield = Vec3(-4.0, 5.0, -6.0);
    wt_atom.apbs_efg = Mat3::Zero();
    wt_atom.apbs_efg.diagonal() << 1.0, 2.0, -3.0;
    mut_atom.apbs_efg = Mat3::Zero();
    mut_atom.apbs_efg.diagonal() << -4.0, 5.0, -1.0;
    wt_atom.apbs_efield_total_diagnostic = Vec3(101.0, 102.0, 103.0);
    mut_atom.apbs_efield_total_diagnostic = Vec3(-201.0, -202.0, -203.0);
    wt_atom.apbs_efg_total_diagnostic = 100.0 * Mat3::Identity();
    mut_atom.apbs_efg_total_diagnostic = -200.0 * Mat3::Identity();
    wt_atom.apbs_efg_total_diagnostic_spherical =
        SphericalTensor::Decompose(wt_atom.apbs_efg_total_diagnostic);
    mut_atom.apbs_efg_total_diagnostic_spherical =
        SphericalTensor::Decompose(mut_atom.apbs_efg_total_diagnostic);

    auto delta = MutationDeltaResult::Compute(wt_conf, ala_conf);
    ASSERT_NE(delta, nullptr);

    int checked = 0;
    double max_abs_t0_t1 = 0.0;
    for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
        if (!delta->HasMatch(ai)) continue;
        const auto& m = delta->MatchedDataAt(ai);
        if (!m.has_apbs_delta) continue;
        const auto& st = m.delta_efg_spherical;
        max_abs_t0_t1 = std::max(max_abs_t0_t1, std::abs(st.T0));
        for (int k = 0; k < 3; ++k)
            max_abs_t0_t1 = std::max(max_abs_t0_t1, std::abs(st.T1[k]));
        ++checked;
    }
    EXPECT_GT(checked, 400);
    EXPECT_LT(max_abs_t0_t1, 1e-10)
        << "delta_apbs full9 compatibility T0/T1 slots must stay structural zeros";

    const fs::path out_dir = fs::temp_directory_path() /
        ("mutation_delta_apbs_" + std::to_string(::getpid()));
    fs::create_directories(out_dir);
    delta->WriteFeatures(wt_conf, out_dir.string());
    const auto payload = ReadFloat64Npy(
        out_dir / "delta_apbs.npy", wt_conf.AtomCount() * 12);
    ASSERT_EQ(payload.size(), wt_conf.AtomCount() * 12);
    const size_t base = pinned_wt * 12;
    const std::array<double, 12> expected = {
        5.0, -7.0, 9.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -std::sqrt(6.0), 0.0, 4.0 * std::sqrt(2.0)};
    for (size_t column = 0; column < expected.size(); ++column) {
        EXPECT_NEAR(payload[base + column], expected[column], 1e-12)
            << "delta_apbs column " << column;
    }
    for (const char* name : {
            "delta_shielding.npy",
            "wt_shielding_diamagnetic.npy",
            "wt_shielding_paramagnetic.npy",
            "mut_shielding_diamagnetic.npy",
            "mut_shielding_paramagnetic.npy",
            "delta_shielding_diamagnetic.npy",
            "delta_shielding_paramagnetic.npy",
            "delta_scalars.npy", "delta_graph.npy", "delta_apbs.npy",
            "delta_ring_proximity.npy"}) {
        std::error_code ec;
        fs::remove(out_dir / name, ec);
    }
    std::error_code ec;
    fs::remove(out_dir, ec);
}


TEST_F(MutationDeltaTest, RingProximityComputed) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    // Count total removed rings
    size_t total_rings = 0;
    for (const auto& s : delta->MutationSites())
        total_rings += s.wt_ring_indices.size();

    int has_proximity = 0;
    for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
        if (!delta->HasMatch(ai)) continue;
        const auto& m = delta->MatchedDataAt(ai);
        if (!m.removed_ring_proximity.empty()) {
            // Each matched atom should have proximity to ALL removed rings
            EXPECT_EQ(m.removed_ring_proximity.size(), total_rings);

            // Nearest removed ring distance should be non-negative
            EXPECT_GE(m.nearest_removed_ring_dist, 0.0);
            EXPECT_LT(m.nearest_removed_ring_dist, 50.0);

            // Cylindrical coords should be consistent
            for (const auto& rp : m.removed_ring_proximity) {
                EXPECT_GE(rp.distance, 0.0);
                double r_from_cyl = std::sqrt(rp.z * rp.z + rp.rho * rp.rho);
                EXPECT_NEAR(r_from_cyl, rp.distance, 0.01)
                    << "Cylindrical coords inconsistent at atom " << ai;
            }
            has_proximity++;
        }
    }
    EXPECT_GT(has_proximity, 400) << "Most matched atoms should have ring proximity";
}


TEST_F(MutationDeltaTest, DistanceDecayCurve) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    const auto& summary = delta->Summary();

    std::cout << "  Distance decay of |delta T0|:\n";
    for (const auto& bin : summary.by_distance) {
        if (bin.count > 0) {
            std::cout << "    " << (int)bin.bin_start << "-"
                      << (int)bin.bin_end << " A: n="
                      << bin.count << " mean_|dT0|="
                      << bin.mean_abs_delta_t0
                      << " mean_|T2|=" << bin.mean_t2_magnitude << "\n";
        }
    }

    // Signal should be stronger near rings than far away
    // Find the 3-4A bin and the 10-11A bin
    double near_signal = 0, far_signal = 0;
    for (const auto& bin : summary.by_distance) {
        if (bin.bin_start == 3.0 && bin.count > 0) near_signal = bin.mean_abs_delta_t0;
        if (bin.bin_start == 10.0 && bin.count > 0) far_signal = bin.mean_abs_delta_t0;
    }
    if (near_signal > 0 && far_signal > 0) {
        EXPECT_GT(near_signal, far_signal)
            << "Signal should decay with distance from removed rings";
    }
}


TEST_F(MutationDeltaTest, ElementStratification) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    const auto& summary = delta->Summary();

    std::cout << "  Delta by element:\n";
    for (const auto& bin : summary.by_element) {
        std::cout << "    " << SymbolForElement(bin.element)
                  << ": n=" << bin.count
                  << " mean_dT0=" << bin.mean_delta_t0
                  << " mean_|dT0|=" << bin.mean_abs_delta_t0
                  << " max_|dT0|=" << bin.max_abs_delta_t0
                  << " mean_|T2|=" << bin.mean_t2_magnitude << "\n";
    }

    EXPECT_FALSE(summary.by_element.empty());

    // Backbone vs sidechain
    std::cout << "  Backbone: n=" << summary.backbone_count
              << " mean_|dT0|=" << summary.backbone_mean_abs_t0 << "\n";
    std::cout << "  Sidechain: n=" << summary.sidechain_count
              << " mean_|dT0|=" << summary.sidechain_mean_abs_t0 << "\n";
}


TEST_F(MutationDeltaTest, DsspDeltaAvailable) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    EXPECT_TRUE(delta->HasDsspDelta());

    if (delta->HasDsspDelta()) {
        // SASA should change at mutation sites (exposed surface changes)
        double max_sasa_delta = 0;
        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            if (!delta->HasMatch(ai)) continue;
            max_sasa_delta = std::max(max_sasa_delta,
                std::abs(delta->MatchedDataAt(ai).delta_sasa));
        }
        std::cout << "  Max |SASA delta|: " << max_sasa_delta << " A^2\n";
        EXPECT_GT(max_sasa_delta, 0.0) << "Some SASA should change";
    }
}


TEST_F(MutationDeltaTest, GraphDeltaAvailable) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    EXPECT_TRUE(delta->HasGraphDelta());

    if (delta->HasGraphDelta()) {
        // Ring atoms in WT had graph_dist_ring=0; in mutant those rings
        // don't exist, so graph_dist_ring is large. Delta should be negative
        // (WT was closer to rings).
        int large_delta = 0;
        int nonzero_graph_delta = 0;
        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            if (!delta->HasMatch(ai)) continue;
            const auto& production = delta->MatchedDataAt(ai);
            int dd = production.delta_graph_dist_ring;
            if (std::abs(dd) > 3) large_delta++;
            if (dd != 0 || std::abs(production.delta_bfs_decay) > 1e-15 ||
                production.delta_is_conjugated != 0) {
                ++nonzero_graph_delta;
            }
        }
        std::cout << "  Atoms with |graph delta| > 3: " << large_delta << "\n";
        EXPECT_GT(large_delta, 0)
            << "aromatic-to-ALA production graphs must change long-range "
               "ring distance for at least one matched atom";
        EXPECT_GT(nonzero_graph_delta, 0);

        const fs::path out_dir = fs::temp_directory_path() /
            ("mutation_delta_graph_production_" +
             std::to_string(::getpid()));
        fs::create_directories(out_dir);
        ASSERT_GT(delta->WriteFeatures(wt_conf, out_dir.string()), 0);
        const auto emitted = ReadFloat64Npy(
            out_dir / "delta_graph.npy", wt_conf.AtomCount() * 5);
        ASSERT_EQ(emitted.size(), wt_conf.AtomCount() * 5);

        for (size_t ai = 0; ai < wt_conf.AtomCount(); ++ai) {
            const size_t row = ai * 5;
            if (!delta->HasMatch(ai)) {
                for (size_t column = 0; column < 5; ++column) {
                    EXPECT_DOUBLE_EQ(emitted[row + column], 0.0);
                }
                continue;
            }
            const auto& production = delta->MatchedDataAt(ai);
            EXPECT_DOUBLE_EQ(emitted[row + 0], 1.0);
            EXPECT_DOUBLE_EQ(emitted[row + 1],
                             production.has_graph_delta ? 1.0 : 0.0);
            EXPECT_DOUBLE_EQ(emitted[row + 2],
                             static_cast<double>(
                                 production.delta_graph_dist_ring));
            EXPECT_DOUBLE_EQ(emitted[row + 3],
                             production.delta_bfs_decay);
            EXPECT_DOUBLE_EQ(emitted[row + 4],
                             static_cast<double>(
                                 production.delta_is_conjugated));
        }
        for (const char* name : {
                "delta_shielding.npy",
                "wt_shielding_diamagnetic.npy",
                "wt_shielding_paramagnetic.npy",
                "mut_shielding_diamagnetic.npy",
                "mut_shielding_paramagnetic.npy",
                "delta_shielding_diamagnetic.npy",
                "delta_shielding_paramagnetic.npy",
                "delta_scalars.npy", "delta_graph.npy", "delta_apbs.npy",
                "delta_ring_proximity.npy"}) {
            std::remove((out_dir / name).string().c_str());
        }
        ::rmdir(out_dir.string().c_str());
    }
}


TEST_F(MutationDeltaTest, AccessViaTemplateAfterAttach) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);
    ASSERT_TRUE(wt_conf.AttachResult(std::move(delta)));

    ASSERT_TRUE(wt_conf.HasResult<MutationDeltaResult>());
    const auto& d = wt_conf.Result<MutationDeltaResult>();
    EXPECT_GT(d.MatchedAtomCount(), 0u);
    EXPECT_EQ(d.MutationSites().size(), 4u);
    EXPECT_FALSE(d.Summary().by_element.empty());
}


// ============================================================================
// Typed-identity matcher invariants (rewrite landed 2026-05-08).
//
// These tests exercise the new substrate-typed-identity matching block
// and document the methodology guarantees the rewrite delivers:
//   * Mechanical-swap fixtures produce zero-drift bound matches.
//   * Mutation residues are wholly excluded from binding.
//   * Non-mutation residues bijection: every WT atom finds a unique
//     mut atom (for this fixture there are no variant differences).
//
// CategoryInfo identity drives the matchup; this is not position-only.
// ============================================================================

TEST_F(MutationDeltaTest, MechanicalSwapDriftIsSubAngstrom) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    // For the A0A7C5FAR6 fixture, every non-mutation atom should be in
    // (essentially) the same coordinates between WT and mut — that's
    // the structural assumption of the mechanical-mutation pipeline.
    // The typed matcher binds by identity; this test documents that
    // those bindings DO agree with spatial proximity at sub-Angstrom
    // distances on real fixture data.
    double max_drift = 0.0;
    for (size_t wi = 0; wi < wt_conf.AtomCount(); ++wi) {
        if (!delta->HasMatch(wi)) continue;
        max_drift = std::max(max_drift, delta->MatchDistanceAt(wi));
    }
    EXPECT_LT(max_drift, 0.5)
        << "max bound-atom drift " << max_drift << " Å — mechanical-swap "
        "fixture should produce near-zero drift across all bound atoms. "
        "Drift > 0.5 Å indicates the mut-generation pipeline changed the "
        "non-mutation coordinates, the substrate is mis-tagging atoms, "
        "or rotamers flipped (in which case typed match is correct and "
        "this assertion's threshold should be relaxed to characterise "
        "the new fixture).";
}


TEST_F(MutationDeltaTest, MutationSiteAtomsAreNotBound) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    const Protein& p = wt_conf.ProteinRef();

    // Build the set of mutation residue indices.
    std::set<size_t> mutation_residues;
    for (const auto& site : delta->MutationSites()) {
        mutation_residues.insert(site.residue_index);
    }
    ASSERT_FALSE(mutation_residues.empty())
        << "Test fixture should have at least one mutation site.";

    // Every WT atom in a mutation residue must be unmatched. The pre-
    // 2026-05-08 spatial matcher silently bound across this boundary
    // (e.g., WT PHE-CG to mut ALA-CB by proximity) — the typed matcher
    // refuses, surfacing the chemistry difference cleanly.
    size_t mutation_atoms_total = 0;
    size_t mutation_atoms_bound_in_error = 0;
    for (size_t wi = 0; wi < p.AtomCount(); ++wi) {
        if (!mutation_residues.count(p.AtomAt(wi).residue_index)) continue;
        ++mutation_atoms_total;
        if (delta->HasMatch(wi)) ++mutation_atoms_bound_in_error;
    }
    EXPECT_GT(mutation_atoms_total, 0u)
        << "Mutation residues should contain WT atoms.";
    EXPECT_EQ(0u, mutation_atoms_bound_in_error)
        << mutation_atoms_bound_in_error << " WT atoms in mutation residues "
        "were bound to mut atoms — the typed matcher should refuse all "
        "bindings across mutation chemistry boundaries.";
}


TEST_F(MutationDeltaTest, NonMutationResiduesAchieveBijection) {
    auto& wt_conf = wt_.protein->Conformation();
    auto delta = MutationDeltaResult::Compute(wt_conf, ala_.protein->Conformation());
    ASSERT_NE(delta, nullptr);

    const Protein& p = wt_conf.ProteinRef();

    // For mutation-free residues with matching variant index, every
    // WT atom must find a unique mut atom (the typed identity tuple
    // is unique-within-residue except for equivalent-H sets, which
    // bind by consume-in-order). Bijection is the strict guarantee.
    std::set<size_t> mutation_residues;
    for (const auto& site : delta->MutationSites()) {
        mutation_residues.insert(site.residue_index);
    }

    std::set<size_t> mut_atoms_used;
    size_t non_mut_atoms = 0;
    size_t non_mut_unmatched = 0;
    for (size_t wi = 0; wi < p.AtomCount(); ++wi) {
        const size_t ri = p.AtomAt(wi).residue_index;
        if (mutation_residues.count(ri)) continue;
        ++non_mut_atoms;
        if (!delta->HasMatch(wi)) {
            ++non_mut_unmatched;
            continue;
        }
        const size_t mi = delta->MutantAtomFor(wi);
        EXPECT_EQ(0u, mut_atoms_used.count(mi))
            << "mut atom " << mi << " is bound to two WT atoms — typed "
            "matcher's consume-in-order should produce 1:1 binding.";
        mut_atoms_used.insert(mi);
    }
    EXPECT_GT(non_mut_atoms, 0u);
    // For this fixture (no variant differences), every non-mutation
    // WT atom should match. The new matcher's `no_id_match` and
    // `variant_unmatched` counters surface deviations elsewhere.
    EXPECT_EQ(0u, non_mut_unmatched)
        << non_mut_unmatched << " WT atoms in non-mutation residues "
        "were not matched. Either the substrate is incomplete on this "
        "fixture (regression in ComposeAtomSemantic), or the fixture "
        "contains a variant difference the test wasn't expecting.";
}
