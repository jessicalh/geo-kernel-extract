#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include "PdbFileReader.h"
#include "GeometryResult.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "ApbsFieldResult.h"
#include "CalculatorConfig.h"
#include "PhysicalConstants.h"
#include <algorithm>
#include <filesystem>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <cstring>
#include <unistd.h>

using namespace nmr;

namespace {

void WriteApbsTestConfig(const std::filesystem::path& path,
                         double pdie,
                         double sdie,
                         double ionic_strength,
                         double clamp) {
    std::ofstream out(path);
    // 161 is the production grid and avoids APBS's coarse-grid direct-solver
    // fallback, whose legacy teardown is not safe at 65^3 for this protein.
    out << "apbs_grid_dim = 161\n"
        << "apbs_manual_grid_padding_A = 70.0\n"
        << "apbs_manual_grid_min_dim_A = 70.0\n"
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

void RunHomogeneousReferenceForcingFunctionInChild() {
    const auto config_path = std::filesystem::temp_directory_path() /
        ("apbs_homogeneous_reference_" + std::to_string(::getpid()) +
         ".toml");
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
        }
    }

    std::filesystem::remove(config_path);
    _exit(ok ? 0 : 1);
}

void RunAnalyticSignAndClampForcingFunctionInChild() {
    constexpr double kClampVPerA = 100.0;
    constexpr double kChargeE = 500.0;
    constexpr double kDistanceA = 3.5;
    constexpr double kDielectric = 4.0;
    const auto config_path = std::filesystem::temp_directory_path() /
        ("apbs_analytic_sign_clamp_" + std::to_string(::getpid()) +
         ".toml");
    const auto output_dir = std::filesystem::temp_directory_path() /
        ("apbs_analytic_sign_clamp_" + std::to_string(::getpid()));
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

            std::filesystem::create_directories(output_dir);
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
            ok = ok && emitted_E.size() == 6
                    && emitted_efg.size() == 10
                    && emitted_phi.size() == 2
                    && emitted_total_E.size() == 6
                    && emitted_total_efg.size() == 10;
            if (ok) {
                const size_t e = kProbe * 3;
                const size_t t2 = kProbe * 5;
                // These are direct pins on the frozen payloads, derived from
                // the continuum configuration above rather than from the
                // in-memory writer source fields.
                ok = emitted_E[e] < 0.0
                    && std::abs(emitted_E[e] + kClampVPerA) < 1e-10
                    && emitted_phi[kProbe] < 0.0
                    && emitted_efg[t2 + 2] > 0.0
                    && emitted_efg[t2 + 4] < 0.0
                    && emitted_total_E[e] > 0.0
                    && emitted_total_efg[t2 + 2] < 0.0
                    && emitted_total_efg[t2 + 4] > 0.0;
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
