#include <gtest/gtest.h>

#include "CalculatorConfig.h"
#include "ChargeAssignmentResult.h"
#include "CoulombResult.h"
#include "EeqCoulombResult.h"
#include "EeqResult.h"
#include "GeometryResult.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "SpatialIndexResult.h"
#include "TestEnvironment.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct NpyArray {
    std::string descr;
    std::vector<size_t> shape;
    std::vector<char> bytes;
};

std::string Trim(std::string value) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(),
                             [&](char c) { return !is_space(c); }));
    value.erase(std::find_if(value.rbegin(), value.rend(),
                             [&](char c) { return !is_space(c); }).base(),
                value.end());
    return value;
}

NpyArray ReadNpy(const fs::path& path) {
    std::ifstream input(path, std::ios::binary);
    EXPECT_TRUE(input.is_open()) << path;
    NpyArray array;
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
    input.read(reinterpret_cast<char*>(&header_length),
               sizeof(header_length));
    std::string header(header_length, '\0');
    input.read(header.data(), static_cast<std::streamsize>(header.size()));

    const std::string descr_key = "'descr': '";
    auto descr_begin = header.find(descr_key);
    if (descr_begin == std::string::npos) {
        ADD_FAILURE() << header;
        return array;
    }
    descr_begin += descr_key.size();
    const auto descr_end = header.find('\'', descr_begin);
    if (descr_end == std::string::npos) {
        ADD_FAILURE() << header;
        return array;
    }
    array.descr = header.substr(descr_begin, descr_end - descr_begin);

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
    while (std::getline(shape_stream, token, ',')) {
        token = Trim(token);
        if (!token.empty()) {
            array.shape.push_back(
                static_cast<size_t>(std::stoull(token)));
        }
    }

    array.bytes.assign(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
    return array;
}

template <typename T>
const T* DataAs(const NpyArray& array) {
    return reinterpret_cast<const T*>(array.bytes.data());
}

void RemoveDirectoryContents(const fs::path& directory) {
    std::error_code ec;
    if (fs::exists(directory, ec)) {
        for (const auto& entry : fs::directory_iterator(directory, ec)) {
            fs::remove(entry.path(), ec);
        }
        fs::remove(directory, ec);
    }
}

std::unique_ptr<nmr::Protein> MakeEeqFixture() {
    auto protein = std::make_unique<nmr::Protein>();
    auto carbon = nmr::Atom::Create(nmr::Element::C);
    auto nitrogen = nmr::Atom::Create(nmr::Element::N);
    auto hydrogen = nmr::Atom::Create(nmr::Element::H);
    carbon->residue_index = 0;
    nitrogen->residue_index = 0;
    hydrogen->residue_index = 0;
    hydrogen->parent_atom_index = 1;
    protein->AddAtom(std::move(carbon));
    protein->AddAtom(std::move(nitrogen));
    protein->AddAtom(std::move(hydrogen));

    nmr::Residue residue;
    residue.type = nmr::AminoAcid::ALA;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    residue.atom_indices = {0, 1, 2};
    residue.CA = 0;
    residue.N = 1;
    residue.H = 2;
    protein->AddResidue(residue);

    const std::vector<nmr::Vec3> positions = {
        nmr::Vec3(0.0, 0.0, 0.0),
        nmr::Vec3(1.4, 0.2, 0.0),
        nmr::Vec3(2.1, 0.9, 0.3),
    };
    protein->FinalizeConstruction(positions);
    protein->MutableAtomAt(2).parent_atom_index = 1;
    protein->AddCrystalConformation(positions, 0.0, 0.0, 0.0, "eeq-test");
    return protein;
}

// Deliberately non-degenerate equivalence fixture.  It contains typed
// backbone sources, a substrate-recognised PHE aromatic ring, a non-ring
// sidechain source, and two geometrically remote sources.  The remote pair
// makes cutoff provenance observable even under the historical off-by-one;
// the large test charges below force the shared 100 V/A clamp.
std::unique_ptr<nmr::Protein> MakeNondegenerateCoulombFixture() {
    auto protein = std::make_unique<nmr::Protein>();

    auto add_residue = [&](nmr::AminoAcid type, int sequence_number) {
        nmr::Residue residue;
        residue.type = type;
        residue.sequence_number = sequence_number;
        residue.chain_id = "A";
        return protein->AddResidue(std::move(residue));
    };
    const size_t ala = add_residue(nmr::AminoAcid::ALA, 1);
    const size_t phe = add_residue(nmr::AminoAcid::PHE, 2);
    const size_t cys_near = add_residue(nmr::AminoAcid::CYS, 3);
    const size_t cys_far_1 = add_residue(nmr::AminoAcid::CYS, 4);
    const size_t cys_far_2 = add_residue(nmr::AminoAcid::CYS, 5);

    auto add_atom = [&](size_t residue_index, nmr::Element element,
                        const char* name) {
        auto atom = nmr::Atom::Create(element);
        atom->residue_index = residue_index;
        atom->pdb_atom_name = name;
        const size_t index = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(residue_index).atom_indices.push_back(index);
        return index;
    };

    const size_t ala_n = add_atom(ala, nmr::Element::N, "N");
    add_atom(ala, nmr::Element::C, "CA");
    add_atom(ala, nmr::Element::C, "C");
    add_atom(ala, nmr::Element::O, "O");
    const size_t ala_h = add_atom(ala, nmr::Element::H, "H");

    for (const char* name : {"CG", "CD1", "CE1", "CZ", "CE2", "CD2"}) {
        add_atom(phe, nmr::Element::C, name);
    }
    add_atom(cys_near, nmr::Element::S, "SG");
    add_atom(cys_far_1, nmr::Element::S, "SG");
    add_atom(cys_far_2, nmr::Element::S, "SG");

    constexpr double ring_r = nmr::BENZENE_CC_BOND_LENGTH;
    constexpr double ring_half = ring_r * 0.5;
    constexpr double ring_sin60 = ring_r * 0.8660254037844386;
    const std::vector<nmr::Vec3> positions = {
        // ALA backbone cluster.
        nmr::Vec3(10.000, 0.000, 0.000),
        nmr::Vec3(11.458, 0.000, 0.000),
        nmr::Vec3(12.983, 0.000, 0.000),
        nmr::Vec3(12.983, 1.231, 0.000),
        nmr::Vec3(9.000, 0.000, 0.000),
        // Regular PHE ring: substrate topology recognises these six atoms.
        nmr::Vec3(0.000, 0.000, 0.000),
        nmr::Vec3(ring_r, 0.000, 0.000),
        nmr::Vec3(ring_r + ring_half, ring_sin60, 0.000),
        nmr::Vec3(ring_r, 2.0 * ring_sin60, 0.000),
        nmr::Vec3(0.000, 2.0 * ring_sin60, 0.000),
        nmr::Vec3(-ring_half, ring_sin60, 0.000),
        // One nearby non-ring sidechain and two sources beyond 20 A from
        // the main cluster (but 10 A apart from each other).
        nmr::Vec3(5.000, 3.000, 1.000),
        nmr::Vec3(40.000, 0.000, 0.000),
        nmr::Vec3(50.000, 0.000, 0.000),
    };

    protein->FinalizeConstruction(positions);
    protein->MutableAtomAt(ala_h).parent_atom_index = ala_n;
    protein->AddCrystalConformation(
        positions, 0.0, 0.0, 0.0, "eeq-coulomb-nondegenerate");
    return protein;
}

const nmr::NamedNumber* FindNumber(const nmr::GeometryChoice& choice,
                                   const char* name) {
    for (const auto& number : choice.Numbers()) {
        if (number.name == name) return &number;
    }
    return nullptr;
}

}  // namespace

TEST(EeqDiagnostics, StoredValuesMatchIndependentD4Formulas) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = MakeEeqFixture();
    auto& conf = protein->Conformation();

    auto eeq = nmr::EeqResult::Compute(conf, 0);
    ASSERT_NE(eeq, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(eeq)));

    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto params = nmr::D4EeqParamsFor(protein->AtomAt(i).element);
        const auto& atom = conf.AtomAt(i);
        const double expected_chi_eff = params.chi + params.kappa *
            std::sqrt(atom.eeq_cn + 1e-14);
        const double expected_diag = params.gam +
            nmr::SQRT_2_OVER_PI / params.rad;
        EXPECT_NEAR(atom.eeq_chi_eff, expected_chi_eff, 1e-14);
        EXPECT_DOUBLE_EQ(atom.eeq_eta, params.gam);
        EXPECT_DOUBLE_EQ(atom.eeq_self_hardness_diag, expected_diag);
        EXPECT_TRUE(std::isfinite(atom.eeq_chi_eff));
        EXPECT_GT(atom.eeq_eta, 0.0);
        EXPECT_GT(atom.eeq_self_hardness_diag, 0.0);
    }

    // Output-level composition pin: the NPY mapping must preserve the
    // independently formula-checked fields above, including hardness column
    // order [eta, A_diag].
    const fs::path output_dir = fs::temp_directory_path() /
        ("eeq_diagnostics_" + std::to_string(::getpid()));
    RemoveDirectoryContents(output_dir);
    ASSERT_TRUE(fs::create_directories(output_dir));
    ASSERT_EQ(conf.Result<nmr::EeqResult>().WriteFeatures(
                  conf, output_dir.string()), 4);

    const auto charges = ReadNpy(output_dir / "eeq_charges.npy");
    const auto cn = ReadNpy(output_dir / "eeq_cn.npy");
    const auto chi = ReadNpy(output_dir / "eeq_chi_eff.npy");
    const auto hardness = ReadNpy(output_dir / "eeq_hardness.npy");
    const size_t N = conf.AtomCount();
    ASSERT_EQ(charges.descr, "<f8");
    ASSERT_EQ(cn.descr, "<f8");
    ASSERT_EQ(chi.descr, "<f8");
    ASSERT_EQ(hardness.descr, "<f8");
    ASSERT_EQ(charges.shape, (std::vector<size_t>{N}));
    ASSERT_EQ(cn.shape, (std::vector<size_t>{N}));
    ASSERT_EQ(chi.shape, (std::vector<size_t>{N}));
    ASSERT_EQ(hardness.shape, (std::vector<size_t>{N, 2u}));
    ASSERT_EQ(charges.bytes.size(), N * sizeof(double));
    ASSERT_EQ(cn.bytes.size(), N * sizeof(double));
    ASSERT_EQ(chi.bytes.size(), N * sizeof(double));
    ASSERT_EQ(hardness.bytes.size(), N * 2 * sizeof(double));

    const double* charge_data = DataAs<double>(charges);
    const double* cn_data = DataAs<double>(cn);
    const double* chi_data = DataAs<double>(chi);
    const double* hardness_data = DataAs<double>(hardness);
    for (size_t i = 0; i < N; ++i) {
        const auto params = nmr::D4EeqParamsFor(protein->AtomAt(i).element);
        const double expected_chi = params.chi + params.kappa *
            std::sqrt(conf.AtomAt(i).eeq_cn + 1e-14);
        const double expected_diag = params.gam +
            nmr::SQRT_2_OVER_PI / params.rad;
        EXPECT_DOUBLE_EQ(charge_data[i], conf.AtomAt(i).eeq_charge);
        EXPECT_DOUBLE_EQ(cn_data[i], conf.AtomAt(i).eeq_cn);
        EXPECT_NEAR(chi_data[i], expected_chi, 1e-14);
        EXPECT_DOUBLE_EQ(hardness_data[i*2], params.gam);
        EXPECT_DOUBLE_EQ(hardness_data[i*2 + 1], expected_diag);
    }
    RemoveDirectoryContents(output_dir);
}

TEST(EeqCoulomb, EqualsForceFieldCoulombWhenSourceChargesAreEqual) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = MakeNondegenerateCoulombFixture();
    auto& conf = protein->Conformation();
    ASSERT_EQ(protein->RingCount(), 1u)
        << "fixture must expose a typed PHE aromatic ring";
    ASSERT_TRUE(conf.AttachResult(nmr::GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::SpatialIndexResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::EeqResult::Compute(conf, 0)));

    // Independent forcing function: give both calculators the same explicit,
    // large, alternating source charges.  These are intentionally not the
    // small native EEQ solution: they make backbone, non-ring sidechain and
    // aromatic channels all nonzero and force the shared E-vector clamp.
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const double magnitude = 40.0 + 3.0 * static_cast<double>(i);
        const double charge = (i % 2 == 0) ? magnitude : -magnitude;
        conf.MutableAtomAt(i).eeq_charge = charge;
        conf.MutableAtomAt(i).partial_charge = charge;
    }
    conf.ForceAttachResultForTesting(
        std::make_unique<nmr::ChargeAssignmentResult>());

    auto ff = nmr::CoulombResult::Compute(conf);
    auto eeq_coulomb = nmr::EeqCoulombResult::Compute(conf);
    ASSERT_NE(ff, nullptr);
    ASSERT_NE(eeq_coulomb, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(ff)));
    ASSERT_TRUE(conf.AttachResult(std::move(eeq_coulomb)));

    double max_backbone = 0.0;
    double max_sidechain = 0.0;
    double max_aromatic = 0.0;
    int max_aromatic_source_count = 0;
    int clamped_atoms = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& atom = conf.AtomAt(i);
        max_backbone = std::max(max_backbone,
                                atom.eeq_coulomb_E_backbone.norm());
        max_sidechain = std::max(max_sidechain,
                                 atom.eeq_coulomb_E_sidechain.norm());
        max_aromatic = std::max(max_aromatic,
                                atom.eeq_coulomb_E_aromatic.norm());
        max_aromatic_source_count = std::max(
            max_aromatic_source_count,
            atom.eeq_coulomb_aromatic_n_sidechain_atoms);
        if (std::abs(atom.eeq_coulomb_E_magnitude -
                     nmr::CalculatorConfig::Get(
                         "efield_magnitude_sanity_clamp")) < 1e-10) {
            ++clamped_atoms;
        }
        for (int d = 0; d < 3; ++d) {
            EXPECT_DOUBLE_EQ(atom.eeq_coulomb_E_total(d),
                             atom.coulomb_E_total(d));
            EXPECT_DOUBLE_EQ(atom.eeq_coulomb_E_backbone(d),
                             atom.coulomb_E_backbone(d));
            EXPECT_DOUBLE_EQ(atom.eeq_coulomb_E_sidechain(d),
                             atom.coulomb_E_sidechain(d));
            EXPECT_DOUBLE_EQ(atom.eeq_coulomb_E_aromatic(d),
                             atom.coulomb_E_aromatic(d));
        }
        for (int row = 0; row < 3; ++row) {
            for (int col = 0; col < 3; ++col) {
                EXPECT_DOUBLE_EQ(atom.eeq_coulomb_EFG_total(row, col),
                                 atom.coulomb_EFG_total(row, col));
                EXPECT_DOUBLE_EQ(atom.eeq_coulomb_EFG_backbone(row, col),
                                 atom.coulomb_EFG_backbone(row, col));
                EXPECT_DOUBLE_EQ(atom.eeq_coulomb_EFG_sidechain(row, col),
                                 atom.coulomb_EFG_sidechain(row, col));
                EXPECT_DOUBLE_EQ(atom.eeq_coulomb_EFG_aromatic(row, col),
                                 atom.coulomb_EFG_aromatic(row, col));
            }
        }
        EXPECT_DOUBLE_EQ(atom.eeq_coulomb_E_magnitude,
                         atom.coulomb_E_magnitude);
        EXPECT_DOUBLE_EQ(atom.eeq_coulomb_E_backbone_frac,
                         atom.coulomb_E_backbone_frac);
        EXPECT_DOUBLE_EQ(atom.eeq_coulomb_aromatic_E_magnitude,
                         atom.aromatic_E_magnitude);
        EXPECT_EQ(atom.eeq_coulomb_aromatic_n_sidechain_atoms,
                  atom.aromatic_n_sidechain_atoms);
        if (protein->AtomAt(i).element == nmr::Element::H) {
            EXPECT_DOUBLE_EQ(atom.eeq_coulomb_E_bond_proj,
                             atom.coulomb_E_bond_proj);
            EXPECT_DOUBLE_EQ(atom.eeq_coulomb_aromatic_E_bond_proj,
                             atom.aromatic_E_bond_proj);
        } else {
            EXPECT_TRUE(std::isnan(atom.eeq_coulomb_E_bond_proj));
            EXPECT_TRUE(std::isnan(atom.coulomb_E_bond_proj));
            EXPECT_TRUE(std::isnan(atom.eeq_coulomb_aromatic_E_bond_proj));
            EXPECT_TRUE(std::isnan(atom.aromatic_E_bond_proj));
        }
    }
    EXPECT_GT(max_backbone, 0.0);
    EXPECT_GT(max_sidechain, 0.0);
    EXPECT_GT(max_aromatic, 0.0);
    EXPECT_GT(max_aromatic_source_count, 0);
    EXPECT_GT(clamped_atoms, 0)
        << "fixture must exercise the shared E-vector clamp";

    int ff_clamp_choices = 0;
    int eeq_clamp_choices = 0;
    int eeq_cutoff_choices = 0;
    const double cutoff =
        nmr::CalculatorConfig::Get("coulomb_efield_cutoff");
    for (const auto& choice : conf.geometry_choices) {
        if (choice.Calculator() == nmr::CalculatorId::Coulomb &&
            choice.Label() == "E-field clamp") {
            ++ff_clamp_choices;
        }
        if (choice.Calculator() == nmr::CalculatorId::EEQCoulomb &&
            choice.Label() == "EEQ E-field clamp") {
            ++eeq_clamp_choices;
        }
        if (choice.Calculator() != nmr::CalculatorId::EEQCoulomb ||
            choice.Label() != "EEQ Coulomb cutoff") {
            continue;
        }

        ++eeq_cutoff_choices;
        const size_t target = choice.GroupKey();
        ASSERT_LT(target, conf.AtomCount());
        int expected_within = 0;
        for (size_t source = 0; source < conf.AtomCount(); ++source) {
            if (source == target) continue;
            if ((conf.PositionAt(target) - conf.PositionAt(source)).norm() <=
                cutoff) {
                ++expected_within;
            }
        }
        const int expected_beyond =
            static_cast<int>(conf.AtomCount() - 1) - expected_within;
        ASSERT_GT(expected_beyond, 0);
        const auto* within = FindNumber(choice, "sources_within_cutoff");
        const auto* beyond = FindNumber(choice, "sources_beyond_cutoff");
        ASSERT_NE(within, nullptr);
        ASSERT_NE(beyond, nullptr);
        EXPECT_DOUBLE_EQ(within->value,
                         static_cast<double>(expected_within));
        EXPECT_DOUBLE_EQ(beyond->value,
                         static_cast<double>(expected_beyond));
    }
    EXPECT_GT(ff_clamp_choices, 0);
    EXPECT_GT(eeq_clamp_choices, 0);
    EXPECT_GT(eeq_cutoff_choices, 0)
        << "two far sources must exercise exact cutoff provenance";

    const fs::path out = fs::temp_directory_path() /
        ("eeq_coulomb_features_" + std::to_string(::getpid()));
    RemoveDirectoryContents(out);
    ASSERT_TRUE(fs::create_directories(out));
    ASSERT_EQ(conf.Result<nmr::EeqCoulombResult>().WriteFeatures(
                  conf, out.string()), 11);
    ASSERT_EQ(conf.Result<nmr::CoulombResult>().WriteFeatures(
                  conf, out.string()), 12);
    for (const char* name : {
            "eeq_coulomb_efg.npy", "eeq_coulomb_E.npy",
            "eeq_coulomb_E_backbone.npy", "eeq_coulomb_E_sidechain.npy",
            "eeq_coulomb_E_aromatic.npy", "eeq_coulomb_efg_backbone.npy",
            "eeq_coulomb_efg_sidechain.npy", "eeq_coulomb_efg_aromatic.npy",
            "eeq_coulomb_scalars.npy", "eeq_coulomb_aromatic_E_proj.npy",
            "eeq_coulomb_aromatic_n_src.npy"}) {
        EXPECT_TRUE(fs::exists(out / name)) << name;
    }
    EXPECT_FALSE(fs::exists(out / "eeq_coulomb_efg_t2.npy"));
    EXPECT_FALSE(fs::exists(out / "eeq_coulomb_E_solvent.npy"));

    // Output-level equivalence pin.  Equal source charges already proved the
    // independent FF calculator is a numeric authority above; comparing the
    // complete NPY payloads now catches every EEQ writer column/index/dtype
    // mapping, including NaN projection sentinels and int32 source counts.
    const std::pair<const char*, const char*> output_pairs[] = {
        {"eeq_coulomb_efg.npy", "coulomb_efg.npy"},
        {"eeq_coulomb_E.npy", "coulomb_E.npy"},
        {"eeq_coulomb_E_backbone.npy", "coulomb_E_backbone.npy"},
        {"eeq_coulomb_E_sidechain.npy", "coulomb_E_sidechain.npy"},
        {"eeq_coulomb_E_aromatic.npy", "coulomb_E_aromatic.npy"},
        {"eeq_coulomb_efg_backbone.npy", "coulomb_efg_backbone.npy"},
        {"eeq_coulomb_efg_sidechain.npy", "coulomb_efg_sidechain.npy"},
        {"eeq_coulomb_efg_aromatic.npy", "coulomb_efg_aromatic.npy"},
        {"eeq_coulomb_scalars.npy", "coulomb_scalars.npy"},
        {"eeq_coulomb_aromatic_E_proj.npy",
         "coulomb_aromatic_E_proj.npy"},
        {"eeq_coulomb_aromatic_n_src.npy",
         "coulomb_aromatic_n_src.npy"},
    };
    for (const auto& [eeq_name, ff_name] : output_pairs) {
        const auto eeq_array = ReadNpy(out / eeq_name);
        const auto ff_array = ReadNpy(out / ff_name);
        EXPECT_EQ(eeq_array.descr, ff_array.descr) << eeq_name;
        EXPECT_EQ(eeq_array.shape, ff_array.shape) << eeq_name;
        EXPECT_EQ(eeq_array.bytes, ff_array.bytes) << eeq_name;
    }
    RemoveDirectoryContents(out);
}

TEST(EeqCoulomb, NonfiniteSourceChargeFailsInsteadOfZeroSanitizing) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = MakeEeqFixture();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(nmr::GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::SpatialIndexResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::EeqResult::Compute(conf, 0)));

    conf.MutableAtomAt(0).eeq_charge =
        std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(nmr::EeqCoulombResult::Compute(conf), nullptr);
}

TEST(EeqWriters, FailedNpyOpensReturnZeroInsteadOfOptimisticCounts) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = MakeEeqFixture();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(nmr::GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::SpatialIndexResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::EeqResult::Compute(conf, 0)));
    ASSERT_TRUE(conf.AttachResult(nmr::EeqCoulombResult::Compute(conf)));

    // A regular file used as the would-be parent makes every
    // <parent>/<feature>.npy open fail with ENOTDIR.  Both writers must count
    // actual successful writes and log each full attempted path.
    const fs::path blocker = fs::temp_directory_path() /
        ("eeq_writer_blocker_" + std::to_string(::getpid()));
    std::error_code ec;
    fs::remove(blocker, ec);
    {
        std::ofstream output(blocker, std::ios::binary);
        ASSERT_TRUE(output.is_open());
        output << "not a directory";
    }
    ASSERT_TRUE(fs::is_regular_file(blocker));
    ::testing::internal::CaptureStderr();
    const int eeq_written = conf.Result<nmr::EeqResult>().WriteFeatures(
        conf, blocker.string());
    const int coulomb_written =
        conf.Result<nmr::EeqCoulombResult>().WriteFeatures(
            conf, blocker.string());
    const std::string error_log =
        ::testing::internal::GetCapturedStderr();
    EXPECT_EQ(eeq_written, 0);
    EXPECT_EQ(coulomb_written, 0);
    for (const char* filename : {
            "eeq_charges.npy", "eeq_cn.npy", "eeq_chi_eff.npy",
            "eeq_hardness.npy", "eeq_coulomb_efg.npy",
            "eeq_coulomb_E.npy", "eeq_coulomb_E_backbone.npy",
            "eeq_coulomb_E_sidechain.npy", "eeq_coulomb_E_aromatic.npy",
            "eeq_coulomb_efg_backbone.npy",
            "eeq_coulomb_efg_sidechain.npy",
            "eeq_coulomb_efg_aromatic.npy", "eeq_coulomb_scalars.npy",
            "eeq_coulomb_aromatic_E_proj.npy",
            "eeq_coulomb_aromatic_n_src.npy"}) {
        const std::string attempted_path =
            (blocker / filename).string();
        EXPECT_NE(error_log.find(attempted_path), std::string::npos)
            << "missing full failed-write path in error log: "
            << attempted_path;
    }
    fs::remove(blocker, ec);
}
