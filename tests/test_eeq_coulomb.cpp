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

#include <cmath>
#include <filesystem>
#include <limits>
#include <memory>
#include <string>
#include <unistd.h>
#include <vector>

namespace fs = std::filesystem;

namespace {

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
}

TEST(EeqCoulomb, EqualsForceFieldCoulombWhenSourceChargesAreEqual) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto protein = MakeEeqFixture();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(nmr::GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::SpatialIndexResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(nmr::EeqResult::Compute(conf, 0)));

    // Independent forcing function: make the FF source column exactly the
    // already-computed EEQ source column, then compare the two calculators.
    for (size_t i = 0; i < conf.AtomCount(); ++i)
        conf.MutableAtomAt(i).partial_charge = conf.AtomAt(i).eeq_charge;
    conf.ForceAttachResultForTesting(
        std::make_unique<nmr::ChargeAssignmentResult>());

    auto ff = nmr::CoulombResult::Compute(conf);
    auto eeq_coulomb = nmr::EeqCoulombResult::Compute(conf);
    ASSERT_NE(ff, nullptr);
    ASSERT_NE(eeq_coulomb, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(ff)));
    ASSERT_TRUE(conf.AttachResult(std::move(eeq_coulomb)));

    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        const auto& atom = conf.AtomAt(i);
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

    const fs::path out = fs::temp_directory_path() /
        ("eeq_coulomb_features_" + std::to_string(::getpid()));
    fs::create_directories(out);
    EXPECT_EQ(conf.Result<nmr::EeqCoulombResult>().WriteFeatures(
                  conf, out.string()), 11);
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
    for (const char* name : {
            "eeq_coulomb_efg.npy", "eeq_coulomb_E.npy",
            "eeq_coulomb_E_backbone.npy", "eeq_coulomb_E_sidechain.npy",
            "eeq_coulomb_E_aromatic.npy", "eeq_coulomb_efg_backbone.npy",
            "eeq_coulomb_efg_sidechain.npy", "eeq_coulomb_efg_aromatic.npy",
            "eeq_coulomb_scalars.npy", "eeq_coulomb_aromatic_E_proj.npy",
            "eeq_coulomb_aromatic_n_src.npy"}) {
        fs::remove(out / name);
    }
    fs::remove(out);
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
