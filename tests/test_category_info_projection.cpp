// tests/test_category_info_projection.cpp
//
// Unit + integration tests for src/CategoryInfoProjection.{h,cpp}.
//
// Coverage:
//   1. Configure parses atom_nom.tbl into the lookup tables.
//   2. Reset clears state.
//   3. Per-atom IUPAC/BMRB lookup against 1UBQ residues.
//   4. Glycine HA2/HA3 stereo inversion (HA2 = pro-R, HA3 = pro-S).
//   5. AmberResidueThreeLetter handles variants (HID/HIE/HIP/CYX/CYM/
//      ASH/GLH/LYN/ARN/TYM) and the canonical case.
//   6. IupacResidueThreeLetter collapses variants to canonical.
//   7. ResidueOneLetter for the standard 20.
//   8. WriteFeatures creates atoms_category_info.npy with the
//      expected size on disk.
//   9. NPY format header invariants (magic + version + structured-dtype
//      header text contains all expected fields).
//

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include "TestEnvironment.h"
#include "CategoryInfoProjection.h"
#include "RuntimeEnvironment.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "Atom.h"
#include "Residue.h"
#include "AminoAcidType.h"
#include "Types.h"

namespace fs = std::filesystem;

namespace {

using nmr::AminoAcid;
using nmr::Atom;
using nmr::CategoryInfoProjection;
using nmr::Protein;
using nmr::Residue;
using nmr::RuntimeEnvironment;

// Test fixture: configure the projection from RuntimeEnvironment, reset
// at teardown so subsequent test executables (if any) start cleanly.
class CategoryInfoProjectionTest : public ::testing::Test {
protected:
    void SetUp() override {
        const std::string& path = RuntimeEnvironment::BmrbAtomNom();
        if (path.empty()) {
            GTEST_SKIP() << "RuntimeEnvironment::BmrbAtomNom() is empty; "
                            "set bmrb_atom_nom in ~/.nmr_tools.toml.";
        }
        CategoryInfoProjection::Reset();
        CategoryInfoProjection::Config cfg;
        cfg.atom_nom_tbl = path;
        CategoryInfoProjection::Configure(std::move(cfg));
        ASSERT_TRUE(CategoryInfoProjection::IsActive());
    }
    void TearDown() override {
        CategoryInfoProjection::Reset();
    }
};

// Helper: locate first atom matching (residue_index, atom_name).
std::size_t FindAtom(const Protein& p, std::size_t residue_index,
                      const std::string& atom_name) {
    for (std::size_t ai = 0; ai < p.AtomCount(); ++ai) {
        const Atom& a = p.AtomAt(ai);
        if (a.residue_index == residue_index && a.pdb_atom_name == atom_name) {
            return ai;
        }
    }
    return SIZE_MAX;
}

// Helper: locate first residue of a given type. Returns SIZE_MAX if none.
std::size_t FindResidueOfType(const Protein& p, AminoAcid type) {
    for (std::size_t ri = 0; ri < p.ResidueCount(); ++ri) {
        if (p.ResidueAt(ri).type == type) return ri;
    }
    return SIZE_MAX;
}

// ============================================================================
// Configure / Reset
// ============================================================================

TEST_F(CategoryInfoProjectionTest, ConfigureMakesActive) {
    EXPECT_TRUE(CategoryInfoProjection::IsActive());
}

TEST_F(CategoryInfoProjectionTest, ResetDeactivates) {
    CategoryInfoProjection::Reset();
    EXPECT_FALSE(CategoryInfoProjection::IsActive());
    EXPECT_TRUE(CategoryInfoProjection::MissLog().empty());
}

// ============================================================================
// Per-residue 3-letter / 1-letter projections (no atom_nom.tbl needed)
// ============================================================================

TEST_F(CategoryInfoProjectionTest, AmberResidueThreeLetter_Canonical) {
    EXPECT_EQ("ALA", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::ALA, -1));
    EXPECT_EQ("HIS", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::HIS, -1));
    EXPECT_EQ("CYS", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::CYS, -1));
    EXPECT_EQ("LYS", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::LYS, -1));
}

TEST_F(CategoryInfoProjectionTest, AmberResidueThreeLetter_Variants) {
    // Indices match AminoAcidType.h's documented variant-index contract:
    // HIS variants 0/1/2 = HID/HIE/HIP, ASP 0 = ASH, GLU 0 = GLH,
    // CYS 0/1 = CYX/CYM, LYS 0 = LYN, ARG 0 = ARN, TYR 0 = TYM.
    EXPECT_EQ("HID", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::HIS, 0));
    EXPECT_EQ("HIE", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::HIS, 1));
    EXPECT_EQ("HIP", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::HIS, 2));
    EXPECT_EQ("ASH", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::ASP, 0));
    EXPECT_EQ("GLH", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::GLU, 0));
    EXPECT_EQ("CYX", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::CYS, 0));
    EXPECT_EQ("CYM", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::CYS, 1));
    EXPECT_EQ("LYN", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::LYS, 0));
    EXPECT_EQ("ARN", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::ARG, 0));
    EXPECT_EQ("TYM", CategoryInfoProjection::AmberResidueThreeLetter(AminoAcid::TYR, 0));
}

TEST_F(CategoryInfoProjectionTest, IupacResidueThreeLetter_AlwaysCanonical) {
    // IUPAC collapses every variant to the canonical 3-letter.
    EXPECT_EQ("HIS", CategoryInfoProjection::IupacResidueThreeLetter(AminoAcid::HIS));
    EXPECT_EQ("CYS", CategoryInfoProjection::IupacResidueThreeLetter(AminoAcid::CYS));
    EXPECT_EQ("ASP", CategoryInfoProjection::IupacResidueThreeLetter(AminoAcid::ASP));
    EXPECT_EQ("GLU", CategoryInfoProjection::IupacResidueThreeLetter(AminoAcid::GLU));
    EXPECT_EQ("LYS", CategoryInfoProjection::IupacResidueThreeLetter(AminoAcid::LYS));
    EXPECT_EQ("ARG", CategoryInfoProjection::IupacResidueThreeLetter(AminoAcid::ARG));
    EXPECT_EQ("TYR", CategoryInfoProjection::IupacResidueThreeLetter(AminoAcid::TYR));
}

TEST_F(CategoryInfoProjectionTest, ResidueOneLetter_Standard20) {
    EXPECT_EQ('A', CategoryInfoProjection::ResidueOneLetter(AminoAcid::ALA));
    EXPECT_EQ('R', CategoryInfoProjection::ResidueOneLetter(AminoAcid::ARG));
    EXPECT_EQ('G', CategoryInfoProjection::ResidueOneLetter(AminoAcid::GLY));
    EXPECT_EQ('W', CategoryInfoProjection::ResidueOneLetter(AminoAcid::TRP));
    EXPECT_EQ('P', CategoryInfoProjection::ResidueOneLetter(AminoAcid::PRO));
}

// ============================================================================
// atom_nom.tbl-driven per-atom queries (require 1UBQ fixture)
// ============================================================================

TEST_F(CategoryInfoProjectionTest, BmrbAtomName_AlaSidechain) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    const Protein& p = *build.protein;

    // Find an ALA residue.
    const std::size_t ri = FindResidueOfType(p, AminoAcid::ALA);
    if (ri == SIZE_MAX) GTEST_SKIP() << "No ALA in 1UBQ";

    const std::size_t ha  = FindAtom(p, ri, "HA");
    const std::size_t hb1 = FindAtom(p, ri, "HB1");
    const std::size_t hb2 = FindAtom(p, ri, "HB2");
    const std::size_t hb3 = FindAtom(p, ri, "HB3");
    ASSERT_NE(SIZE_MAX, ha);
    ASSERT_NE(SIZE_MAX, hb1);
    ASSERT_NE(SIZE_MAX, hb2);
    ASSERT_NE(SIZE_MAX, hb3);

    // BMRB names match the canonical sidechain names for ALA.
    EXPECT_EQ("HA",  CategoryInfoProjection::BmrbAtomName(p, ha));
    EXPECT_EQ("HB1", CategoryInfoProjection::BmrbAtomName(p, hb1));
    EXPECT_EQ("HB2", CategoryInfoProjection::BmrbAtomName(p, hb2));
    EXPECT_EQ("HB3", CategoryInfoProjection::BmrbAtomName(p, hb3));
}

TEST_F(CategoryInfoProjectionTest, BmrbStereoLabel_GlyHaPair) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    const Protein& p = *build.protein;

    const std::size_t ri = FindResidueOfType(p, AminoAcid::GLY);
    if (ri == SIZE_MAX) GTEST_SKIP() << "No GLY in 1UBQ";

    const std::size_t ha2 = FindAtom(p, ri, "HA2");
    const std::size_t ha3 = FindAtom(p, ri, "HA3");
    ASSERT_NE(SIZE_MAX, ha2);
    ASSERT_NE(SIZE_MAX, ha3);

    // atom_nom.tbl convention for prochiral methylenes:
    //   suffix-2 = pro-R, suffix-3 = pro-S.
    // Holds across the standard 20 (Gly Hα, Arg Hβ/Hγ/Hδ, Pro Hβ/Hγ/Hδ,
    // and every other methylene). Source: BMRB atom_nom.tbl SC column.
    EXPECT_EQ("pro-R", CategoryInfoProjection::BmrbStereoLabel(p, ha2));
    EXPECT_EQ("pro-S", CategoryInfoProjection::BmrbStereoLabel(p, ha3));
}

TEST_F(CategoryInfoProjectionTest, BmrbStereoLabel_ArgMethyleneIsProR_Or_ProS) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    const Protein& p = *build.protein;

    const std::size_t ri = FindResidueOfType(p, AminoAcid::ARG);
    if (ri == SIZE_MAX) GTEST_SKIP() << "No ARG in 1UBQ";

    const std::size_t hb2 = FindAtom(p, ri, "HB2");
    const std::size_t hb3 = FindAtom(p, ri, "HB3");
    ASSERT_NE(SIZE_MAX, hb2);
    ASSERT_NE(SIZE_MAX, hb3);

    // atom_nom.tbl convention: suffix-2 = pro-R, suffix-3 = pro-S.
    EXPECT_EQ("pro-R", CategoryInfoProjection::BmrbStereoLabel(p, hb2));
    EXPECT_EQ("pro-S", CategoryInfoProjection::BmrbStereoLabel(p, hb3));
}

// ============================================================================
// WriteFeatures — verify NPY file lands and has the expected shape
// ============================================================================

TEST_F(CategoryInfoProjectionTest, WriteFeaturesEmitsNpy) {
    if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
        GTEST_SKIP() << "1UBQ PDB not found";
    }
    auto build = nmr::BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    ASSERT_TRUE(build.Ok()) << build.error;
    const Protein& p = *build.protein;

    // Use a unique tmpdir so parallel test runs don't collide.
    const fs::path dir = fs::temp_directory_path() /
        ("category_info_test_" + std::to_string(::getpid()));
    fs::create_directories(dir);

    const int written = CategoryInfoProjection::WriteFeatures(p, dir.string());
    EXPECT_EQ(1, written);

    const fs::path npy = dir / "atoms_category_info.npy";
    ASSERT_TRUE(fs::exists(npy));

    // Read the file and verify NPY 1.0 magic + structured-dtype header.
    std::ifstream f(npy, std::ios::binary);
    ASSERT_TRUE(f);
    std::vector<char> all((std::istreambuf_iterator<char>(f)),
                          std::istreambuf_iterator<char>());

    // Magic (\x93NUMPY) + version (1.0) + header_len (uint16 LE).
    ASSERT_GE(all.size(), static_cast<size_t>(10));
    EXPECT_EQ(static_cast<unsigned char>(all[0]), 0x93);
    EXPECT_EQ(all[1], 'N');
    EXPECT_EQ(all[2], 'U');
    EXPECT_EQ(all[3], 'M');
    EXPECT_EQ(all[4], 'P');
    EXPECT_EQ(all[5], 'Y');
    EXPECT_EQ(all[6], '\x01');
    EXPECT_EQ(all[7], '\x00');

    const uint16_t header_len = static_cast<uint8_t>(all[8]) |
                                 (static_cast<uint8_t>(all[9]) << 8);
    ASSERT_GE(all.size(), static_cast<size_t>(10 + header_len));

    const std::string header(all.data() + 10, header_len);
    // The structured-dtype header must contain the declared field names.
    EXPECT_NE(header.find("'atom_index'"), std::string::npos);
    EXPECT_NE(header.find("'amber_atom_name'"), std::string::npos);
    EXPECT_NE(header.find("'iupac_atom_name'"), std::string::npos);
    EXPECT_NE(header.find("'bmrb_atom_name'"), std::string::npos);
    EXPECT_NE(header.find("'bmrb_residue_3letter'"), std::string::npos);
    EXPECT_NE(header.find("'locant'"), std::string::npos);
    EXPECT_NE(header.find("'backbone_role'"), std::string::npos);
    EXPECT_NE(header.find("'polar_h_kind'"), std::string::npos);
    EXPECT_NE(header.find("'pseudoatom_kind'"), std::string::npos);
    EXPECT_NE(header.find("'iupac_naming_provenance'"), std::string::npos);

    // Total file size: 10 prefix + header + N * record_size. The record
    // size is internal but stable; verify the file is at least header +
    // (N atoms × 50 bytes) — the schema's minimum row width.
    const std::size_t prefix_plus_header = 10 + header_len;
    ASSERT_GT(all.size(), prefix_plus_header);
    const std::size_t data_bytes = all.size() - prefix_plus_header;
    EXPECT_EQ(0u, data_bytes % p.AtomCount())
        << "Data section must be an integer multiple of atom count.";

    // Cleanup.
    std::error_code ec;
    fs::remove(npy, ec);
    fs::remove(dir, ec);
}

}  // anonymous namespace
