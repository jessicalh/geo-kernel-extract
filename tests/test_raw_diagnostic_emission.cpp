#include "TestEnvironment.h"
#include <gtest/gtest.h>

#include "CalculatorConfig.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "LegacyAmberTopology.h"
#include "MolecularGraphResult.h"
#include "MutationDeltaResult.h"
#include "NpyWriter.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"
#include "Protein.h"
#include "Ring.h"
#include "SasaResult.h"
#include "SpatialIndexResult.h"
#include "WaterHBondGeometryResult.h"

#include <cmath>
#include <array>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <limits>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include <unistd.h>

namespace fs = std::filesystem;

namespace {

using namespace nmr;

template <typename T>
std::vector<T> ReadNpyPayload(const fs::path& path,
                              const std::string& dtype,
                              const std::string& shape_fragment) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    std::vector<char> bytes((std::istreambuf_iterator<char>(in)),
                            std::istreambuf_iterator<char>());
    if (bytes.size() < 10) return {};
    const uint16_t header_len =
        static_cast<uint8_t>(bytes[8]) |
        (static_cast<uint16_t>(static_cast<uint8_t>(bytes[9])) << 8);
    if (bytes.size() < 10u + header_len) return {};
    const std::string header(bytes.data() + 10, header_len);
    EXPECT_NE(header.find("'descr': '" + dtype + "'"), std::string::npos);
    EXPECT_NE(header.find(shape_fragment), std::string::npos) << header;
    const size_t offset = 10u + header_len;
    EXPECT_EQ((bytes.size() - offset) % sizeof(T), 0u);
    std::vector<T> values((bytes.size() - offset) / sizeof(T));
    if (!values.empty())
        std::memcpy(values.data(), bytes.data() + offset,
                    values.size() * sizeof(T));
    return values;
}

fs::path TempDir(const char* stem) {
    const auto dir = fs::temp_directory_path() /
        (std::string(stem) + "_" + std::to_string(::getpid()));
    fs::create_directories(dir);
    return dir;
}

void RemoveFilesAndDir(const fs::path& dir,
                       std::initializer_list<const char*> files) {
    std::error_code ec;
    for (const char* file : files) fs::remove(dir / file, ec);
    fs::remove(dir, ec);
}

std::unique_ptr<Protein> BuildUnnamedProtein(
        const std::vector<Element>& elements,
        const std::vector<Vec3>& positions) {
    auto protein = std::make_unique<Protein>();
    Residue residue;
    residue.type = AminoAcid::ALA;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const size_t ri = protein->AddResidue(residue);
    for (Element element : elements) {
        auto atom = Atom::Create(element);
        atom->residue_index = ri;
        const size_t ai = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(ai);
    }
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "raw-emission-test");
    return protein;
}

std::unique_ptr<Protein> BuildAspCarboxylateProtein() {
    auto protein = std::make_unique<Protein>();
    Residue residue;
    residue.type = AminoAcid::ASP;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const size_t ri = protein->AddResidue(residue);

    const std::array<std::pair<const char*, Element>, 3> atoms{{
        {"CG", Element::C}, {"OD1", Element::O}, {"OD2", Element::O}}};
    for (const auto& [name, element] : atoms) {
        auto atom = Atom::Create(element);
        atom->pdb_atom_name = name;
        atom->residue_index = ri;
        const size_t ai = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(ai);
    }
    const std::vector<Vec3> positions{
        Vec3(0.0, 0.0, 0.0), Vec3(1.25, 0.0, 0.0),
        Vec3(-1.25, 0.0, 0.0)};
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "typed-asp-carboxylate");
    return protein;
}

std::unique_ptr<Protein> BuildAlaHydrogenProtein() {
    auto protein = std::make_unique<Protein>();
    Residue residue;
    residue.type = AminoAcid::ALA;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    const size_t ri = protein->AddResidue(residue);

    const std::array<std::tuple<const char*, Element, Vec3>, 8> atoms{{
        {"N",   Element::N, Vec3(0.0, 0.0, 0.0)},
        {"H",   Element::H, Vec3(1.0, 0.0, 0.0)},
        {"CA",  Element::C, Vec3(0.0, -1.45, 0.0)},
        {"HA",  Element::H, Vec3(0.0, -1.45, 1.09)},
        {"C",   Element::C, Vec3(-1.30, -2.10, 0.0)},
        {"O",   Element::O, Vec3(-2.40, -1.70, 0.0)},
        {"CB",  Element::C, Vec3(1.20, -2.20, 0.0)},
        {"HB1", Element::H, Vec3(2.10, -1.60, 0.0)},
    }};
    std::vector<Vec3> positions;
    positions.reserve(atoms.size());
    for (const auto& [name, element, position] : atoms) {
        auto atom = Atom::Create(element);
        atom->pdb_atom_name = name;
        atom->residue_index = ri;
        const size_t ai = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(ai);
        positions.push_back(position);
    }
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "typed-ala-hydrogens");
    return protein;
}

TEST(RawDiagnosticEmission, SasaFractionUsesProductionWriter) {
    auto protein = BuildUnnamedProtein({Element::C}, {Vec3::Zero()});
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));
    auto sasa = SasaResult::Compute(conf);
    ASSERT_NE(sasa, nullptr);
    const double probe = CalculatorConfig::Get("sasa_probe_radius");
    const double radius = BondiVdwRadius(Element::C) + probe;
    EXPECT_NEAR(sasa->AtomSASA(0), 4.0 * PI * radius * radius, 1e-12);

    const fs::path dir = TempDir("sasa_fraction_emit");
    ASSERT_EQ(sasa->WriteFeatures(conf, dir.string()), 3);
    const auto fraction = ReadNpyPayload<double>(
        dir / "atom_sasa_fraction.npy", "<f8", "'shape': (1,)");
    ASSERT_EQ(fraction.size(), 1u);
    EXPECT_DOUBLE_EQ(fraction[0], 1.0);
    EXPECT_DOUBLE_EQ(sasa_result_detail::NormalizedSasaFraction(
                         Element::C, conf.AtomAt(0).atom_sasa, probe),
                     fraction[0]);

    RemoveFilesAndDir(dir, {"atom_sasa.npy", "sasa_normal.npy",
                            "atom_sasa_fraction.npy"});
}

TEST(RawDiagnosticEmission, RingPairGeometryUsesProductionRowPacker) {
    RingGeometry a;
    a.center = Vec3::Zero();
    a.normal = Vec3(0.0, 0.0, 1.0);
    RingGeometry b;
    b.center = Vec3(3.0, 4.0, 12.0);
    b.normal = Vec3(0.0, 1.0, 0.0);

    const auto row =
        conformation_result_detail::PackRingPairGeometryRow(
            2, 5, 7, 11,
            static_cast<int>(RingTypeIndex::TrpBenzene),
            static_cast<int>(RingTypeIndex::TrpPyrrole),
            a, b, true);
    EXPECT_EQ(row[0], 2.0);
    EXPECT_EQ(row[1], 5.0);
    EXPECT_EQ(row[2], 7.0);
    EXPECT_EQ(row[3], 11.0);
    EXPECT_EQ(row[4], static_cast<double>(RingTypeIndex::TrpBenzene));
    EXPECT_EQ(row[5], static_cast<double>(RingTypeIndex::TrpPyrrole));
    EXPECT_DOUBLE_EQ(row[6], 13.0);
    EXPECT_DOUBLE_EQ(row[7], 0.0);
    EXPECT_DOUBLE_EQ(row[8], 1.0);
    EXPECT_DOUBLE_EQ(row[9], 12.0);
    EXPECT_DOUBLE_EQ(row[10], 4.0);
    EXPECT_DOUBLE_EQ(row[11], 5.0);
    EXPECT_EQ(row[12], 1.0);
}

TEST(RawDiagnosticEmission, SpatialRowsReadBackExactProductionNeighbours) {
    auto protein = BuildUnnamedProtein(
        {Element::C, Element::C, Element::C},
        {Vec3(0.0, 0.0, 0.0), Vec3(3.0, 0.0, 0.0),
         Vec3(20.0, 0.0, 0.0)});
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    auto spatial = SpatialIndexResult::Compute(conf);
    ASSERT_NE(spatial, nullptr);

    const fs::path dir = TempDir("spatial_neighbour_emit");
    ASSERT_EQ(spatial->WriteFeatures(conf, dir.string()), 1);
    const auto rows = ReadNpyPayload<double>(
        dir / "spatial_neighbors.npy", "<f8", "'shape': (2,6)");
    ASSERT_EQ(rows.size(), 12u);

    // Production storage is directed: 0->1 and 1->0, with no self or
    // out-of-cutoff atom-2 row.
    EXPECT_EQ(rows[0], 0.0);
    EXPECT_EQ(rows[1], 1.0);
    EXPECT_DOUBLE_EQ(rows[2], 1.0);
    EXPECT_DOUBLE_EQ(rows[5], 3.0);
    EXPECT_EQ(rows[6], 1.0);
    EXPECT_EQ(rows[7], 0.0);
    EXPECT_DOUBLE_EQ(rows[8], -1.0);
    EXPECT_DOUBLE_EQ(rows[11], 3.0);

    RemoveFilesAndDir(dir, {"spatial_neighbors.npy"});
}

TEST(RawDiagnosticEmission, MolecularGraphWriterUsesProductionCompute) {
    const fs::path fixture = nmr::test::TestEnvironment::UbqProtonated();
    ASSERT_TRUE(fs::is_regular_file(fixture)) << fixture;
    auto loaded = BuildFromProtonatedPdb(fixture.string());
    ASSERT_TRUE(loaded.Ok()) << loaded.error;
    auto& protein = *loaded.protein;
    auto& conf = protein.Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));
    auto graph = MolecularGraphResult::Compute(conf);
    ASSERT_NE(graph, nullptr);
    ASSERT_GT(protein.RingCount(), 0u);

    std::set<size_t> ring_atoms;
    for (size_t ri = 0; ri < protein.RingCount(); ++ri) {
        for (size_t ai : protein.RingAt(ri).atom_indices) {
            ring_atoms.insert(ai);
        }
    }
    const size_t ring_atom = *ring_atoms.begin();
    size_t bonded_non_ring = SIZE_MAX;
    for (size_t bond_index : protein.AtomAt(ring_atom).bond_indices) {
        const Bond& bond = protein.BondAt(bond_index);
        const size_t other = bond.atom_index_a == ring_atom
            ? bond.atom_index_b : bond.atom_index_a;
        if (!ring_atoms.count(other)) {
            bonded_non_ring = other;
            break;
        }
    }
    ASSERT_NE(bonded_non_ring, SIZE_MAX);

    const ConformationAtom& source = conf.AtomAt(ring_atom);
    EXPECT_EQ(source.graph_dist_ring, 0);
    EXPECT_EQ(source.nearest_ring_atom_index,
              static_cast<int>(ring_atom));
    EXPECT_DOUBLE_EQ(source.bfs_decay, 1.0);
    EXPECT_TRUE(source.is_conjugated);
    EXPECT_GT(source.n_pi_bonds_3, 0);

    const ConformationAtom& one_bond = conf.AtomAt(bonded_non_ring);
    EXPECT_EQ(one_bond.graph_dist_ring, 1);
    EXPECT_NEAR(one_bond.bfs_decay,
                std::exp(-1.0 / BFS_DECAY_LENGTH), 1e-15);

    double independent_eneg_sum_1 = 0.0;
    for (size_t bond_index : protein.AtomAt(ring_atom).bond_indices) {
        const Bond& bond = protein.BondAt(bond_index);
        const size_t other = bond.atom_index_a == ring_atom
            ? bond.atom_index_b : bond.atom_index_a;
        independent_eneg_sum_1 += protein.AtomAt(other).Electronegativity();
    }
    EXPECT_DOUBLE_EQ(source.eneg_sum_1, independent_eneg_sum_1);

    const fs::path dir = TempDir("molecular_graph_emit");
    ASSERT_EQ(graph->WriteFeatures(conf, dir.string()), 3);
    const std::string n_rows = std::to_string(conf.AtomCount());
    const auto integer = ReadNpyPayload<int32_t>(
        dir / "molecular_graph_int.npy", "<i4",
        "'shape': (" + n_rows + ",6)");
    const auto floating = ReadNpyPayload<double>(
        dir / "molecular_graph_float.npy", "<f8",
        "'shape': (" + n_rows + ",3)");
    const auto compatibility = ReadNpyPayload<double>(
        dir / "molecular_graph.npy", "<f8",
        "'shape': (" + n_rows + ",9)");
    const size_t N = conf.AtomCount();
    ASSERT_EQ(integer.size(), N * 6);
    ASSERT_EQ(floating.size(), N * 3);
    ASSERT_EQ(compatibility.size(), N * 9);
    EXPECT_EQ(integer[ring_atom*6 + 0], 0);
    EXPECT_EQ(integer[ring_atom*6 + 3], source.n_pi_bonds_3);
    EXPECT_EQ(integer[ring_atom*6 + 4], 1);
    EXPECT_EQ(integer[ring_atom*6 + 5], static_cast<int>(ring_atom));
    EXPECT_DOUBLE_EQ(floating[ring_atom*3 + 0], independent_eneg_sum_1);
    EXPECT_DOUBLE_EQ(floating[ring_atom*3 + 2], 1.0);
    EXPECT_EQ(integer[bonded_non_ring*6 + 0], 1);
    EXPECT_NEAR(floating[bonded_non_ring*3 + 2],
                std::exp(-1.0 / BFS_DECAY_LENGTH), 1e-15);
    EXPECT_DOUBLE_EQ(compatibility[ring_atom*9 + 7],
                     static_cast<double>(ring_atom));

    RemoveFilesAndDir(dir, {"molecular_graph_int.npy",
                            "molecular_graph_float.npy",
                            "molecular_graph.npy"});
}

TEST(RawDiagnosticEmission, EnrichmentRawSemanticsReadBack) {
    auto protein = BuildAspCarboxylateProtein();
    ASSERT_TRUE(protein->LegacyAmber().HasAtomSemantic());
    auto& conf = protein->Conformation();
    auto enrichment = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrichment, nullptr);

    const fs::path dir = TempDir("enrichment_semantic_emit");
    ASSERT_EQ(enrichment->WriteFeatures(conf, dir.string()), 12);

    const auto parent = ReadNpyPayload<uint8_t>(
        dir / "enrichment_parent_is_sp2.npy", "|u1", "'shape': (3,)");
    const auto polar = ReadNpyPayload<uint8_t>(
        dir / "semantic_polar_h_kind.npy", "|u1", "'shape': (3,)");
    const auto planar = ReadNpyPayload<uint8_t>(
        dir / "semantic_planar_group_kind.npy", "|u1", "'shape': (3,)");
    const auto charge = ReadNpyPayload<int8_t>(
        dir / "semantic_formal_charge.npy", "|i1", "'shape': (3,)");
    const auto ring = ReadNpyPayload<uint8_t>(
        dir / "semantic_ring_position.npy", "|u1", "'shape': (3,)");
    const auto locant = ReadNpyPayload<uint8_t>(
        dir / "semantic_locant.npy", "|u1", "'shape': (3,)");
    const auto acceptor = ReadNpyPayload<uint8_t>(
        dir / "enrichment_acceptor_class.npy", "|u1", "'shape': (3,)");
    const auto hybrid = ReadNpyPayload<uint8_t>(
        dir / "enrichment_hybridisation_class.npy", "|u1", "'shape': (3,)");
    ASSERT_EQ(parent.size(), 3u);
    EXPECT_EQ(parent[0], 0u);

    for (size_t i = 0; i < 3; ++i) {
        const AtomSemanticTable& sem = protein->LegacyAmber().SemanticAt(i);
        EXPECT_EQ(polar[i], static_cast<uint8_t>(sem.polar_h));
        EXPECT_EQ(planar[i], static_cast<uint8_t>(sem.planar_group));
        EXPECT_EQ(charge[i], sem.formal_charge);
        EXPECT_EQ(ring[i], static_cast<uint8_t>(
                               sem.ring_position.primary.position));
        EXPECT_EQ(locant[i], static_cast<uint8_t>(sem.locant));
    }
    // OD1/OD2 are typed sidechain carboxylate acceptors and every member of
    // the carboxylate planar group projects to sp2.
    EXPECT_EQ(acceptor[1], 3u);
    EXPECT_EQ(acceptor[2], 3u);
    EXPECT_EQ(hybrid[0], static_cast<uint8_t>(Hybridisation::sp2));
    EXPECT_EQ(hybrid[1], static_cast<uint8_t>(Hybridisation::sp2));
    EXPECT_EQ(hybrid[2], static_cast<uint8_t>(Hybridisation::sp2));

    RemoveFilesAndDir(dir, {
        "enrichment_role.npy", "enrichment_hybridisation.npy",
        "enrichment_flags.npy", "enrichment_parent_is_sp2.npy",
        "semantic_polar_h_kind.npy", "semantic_planar_group_kind.npy",
        "semantic_formal_charge.npy", "semantic_ring_position.npy",
        "semantic_locant.npy", "enrichment_donor_class.npy",
        "enrichment_acceptor_class.npy",
        "enrichment_hybridisation_class.npy"});
}

TEST(RawDiagnosticEmission, EnrichmentParentSp2CallsProductionCompute) {
    auto protein = BuildAlaHydrogenProtein();
    auto& conf = protein->Conformation();
    auto enrichment = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrichment, nullptr);

    // The named fixture resolves parents through the covalent topology:
    // backbone H -> peptide N (sp2), while HA/HB1 -> aliphatic C (sp3).
    ASSERT_EQ(protein->AtomAt(1).parent_atom_index, 0u);
    ASSERT_EQ(protein->AtomAt(3).parent_atom_index, 2u);
    ASSERT_EQ(protein->AtomAt(7).parent_atom_index, 6u);
    EXPECT_TRUE(conf.AtomAt(1).parent_is_sp2);
    EXPECT_FALSE(conf.AtomAt(3).parent_is_sp2);
    EXPECT_FALSE(conf.AtomAt(7).parent_is_sp2);

    const fs::path dir = TempDir("enrichment_parent_sp2_production");
    ASSERT_EQ(enrichment->WriteFeatures(conf, dir.string()), 12);
    const auto parent = ReadNpyPayload<uint8_t>(
        dir / "enrichment_parent_is_sp2.npy", "|u1", "'shape': (8,)");
    ASSERT_EQ(parent.size(), 8u);
    EXPECT_EQ(parent[1], 1u);
    EXPECT_EQ(parent[3], 0u);
    EXPECT_EQ(parent[7], 0u);

    RemoveFilesAndDir(dir, {
        "enrichment_role.npy", "enrichment_hybridisation.npy",
        "enrichment_flags.npy", "enrichment_parent_is_sp2.npy",
        "semantic_polar_h_kind.npy", "semantic_planar_group_kind.npy",
        "semantic_formal_charge.npy", "semantic_ring_position.npy",
        "semantic_locant.npy", "enrichment_donor_class.npy",
        "enrichment_acceptor_class.npy",
        "enrichment_hybridisation_class.npy"});
}

TEST(RawDiagnosticEmission, WaterHBondGeometryUsesProductionKernelAndWriter) {
    static_assert(
        water_hbond_geometry_detail::kCandidateHeavyDistance_A == 3.5,
        "A27 raw-candidate shell is a frozen producer contract");
    const auto ideal = water_hbond_geometry_detail::EvaluateGeometry(
        Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0),
        Vec3(2.0, 0.0, 0.0));
    EXPECT_DOUBLE_EQ(ideal.donor_heavy_acceptor_distance_A, 2.0);
    EXPECT_DOUBLE_EQ(ideal.h_acceptor_distance_A, 1.0);
    EXPECT_DOUBLE_EQ(ideal.angle_deg, 180.0);
    EXPECT_TRUE(ideal.passes_geometry);
    const auto bent = water_hbond_geometry_detail::EvaluateGeometry(
        Vec3(0.0, 0.0, 0.0), Vec3(1.0, 0.0, 0.0),
        Vec3(1.0, 1.0, 0.0));
    EXPECT_DOUBLE_EQ(bent.angle_deg, 90.0);
    EXPECT_FALSE(bent.passes_geometry);

    auto protein = BuildAspCarboxylateProtein();
    auto& conf = protein->Conformation();
    SolventEnvironment solvent;
    WaterMolecule water;
    water.O_pos = Vec3(2.80, 0.0, 0.0);
    water.H1_pos = Vec3(1.80, 0.0, 0.0);
    water.H2_pos = Vec3(3.10, 0.8, 0.0);
    water.O_charge = -0.834;
    water.H_charge = 0.417;
    solvent.waters.push_back(water);
    solvent.water_O_positions.push_back(water.O_pos);

    auto result = WaterHBondGeometryResult::Compute(conf, solvent);
    ASSERT_NE(result, nullptr);
    const fs::path dir = TempDir("water_hbond_geometry_emit");
    ASSERT_EQ(result->WriteFeatures(conf, dir.string()), 3);
    const auto candidates = ReadNpyPayload<double>(
        dir / "water_hbond_candidates.npy", "<f8", "'shape': (1,16)");
    ASSERT_EQ(candidates.size(), 16u);
    EXPECT_EQ(candidates[0], 1.0);  // OD1 atom
    EXPECT_EQ(candidates[3], 1.0);  // water donates
    EXPECT_EQ(candidates[4], 3.0);  // typed carboxylate acceptor
    EXPECT_NEAR(candidates[5], 1.55, 1e-12);
    EXPECT_NEAR(candidates[6], 0.55, 1e-12);
    EXPECT_DOUBLE_EQ(candidates[8], 180.0);
    EXPECT_EQ(candidates[15], 1.0);
    const auto counts = ReadNpyPayload<int32_t>(
        dir / "water_hbond_counts.npy", "<i4", "'shape': (3,6)");
    ASSERT_EQ(counts.size(), 18u);
    EXPECT_EQ(counts[1*6 + 0], 1);
    EXPECT_EQ(counts[1*6 + 1], 1);

    RemoveFilesAndDir(dir, {"water_hbond_candidates.npy",
                            "water_hbond_counts.npy",
                            "water_hbond_nearest.npy"});
}

TEST(RawDiagnosticEmission,
     WaterHBondNearestPassModeUsesFartherPassingCandidate) {
    auto protein = BuildAspCarboxylateProtein();
    auto& conf = protein->Conformation();

    SolventEnvironment solvent;
    WaterMolecule nearer_failing;
    nearer_failing.O_pos = Vec3(3.0, 0.0, 0.0);
    nearer_failing.H1_pos = Vec3(3.0, 1.0, 0.0);
    nearer_failing.H2_pos = Vec3(3.8, -0.6, 0.0);
    nearer_failing.O_charge = -0.834;
    nearer_failing.H_charge = 0.417;
    solvent.waters.push_back(nearer_failing);
    solvent.water_O_positions.push_back(nearer_failing.O_pos);

    WaterMolecule farther_passing;
    farther_passing.O_pos = Vec3(4.0, 0.0, 0.0);
    farther_passing.H1_pos = Vec3(3.0, 0.0, 0.0);
    farther_passing.H2_pos = Vec3(4.5, 0.8, 0.0);
    farther_passing.O_charge = -0.834;
    farther_passing.H_charge = 0.417;
    solvent.waters.push_back(farther_passing);
    solvent.water_O_positions.push_back(farther_passing.O_pos);

    auto result = WaterHBondGeometryResult::Compute(conf, solvent);
    ASSERT_NE(result, nullptr);
    const fs::path dir = TempDir("water_hbond_nearest_pass_mode");
    ASSERT_EQ(result->WriteFeatures(conf, dir.string()), 3);

    const auto candidates = ReadNpyPayload<double>(
        dir / "water_hbond_candidates.npy", "<f8", "'shape': (2,16)");
    ASSERT_EQ(candidates.size(), 32u);
    EXPECT_EQ(candidates[0*16 + 0], 1.0);  // typed OD1 acceptor
    EXPECT_EQ(candidates[0*16 + 2], 0.0);  // nearer water
    EXPECT_DOUBLE_EQ(candidates[0*16 + 5], 1.75);
    EXPECT_LT(candidates[0*16 + 8], 120.0);
    EXPECT_EQ(candidates[0*16 + 15], 0.0);  // fails angle
    EXPECT_EQ(candidates[1*16 + 0], 1.0);
    EXPECT_EQ(candidates[1*16 + 2], 1.0);  // farther water
    EXPECT_DOUBLE_EQ(candidates[1*16 + 5], 2.75);
    EXPECT_DOUBLE_EQ(candidates[1*16 + 8], 180.0);
    EXPECT_EQ(candidates[1*16 + 15], 1.0);  // passes geometry

    const auto counts = ReadNpyPayload<int32_t>(
        dir / "water_hbond_counts.npy", "<i4", "'shape': (3,6)");
    ASSERT_EQ(counts.size(), 18u);
    EXPECT_EQ(counts[1*6 + 0], 2);  // water-donor candidates
    EXPECT_EQ(counts[1*6 + 1], 1);  // water-donor passes
    EXPECT_EQ(counts[1*6 + 4], 0);  // overall-nearest water still index 0
    EXPECT_EQ(counts[1*6 + 5], 1);  // nearest passing mode: water donates

    const auto nearest = ReadNpyPayload<double>(
        dir / "water_hbond_nearest.npy", "<f8", "'shape': (3,8)");
    ASSERT_EQ(nearest.size(), 24u);
    EXPECT_DOUBLE_EQ(nearest[1*8 + 0], 1.75);
    EXPECT_EQ(nearest[1*8 + 3], 1.0);
    EXPECT_EQ(nearest[1*8 + 4], 0.0);
    EXPECT_EQ(nearest[1*8 + 5], 0.0);  // overall nearest still fails
    EXPECT_EQ(nearest[1*8 + 6], 2.0);
    EXPECT_EQ(nearest[1*8 + 7], 1.0);

    RemoveFilesAndDir(dir, {"water_hbond_candidates.npy",
                            "water_hbond_counts.npy",
                            "water_hbond_nearest.npy"});
}

TEST(RawDiagnosticEmission,
     WaterHBondGeometryProteinDonorUsesFrozenBoundaryAndZeroRows) {
    auto protein = BuildAlaHydrogenProtein();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(protein->LegacyAmber().SemanticAt(1).IsPolarH());
    ASSERT_EQ(protein->AtomAt(1).parent_atom_index, 0u);

    SolventEnvironment solvent;
    WaterMolecule boundary;
    boundary.O_pos = Vec3(3.5, 0.0, 0.0);  // exactly 3.5 A from donor N
    boundary.H1_pos = boundary.O_pos + Vec3(0.7, 0.4, 0.0);
    boundary.H2_pos = boundary.O_pos + Vec3(-0.7, 0.4, 0.0);
    solvent.waters.push_back(boundary);
    solvent.water_O_positions.push_back(boundary.O_pos);

    WaterMolecule outside = boundary;
    outside.O_pos = Vec3(
        std::nextafter(
            water_hbond_geometry_detail::kCandidateHeavyDistance_A,
            std::numeric_limits<double>::infinity()),
        0.0, 0.0);
    outside.H1_pos = outside.O_pos + Vec3(0.7, 0.4, 0.0);
    outside.H2_pos = outside.O_pos + Vec3(-0.7, 0.4, 0.0);
    solvent.waters.push_back(outside);
    solvent.water_O_positions.push_back(outside.O_pos);

    auto result = WaterHBondGeometryResult::Compute(conf, solvent);
    ASSERT_NE(result, nullptr);
    const fs::path dir = TempDir("water_hbond_protein_donor_boundary");
    ASSERT_EQ(result->WriteFeatures(conf, dir.string()), 3);

    const auto candidates = ReadNpyPayload<double>(
        dir / "water_hbond_candidates.npy", "<f8", "'shape': (1,16)");
    ASSERT_EQ(candidates.size(), 16u);
    EXPECT_EQ(candidates[0], 1.0);  // backbone amide H
    EXPECT_EQ(candidates[2], 0.0);  // exact-boundary water; outside omitted
    EXPECT_EQ(candidates[3], 2.0);  // protein donates to water
    EXPECT_GT(candidates[4], 0.0);  // typed PolarHKind, not a name rule
    EXPECT_DOUBLE_EQ(candidates[5], 2.5);  // H...water-O
    EXPECT_DOUBLE_EQ(candidates[6], 2.5);
    EXPECT_DOUBLE_EQ(candidates[7], 3.5);  // donor-heavy...water-O
    EXPECT_DOUBLE_EQ(candidates[8], 180.0);
    EXPECT_TRUE(std::isnan(candidates[12]));
    EXPECT_TRUE(std::isnan(candidates[13]));
    EXPECT_TRUE(std::isnan(candidates[14]));
    EXPECT_EQ(candidates[15], 1.0);

    const auto counts = ReadNpyPayload<int32_t>(
        dir / "water_hbond_counts.npy", "<i4", "'shape': (8,6)");
    ASSERT_EQ(counts.size(), 48u);
    EXPECT_EQ(counts[1*6 + 0], 0);
    EXPECT_EQ(counts[1*6 + 1], 0);
    EXPECT_EQ(counts[1*6 + 2], 1);
    EXPECT_EQ(counts[1*6 + 3], 1);
    EXPECT_EQ(counts[1*6 + 4], 0);
    EXPECT_EQ(counts[1*6 + 5], 2);

    RemoveFilesAndDir(dir, {"water_hbond_candidates.npy",
                            "water_hbond_counts.npy",
                            "water_hbond_nearest.npy"});

    SolventEnvironment empty_solvent;
    auto empty = WaterHBondGeometryResult::Compute(conf, empty_solvent);
    ASSERT_NE(empty, nullptr);
    const fs::path empty_dir = TempDir("water_hbond_zero_candidates");
    ASSERT_EQ(empty->WriteFeatures(conf, empty_dir.string()), 3);
    const auto empty_candidates = ReadNpyPayload<double>(
        empty_dir / "water_hbond_candidates.npy", "<f8",
        "'shape': (0,16)");
    EXPECT_TRUE(empty_candidates.empty());
    const auto empty_counts = ReadNpyPayload<int32_t>(
        empty_dir / "water_hbond_counts.npy", "<i4", "'shape': (8,6)");
    ASSERT_EQ(empty_counts.size(), 48u);
    for (size_t ai = 0; ai < 8; ++ai) {
        EXPECT_EQ(empty_counts[ai*6 + 0], 0);
        EXPECT_EQ(empty_counts[ai*6 + 1], 0);
        EXPECT_EQ(empty_counts[ai*6 + 2], 0);
        EXPECT_EQ(empty_counts[ai*6 + 3], 0);
        EXPECT_EQ(empty_counts[ai*6 + 4], -1);
        EXPECT_EQ(empty_counts[ai*6 + 5], 0);
    }
    const auto empty_nearest = ReadNpyPayload<double>(
        empty_dir / "water_hbond_nearest.npy", "<f8", "'shape': (8,8)");
    ASSERT_EQ(empty_nearest.size(), 64u);
    for (size_t ai = 0; ai < 8; ++ai) {
        EXPECT_TRUE(std::isnan(empty_nearest[ai*8 + 0]));
        EXPECT_EQ(empty_nearest[ai*8 + 6], 0.0);
        EXPECT_EQ(empty_nearest[ai*8 + 7], 0.0);
    }
    RemoveFilesAndDir(empty_dir, {"water_hbond_candidates.npy",
                                  "water_hbond_counts.npy",
                                  "water_hbond_nearest.npy"});
}

}  // namespace
