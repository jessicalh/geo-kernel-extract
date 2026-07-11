#include <gtest/gtest.h>

#include "CalculatorConfig.h"
#include "EnrichmentResult.h"
#include "GeometryResult.h"
#include "LegacyAmberTopology.h"
#include "MolecularGraphResult.h"
#include "MutationDeltaResult.h"
#include "NpyWriter.h"
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
#include <string>
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

TEST(RawDiagnosticEmission, SasaFractionUsesProductionWriter) {
    auto protein = BuildUnnamedProtein({Element::C}, {Vec3::Zero()});
    auto& conf = protein->Conformation();
    const double probe = CalculatorConfig::Get("sasa_probe_radius");
    const double radius = BondiVdwRadius(Element::C) + probe;
    conf.MutableAtomAt(0).atom_sasa = 4.0 * PI * radius * radius;

    const fs::path dir = TempDir("sasa_fraction_emit");
    SasaResult production_writer;
    ASSERT_EQ(production_writer.WriteFeatures(conf, dir.string()), 3);
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

TEST(RawDiagnosticEmission, MolecularGraphWriterPreservesRawFields) {
    auto protein = BuildUnnamedProtein(
        {Element::C, Element::N},
        {Vec3(0.0, 0.0, 0.0), Vec3(3.0, 0.0, 0.0)});
    auto& conf = protein->Conformation();
    auto& a = conf.MutableAtomAt(0);
    a.graph_dist_ring = 2;
    a.graph_dist_N = 1;
    a.graph_dist_O = -1;
    a.n_pi_bonds_3 = 4;
    a.is_conjugated = true;
    a.nearest_ring_atom_index = 17;
    a.eneg_sum_1 = 3.04;
    a.eneg_sum_2 = 5.55;
    a.bfs_decay = 0.6065306597126334;

    const fs::path dir = TempDir("molecular_graph_emit");
    MolecularGraphResult production_writer;
    ASSERT_EQ(production_writer.WriteFeatures(conf, dir.string()), 3);
    const auto integer = ReadNpyPayload<int32_t>(
        dir / "molecular_graph_int.npy", "<i4", "'shape': (2,6)");
    const auto floating = ReadNpyPayload<double>(
        dir / "molecular_graph_float.npy", "<f8", "'shape': (2,3)");
    const auto compatibility = ReadNpyPayload<double>(
        dir / "molecular_graph.npy", "<f8", "'shape': (2,9)");
    ASSERT_EQ(integer.size(), 12u);
    ASSERT_EQ(floating.size(), 6u);
    ASSERT_EQ(compatibility.size(), 18u);
    EXPECT_EQ(integer[0], 2);
    EXPECT_EQ(integer[1], 1);
    EXPECT_EQ(integer[2], -1);
    EXPECT_EQ(integer[3], 4);
    EXPECT_EQ(integer[4], 1);
    EXPECT_EQ(integer[5], 17);
    EXPECT_DOUBLE_EQ(floating[0], 3.04);
    EXPECT_DOUBLE_EQ(floating[1], 5.55);
    EXPECT_DOUBLE_EQ(floating[2], 0.6065306597126334);
    EXPECT_DOUBLE_EQ(compatibility[7], 17.0);

    RemoveFilesAndDir(dir, {"molecular_graph_int.npy",
                            "molecular_graph_float.npy",
                            "molecular_graph.npy"});
}

TEST(RawDiagnosticEmission, MutationGraphRowUsesProductionSchema) {
    MatchedAtomData data;
    data.has_graph_delta = true;
    data.delta_graph_dist_ring = -4;
    data.delta_bfs_decay = 0.375;
    data.delta_is_conjugated = 1;
    const auto row = mutation_delta_detail::PackDeltaGraphRow(true, &data);
    EXPECT_EQ(row[0], 1.0);
    EXPECT_EQ(row[1], 1.0);
    EXPECT_EQ(row[2], -4.0);
    EXPECT_DOUBLE_EQ(row[3], 0.375);
    EXPECT_EQ(row[4], 1.0);

    const auto unmatched =
        mutation_delta_detail::PackDeltaGraphRow(false, nullptr);
    EXPECT_EQ(unmatched, (std::array<double, 5>{0, 0, 0, 0, 0}));
}

TEST(RawDiagnosticEmission, EnrichmentRawSemanticsAndSp2FlagReadBack) {
    auto protein = BuildAspCarboxylateProtein();
    ASSERT_TRUE(protein->LegacyAmber().HasAtomSemantic());
    auto& conf = protein->Conformation();
    auto enrichment = EnrichmentResult::Compute(conf);
    ASSERT_NE(enrichment, nullptr);

    // A31 is an emit-only finding: pin the stored bit independently of the
    // chemistry classifier, then ensure the production writer preserves it.
    conf.MutableAtomAt(0).parent_is_sp2 = true;
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
    EXPECT_EQ(parent[0], 1u);

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

TEST(RawDiagnosticEmission, WaterHBondGeometryUsesProductionKernelAndWriter) {
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

}  // namespace
