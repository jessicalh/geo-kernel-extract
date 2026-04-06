#include <gtest/gtest.h>
#include "Protein.h"
#include "ConformationResult.h"
#include "ProteinConformation.h"

using namespace nmr;

// Build a minimal protein for testing
static std::unique_ptr<Protein> MakeTinyProtein() {
    auto pp = std::make_unique<Protein>();
    Protein& p = *pp;

    // Two residues: ALA and PHE
    Residue ala;
    ala.type = AminoAcid::ALA;
    ala.sequence_number = 1;
    ala.chain_id = "A";

    Residue phe;
    phe.type = AminoAcid::PHE;
    phe.sequence_number = 2;
    phe.chain_id = "A";

    size_t ri_ala = p.AddResidue(ala);
    size_t ri_phe = p.AddResidue(phe);

    // Add backbone atoms for ALA
    auto n = Atom::Create(Element::N); n->pdb_atom_name = "N"; n->residue_index = ri_ala;
    auto ca = Atom::Create(Element::C); ca->pdb_atom_name = "CA"; ca->residue_index = ri_ala;
    auto c = Atom::Create(Element::C); c->pdb_atom_name = "C"; c->residue_index = ri_ala;
    auto o = Atom::Create(Element::O); o->pdb_atom_name = "O"; o->residue_index = ri_ala;

    size_t n_idx = p.AddAtom(std::move(n));
    size_t ca_idx = p.AddAtom(std::move(ca));
    size_t c_idx = p.AddAtom(std::move(c));
    size_t o_idx = p.AddAtom(std::move(o));

    p.MutableResidueAt(ri_ala).atom_indices = {n_idx, ca_idx, c_idx, o_idx};

    // Add ring atoms for PHE
    const char* ring_names[] = {"CG", "CD1", "CE1", "CZ", "CE2", "CD2"};
    std::vector<size_t> phe_indices;
    for (auto name : ring_names) {
        auto a = Atom::Create(Element::C); a->pdb_atom_name = name; a->residue_index = ri_phe;
        phe_indices.push_back(p.AddAtom(std::move(a)));
    }
    p.MutableResidueAt(ri_phe).atom_indices = phe_indices;

    p.CacheResidueBackboneIndices();
    p.DetectAromaticRings();

    return pp;
}


TEST(ObjectModel, ProteinOwnsResidues) {
    auto p = MakeTinyProtein();
    EXPECT_EQ(p->ResidueCount(), 2);
    EXPECT_EQ(p->ResidueAt(0).type, AminoAcid::ALA);
    EXPECT_EQ(p->ResidueAt(1).type, AminoAcid::PHE);
}

TEST(ObjectModel, ProteinOwnsAtoms) {
    auto p = MakeTinyProtein();
    EXPECT_EQ(p->AtomCount(), 10);  // 4 ALA + 6 PHE
    EXPECT_EQ(p->AtomAt(0).element, Element::N);
    EXPECT_EQ(p->AtomAt(4).element, Element::C);  // PHE CG
}

TEST(ObjectModel, RingDetected) {
    auto p = MakeTinyProtein();
    EXPECT_EQ(p->RingCount(), 1);
    EXPECT_EQ(p->RingAt(0).type_index, RingTypeIndex::PheBenzene);
    EXPECT_DOUBLE_EQ(p->RingAt(0).Intensity(), -12.0);
}

TEST(ObjectModel, BackboneIndicesCached) {
    auto p = MakeTinyProtein();
    const Residue& ala = p->ResidueAt(0);
    EXPECT_NE(ala.N, Residue::NONE);
    EXPECT_NE(ala.CA, Residue::NONE);
    EXPECT_NE(ala.C, Residue::NONE);
    EXPECT_NE(ala.O, Residue::NONE);
}

TEST(ObjectModel, ResidueQueryMethods) {
    auto p = MakeTinyProtein();
    EXPECT_FALSE(p->ResidueAt(0).IsAromatic());
    EXPECT_TRUE(p->ResidueAt(1).IsAromatic());
}

TEST(ObjectModel, ConformationCreation) {
    auto p = MakeTinyProtein();
    std::vector<Vec3> positions(p->AtomCount(), Vec3::Zero());
    for (size_t i = 0; i < positions.size(); ++i)
        positions[i] = Vec3(i * 1.5, 0.0, 0.0);

    p->AddCrystalConformation(positions, 1.5, 0.2, 293.0, "TEST");
    EXPECT_TRUE(p->HasCrystalConformation());
    EXPECT_EQ(p->ConformationCount(), 1);

    auto& conf = p->CrystalConf();
    EXPECT_EQ(conf.AtomCount(), p->AtomCount());
    EXPECT_EQ(&conf.ProteinRef(), p.get());
}

TEST(ObjectModel, ConformationAtomPosition) {
    auto p = MakeTinyProtein();
    std::vector<Vec3> positions(p->AtomCount());
    positions[0] = Vec3(1.0, 2.0, 3.0);
    positions[1] = Vec3(4.0, 5.0, 6.0);
    for (size_t i = 2; i < positions.size(); ++i)
        positions[i] = Vec3::Zero();

    p->AddCrystalConformation(positions, 0.0, 0.0, 0.0, "T");
    auto& conf = p->CrystalConf();

    EXPECT_NEAR(conf.AtomAt(0).Position().x(), 1.0, 1e-15);
    EXPECT_NEAR(conf.AtomAt(1).Position().y(), 5.0, 1e-15);
}

TEST(ObjectModel, ConformationAtomDefaultFields) {
    auto p = MakeTinyProtein();
    std::vector<Vec3> positions(p->AtomCount(), Vec3::Zero());
    p->AddCrystalConformation(positions, 0.0, 0.0, 0.0, "T");
    auto& conf = p->CrystalConf();

    const auto& ca = conf.AtomAt(0);
    EXPECT_EQ(ca.role, AtomRole::Unknown);
    EXPECT_DOUBLE_EQ(ca.partial_charge, 0.0);
    EXPECT_DOUBLE_EQ(ca.total_B_field.norm(), 0.0);
    EXPECT_DOUBLE_EQ(ca.bs_shielding_contribution.T0, 0.0);
    EXPECT_DOUBLE_EQ(ca.demo_nearest_ring_distance, 0.0);
}


// A trivial ConformationResult for testing the attach mechanism
class TrivialResult : public ConformationResult {
public:
    std::string Name() const override { return "TrivialResult"; }
    std::vector<std::type_index> Dependencies() const override { return {}; }
};

class DependentResult : public ConformationResult {
public:
    std::string Name() const override { return "DependentResult"; }
    std::vector<std::type_index> Dependencies() const override {
        return { std::type_index(typeid(TrivialResult)) };
    }
};


TEST(ObjectModel, AttachAndAccessResult) {
    auto p = MakeTinyProtein();
    std::vector<Vec3> positions(p->AtomCount(), Vec3::Zero());
    p->AddCrystalConformation(positions, 0.0, 0.0, 0.0, "T");
    auto& conf = p->CrystalConf();

    EXPECT_FALSE(conf.HasResult<TrivialResult>());
    auto r = std::make_unique<TrivialResult>();
    EXPECT_TRUE(conf.AttachResult(std::move(r)));
    EXPECT_TRUE(conf.HasResult<TrivialResult>());
}

TEST(ObjectModel, SingletonGuarantee) {
    auto p = MakeTinyProtein();
    std::vector<Vec3> positions(p->AtomCount(), Vec3::Zero());
    p->AddCrystalConformation(positions, 0.0, 0.0, 0.0, "T");
    auto& conf = p->CrystalConf();

    conf.AttachResult(std::make_unique<TrivialResult>());
    // Second attach of same type should fail
    EXPECT_FALSE(conf.AttachResult(std::make_unique<TrivialResult>()));
}

TEST(ObjectModel, DependencyCheckPass) {
    auto p = MakeTinyProtein();
    std::vector<Vec3> positions(p->AtomCount(), Vec3::Zero());
    p->AddCrystalConformation(positions, 0.0, 0.0, 0.0, "T");
    auto& conf = p->CrystalConf();

    conf.AttachResult(std::make_unique<TrivialResult>());
    EXPECT_TRUE(conf.AttachResult(std::make_unique<DependentResult>()));
}

TEST(ObjectModel, DependencyCheckFail) {
    auto p = MakeTinyProtein();
    std::vector<Vec3> positions(p->AtomCount(), Vec3::Zero());
    p->AddCrystalConformation(positions, 0.0, 0.0, 0.0, "T");
    auto& conf = p->CrystalConf();

    // DependentResult requires TrivialResult, which is not attached
    EXPECT_FALSE(conf.AttachResult(std::make_unique<DependentResult>()));
}

TEST(ObjectModel, BuildContext) {
    auto p = MakeTinyProtein();
    auto ctx = std::make_unique<ProteinBuildContext>();
    ctx->pdb_source = "1ubq.pdb";
    ctx->force_field = "ff14SB";
    p->SetBuildContext(std::move(ctx));

    EXPECT_EQ(p->BuildContext().pdb_source, "1ubq.pdb");
    EXPECT_EQ(p->BuildContext().force_field, "ff14SB");
}
