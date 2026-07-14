#include <gtest/gtest.h>

#include "Atom.h"
#include "ChargeAssignmentResult.h"
#include "ChargeSource.h"
#include "MopacResult.h"
#include "OperationRunner.h"
#include "Protein.h"
#include "Residue.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <string>
#include <unistd.h>
#include <vector>

using namespace nmr;
namespace fs = std::filesystem;

namespace {

std::unique_ptr<Protein> BuildUnnamedMolecule(
        const std::vector<Element>& elements,
        const std::vector<Vec3>& positions,
        const std::string& description) {
    auto protein = std::make_unique<Protein>();

    Residue residue;
    residue.type = AminoAcid::Unknown;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    protein->AddResidue(std::move(residue));

    for (Element element : elements) {
        auto atom = Atom::Create(element);
        atom->residue_index = 0;
        const size_t atom_index = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(0).atom_indices.push_back(atom_index);
    }

    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, description);
    return protein;
}
std::unique_ptr<Protein> BuildWater() {
    return BuildUnnamedMolecule(
        {Element::O, Element::H, Element::H},
        {
            Vec3(0.0, 0.0, 0.0),
            Vec3(0.9572, 0.0, 0.0),
            Vec3(-0.2399872, 0.927297, 0.0),
        },
        "direct MOPAC water probe");
}

std::unique_ptr<Protein> BuildDissociatedHydrogen() {
    return BuildUnnamedMolecule(
        {Element::H, Element::H},
        {Vec3(0.0, 0.0, 0.0), Vec3(100.0, 0.0, 0.0)},
        "libmopac crash-containment probe");
}

std::unique_ptr<Protein> BuildRunnerFixture() {
    auto protein = std::make_unique<Protein>();

    auto atom_n = Atom::Create(Element::N);
    auto atom_c = Atom::Create(Element::C);
    auto atom_o = Atom::Create(Element::O);
    auto atom_h = Atom::Create(Element::H);
    atom_n->pdb_atom_name = "N";
    atom_c->pdb_atom_name = "C";
    atom_o->pdb_atom_name = "O";
    atom_h->pdb_atom_name = "H";
    atom_n->residue_index = 0;
    atom_c->residue_index = 0;
    atom_o->residue_index = 0;
    atom_h->residue_index = 0;
    protein->AddAtom(std::move(atom_n));
    protein->AddAtom(std::move(atom_c));
    protein->AddAtom(std::move(atom_o));
    protein->AddAtom(std::move(atom_h));

    Residue residue;
    residue.type = AminoAcid::ALA;
    residue.sequence_number = 1;
    residue.chain_id = "A";
    residue.atom_indices = {0, 1, 2, 3};
    residue.N = 0;
    residue.C = 1;
    residue.O = 2;
    residue.H = 3;
    protein->AddResidue(std::move(residue));

    const std::vector<Vec3> positions = {
        Vec3(0.0, 1.45, 0.0),
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.2, 0.0, 0.0),
        Vec3(0.6, 3.0, 0.0),
    };
    protein->FinalizeConstruction(positions);
    protein->AddConformation(positions, "requested MOPAC hard-abort probe");
    return protein;
}

bool ContainsAttached(const RunResult& result, const std::string& name) {
    return std::find(result.attached.begin(), result.attached.end(), name)
        != result.attached.end();
}

std::string ShellQuote(const std::string& value) {
    std::string quoted = "'";
    for (char c : value) {
        if (c == '\'') quoted += "'\\''";
        else quoted += c;
    }
    quoted += "'";
    return quoted;
}

struct TemporaryDirectory {
    fs::path path;

    static void Remove(const fs::path& directory) {
        std::error_code error;
        if (fs::exists(directory, error)) {
            for (const auto& entry : fs::directory_iterator(directory)) {
                fs::remove(entry.path(), error);
                error.clear();
            }
            fs::remove(directory, error);
        }
    }

    explicit TemporaryDirectory(const std::string& label) {
        path = fs::temp_directory_path() /
            (label + "_" + std::to_string(::getpid()));
        Remove(path);
        std::error_code error;
        fs::create_directories(path, error);
        EXPECT_FALSE(error) << error.message();
    }

    ~TemporaryDirectory() {
        Remove(path);
    }
};

}  // namespace

TEST(MopacDirectApi, WaterProbeMatchesProductionAndConsumedSurface) {
    auto protein = BuildWater();
    auto& conf = protein->Conformation();
    ASSERT_EQ(protein->BondCount(), 2u);

    std::string error;
    auto result = MopacResult::Compute(conf, 0, 1, &error);
    ASSERT_NE(result, nullptr) << error;
    EXPECT_TRUE(error.empty());
    EXPECT_TRUE(result->Dependencies().empty());

    // The consumed surface preserves the numeric values recovered from the
    // legacy fixed-precision text fields. The exact API values are emitted in
    // explicit full-precision arrays and checked against the independent
    // mopac_feature_probe.py archive in the schema test below.
    constexpr std::array<double, 3> expected_charges = {
        -0.648152,
        0.324018,
        0.324134,
    };
    constexpr std::array<double, 3> expected_s = {
        1.81616,
        0.67598,
        0.67587,
    };
    constexpr std::array<double, 3> expected_p = {
        4.83199,
        0.0,
        0.0,
    };
    constexpr std::array<double, 3> expected_api_valencies = {
        1.788263447550035,
        0.8950121991514495,
        0.8949370367721683,
    };
    constexpr std::array<double, 3> expected_consumed_valencies = {
        1.788,
        0.894,
        0.894,
    };

    for (size_t atom = 0; atom < conf.AtomCount(); ++atom) {
        EXPECT_NEAR(result->ChargeAt(atom), expected_charges[atom], 1e-12);
        EXPECT_NEAR(result->SPopAt(atom), expected_s[atom], 1e-12);
        EXPECT_NEAR(result->PPopAt(atom), expected_p[atom], 1e-12);
        EXPECT_NEAR(result->ValencyAt(atom),
                    expected_consumed_valencies[atom], 1e-12);

        const auto& stored = conf.AtomAt(atom);
        EXPECT_DOUBLE_EQ(stored.mopac_charge, result->ChargeAt(atom));
        EXPECT_DOUBLE_EQ(stored.mopac_s_pop, result->SPopAt(atom));
        EXPECT_DOUBLE_EQ(stored.mopac_p_pop, result->PPopAt(atom));
        EXPECT_DOUBLE_EQ(stored.mopac_valency, result->ValencyAt(atom));
        for (size_t i = 1; i < stored.mopac_bond_neighbours.size(); ++i) {
            EXPECT_GE(stored.mopac_bond_neighbours[i - 1].wiberg_order,
                      stored.mopac_bond_neighbours[i].wiberg_order);
        }
    }

    // The distinct libmopac CSC-diagonal valencies are retained separately
    // and checked by WaterEmitsCompleteLoadableSchema against the probe.
    EXPECT_NE(expected_api_valencies[0], expected_consumed_valencies[0]);

    EXPECT_DOUBLE_EQ(result->BondOrder(0, 1), 0.894);
    EXPECT_DOUBLE_EQ(result->BondOrder(0, 2), 0.894);
    EXPECT_DOUBLE_EQ(result->BondOrder(1, 2), 0.0);
    EXPECT_EQ(result->AllBondOrders().size(), 4u);
    for (size_t bond = 0; bond < protein->BondCount(); ++bond) {
        const Bond& topology_bond = protein->BondAt(bond);
        EXPECT_DOUBLE_EQ(
            result->TopologyBondOrder(bond),
            result->BondOrder(topology_bond.atom_index_a,
                              topology_bond.atom_index_b));
    }

    EXPECT_DOUBLE_EQ(result->HeatOfFormation(), -57.79012);
    const Vec3 dipole = result->Dipole();
    EXPECT_DOUBLE_EQ(dipole.x(), 1.314);
    EXPECT_DOUBLE_EQ(dipole.y(), 1.699);
    EXPECT_DOUBLE_EQ(dipole.z(), 0.0);
}

TEST(MopacDirectApi, WaterEmitsCompleteLoadableSchema) {
    auto protein = BuildWater();
    auto& conf = protein->Conformation();
    std::string error;
    auto result = MopacResult::Compute(conf, 0, 1, &error);
    ASSERT_NE(result, nullptr) << error;

    TemporaryDirectory output("mopac_direct_npy");
    ASSERT_EQ(result->WriteFeatures(conf, output.path.string()), 50);

    const fs::path probe_archive = output.path / "water_probe_reference.npz";
    const std::string probe_command =
        "OMP_NUM_THREADS=1 OMP_STACKSIZE=2G " +
        ShellQuote(NMR_MOPAC_PROBE_ENV) + " " +
        ShellQuote(NMR_TEST_PYTHON_EXECUTABLE) + " " +
        ShellQuote(NMR_MOPAC_PROBE_SCRIPT) + " " +
        ShellQuote(NMR_MOPAC_PROBE_INPUT) +
        " --library " + ShellQuote(NMR_MOPAC_PROBE_LIBRARY) +
        " --electronic-features --output " +
        ShellQuote(probe_archive.string());
    ASSERT_EQ(std::system(probe_command.c_str()), 0) << probe_command;

    const std::string validation_command =
        "PYTHONPATH=" + ShellQuote(NMR_TEST_PYTHONPATH) + " " +
        ShellQuote(NMR_TEST_PYTHON_EXECUTABLE) + " " +
        ShellQuote(NMR_MOPAC_NPY_VALIDATOR) + " " +
        ShellQuote(output.path.string()) + " " +
        ShellQuote(probe_archive.string()) + " " +
        ShellQuote(NMR_MOPAC_PROBE_INPUT);
    EXPECT_EQ(std::system(validation_command.c_str()), 0)
        << validation_command;
}

TEST(MopacDirectApi, DissociatedHydrogenCrashIsContained) {
    auto protein = BuildDissociatedHydrogen();
    auto& conf = protein->Conformation();

    std::string error;
    auto result = MopacResult::Compute(conf, 0, 1, &error);
    EXPECT_EQ(result, nullptr);
    EXPECT_NE(error.find("MOPAC worker terminated by signal 6"),
              std::string::npos) << error;

    // The extractor process is still alive and no partial calculator-owned
    // state was committed after the worker's real SIGABRT.
    for (size_t atom = 0; atom < conf.AtomCount(); ++atom) {
        const auto& stored = conf.AtomAt(atom);
        EXPECT_DOUBLE_EQ(stored.mopac_charge, 0.0);
        EXPECT_DOUBLE_EQ(stored.mopac_s_pop, 0.0);
        EXPECT_DOUBLE_EQ(stored.mopac_p_pop, 0.0);
        EXPECT_DOUBLE_EQ(stored.mopac_valency, 0.0);
        EXPECT_TRUE(stored.mopac_bond_neighbours.empty());
    }
}

TEST(MopacDirectApi, RequestedFailureHardAbortsOperationRunner) {
    auto protein = BuildRunnerFixture();
    auto& conf = protein->Conformation();

    std::vector<AtomChargeRadius> charges(4);
    for (auto& charge : charges) {
        charge.partial_charge = 25.0;
        charge.pb_radius = 1.5;
        charge.status = ChargeAssignmentStatus::Matched;
    }
    PreloadedChargeSource source(
        std::move(charges), ForceField::Amber_ff14SB);

    RunOptions options;
    options.charge_source = &source;
    options.net_charge = 100;
    options.skip_dssp = true;
    options.skip_apbs = true;
    options.skip_coulomb = true;

    const RunResult run = OperationRunner::Run(conf, options);
    EXPECT_FALSE(run.Ok());
    EXPECT_NE(run.error.find("PM7/MOZYME/1SCF failed"), std::string::npos)
        << run.error;
    EXPECT_NE(run.error.find("CHECK"), std::string::npos) << run.error;
    EXPECT_TRUE(ContainsAttached(run, "ChargeAssignmentResult"));
    EXPECT_FALSE(ContainsAttached(run, "MopacResult"));
    EXPECT_FALSE(ContainsAttached(run, "BiotSavartResult"));
    EXPECT_TRUE(conf.HasResult<ChargeAssignmentResult>());
    EXPECT_FALSE(conf.HasResult<MopacResult>());

    for (size_t atom = 0; atom < conf.AtomCount(); ++atom) {
        const auto& stored = conf.AtomAt(atom);
        EXPECT_DOUBLE_EQ(stored.mopac_charge, 0.0);
        EXPECT_TRUE(stored.mopac_bond_neighbours.empty());
    }
}
