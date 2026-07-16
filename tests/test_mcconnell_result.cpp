#include "TestEnvironment.h"
#include <gtest/gtest.h>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <fstream>
#include <iterator>
#include <cstdio>
#include <sstream>

#include "McConnellResult.h"
#include "CalculatorConfig.h"
#include "MopacResult.h"
#include "MopacMcConnellResult.h"
#include "MopacMcConnellShieldingTimeSeriesTrajectoryResult.h"
#include "SidechainCarbonylAnisotropyResult.h"
#include "GeometryResult.h"
#include "SpatialIndexResult.h"
#include "Protein.h"
#include "PdbFileReader.h"
#include "OrcaRunLoader.h"
#include "PhysicalConstants.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "DirectionalTestHelpers.h"

#include <filesystem>
#include <Eigen/Geometry>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <unistd.h>
namespace fs = std::filesystem;
using namespace nmr;

namespace {

double MaxAbs(const Mat3& m) {
    double v = 0.0;
    for (int r = 0; r < 3; ++r)
        for (int c = 0; c < 3; ++c)
            v = std::max(v, std::abs(m(r, c)));
    return v;
}

bool SphericalAllZero(const SphericalTensor& st) {
    if (st.T0 != 0.0) return false;
    for (double v : st.T1) if (v != 0.0) return false;
    for (double v : st.T2) if (v != 0.0) return false;
    return true;
}

std::string ReadText(const fs::path& path) {
    std::ifstream in(path);
    return std::string(std::istreambuf_iterator<char>(in),
                       std::istreambuf_iterator<char>());
}

std::unique_ptr<Protein> BuildSyntheticPeptideCOProtein(
        bool degenerate = false,
        const Mat3& rotation = Mat3::Identity()) {
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

    protein->AddAtom(std::move(atom_n));  // 0: amide N
    protein->AddAtom(std::move(atom_c));  // 1: carbonyl C
    protein->AddAtom(std::move(atom_o));  // 2: carbonyl O
    protein->AddAtom(std::move(atom_h));  // 3: target probe

    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.atom_indices = {0, 1, 2, 3};
    res.N = 0;
    res.C = 1;
    res.O = 2;
    res.H = 3;
    protein->AddResidue(res);

    std::vector<Vec3> positions = {
        degenerate ? Vec3(-1.45, 0.0, 0.0) : Vec3(0.0, 1.45, 0.0),
        Vec3(0.0, 0.0, 0.0),
        Vec3(1.2, 0.0, 0.0),
        Vec3(0.6, 3.0, 0.0)
    };
    for (Vec3& pos : positions) pos = rotation * pos;

    protein->FinalizeConstruction(positions);
    protein->AddCrystalConformation(positions, 0, 0, 0, "synthetic_peptide_co");
    return protein;
}

std::unique_ptr<Protein> BuildSyntheticXHCategoryProtein() {
    auto protein = std::make_unique<Protein>();
    for (int ri = 0; ri < 3; ++ri) {
        Residue res;
        res.type = AminoAcid::Unknown;
        res.sequence_number = ri + 1;
        res.chain_id = "A";
        protein->AddResidue(std::move(res));
    }

    std::vector<Vec3> positions;
    auto add = [&](size_t ri, Element element, Vec3 pos) {
        auto atom = Atom::Create(element);
        atom->residue_index = ri;
        const size_t ai = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(ri).atom_indices.push_back(ai);
        positions.push_back(pos);
        return ai;
    };

    const size_t backbone_n = add(0, Element::N, Vec3(0.0, 0.0, 0.0));
    add(0, Element::H, Vec3(1.0, 0.0, 0.0));
    add(1, Element::C, Vec3(10.0, 0.0, 0.0));
    add(1, Element::H, Vec3(11.0, 0.0, 0.0));
    add(2, Element::S, Vec3(20.0, 0.0, 0.0));
    add(2, Element::H, Vec3(21.3, 0.0, 0.0));
    protein->MutableResidueAt(0).N = backbone_n;

    protein->FinalizeConstruction(positions);
    protein->AddCrystalConformation(
        positions, 0.0, 0.0, 0.0, "synthetic_xh_categories");
    return protein;
}

std::unique_ptr<Protein> BuildSyntheticSidechainCOProtein() {
    auto protein = std::make_unique<Protein>();
    std::vector<Vec3> positions;

    auto add_residue = [&](AminoAcid type, int sequence_number) {
        Residue residue;
        residue.type = type;
        residue.sequence_number = sequence_number;
        residue.chain_id = "A";
        return protein->AddResidue(std::move(residue));
    };
    auto add_atom = [&](std::size_t residue_index, const char* name,
                        Element element, const Vec3& position) {
        auto atom = Atom::Create(element);
        atom->pdb_atom_name = name;
        atom->residue_index = residue_index;
        const std::size_t atom_index = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(residue_index).atom_indices.push_back(
            atom_index);
        positions.push_back(position);
        return atom_index;
    };

    // ASN primary amide: one typed SidechainCO source.  The amide N is the
    // second same-planar-group atom that fixes the local plane.
    const std::size_t asn = add_residue(AminoAcid::ASN, 1);
    add_atom(asn, "CG", Element::C, Vec3(0.0, 0.0, 0.0));
    add_atom(asn, "OD1", Element::O, Vec3(1.23, 0.0, 0.0));
    add_atom(asn, "ND2", Element::N, Vec3(-0.60, 1.05, 0.0));

    // ASP carboxylate: two typed SidechainCO sources, each using the other
    // oxygen as its deterministic in-plane reference.
    const std::size_t asp = add_residue(AminoAcid::ASP, 2);
    add_atom(asp, "CG", Element::C, Vec3(5.0, 0.0, 0.0));
    add_atom(asp, "OD1", Element::O, Vec3(6.25, 0.0, 0.0));
    add_atom(asp, "OD2", Element::O, Vec3(4.375, 1.08253175, 0.0));

    // A real peptide C=O in the same synthetic topology is a negative
    // control: source enumeration must exclude it by typed BondCategory.
    const std::size_t ala = add_residue(AminoAcid::ALA, 3);
    const std::size_t ala_n =
        add_atom(ala, "N", Element::N, Vec3(12.0, 1.35, 0.0));
    const std::size_t ala_c =
        add_atom(ala, "C", Element::C, Vec3(12.0, 0.0, 0.0));
    const std::size_t ala_o =
        add_atom(ala, "O", Element::O, Vec3(13.20, 0.0, 0.0));
    protein->MutableResidueAt(ala).N = ala_n;
    protein->MutableResidueAt(ala).C = ala_c;
    protein->MutableResidueAt(ala).O = ala_o;

    protein->FinalizeConstruction(positions);
    protein->AddCrystalConformation(
        positions, 0.0, 0.0, 0.0, "synthetic_sidechain_co");
    return protein;
}

// A chemically much less pathological external-MOPAC forcing fixture than
// the disconnected nine-atom classification fixture above.  Clone the
// neutral TYR59-ASN60-ILE61 window from the checked-in protonated 1UBQ input;
// it retains a typed ASN sidechain carbonyl and two real peptide links while
// keeping each PM7+MOZYME rerun small.
std::unique_ptr<Protein> BuildProtonatedAsnWindowProtein() {
    auto loaded = BuildFromProtonatedPdb(
        nmr::test::TestEnvironment::UbqProtonated());
    if (!loaded.Ok()) return nullptr;

    const Protein& source = *loaded.protein;
    const ProteinConformation& source_conf = source.Conformation();
    auto protein = std::make_unique<Protein>();
    std::vector<Vec3> positions;

    for (std::size_t source_residue_index = 0;
         source_residue_index < source.ResidueCount();
         ++source_residue_index) {
        const Residue& source_residue =
            source.ResidueAt(source_residue_index);
        if (source_residue.chain_id != "A" ||
            source_residue.sequence_number < 59 ||
            source_residue.sequence_number > 61) {
            continue;
        }

        Residue residue;
        residue.type = source_residue.type;
        residue.sequence_number = source_residue.sequence_number;
        residue.chain_id = source_residue.chain_id;
        residue.insertion_code = source_residue.insertion_code;
        residue.protonation_variant_index =
            source_residue.protonation_variant_index;
        residue.protonation_state_resolved =
            source_residue.protonation_state_resolved;
        residue.terminal_state = ResidueTerminalState::Internal;
        const std::size_t residue_index =
            protein->AddResidue(std::move(residue));

        for (std::size_t source_atom_index :
             source_residue.atom_indices) {
            const Atom& source_atom = source.AtomAt(source_atom_index);
            auto atom = Atom::Create(source_atom.element);
            atom->pdb_atom_name = source_atom.pdb_atom_name;
            atom->residue_index = residue_index;
            const std::size_t atom_index =
                protein->AddAtom(std::move(atom));
            Residue& added_residue =
                protein->MutableResidueAt(residue_index);
            added_residue.atom_indices.push_back(atom_index);
            if (source_atom.pdb_atom_name == "N") added_residue.N = atom_index;
            else if (source_atom.pdb_atom_name == "CA") added_residue.CA = atom_index;
            else if (source_atom.pdb_atom_name == "C") added_residue.C = atom_index;
            else if (source_atom.pdb_atom_name == "O") added_residue.O = atom_index;
            else if (source_atom.pdb_atom_name == "H") added_residue.H = atom_index;
            else if (source_atom.pdb_atom_name == "HA" ||
                     source_atom.pdb_atom_name == "HA2") {
                added_residue.HA = atom_index;
            } else if (source_atom.pdb_atom_name == "CB") {
                added_residue.CB = atom_index;
            }
            positions.push_back(source_conf.PositionAt(source_atom_index));
        }
    }

    if (protein->ResidueCount() != 3 || positions.empty()) return nullptr;
    protein->FinalizeConstruction(positions);
    protein->AddCrystalConformation(
        positions, 0.0, 0.0, 0.0, "1UBQ_TYR59_ASN60_ILE61_window");
    return protein;
}

template<typename T>
struct NpyArray {
    std::vector<std::size_t> shape;
    std::vector<T> values;
};

std::string TrimNpyToken(std::string token) {
    auto is_space = [](unsigned char c) { return std::isspace(c) != 0; };
    token.erase(token.begin(),
                std::find_if(token.begin(), token.end(),
                             [&](char c) { return !is_space(c); }));
    token.erase(std::find_if(token.rbegin(), token.rend(),
                             [&](char c) { return !is_space(c); }).base(),
                token.end());
    return token;
}

template<typename T>
NpyArray<T> ReadNpy(const fs::path& path, const char* dtype) {
    std::ifstream in(path, std::ios::binary);
    EXPECT_TRUE(in.is_open()) << path;
    NpyArray<T> result;
    if (!in.is_open()) return result;

    char magic[6] = {};
    in.read(magic, 6);
    EXPECT_EQ(std::string(magic, 6), std::string("\x93NUMPY", 6));
    char version[2] = {};
    in.read(version, 2);
    EXPECT_EQ(version[0], 1);
    EXPECT_EQ(version[1], 0);
    std::uint16_t header_length = 0;
    in.read(reinterpret_cast<char*>(&header_length), sizeof(header_length));
    std::string header(header_length, '\0');
    in.read(header.data(), header_length);
    EXPECT_NE(header.find(std::string("'descr': '") + dtype + "'"),
              std::string::npos);

    const std::size_t shape_begin = header.find('(');
    const std::size_t shape_end = header.find(')', shape_begin);
    EXPECT_NE(shape_begin, std::string::npos);
    EXPECT_NE(shape_end, std::string::npos);
    if (shape_begin == std::string::npos || shape_end == std::string::npos)
        return result;

    std::stringstream shape_stream(
        header.substr(shape_begin + 1, shape_end - shape_begin - 1));
    std::string token;
    std::size_t value_count = 1;
    while (std::getline(shape_stream, token, ',')) {
        token = TrimNpyToken(std::move(token));
        if (token.empty()) continue;
        const std::size_t extent = static_cast<std::size_t>(
            std::stoull(token));
        result.shape.push_back(extent);
        value_count *= extent;
    }
    if (result.shape.empty()) value_count = 0;
    result.values.resize(value_count);
    if (value_count > 0) {
        in.read(reinterpret_cast<char*>(result.values.data()),
                static_cast<std::streamsize>(value_count * sizeof(T)));
    }
    return result;
}

template<typename T>
std::vector<T> ReadH5Flat(const fs::path& path,
                          const std::string& dataset,
                          std::vector<std::size_t>* dimensions = nullptr) {
    HighFive::File file(path.string(), HighFive::File::ReadOnly);
    auto data_set = file.getDataSet(dataset);
    const std::vector<std::size_t> dims =
        data_set.getSpace().getDimensions();
    if (dimensions) *dimensions = dims;
    std::size_t count = 1;
    for (std::size_t extent : dims) count *= extent;
    std::vector<T> values(count);
    if (!values.empty()) data_set.read(values.data());
    return values;
}

SphericalTensor UnpackFull9(const double* values) {
    SphericalTensor tensor;
    tensor.T0 = values[0];
    for (std::size_t component = 0; component < 3; ++component)
        tensor.T1[component] = values[component + 1];
    for (std::size_t component = 0; component < 5; ++component)
        tensor.T2[component] = values[component + 4];
    return tensor;
}

double MaxFiniteDifference(const std::vector<double>& lhs,
                           const std::vector<double>& rhs) {
    EXPECT_EQ(lhs.size(), rhs.size());
    double maximum = 0.0;
    for (std::size_t i = 0; i < std::min(lhs.size(), rhs.size()); ++i) {
        if (std::isnan(lhs[i]) && std::isnan(rhs[i])) continue;
        EXPECT_TRUE(std::isfinite(lhs[i])) << "lhs[" << i << "]";
        EXPECT_TRUE(std::isfinite(rhs[i])) << "rhs[" << i << "]";
        if (std::isfinite(lhs[i]) && std::isfinite(rhs[i]))
            maximum = std::max(maximum, std::abs(lhs[i] - rhs[i]));
    }
    return maximum;
}

void RemoveTempTree(const fs::path& root) {
    if (!fs::exists(root)) return;
    std::vector<fs::path> directories;
    for (const auto& entry : fs::recursive_directory_iterator(root)) {
        if (entry.is_directory()) directories.push_back(entry.path());
        else std::remove(entry.path().string().c_str());
    }
    std::sort(directories.begin(), directories.end(),
              [](const fs::path& lhs, const fs::path& rhs) {
                  return lhs.native().size() > rhs.native().size();
              });
    for (const fs::path& directory : directories)
        ::rmdir(directory.string().c_str());
    ::rmdir(root.string().c_str());
}

fs::path SidechainCOTempDir(const char* stem) {
    const fs::path dir = fs::temp_directory_path() /
        (std::string(stem) + "_" + std::to_string(::getpid()));
    fs::create_directories(dir);
    return dir;
}

void RemoveSidechainCOOutputs(const fs::path& dir) {
    constexpr std::array<const char*, 6> files = {
        "sidechain_co_source_bonds.npy",
        "sidechain_co_frame.npy",
        "sidechain_co_frame_quality.npy",
        "sidechain_co_fixed_T2.npy",
        "sidechain_co_bo_T2.npy",
        "sidechain_co_scalar_audit.npy",
    };
    for (const char* file : files) {
        std::remove((dir / file).string().c_str());
    }
    ::rmdir(dir.string().c_str());
}

void RemoveMcConnellOutputs(const fs::path& dir) {
    for (size_t category = 0;
         category < kMcConnellSourceCategoryCount; ++category) {
        const auto cat = static_cast<McConnellSourceCategory>(category);
        for (size_t channel = 0; channel < kMcConnellChannelCount;
             ++channel) {
            const auto ch = static_cast<McConnellChannel>(channel);
            const std::string filename = std::string("mc_") +
                McConnellSourceCategoryStem(cat) + "_" +
                McConnellChannelStem(ch) + ".npy";
            std::remove((dir / filename).string().c_str());
        }
    }
    for (const char* filename : {
            "mc_peptide_co_rhombic.npy",
            "mc_nearfield_counts.npy",
            "mc_nearest_co_dir.npy",
            "mc_nearest_co_midpoint.npy",
            "mc_nearest_co_T2.npy",
            "mc_nearest_cn_T2.npy",
            "mc_nearest_co_bond_index.npy",
            "mc_nearest_cn_bond_index.npy",
            "mopac_mc_co_sum.npy",
            "mopac_mc_cn_sum.npy",
            "mopac_mc_sidechain_sum.npy",
            "mopac_mc_aromatic_sum.npy",
            "mopac_mc_co_nearest.npy",
            "mopac_mc_nearest_co_dist.npy",
            "mopac_mc_nearest_cn_dist.npy",
            "mopac_mc_nearest_co_T2.npy",
            "mopac_mc_nearest_cn_T2.npy",
            "mopac_mc_backbone_total.npy",
            "mopac_mc_sidechain_total.npy",
            "mopac_mc_aromatic_total.npy",
            "mopac_mc_shielding.npy",
            "extraction_manifest.json"}) {
        std::remove((dir / filename).string().c_str());
    }
    ::rmdir(dir.string().c_str());
}

size_t FirstBondWithCategory(const Protein& protein, BondCategory category) {
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        if (protein.BondAt(bi).category == category)
            return bi;
    }
    return SIZE_MAX;
}

size_t FirstBondBetween(const Protein& protein, size_t a, size_t b) {
    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const Bond& bond = protein.BondAt(bi);
        if ((bond.atom_index_a == a && bond.atom_index_b == b) ||
            (bond.atom_index_a == b && bond.atom_index_b == a))
            return bi;
    }
    return SIZE_MAX;
}

struct IndependentAromaticOracle {
    size_t aromatic_bond_count = 0;
    size_t retained_bo_bond_count = 0;
    std::vector<Mat3> fixed;
    std::vector<Mat3> bond_order;
};

IndependentAromaticOracle ComputeIndependentAromaticOracle(
        const ProteinConformation& conf,
        const std::vector<double>* topology_bond_orders = nullptr) {
    const Protein& protein = conf.ProteinRef();
    IndependentAromaticOracle oracle;
    oracle.fixed.assign(conf.AtomCount(), Mat3::Zero());
    oracle.bond_order.assign(conf.AtomCount(), Mat3::Zero());

    const double cutoff =
        CalculatorConfig::Get("mcconnell_bond_anisotropy_cutoff");
    const double singularity_guard =
        CalculatorConfig::Get("singularity_guard_distance");
    const double near_field_ratio =
        CalculatorConfig::Get("near_field_exclusion_ratio");
    const double bo_floor =
        CalculatorConfig::Get("mopac_bond_order_noise_floor");

    for (size_t bi = 0; bi < protein.BondCount(); ++bi) {
        const Bond& bond = protein.BondAt(bi);
        if (bond.category != BondCategory::Aromatic) continue;
        ++oracle.aromatic_bond_count;

        const Vec3 endpoint_a = conf.PositionAt(bond.atom_index_a);
        const Vec3 endpoint_b = conf.PositionAt(bond.atom_index_b);
        const Vec3 bond_vector = endpoint_b - endpoint_a;
        const double bond_length = bond_vector.norm();
        if (bond_length < 1e-15) continue;

        const Vec3 midpoint = 0.5 * (endpoint_a + endpoint_b);
        const Vec3 u = bond_vector / bond_length;
        const Mat3 source_shape =
            u * u.transpose() - Mat3::Identity() / 3.0;
        const double bo = topology_bond_orders &&
                          (*topology_bond_orders)[bi] >= bo_floor
            ? (*topology_bond_orders)[bi]
            : 0.0;
        if (bo > 0.0) ++oracle.retained_bo_bond_count;

        for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
            if (ai == bond.atom_index_a || ai == bond.atom_index_b) continue;

            const Vec3 displacement = conf.PositionAt(ai) - midpoint;
            const double distance = displacement.norm();
            if (distance > cutoff || distance < singularity_guard ||
                distance <= near_field_ratio * bond_length) {
                continue;
            }

            const Vec3 n = displacement / distance;
            const double r3 = distance * distance * distance;
            const Mat3 dipolar =
                (3.0 * n * n.transpose() - Mat3::Identity()) / r3;
            const Mat3 response = dipolar * source_shape;
            oracle.fixed[ai] += response;
            oracle.bond_order[ai] += bo * response;
        }
    }
    return oracle;
}

Mat3 SumMcConnellCategoriesExcept(const ConformationAtom& atom,
                                  size_t excluded_category,
                                  McConnellChannel channel) {
    Mat3 sum = Mat3::Zero();
    const size_t channel_index = static_cast<size_t>(channel);
    for (size_t category = 0;
         category < kMcConnellSourceCategoryCount; ++category) {
        if (category == excluded_category) continue;
        sum += atom.mcconnell_source_tensors[category][channel_index]
                   .Reconstruct();
    }
    return sum;
}

}  // namespace



// ============================================================================
// Analytical test: known geometry, hand-calculable result
// ============================================================================

TEST(McConnellAnalytical, DipolarKernelAtKnownGeometry) {
    // Atom at (3, 0, 0). Bond midpoint at origin. Bond direction (0, 0, 1).
    //
    // d = (3, 0, 0), r = 3, d_hat = (1, 0, 0)
    // b_hat = (0, 0, 1)
    // cos_theta = d_hat . b_hat = 0
    //
    // Unit axial source shape:
    //   Qhat = u u^T - I/3 = diag(-1/3, -1/3, 2/3)
    // Clean McConnell response:
    //   A = D(r) Qhat
    // PCS scalar branch:
    //   T0(A) = trace(A)/3 = n^T Qhat n / r^3 = (-1/3)/27 = -1/81
    //
    // Symmetric dipolar kernel K:
    //   K_11 = (3*1*1 - 1) / 27 =  2/27 =  0.07407
    //   K_22 = (3*0*0 - 1) / 27 = -1/27 = -0.03704
    //   K_33 = (3*0*0 - 1) / 27 = -1/27 = -0.03704
    //   Off-diagonal = 0 (d_hat only has x component)
    //   Trace = 0 (traceless)
    //
    // Response A = D Qhat:
    //   diag(-2/81, 1/81, -2/81), so trace(A)/3 = -1/81.

    // Build a minimal protein with one bond
    auto protein = std::make_unique<Protein>();

    auto atom_a = Atom::Create(Element::C);
    auto atom_b = Atom::Create(Element::O);
    auto atom_field = Atom::Create(Element::H);
    atom_a->residue_index = 0;
    atom_b->residue_index = 0;
    atom_field->residue_index = 0;

    protein->AddAtom(std::move(atom_a));  // 0: bond endpoint C at origin
    protein->AddAtom(std::move(atom_b));  // 1: bond endpoint O at (0,0,1)
    protein->AddAtom(std::move(atom_field));  // 2: field point H at (3,0,0)

    Residue res;
    res.type = AminoAcid::ALA;
    res.sequence_number = 1;
    res.chain_id = "A";
    res.atom_indices = {0, 1, 2};
    protein->AddResidue(res);

    // Positions: C at (-0.5,0,0), O at (0.5,0,0) gives midpoint at origin,
    // but we want bond direction (0,0,1). So:
    // C at (0,0,-0.5), O at (0,0,0.5), midpoint at origin, direction = (0,0,1)
    std::vector<Vec3> positions = {
        Vec3(0.0, 0.0, -0.5),   // C
        Vec3(0.0, 0.0,  0.5),   // O
        Vec3(3.0, 0.0,  0.0)    // H (field point)
    };

    protein->FinalizeConstruction(positions);

    // Manually set bond category for testing
    // (FinalizeConstruction may not classify correctly for this toy geometry)

    auto& conf = protein->AddCrystalConformation(positions, 0, 0, 0, "test");

    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));

    auto mc = McConnellResult::Compute(conf);
    ASSERT_NE(mc, nullptr);
    conf.AttachResult(std::move(mc));

    // Check the field point atom (index 2)
    const auto& ca = conf.AtomAt(2);

    // Should have at least one bond neighbour
    ASSERT_GT(ca.bond_neighbours.size(), 0u)
        << "Field atom should see at least one bond";

    // Find the bond to the C-O pair
    const BondNeighbourhood& bn = ca.bond_neighbours[0];

    // Distance should be 3.0 A (from midpoint at origin to (3,0,0))
    EXPECT_NEAR(bn.distance_to_midpoint, 3.0, 0.01)
        << "Distance to bond midpoint should be 3.0 A";

    // PCS scalar: n^T Qhat n / r^3 = -1/81.
    EXPECT_NEAR(bn.mcconnell_scalar, -1.0/81.0, 1e-6)
        << "McConnell PCS scalar should be -1/81";

    // Dipolar kernel K should be traceless
    double trace_K = bn.dipolar_tensor.trace();
    EXPECT_NEAR(trace_K, 0.0, 1e-10) << "Dipolar kernel must be traceless";

    // K diagonal: (2/27, -1/27, -1/27)
    EXPECT_NEAR(bn.dipolar_tensor(0,0),  2.0/27.0, 1e-6);
    EXPECT_NEAR(bn.dipolar_tensor(1,1), -1.0/27.0, 1e-6);
    EXPECT_NEAR(bn.dipolar_tensor(2,2), -1.0/27.0, 1e-6);

    // K off-diagonal: all zero (d_hat = (1,0,0))
    EXPECT_NEAR(bn.dipolar_tensor(0,1), 0.0, 1e-10);
    EXPECT_NEAR(bn.dipolar_tensor(0,2), 0.0, 1e-10);
    EXPECT_NEAR(bn.dipolar_tensor(1,2), 0.0, 1e-10);

    EXPECT_NEAR(ca.mc_shielding_contribution.T0, -1.0/81.0, 1e-6)
        << "T0 of D(r)Qhat should equal the PCS scalar";

    std::cout << "  Analytical test passed:\n"
              << "    PCS scalar = " << bn.mcconnell_scalar
              << " (expected -1/81 = " << -1.0/81.0 << ")\n"
              << "    T0 = " << ca.mc_shielding_contribution.T0 << "\n"
              << "    K trace = " << trace_K << "\n"
              << "    K diag = (" << bn.dipolar_tensor(0,0) << ", "
              << bn.dipolar_tensor(1,1) << ", "
              << bn.dipolar_tensor(2,2) << ")\n";
}


// ============================================================================
// Full protein test
// ============================================================================

class McConnellProteinTest : public ::testing::Test {
protected:
    void SetUp() override {
        if (!fs::exists(nmr::test::TestEnvironment::UbqProtonated())) {
            GTEST_SKIP() << "1UBQ.pdb not found";
        }
        auto r = BuildFromProtonatedPdb(nmr::test::TestEnvironment::UbqProtonated());
        if (!r.Ok()) GTEST_SKIP() << "Failed to load 1UBQ";
        protein = std::move(r.protein);

        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        conf.AttachResult(SpatialIndexResult::Compute(conf));
    }

    std::unique_ptr<Protein> protein;
};


TEST_F(McConnellProteinTest, ComputeAndAttach) {
    auto& conf = protein->Conformation();
    auto mc = McConnellResult::Compute(conf);
    ASSERT_NE(mc, nullptr);
    ASSERT_TRUE(conf.AttachResult(std::move(mc)));
    ASSERT_TRUE(conf.HasResult<McConnellResult>());
}


TEST_F(McConnellProteinTest, EveryAtomHasBondNeighbours) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int with_neighbours = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        if (!conf.AtomAt(ai).bond_neighbours.empty())
            with_neighbours++;
    }

    // Every atom in a protein should see nearby bonds
    EXPECT_GT(with_neighbours, static_cast<int>(conf.AtomCount()) - 5)
        << "Almost every atom should have bond neighbours";

    std::cout << "  Atoms with bond neighbours: " << with_neighbours
              << " / " << conf.AtomCount() << "\n";
}


TEST_F(McConnellProteinTest, DipolarKernelsAreTraceless) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    double max_trace = 0.0;
    int checked = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& bn : conf.AtomAt(ai).bond_neighbours) {
            double trace = bn.dipolar_tensor.trace();
            max_trace = std::max(max_trace, std::abs(trace));
            checked++;
        }
    }

    EXPECT_LT(max_trace, 1e-10)
        << "All dipolar kernels K must be traceless";

    std::cout << "  Checked " << checked << " kernels, max |trace| = "
              << max_trace << "\n";
}


TEST_F(McConnellProteinTest, McConnellScalarNonZero) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    // Most atoms should have non-zero CO sum
    int with_co = 0;
    double max_co = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        double co = std::abs(conf.AtomAt(ai).mcconnell_co_sum);
        if (co > 1e-6) with_co++;
        max_co = std::max(max_co, co);
    }

    EXPECT_GT(with_co, static_cast<int>(conf.AtomCount() / 2))
        << "Most atoms should have non-zero CO McConnell sum";

    std::cout << "  Atoms with |CO sum| > 1e-6: " << with_co
              << " / " << conf.AtomCount()
              << ", max |CO sum| = " << max_co << "\n";
}


TEST_F(McConnellProteinTest, ShieldingContributionHasT0AndT2) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int nonzero_t0 = 0, nonzero_t2 = 0;
    double max_t0 = 0, max_t2 = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& sc = conf.AtomAt(ai).mc_shielding_contribution;
        if (std::abs(sc.T0) > 1e-8) nonzero_t0++;
        double t2mag = sc.T2Magnitude();
        if (t2mag > 1e-8) nonzero_t2++;
        max_t0 = std::max(max_t0, std::abs(sc.T0));
        max_t2 = std::max(max_t2, t2mag);
    }

    EXPECT_GT(nonzero_t0, 0) << "Full McConnell tensor must have non-zero T0";
    EXPECT_GT(nonzero_t2, 0) << "Full McConnell tensor must have non-zero T2";

    std::cout << "  T0 nonzero: " << nonzero_t0 << ", max |T0| = " << max_t0 << "\n";
    std::cout << "  T2 nonzero: " << nonzero_t2 << ", max |T2| = " << max_t2 << "\n";
}


TEST_F(McConnellProteinTest, PCSScalarMatchesTraceBranchForSingleBond) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int checked = 0;
    double max_diff = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        for (const auto& bn : conf.AtomAt(ai).bond_neighbours) {
            Vec3 atom_pos = conf.PositionAt(ai);
            Vec3 midpoint = conf.bond_midpoints[bn.bond_index];
            Vec3 axis = conf.bond_directions[bn.bond_index];
            auto kernel = McConnellResult::ComputePairKernel(
                atom_pos, midpoint, axis);
            double r = kernel.distance;
            if (r < MIN_DISTANCE) continue;

            double r3 = r * r * r;
            double pcs = kernel.direction.dot(
                kernel.source_shape * kernel.direction) / r3;
            double t0 = SphericalTensor::Decompose(kernel.response).T0;
            max_diff = std::max(max_diff, std::abs(t0 - bn.mcconnell_scalar));
            max_diff = std::max(max_diff, std::abs(t0 - pcs));
            checked++;
        }
    }

    EXPECT_LT(max_diff, 1e-10)
        << "T0(DQhat) must equal the PCS scalar n^T Qhat n/r^3";

    std::cout << "  Checked " << checked
              << " pairs, max PCS scalar branch diff = "
              << max_diff << "\n";
}


TEST(McConnellImplementationChecks, RotationEquivariance) {
    Vec3 source_center(0.3, -0.7, 1.1);
    Vec3 source_axis(0.4, 0.8, -0.2);
    Vec3 target(3.2, 1.4, -0.6);

    const auto k = McConnellResult::ComputePairKernel(
        target, source_center, source_axis);
    const SphericalTensor st = SphericalTensor::Decompose(k.response);

    const Mat3 R = Eigen::AngleAxisd(
        0.73, Vec3(0.2, -0.4, 0.9).normalized()).toRotationMatrix();
    const auto kr = McConnellResult::ComputePairKernel(
        R * target, R * source_center, R * source_axis);
    const SphericalTensor sr = SphericalTensor::Decompose(kr.response);

    EXPECT_NEAR(sr.T0, st.T0, 1e-12);

    Vec3 t1(st.T1[0], st.T1[1], st.T1[2]);
    Vec3 t1_rot(sr.T1[0], sr.T1[1], sr.T1[2]);
    EXPECT_LT((t1_rot - R * t1).norm(), 1e-12);

    SphericalTensor st_t2 = st;
    st_t2.T0 = 0.0;
    st_t2.T1 = {0.0, 0.0, 0.0};
    SphericalTensor sr_t2 = sr;
    sr_t2.T0 = 0.0;
    sr_t2.T1 = {0.0, 0.0, 0.0};
    const Mat3 expected_t2 = R * st_t2.Reconstruct() * R.transpose();
    EXPECT_LT(MaxAbs(sr_t2.Reconstruct() - expected_t2), 1e-12);

    EXPECT_LT(MaxAbs(kr.response - R * k.response * R.transpose()), 1e-12);
}


TEST(McConnellImplementationChecks, PairPCSScalarIdentity) {
    const auto k = McConnellResult::ComputePairKernel(
        Vec3(2.0, -1.0, 3.0),
        Vec3(-0.5, 0.25, 0.75),
        Vec3(0.3, 0.7, 0.2));
    ASSERT_GT(k.distance, 0.0);
    const double r3 = k.distance * k.distance * k.distance;
    const double pcs = k.direction.dot(k.source_shape * k.direction) / r3;
    const double t0 = SphericalTensor::Decompose(k.response).T0;
    EXPECT_NEAR(t0, k.response.trace() / 3.0, 1e-15);
    EXPECT_NEAR(t0, pcs, 1e-15);
}


TEST(McConnellImplementationChecks, PeptideCORhombicSourceShapePinnedAnalytic) {
    auto protein = BuildSyntheticPeptideCOProtein();
    auto& conf = protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);

    const auto k = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));

    // C->O = +x, plane normal = +z, e_in = +y.  From
    // (-5.4, +4.0, -14) the mean-subtracted principal values are
    // (-4/15, 137/15, -133/15); dividing by axial scale 137/10 gives
    // (out, para, in) = (-8/411, 2/3, -266/411).
    Mat3 expected = Mat3::Zero();
    expected(0, 0) = 2.0 / 3.0;
    expected(1, 1) = -266.0 / 411.0;
    expected(2, 2) = -8.0 / 411.0;

    EXPECT_LT(MaxAbs(k.source_shape - expected), 1e-12);
    EXPECT_NEAR(k.source_shape.trace(), 0.0, 1e-15);
}


TEST(McConnellImplementationChecks, PeptideCORhombicSourceShapeRotationEquivariance) {
    const Mat3 R = Eigen::AngleAxisd(
        0.91, Vec3(0.3, -0.7, 0.2).normalized()).toRotationMatrix();

    auto protein = BuildSyntheticPeptideCOProtein();
    auto rotated_protein = BuildSyntheticPeptideCOProtein(false, R);
    auto& conf = protein->Conformation();
    auto& rotated_conf = rotated_protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    rotated_conf.AttachResult(GeometryResult::Compute(rotated_conf));
    const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    const size_t rotated_co_bond =
        FirstBondWithCategory(*rotated_protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);
    ASSERT_NE(rotated_co_bond, SIZE_MAX);

    const auto k = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));
    const auto kr = McConnellResult::ComputePeptideCORhombicPairKernel(
        rotated_conf, rotated_co_bond, rotated_conf.PositionAt(3));

    EXPECT_LT(MaxAbs(kr.source_shape -
                     R * k.source_shape * R.transpose()), 1e-12);
    EXPECT_LT(MaxAbs(kr.response -
                     R * k.response * R.transpose()), 1e-12);
}


TEST(McConnellImplementationChecks, PeptideCORhombicContributionNonZeroAndTraceless) {
    auto protein = BuildSyntheticPeptideCOProtein();
    auto& conf = protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    conf.AttachResult(McConnellResult::Compute(conf));

    const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);
    const auto k = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));
    EXPECT_NEAR(k.source_shape.trace(), 0.0, 1e-15);

    // The emitted rhombic array is the additive delta relative to the
    // legacy axial PeptideCO shape.  Here n = e_in = +y, r = 3 A, so
    // D = diag(-1/27, 2/27, -1/27) and
    // Q_delta = diag(0, -43/137, +43/137).
    Mat3 expected_delta = Mat3::Zero();
    expected_delta(1, 1) = -86.0 / 3699.0;
    expected_delta(2, 2) = -43.0 / 3699.0;

    const Mat3 actual =
        conf.AtomAt(3).mcconnell_peptide_co_rhombic.Reconstruct();
    EXPECT_GT(MaxAbs(actual), 1e-6);
    EXPECT_LT(MaxAbs(actual - expected_delta), 1e-12);

    // C33: the canonical fixed PeptideCO channel now contains the full
    // rhombic response.  The separately emitted audit remains exactly the
    // unweighted rhombic-minus-axial delta.
    const auto axial = McConnellResult::ComputePairKernel(
        conf.PositionAt(3), conf.bond_midpoints[co_bond],
        conf.bond_directions[co_bond]);
    const auto rhombic = McConnellResult::ComputePeptideCORhombicPairKernel(
        conf, co_bond, conf.PositionAt(3));
    const size_t cat = static_cast<size_t>(McConnellSourceCategory::PeptideCO);
    const Mat3 fixed = conf.AtomAt(3)
        .mcconnell_source_tensors[cat]
                                   [static_cast<size_t>(McConnellChannel::Fixed)]
        .Reconstruct();
    EXPECT_LT(MaxAbs(fixed - rhombic.response), 1e-12);
    EXPECT_LT(MaxAbs(fixed - (axial.response + actual)), 1e-12);

    // The same production selector owns the BO channel.  Pin a non-unit BO
    // so an axial-only regression cannot hide behind fixed-channel coverage.
    const auto channels = mcconnell_result_detail::SelectChannelResponses(
        McConnellSourceCategory::PeptideCO,
        axial.response, rhombic.response, 1.75);
    EXPECT_LT(MaxAbs(channels.fixed - rhombic.response), 1e-12);
    EXPECT_LT(MaxAbs(channels.bond_order - 1.75 * rhombic.response), 1e-12);
    EXPECT_LT(MaxAbs(channels.rhombic_audit - actual), 1e-12);
}


TEST(McConnellImplementationChecks,
     PeptideCORhombicBondOrderChannelCallsRealMopacAndProductionCompute) {
    auto protein = BuildSyntheticPeptideCOProtein();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    const size_t co_bond =
        FirstBondWithCategory(*protein, BondCategory::PeptideCO);
    ASSERT_NE(co_bond, SIZE_MAX);

    // No synthetic bond-order seeding: run the real upstream producer, then
    // feed its topology-parallel BO into the full McConnell Compute path.
    auto mopac = MopacResult::Compute(conf, 0, 1);
    ASSERT_NE(mopac, nullptr) << "real MOPAC calculation failed";
    const double bond_order = mopac->TopologyBondOrder(co_bond);
    EXPECT_GT(bond_order,
              CalculatorConfig::Get("mopac_bond_order_noise_floor"));
    ASSERT_TRUE(conf.AttachResult(std::move(mopac)));
    ASSERT_TRUE(conf.AttachResult(McConnellResult::Compute(conf)));

    constexpr size_t target = 3;
    const Mat3 rhombic =
        McConnellResult::ComputePeptideCORhombicPairKernel(
            conf, co_bond, conf.PositionAt(target)).response;
    const size_t category =
        static_cast<size_t>(McConnellSourceCategory::PeptideCO);
    const size_t channel =
        static_cast<size_t>(McConnellChannel::BondOrder);
    const Mat3 production = conf.AtomAt(target)
        .mcconnell_source_tensors[category][channel].Reconstruct();
    EXPECT_GT(MaxAbs(production), 1e-8);
    EXPECT_LT(MaxAbs(production - bond_order * rhombic), 1e-11);
    EXPECT_EQ(conf.AtomAt(target).nearest_CO_bond_index, co_bond);
    EXPECT_EQ(conf.AtomAt(target).nearest_CN_bond_index, SIZE_MAX);

    Mat3 bo_category_sum = Mat3::Zero();
    for (size_t source_category = 0;
         source_category < kMcConnellSourceCategoryCount;
         ++source_category) {
        bo_category_sum += conf.AtomAt(target)
            .mcconnell_source_tensors[source_category][channel]
            .Reconstruct();
    }
    EXPECT_GT(MaxAbs(bo_category_sum), 1e-8);
    EXPECT_LT(MaxAbs(
        bo_category_sum - conf.AtomAt(target)
            .mopac_mc_shielding_contribution.Reconstruct()),
        1e-11)
        << "the compatibility/trajectory source must be the full BO sum";
    ASSERT_TRUE(conf.AttachResult(MopacMcConnellResult::Compute(conf)));

    const fs::path out_dir = fs::temp_directory_path() /
        ("mcconnell_real_mopac_bo_" + std::to_string(::getpid()));
    fs::create_directories(out_dir);
    ASSERT_EQ(conf.Result<McConnellResult>().WriteFeatures(
                  conf, out_dir.string()),
              28);
    ASSERT_EQ(conf.Result<MopacMcConnellResult>().WriteFeatures(
                  conf, out_dir.string()),
              13);
    const auto emitted = ReadNpy<double>(
        out_dir / "mc_peptide_co_bo.npy", "<f8");
    ASSERT_EQ(emitted.shape,
              (std::vector<size_t>{conf.AtomCount(), 9}));
    std::array<double, 9> expected{};
    SphericalTensor::Decompose(bond_order * rhombic)
        .PackFull9(expected.data());
    for (size_t component = 0; component < expected.size(); ++component) {
        EXPECT_NEAR(emitted.values[target * 9 + component],
                    expected[component], 1e-11);
    }

    const auto emitted_co_bond_index = ReadNpy<int32_t>(
        out_dir / "mc_nearest_co_bond_index.npy", "<i4");
    const auto emitted_cn_bond_index = ReadNpy<int32_t>(
        out_dir / "mc_nearest_cn_bond_index.npy", "<i4");
    ASSERT_EQ(emitted_co_bond_index.shape,
              (std::vector<size_t>{conf.AtomCount()}));
    ASSERT_EQ(emitted_cn_bond_index.shape,
              (std::vector<size_t>{conf.AtomCount()}));
    for (size_t atom = 0; atom < conf.AtomCount(); ++atom) {
        const auto& ca = conf.AtomAt(atom);
        const int32_t expected_co = ca.nearest_CO_bond_index == SIZE_MAX
            ? -1 : static_cast<int32_t>(ca.nearest_CO_bond_index);
        const int32_t expected_cn = ca.nearest_CN_bond_index == SIZE_MAX
            ? -1 : static_cast<int32_t>(ca.nearest_CN_bond_index);
        EXPECT_EQ(emitted_co_bond_index.values[atom], expected_co);
        EXPECT_EQ(emitted_cn_bond_index.values[atom], expected_cn);
    }
    EXPECT_EQ(emitted_co_bond_index.values[target],
              static_cast<int32_t>(co_bond));
    EXPECT_EQ(emitted_cn_bond_index.values[target], -1);
    EXPECT_NE(std::find(emitted_co_bond_index.values.begin(),
                        emitted_co_bond_index.values.end(), -1),
              emitted_co_bond_index.values.end());
    EXPECT_NE(std::find(emitted_cn_bond_index.values.begin(),
                        emitted_cn_bond_index.values.end(), -1),
              emitted_cn_bond_index.values.end());

    auto expect_scalar_readback = [&](const char* stem, auto field) {
        const auto array = ReadNpy<double>(
            out_dir / (std::string(stem) + ".npy"), "<f8");
        ASSERT_EQ(array.shape,
                  (std::vector<size_t>{conf.AtomCount()}));
        for (size_t atom = 0; atom < conf.AtomCount(); ++atom) {
            const double expected_value = field(conf.AtomAt(atom));
            if (std::isnan(expected_value)) {
                EXPECT_TRUE(std::isnan(array.values[atom])) << stem;
            } else {
                EXPECT_DOUBLE_EQ(array.values[atom], expected_value) << stem;
            }
        }
    };
    expect_scalar_readback("mopac_mc_co_sum",
        [](const ConformationAtom& atom) { return atom.mopac_mc_co_sum; });
    expect_scalar_readback("mopac_mc_cn_sum",
        [](const ConformationAtom& atom) { return atom.mopac_mc_cn_sum; });
    expect_scalar_readback("mopac_mc_sidechain_sum",
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_sidechain_sum;
        });
    expect_scalar_readback("mopac_mc_aromatic_sum",
        [](const ConformationAtom& atom) { return atom.mopac_mc_aromatic_sum; });
    expect_scalar_readback("mopac_mc_co_nearest",
        [](const ConformationAtom& atom) { return atom.mopac_mc_co_nearest; });
    expect_scalar_readback("mopac_mc_nearest_co_dist",
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_nearest_CO_dist;
        });
    expect_scalar_readback("mopac_mc_nearest_cn_dist",
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_nearest_CN_dist;
        });

    auto expect_tensor_readback = [&](const char* stem, auto field,
                                      auto available) {
        const auto array = ReadNpy<double>(
            out_dir / (std::string(stem) + ".npy"), "<f8");
        ASSERT_EQ(array.shape,
                  (std::vector<size_t>{conf.AtomCount(), 9u}));
        for (size_t atom = 0; atom < conf.AtomCount(); ++atom) {
            if (!available(conf.AtomAt(atom))) {
                for (size_t component = 0; component < 9; ++component) {
                    EXPECT_TRUE(std::isnan(
                        array.values[atom * 9 + component])) << stem;
                }
                continue;
            }
            std::array<double, 9> packed{};
            field(conf.AtomAt(atom)).PackFull9(packed.data());
            for (size_t component = 0; component < packed.size(); ++component) {
                EXPECT_DOUBLE_EQ(array.values[atom * 9 + component],
                                 packed[component]) << stem;
            }
        }
    };
    const auto always = [](const ConformationAtom&) { return true; };
    expect_tensor_readback("mopac_mc_nearest_co_T2",
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_T2_CO_nearest;
        },
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_nearest_CO_dist < NO_DATA_SENTINEL;
        });
    expect_tensor_readback("mopac_mc_nearest_cn_T2",
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_T2_CN_nearest;
        },
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_nearest_CN_dist < NO_DATA_SENTINEL;
        });
    expect_tensor_readback("mopac_mc_backbone_total",
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_T2_backbone_total;
        }, always);
    expect_tensor_readback("mopac_mc_sidechain_total",
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_T2_sidechain_total;
        }, always);
    expect_tensor_readback("mopac_mc_aromatic_total",
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_T2_aromatic_total;
        }, always);
    expect_tensor_readback("mopac_mc_shielding",
        [](const ConformationAtom& atom) {
            return atom.mopac_mc_shielding_contribution;
        }, always);

    // Geometry owns the nearest source; this fixture has one C=O. Pin its
    // distance and the scalar/Full9 MOPAC weighting to the production kernel
    // and real bond order, not to the fields the writer reads.
    const auto nearest_kernel = McConnellResult::ComputePairKernel(
        conf.PositionAt(target), conf.bond_midpoints[co_bond],
        conf.bond_directions[co_bond]);
    EXPECT_NEAR(conf.AtomAt(target).mopac_mc_nearest_CO_dist,
                nearest_kernel.distance, 1e-14);
    EXPECT_NEAR(conf.AtomAt(target).mopac_mc_co_nearest,
                bond_order * nearest_kernel.pcs_scalar, 1e-14);
    const auto emitted_nearest_co = ReadNpy<double>(
        out_dir / "mopac_mc_nearest_co_T2.npy", "<f8");
    std::array<double, 9> expected_nearest{};
    SphericalTensor::Decompose(bond_order * nearest_kernel.response)
        .PackFull9(expected_nearest.data());
    for (size_t component = 0; component < expected_nearest.size();
         ++component) {
        EXPECT_NEAR(emitted_nearest_co.values[target * 9 + component],
                    expected_nearest[component], 1e-11);
    }
    RemoveMcConnellOutputs(out_dir);
}


TEST(McConnellImplementationChecks, XHBondsUseDedicatedProductionCategories) {
    auto protein = BuildSyntheticXHCategoryProtein();
    ASSERT_EQ(protein->BondCount(), 3u);

    std::array<McConnellSourceCategory, 3> categories{};
    for (size_t bi = 0; bi < protein->BondCount(); ++bi) {
        ASSERT_TRUE(mcconnell_result_detail::IsXHBond(
            *protein, protein->BondAt(bi)));
        categories[bi] = mcconnell_result_detail::SourceCategory(
            *protein, protein->BondAt(bi));
    }

    // Bonds are emitted in atom-pair order by CovalentTopology::Resolve.
    EXPECT_EQ(categories[0], McConnellSourceCategory::BackboneXH);
    EXPECT_EQ(categories[1], McConnellSourceCategory::SidechainXH);
    EXPECT_EQ(categories[2], McConnellSourceCategory::SH);
    EXPECT_STREQ(McConnellSourceCategoryStem(categories[0]), "backbone_xh");
    EXPECT_STREQ(McConnellSourceCategoryStem(categories[1]), "sidechain_xh");
    EXPECT_STREQ(McConnellSourceCategoryStem(categories[2]), "s_h");

    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    // Exercise the BO redistribution through the real upstream producer.
    // No bond-order state is seeded in this test.
    // The three disconnected X-H fragments are a deliberately synthetic
    // category probe. Their closed-shell electronic state is -2; neutral
    // spin-0 MOZYME correctly rejects the under-valenced Lewis guess.
    auto mopac = MopacResult::Compute(conf, -2, 1);
    ASSERT_NE(mopac, nullptr) << "real MOPAC calculation failed";
    std::array<double, 3> bond_orders{};
    for (size_t bi = 0; bi < bond_orders.size(); ++bi) {
        bond_orders[bi] = mopac->TopologyBondOrder(bi);
        EXPECT_GT(bond_orders[bi],
                  CalculatorConfig::Get("mopac_bond_order_noise_floor"));
    }
    ASSERT_TRUE(conf.AttachResult(std::move(mopac)));
    ASSERT_TRUE(conf.AttachResult(McConnellResult::Compute(conf)));

    const size_t fixed = static_cast<size_t>(McConnellChannel::Fixed);
    const size_t bo = static_cast<size_t>(McConnellChannel::BondOrder);
    const size_t backbone_xh =
        static_cast<size_t>(McConnellSourceCategory::BackboneXH);
    const size_t sidechain_xh =
        static_cast<size_t>(McConnellSourceCategory::SidechainXH);
    const size_t sh = static_cast<size_t>(McConnellSourceCategory::SH);
    const size_t backbone_other =
        static_cast<size_t>(McConnellSourceCategory::BackboneOther);
    const size_t sidechain_other =
        static_cast<size_t>(McConnellSourceCategory::SidechainOther);

    std::array<bool, 3> dedicated_nonzero{{false, false, false}};
    std::array<bool, 3> dedicated_bo_nonzero{{false, false, false}};
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& atom = conf.AtomAt(ai);
        dedicated_nonzero[0] = dedicated_nonzero[0] ||
            !SphericalAllZero(
                atom.mcconnell_source_tensors[backbone_xh][fixed]);
        dedicated_nonzero[1] = dedicated_nonzero[1] ||
            !SphericalAllZero(
                atom.mcconnell_source_tensors[sidechain_xh][fixed]);
        dedicated_nonzero[2] = dedicated_nonzero[2] ||
            !SphericalAllZero(atom.mcconnell_source_tensors[sh][fixed]);
        dedicated_bo_nonzero[0] = dedicated_bo_nonzero[0] ||
            !SphericalAllZero(
                atom.mcconnell_source_tensors[backbone_xh][bo]);
        dedicated_bo_nonzero[1] = dedicated_bo_nonzero[1] ||
            !SphericalAllZero(
                atom.mcconnell_source_tensors[sidechain_xh][bo]);
        dedicated_bo_nonzero[2] = dedicated_bo_nonzero[2] ||
            !SphericalAllZero(atom.mcconnell_source_tensors[sh][bo]);

        EXPECT_TRUE(SphericalAllZero(
            atom.mcconnell_source_tensors[backbone_other][fixed]));
        EXPECT_TRUE(SphericalAllZero(
            atom.mcconnell_source_tensors[sidechain_other][fixed]));
        EXPECT_TRUE(SphericalAllZero(
            atom.mcconnell_source_tensors[backbone_other][bo]));
        EXPECT_TRUE(SphericalAllZero(
            atom.mcconnell_source_tensors[sidechain_other][bo]));

        Mat3 category_sum = Mat3::Zero();
        for (size_t category = 0;
             category < kMcConnellSourceCategoryCount; ++category) {
            category_sum += atom.mcconnell_source_tensors[category][fixed]
                                .Reconstruct();
        }
        EXPECT_LT(MaxAbs(category_sum -
                         atom.mc_shielding_contribution.Reconstruct()),
                  1e-14);

        Mat3 bo_category_sum = Mat3::Zero();
        for (size_t category = 0;
             category < kMcConnellSourceCategoryCount; ++category) {
            bo_category_sum += atom.mcconnell_source_tensors[category][bo]
                                   .Reconstruct();
        }
        EXPECT_LT(MaxAbs(
            bo_category_sum -
            atom.mopac_mc_shielding_contribution.Reconstruct()),
            1e-12);
    }
    EXPECT_EQ(dedicated_nonzero,
              (std::array<bool, 3>{{true, true, true}}));
    EXPECT_EQ(dedicated_bo_nonzero,
              (std::array<bool, 3>{{true, true, true}}));

    // Independent analytic oracle for atom 2.  It is an endpoint of the
    // sidechain C-H source (self-excluded), the S-H midpoint lies beyond the
    // 10 A cutoff, and the backbone N-H source lies on +x with midpoint 0.5.
    // Therefore its entire fixed response is the single N-H D(r)Qhat term:
    // diag(4/3,1/3,1/3) / 9.5^3.
    const double r3 = 9.5 * 9.5 * 9.5;
    Mat3 expected = Mat3::Zero();
    expected(0, 0) = (4.0 / 3.0) / r3;
    expected(1, 1) = (1.0 / 3.0) / r3;
    expected(2, 2) = (1.0 / 3.0) / r3;
    const Mat3 production = conf.AtomAt(2)
        .mcconnell_source_tensors[backbone_xh][fixed].Reconstruct();
    EXPECT_LT(MaxAbs(production - expected), 1e-14);
    EXPECT_LT(MaxAbs(conf.AtomAt(2).T2_backbone_total.Reconstruct() -
                     expected),
              1e-14);

    // Independent BO-channel redistribution oracle.  Each selected target
    // is a self-excluded endpoint of its own X-H bond and sees exactly one
    // other X-H source inside the 10 A cutoff.  All three sources/targets
    // are collinear, so D(r)Qhat = diag(4/3,1/3,1/3)/r^3, scaled by the
    // real MOPAC Wiberg order of that source bond.
    auto expected_collinear_bo = [](double bond_order, double distance) {
        const double r3 = distance * distance * distance;
        Mat3 result = Mat3::Zero();
        result(0, 0) = bond_order * (4.0 / 3.0) / r3;
        result(1, 1) = bond_order * (1.0 / 3.0) / r3;
        result(2, 2) = bond_order * (1.0 / 3.0) / r3;
        return result;
    };
    const Mat3 expected_backbone_bo =
        expected_collinear_bo(bond_orders[0], 9.5);
    const Mat3 expected_sidechain_bo =
        expected_collinear_bo(bond_orders[1], 9.5);
    const Mat3 expected_sh_bo =
        expected_collinear_bo(bond_orders[2], 9.65);
    EXPECT_LT(MaxAbs(
        conf.AtomAt(2).mcconnell_source_tensors[backbone_xh][bo]
            .Reconstruct() - expected_backbone_bo),
        1e-12);
    EXPECT_LT(MaxAbs(
        conf.AtomAt(1).mcconnell_source_tensors[sidechain_xh][bo]
            .Reconstruct() - expected_sidechain_bo),
        1e-12);
    EXPECT_LT(MaxAbs(
        conf.AtomAt(3).mcconnell_source_tensors[sh][bo]
            .Reconstruct() - expected_sh_bo),
        1e-12);

    const fs::path out_dir = fs::temp_directory_path() /
        ("mcconnell_xh_production_" + std::to_string(::getpid()));
    fs::create_directories(out_dir);
    ASSERT_EQ(conf.Result<McConnellResult>().WriteFeatures(
                  conf, out_dir.string()),
              28);
    const auto emitted = ReadNpy<double>(
        out_dir / "mc_backbone_xh_fixed.npy", "<f8");
    const auto emitted_backbone_bo = ReadNpy<double>(
        out_dir / "mc_backbone_xh_bo.npy", "<f8");
    const auto emitted_sidechain_bo = ReadNpy<double>(
        out_dir / "mc_sidechain_xh_bo.npy", "<f8");
    const auto emitted_sh_bo = ReadNpy<double>(
        out_dir / "mc_s_h_bo.npy", "<f8");
    const auto emitted_old_sidechain_bo = ReadNpy<double>(
        out_dir / "mc_sidechain_other_bo.npy", "<f8");
    ASSERT_EQ(emitted.shape,
              (std::vector<size_t>{conf.AtomCount(), 9}));
    ASSERT_EQ(emitted_backbone_bo.shape, emitted.shape);
    ASSERT_EQ(emitted_sidechain_bo.shape, emitted.shape);
    ASSERT_EQ(emitted_sh_bo.shape, emitted.shape);
    ASSERT_EQ(emitted_old_sidechain_bo.shape, emitted.shape);
    std::array<double, 9> expected_packed{};
    SphericalTensor::Decompose(expected).PackFull9(expected_packed.data());
    for (size_t component = 0; component < expected_packed.size();
         ++component) {
        EXPECT_NEAR(emitted.values[2 * 9 + component],
                    expected_packed[component], 1e-14);
    }
    auto expect_emitted_tensor = [&](const NpyArray<double>& array,
                                     size_t atom_index,
                                     const Mat3& expected_tensor) {
        std::array<double, 9> packed{};
        SphericalTensor::Decompose(expected_tensor).PackFull9(packed.data());
        for (size_t component = 0; component < packed.size(); ++component) {
            EXPECT_NEAR(array.values[atom_index * 9 + component],
                        packed[component], 1e-12);
        }
    };
    expect_emitted_tensor(emitted_backbone_bo, 2, expected_backbone_bo);
    expect_emitted_tensor(emitted_sidechain_bo, 1, expected_sidechain_bo);
    expect_emitted_tensor(emitted_sh_bo, 3, expected_sh_bo);
    for (double value : emitted_old_sidechain_bo.values) {
        EXPECT_DOUBLE_EQ(value, 0.0)
            << "X-H BO response must not remain in SidechainOther";
    }
    RemoveMcConnellOutputs(out_dir);
}


TEST(McConnellImplementationChecks,
     NearestSourceNpyRowsAreNanWhenNoAcceptedCoOrCnExists) {
    auto protein = BuildSyntheticXHCategoryProtein();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(McConnellResult::Compute(conf)));

    for (std::size_t i = 0; i < conf.AtomCount(); ++i) {
        ASSERT_EQ(conf.AtomAt(i).nearest_CO_dist, NO_DATA_SENTINEL);
        ASSERT_EQ(conf.AtomAt(i).nearest_CN_dist, NO_DATA_SENTINEL);
    }

    const fs::path out_dir = fs::temp_directory_path() /
        ("mcconnell_missing_nearest_" + std::to_string(::getpid()));
    fs::create_directories(out_dir);
    ASSERT_EQ(conf.Result<McConnellResult>().WriteFeatures(
                  conf, out_dir.string()),
              28);

    const auto co_dir = ReadNpy<double>(
        out_dir / "mc_nearest_co_dir.npy", "<f8");
    const auto co_midpoint = ReadNpy<double>(
        out_dir / "mc_nearest_co_midpoint.npy", "<f8");
    const auto co_tensor = ReadNpy<double>(
        out_dir / "mc_nearest_co_T2.npy", "<f8");
    const auto cn_tensor = ReadNpy<double>(
        out_dir / "mc_nearest_cn_T2.npy", "<f8");
    const auto co_bond_index = ReadNpy<int32_t>(
        out_dir / "mc_nearest_co_bond_index.npy", "<i4");
    const auto cn_bond_index = ReadNpy<int32_t>(
        out_dir / "mc_nearest_cn_bond_index.npy", "<i4");
    ASSERT_EQ(co_dir.shape,
              (std::vector<std::size_t>{conf.AtomCount(), 3}));
    ASSERT_EQ(co_midpoint.shape,
              (std::vector<std::size_t>{conf.AtomCount(), 3}));
    ASSERT_EQ(co_tensor.shape,
              (std::vector<std::size_t>{conf.AtomCount(), 9}));
    ASSERT_EQ(cn_tensor.shape,
              (std::vector<std::size_t>{conf.AtomCount(), 9}));
    ASSERT_EQ(co_bond_index.shape,
              (std::vector<std::size_t>{conf.AtomCount()}));
    ASSERT_EQ(cn_bond_index.shape,
              (std::vector<std::size_t>{conf.AtomCount()}));
    for (const auto* array : {&co_dir, &co_midpoint, &co_tensor, &cn_tensor}) {
        for (double value : array->values) {
            EXPECT_TRUE(std::isnan(value));
        }
    }
    for (const auto* array : {&co_bond_index, &cn_bond_index}) {
        for (int32_t value : array->values) {
            EXPECT_EQ(value, -1);
        }
    }

    RemoveMcConnellOutputs(out_dir);
}


TEST(SidechainCarbonylAnisotropyProduction,
     TypedSourcesFramesAndNoMopacReadBack) {
    using namespace sidechain_carbonyl_anisotropy_detail;

    auto protein = BuildSyntheticSidechainCOProtein();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    std::vector<SourceBond> expected_sources;
    std::size_t peptide_co_count = 0;
    for (std::size_t bond_index = 0;
         bond_index < protein->BondCount(); ++bond_index) {
        const Bond& bond = protein->BondAt(bond_index);
        if (bond.category == BondCategory::PeptideCO) ++peptide_co_count;
        if (bond.category == BondCategory::SidechainCO) {
            expected_sources.push_back(
                ClassifySourceBond(*protein, bond_index));
        }
    }
    ASSERT_EQ(peptide_co_count, 1u);
    ASSERT_EQ(expected_sources.size(), 3u);

    std::size_t amide_sources = 0;
    std::size_t carboxylate_sources = 0;
    for (const SourceBond& source : expected_sources) {
        EXPECT_TRUE(source.source_valid);
        EXPECT_EQ(source.bond_category, BondCategory::SidechainCO);
        EXPECT_NE(source.plane_reference_atom, SIZE_MAX);
        if (source.oxygen_semantic_class ==
            OxygenSemanticClass::SidechainAmide) {
            ++amide_sources;
            EXPECT_EQ(source.planar_group_kind,
                      PlanarGroupKind::SidechainAmide);
        } else if (source.oxygen_semantic_class ==
                   OxygenSemanticClass::SidechainCarboxylate) {
            ++carboxylate_sources;
            EXPECT_EQ(source.planar_group_kind,
                      PlanarGroupKind::Carboxylate);
        }

        const SourceFrame frame = BuildSourceFrame(conf, source);
        EXPECT_TRUE(frame.frame_valid);
        EXPECT_TRUE(frame.origin.allFinite());
        EXPECT_NEAR(frame.x_axis.norm(), 1.0, 1e-14);
        EXPECT_NEAR(frame.y_axis.norm(), 1.0, 1e-14);
        EXPECT_NEAR(frame.z_axis.norm(), 1.0, 1e-14);
        EXPECT_NEAR(frame.x_axis.dot(frame.y_axis), 0.0, 1e-14);
        EXPECT_NEAR(frame.x_axis.dot(frame.z_axis), 0.0, 1e-14);
        EXPECT_NEAR(frame.y_axis.dot(frame.z_axis), 0.0, 1e-14);
        EXPECT_NEAR(frame.z_axis.cross(frame.x_axis).dot(frame.y_axis),
                    1.0, 1e-14);
        const Vec3 expected_x =
            (conf.PositionAt(source.oxygen_atom) -
             conf.PositionAt(source.carbon_atom)).normalized();
        EXPECT_NEAR(frame.x_axis.dot(expected_x), 1.0, 1e-14);
        EXPECT_LT(frame.orthogonality_error, 1e-14);
        EXPECT_GT(frame.normal_norm, 0.0);
    }
    EXPECT_EQ(amide_sources, 1u);
    EXPECT_EQ(carboxylate_sources, 2u);

    ASSERT_TRUE(conf.AttachResult(McConnellResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(
        SidechainCarbonylAnisotropyResult::Compute(conf)));

    const fs::path output_dir =
        SidechainCOTempDir("sidechain_co_no_mopac");
    ASSERT_EQ(conf.Result<SidechainCarbonylAnisotropyResult>()
                  .WriteFeatures(conf, output_dir.string()),
              6);

    const auto source_rows = ReadNpy<std::int32_t>(
        output_dir / "sidechain_co_source_bonds.npy", "<i4");
    const auto frame_rows = ReadNpy<double>(
        output_dir / "sidechain_co_frame.npy", "<f8");
    const auto quality_rows = ReadNpy<double>(
        output_dir / "sidechain_co_frame_quality.npy", "<f8");
    const auto fixed_rows = ReadNpy<double>(
        output_dir / "sidechain_co_fixed_T2.npy", "<f8");
    const auto bo_rows = ReadNpy<double>(
        output_dir / "sidechain_co_bo_T2.npy", "<f8");
    const auto audit_rows = ReadNpy<double>(
        output_dir / "sidechain_co_scalar_audit.npy", "<f8");

    const std::size_t Q = expected_sources.size();
    const std::size_t N = conf.AtomCount();
    ASSERT_EQ(source_rows.shape, (std::vector<std::size_t>{Q, 8}));
    ASSERT_EQ(frame_rows.shape, (std::vector<std::size_t>{Q, 12}));
    ASSERT_EQ(quality_rows.shape, (std::vector<std::size_t>{Q, 4}));
    ASSERT_EQ(fixed_rows.shape, (std::vector<std::size_t>{N, 9}));
    ASSERT_EQ(bo_rows.shape, (std::vector<std::size_t>{N, 9}));
    ASSERT_EQ(audit_rows.shape, (std::vector<std::size_t>{N, 4}));

    for (std::size_t row = 0; row < Q; ++row) {
        const SourceBond& source = expected_sources[row];
        EXPECT_EQ(source_rows.values[row * 8 + 0],
                  static_cast<std::int32_t>(source.bond_index));
        EXPECT_EQ(source_rows.values[row * 8 + 1],
                  static_cast<std::int32_t>(source.carbon_atom));
        EXPECT_EQ(source_rows.values[row * 8 + 2],
                  static_cast<std::int32_t>(source.oxygen_atom));
        EXPECT_EQ(source_rows.values[row * 8 + 3],
                  static_cast<std::int32_t>(source.residue_index));
        EXPECT_EQ(source_rows.values[row * 8 + 4],
                  static_cast<std::int32_t>(BondCategory::SidechainCO));
        EXPECT_NE(source_rows.values[row * 8 + 4],
                  static_cast<std::int32_t>(BondCategory::PeptideCO));
        EXPECT_EQ(source_rows.values[row * 8 + 5],
                  static_cast<std::int32_t>(source.planar_group_kind));
        EXPECT_EQ(source_rows.values[row * 8 + 6],
                  static_cast<std::int32_t>(
                      source.oxygen_semantic_class));
        EXPECT_EQ(source_rows.values[row * 8 + 7], 1);

        const SourceFrame production_frame =
            BuildSourceFrame(conf, source);
        const std::array<double, 12> packed_frame = {
            production_frame.origin.x(), production_frame.origin.y(),
            production_frame.origin.z(), production_frame.x_axis.x(),
            production_frame.x_axis.y(), production_frame.x_axis.z(),
            production_frame.y_axis.x(), production_frame.y_axis.y(),
            production_frame.y_axis.z(), production_frame.z_axis.x(),
            production_frame.z_axis.y(), production_frame.z_axis.z(),
        };
        for (std::size_t component = 0; component < 12; ++component) {
            EXPECT_DOUBLE_EQ(frame_rows.values[row * 12 + component],
                             packed_frame[component]);
        }
        EXPECT_DOUBLE_EQ(quality_rows.values[row * 4 + 0],
                         production_frame.bond_length);
        EXPECT_DOUBLE_EQ(quality_rows.values[row * 4 + 1],
                         production_frame.orthogonality_error);
        EXPECT_DOUBLE_EQ(quality_rows.values[row * 4 + 2],
                         production_frame.normal_norm);
        EXPECT_DOUBLE_EQ(quality_rows.values[row * 4 + 3], 1.0);
    }

    bool any_nonzero_fixed = false;
    const std::size_t sidechain_category = static_cast<std::size_t>(
        McConnellSourceCategory::SidechainCO);
    const std::size_t fixed_channel = static_cast<std::size_t>(
        McConnellChannel::Fixed);
    for (std::size_t atom_index = 0; atom_index < N; ++atom_index) {
        std::array<double, 9> expected_fixed{};
        const SphericalTensor& tensor = conf.AtomAt(atom_index)
            .mcconnell_source_tensors[sidechain_category][fixed_channel];
        tensor.PackFull9(expected_fixed.data());
        for (std::size_t component = 0; component < 9; ++component) {
            EXPECT_DOUBLE_EQ(fixed_rows.values[atom_index * 9 + component],
                             expected_fixed[component]);
            EXPECT_TRUE(std::isnan(
                bo_rows.values[atom_index * 9 + component]));
            any_nonzero_fixed = any_nonzero_fixed ||
                std::abs(expected_fixed[component]) > 0.0;
        }

        std::size_t expected_count = 0;
        double expected_nearest = std::numeric_limits<double>::infinity();
        for (const BondNeighbourhood& neighbour :
             conf.AtomAt(atom_index).bond_neighbours) {
            if (neighbour.bond_category != BondCategory::SidechainCO)
                continue;
            ++expected_count;
            expected_nearest = std::min(
                expected_nearest, neighbour.distance_to_midpoint);
        }
        EXPECT_DOUBLE_EQ(audit_rows.values[atom_index * 4 + 0],
                         tensor.T2Magnitude());
        EXPECT_TRUE(std::isnan(audit_rows.values[atom_index * 4 + 1]));
        EXPECT_DOUBLE_EQ(audit_rows.values[atom_index * 4 + 2],
                         static_cast<double>(expected_count));
        if (expected_count > 0) {
            EXPECT_DOUBLE_EQ(audit_rows.values[atom_index * 4 + 3],
                             expected_nearest);
        } else {
            EXPECT_TRUE(std::isnan(
                audit_rows.values[atom_index * 4 + 3]));
        }
    }
    EXPECT_TRUE(any_nonzero_fixed);

    RemoveSidechainCOOutputs(output_dir);
}


TEST(SidechainCarbonylAnisotropyProduction,
     RealMopacMakesBondOrderRowsFiniteAndNonzero) {
    auto protein = BuildSyntheticSidechainCOProtein();
    auto& conf = protein->Conformation();
    ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

    auto mopac = MopacResult::Compute(conf, 0, 1);
    ASSERT_NE(mopac, nullptr) << "real MOPAC calculation failed";
    bool any_sidechain_co_bo = false;
    for (size_t bond_index = 0; bond_index < protein->BondCount();
         ++bond_index) {
        if (protein->BondAt(bond_index).category !=
            BondCategory::SidechainCO) {
            continue;
        }
        any_sidechain_co_bo = any_sidechain_co_bo ||
            mopac->TopologyBondOrder(bond_index) >=
                CalculatorConfig::Get("mopac_bond_order_noise_floor");
    }
    ASSERT_TRUE(any_sidechain_co_bo)
        << "real MOPAC must produce at least one retained SidechainCO BO";
    ASSERT_TRUE(conf.AttachResult(std::move(mopac)));
    ASSERT_TRUE(conf.AttachResult(McConnellResult::Compute(conf)));
    ASSERT_TRUE(conf.AttachResult(
        SidechainCarbonylAnisotropyResult::Compute(conf)));

    const fs::path output_dir =
        SidechainCOTempDir("sidechain_co_with_mopac");
    ASSERT_EQ(conf.Result<SidechainCarbonylAnisotropyResult>()
                  .WriteFeatures(conf, output_dir.string()),
              6);
    const auto bo_rows = ReadNpy<double>(
        output_dir / "sidechain_co_bo_T2.npy", "<f8");
    const auto audit_rows = ReadNpy<double>(
        output_dir / "sidechain_co_scalar_audit.npy", "<f8");
    ASSERT_EQ(bo_rows.shape,
              (std::vector<std::size_t>{conf.AtomCount(), 9}));
    ASSERT_EQ(audit_rows.shape,
              (std::vector<std::size_t>{conf.AtomCount(), 4}));
    const size_t category =
        static_cast<size_t>(McConnellSourceCategory::SidechainCO);
    const size_t channel =
        static_cast<size_t>(McConnellChannel::BondOrder);
    bool any_nonzero_emitted = false;
    for (std::size_t atom_index = 0;
         atom_index < conf.AtomCount(); ++atom_index) {
        std::array<double, 9> expected{};
        const SphericalTensor& production = conf.AtomAt(atom_index)
            .mcconnell_source_tensors[category][channel];
        production.PackFull9(expected.data());
        for (size_t component = 0; component < expected.size();
             ++component) {
            const double emitted = bo_rows.values[atom_index * 9 + component];
            EXPECT_TRUE(std::isfinite(emitted));
            EXPECT_DOUBLE_EQ(emitted, expected[component]);
            any_nonzero_emitted = any_nonzero_emitted ||
                std::abs(emitted) > 1e-15;
        }
        EXPECT_TRUE(std::isfinite(audit_rows.values[atom_index * 4 + 1]));
        EXPECT_DOUBLE_EQ(audit_rows.values[atom_index * 4 + 1],
                         production.T2Magnitude());
    }
    EXPECT_TRUE(any_nonzero_emitted);

    RemoveSidechainCOOutputs(output_dir);
}


TEST(MopacExternalDirectionalFreeze,
     RealExecutableRerunO3NpyH5AndNonClosedAoDiagnostics) {
    using nmr::test::directional::Axial;
    using nmr::test::directional::EvenRank2;
    using nmr::test::directional::Near;
    using nmr::test::directional::NearMatrix;
    using nmr::test::directional::NearVector;
    using nmr::test::directional::Polar;
    using nmr::test::directional::Position;
    using nmr::test::directional::Positions;
    using nmr::test::directional::RotateNativeT2;
    using nmr::test::directional::SeededTransform;
    using nmr::test::directional::T1Vector;

    // Recorded forcing seed and tolerances.  The transformed reruns on this
    // protonated 54-atom peptide window establish that production PM7+MOZYME
    // has deterministic axis/origin sensitivity in its localized-orbital
    // SCF, beyond printed decimal quantisation.  The explicit envelopes below
    // name and retain that external-source numerical limitation; exact
    // NPY/H5 serialization and owner geometry still use bitwise/near-roundoff
    // checks.
    constexpr std::uint64_t kSeed = 0x4D4F5041434F3355ULL;
    constexpr double kBondOrderAbsTolerance = 6.0e-3;
    constexpr double kTensorAbsTolerance = 3.0e-3;
    constexpr double kTensorRelTolerance = 1.5e-2;
    constexpr double kDipoleAbsToleranceDebye = 3.0e-1;
    constexpr double kDipoleRelTolerance = 1.0e-2;
    constexpr double kHeatAbsToleranceKcalMol = 1.0e-2;
    constexpr double kAOShellAbsToleranceElectron = 7.0e-3;
    constexpr double kGeometryAbsTolerance = 3.0e-11;

    const auto proper = SeededTransform(kSeed, false);
    const auto improper = SeededTransform(kSeed, true);
    const std::array<nmr::test::directional::OrthogonalTransform, 2>
        transforms{{proper, improper}};

    auto protein = BuildProtonatedAsnWindowProtein();
    ASSERT_NE(protein, nullptr)
        << "could not build checked-in protonated ASN forcing window";
    ASSERT_GT(protein->AtomCount(), 9u);
    ASSERT_NE(FirstBondWithCategory(*protein, BondCategory::SidechainCO),
              SIZE_MAX);
    ASSERT_NE(FirstBondWithCategory(*protein, BondCategory::Aromatic),
              SIZE_MAX);
    const std::vector<Vec3> source_positions =
        protein->Conformation().Positions();
    protein->AddConformation(
        Positions(proper, source_positions),
        "real MOPAC proper transformed-input rerun");
    protein->AddConformation(
        Positions(improper, source_positions),
        "real MOPAC improper transformed-input rerun");

    const fs::path output_root = fs::temp_directory_path() /
        ("mopac_external_directional_freeze_" +
         std::to_string(::getpid()));
    RemoveTempTree(output_root);
    fs::create_directories(output_root);

    std::array<std::vector<double>, 3> topology_bond_orders;
    std::array<Vec3, 3> molecular_dipoles;
    std::array<double, 3> heats{};
    std::array<NpyArray<double>, 3> global_npys;
    std::array<NpyArray<double>, 3> atom_population_npys;
    std::array<NpyArray<double>, 3> ao_npys;
    std::array<NpyArray<double>, 3> ao_total_npys;
    std::array<NpyArray<double>, 3> sidechain_bo_npys;
    std::array<NpyArray<double>, 3> sidechain_frame_npys;
    std::array<NpyArray<double>, 3> sidechain_quality_npys;
    std::array<NpyArray<double>, 3> aromatic_fixed_npys;
    std::array<NpyArray<double>, 3> mopac_aromatic_sum_npys;
    std::array<NpyArray<double>, 3> mopac_aromatic_total_npys;
    std::array<NpyArray<double>, 3> mopac_shielding_npys;
    std::array<IndependentAromaticOracle, 3> aromatic_oracles;

    std::array<std::string, kMcConnellSourceCategoryCount>
        bo_filenames{};
    for (std::size_t category = 0;
         category < kMcConnellSourceCategoryCount; ++category) {
        bo_filenames[category] = std::string("mc_") +
            McConnellSourceCategoryStem(
                static_cast<McConnellSourceCategory>(category)) +
            "_bo.npy";
    }
    std::array<std::array<NpyArray<double>,
                          kMcConnellSourceCategoryCount>, 3>
        category_bo_npys;

    double max_aromatic_expected_fixed = 0.0;
    double max_aromatic_expected_bo = 0.0;
    double max_aromatic_fixed_oracle_error = 0.0;
    double max_aromatic_bo_oracle_error = 0.0;
    double max_aromatic_fixed_aggregate_error = 0.0;
    double max_aromatic_bo_aggregate_error = 0.0;

    for (std::size_t frame = 0; frame < 3; ++frame) {
        ProteinConformation& conf = protein->ConformationAt(frame);
        ASSERT_TRUE(conf.AttachResult(GeometryResult::Compute(conf)));
        ASSERT_TRUE(conf.AttachResult(SpatialIndexResult::Compute(conf)));

        auto mopac = MopacResult::Compute(conf, 0, 1);
        ASSERT_NE(mopac, nullptr)
            << "configured real MOPAC failed on transformed frame " << frame;
        topology_bond_orders[frame].reserve(protein->BondCount());
        for (std::size_t bond = 0; bond < protein->BondCount(); ++bond)
            topology_bond_orders[frame].push_back(
                mopac->TopologyBondOrder(bond));
        molecular_dipoles[frame] = mopac->Dipole();
        heats[frame] = mopac->HeatOfFormation();
        ASSERT_TRUE(conf.AttachResult(std::move(mopac)));
        ASSERT_TRUE(conf.AttachResult(McConnellResult::Compute(conf)));
        ASSERT_TRUE(conf.AttachResult(MopacMcConnellResult::Compute(conf)));
        ASSERT_TRUE(conf.AttachResult(
            SidechainCarbonylAnisotropyResult::Compute(conf)));

        ASSERT_EQ(topology_bond_orders[frame].size(), protein->BondCount());
        aromatic_oracles[frame] = ComputeIndependentAromaticOracle(
            conf, &topology_bond_orders[frame]);
        ASSERT_GT(aromatic_oracles[frame].aromatic_bond_count, 0u);
        ASSERT_GT(aromatic_oracles[frame].retained_bo_bond_count, 0u)
            << "real MOPAC must retain at least one aromatic bond order";

        const size_t aromatic_category = static_cast<size_t>(
            McConnellSourceCategory::Aromatic);
        for (size_t atom_index = 0;
             atom_index < conf.AtomCount(); ++atom_index) {
            const ConformationAtom& atom = conf.AtomAt(atom_index);
            const Mat3& expected_fixed =
                aromatic_oracles[frame].fixed[atom_index];
            const Mat3& expected_bo =
                aromatic_oracles[frame].bond_order[atom_index];
            max_aromatic_expected_fixed = std::max(
                max_aromatic_expected_fixed, MaxAbs(expected_fixed));
            max_aromatic_expected_bo = std::max(
                max_aromatic_expected_bo, MaxAbs(expected_bo));

            const Mat3 fixed = atom.mcconnell_source_tensors
                [aromatic_category]
                [static_cast<size_t>(McConnellChannel::Fixed)]
                .Reconstruct();
            const Mat3 bo = atom.mcconnell_source_tensors
                [aromatic_category]
                [static_cast<size_t>(McConnellChannel::BondOrder)]
                .Reconstruct();
            max_aromatic_fixed_oracle_error = std::max(
                max_aromatic_fixed_oracle_error,
                MaxAbs(fixed - expected_fixed));
            max_aromatic_bo_oracle_error = std::max(
                max_aromatic_bo_oracle_error,
                MaxAbs(bo - expected_bo));

            const SphericalTensor expected_fixed_spherical =
                SphericalTensor::Decompose(expected_fixed);
            const SphericalTensor expected_bo_spherical =
                SphericalTensor::Decompose(expected_bo);
            EXPECT_NEAR(atom.mcconnell_aromatic_sum,
                        expected_fixed_spherical.T0, 1e-12);
            EXPECT_LT(MaxAbs(atom.T2_aromatic_total.Reconstruct() -
                             expected_fixed),
                      1e-12);
            EXPECT_NEAR(atom.mopac_mc_aromatic_sum,
                        expected_bo_spherical.T0, 1e-11);
            EXPECT_LT(MaxAbs(
                atom.mopac_mc_T2_aromatic_total.Reconstruct() - expected_bo),
                1e-11);

            const Mat3 fixed_residual =
                atom.mc_shielding_contribution.Reconstruct() -
                SumMcConnellCategoriesExcept(
                    atom, aromatic_category, McConnellChannel::Fixed);
            const Mat3 bo_residual =
                atom.mopac_mc_shielding_contribution.Reconstruct() -
                SumMcConnellCategoriesExcept(
                    atom, aromatic_category, McConnellChannel::BondOrder);
            max_aromatic_fixed_aggregate_error = std::max(
                max_aromatic_fixed_aggregate_error,
                MaxAbs(fixed_residual - expected_fixed));
            max_aromatic_bo_aggregate_error = std::max(
                max_aromatic_bo_aggregate_error,
                MaxAbs(bo_residual - expected_bo));
        }

        const fs::path output_dir = output_root /
            (frame == 0 ? "original" :
             frame == 1 ? "proper" : "improper");
        fs::create_directories(output_dir);
        EXPECT_GT(conf.Result<MopacResult>().WriteFeatures(
                      conf, output_dir.string()),
                  0);
        EXPECT_EQ(conf.Result<McConnellResult>().WriteFeatures(
                      conf, output_dir.string()),
                  28);
        EXPECT_EQ(conf.Result<MopacMcConnellResult>().WriteFeatures(
                      conf, output_dir.string()),
                  13);
        EXPECT_EQ(conf.Result<SidechainCarbonylAnisotropyResult>()
                      .WriteFeatures(conf, output_dir.string()),
                  6);

        global_npys[frame] = ReadNpy<double>(
            output_dir / "mopac_global.npy", "<f8");
        atom_population_npys[frame] = ReadNpy<double>(
            output_dir / "mopac_atom_populations.npy", "<f8");
        ao_npys[frame] = ReadNpy<double>(
            output_dir / "mopac_atomic_orbital_populations.npy", "<f8");
        ao_total_npys[frame] = ReadNpy<double>(
            output_dir / "mopac_atomic_orbital_population_totals.npy",
            "<f8");
        sidechain_bo_npys[frame] = ReadNpy<double>(
            output_dir / "sidechain_co_bo_T2.npy", "<f8");
        sidechain_frame_npys[frame] = ReadNpy<double>(
            output_dir / "sidechain_co_frame.npy", "<f8");
        sidechain_quality_npys[frame] = ReadNpy<double>(
            output_dir / "sidechain_co_frame_quality.npy", "<f8");
        aromatic_fixed_npys[frame] = ReadNpy<double>(
            output_dir / "mc_aromatic_fixed.npy", "<f8");
        mopac_aromatic_sum_npys[frame] = ReadNpy<double>(
            output_dir / "mopac_mc_aromatic_sum.npy", "<f8");
        mopac_aromatic_total_npys[frame] = ReadNpy<double>(
            output_dir / "mopac_mc_aromatic_total.npy", "<f8");
        mopac_shielding_npys[frame] = ReadNpy<double>(
            output_dir / "mopac_mc_shielding.npy", "<f8");
        for (std::size_t category = 0;
             category < kMcConnellSourceCategoryCount; ++category) {
            category_bo_npys[frame][category] = ReadNpy<double>(
                output_dir / bo_filenames[category], "<f8");
        }

        ASSERT_EQ(global_npys[frame].shape,
                  (std::vector<std::size_t>{4u}));
        EXPECT_DOUBLE_EQ(global_npys[frame].values[0], heats[frame]);
        EXPECT_DOUBLE_EQ(global_npys[frame].values[1],
                         molecular_dipoles[frame].x());
        EXPECT_DOUBLE_EQ(global_npys[frame].values[2],
                         molecular_dipoles[frame].y());
        EXPECT_DOUBLE_EQ(global_npys[frame].values[3],
                         molecular_dipoles[frame].z());

        // The live producer does not populate the intended per-atom dipole
        // contribution columns.  Their serialized NaNs are unavailability,
        // never a manufactured zero or a directional observation.
        ASSERT_EQ(atom_population_npys[frame].shape,
                  (std::vector<std::size_t>{conf.AtomCount(), 12u}));
        for (std::size_t atom = 0; atom < conf.AtomCount(); ++atom) {
            for (std::size_t component = 6; component <= 9; ++component) {
                EXPECT_TRUE(std::isnan(
                    atom_population_npys[frame]
                        .values[atom * 12u + component]));
            }
        }

        // Piece 10 (libmopac) removed the typed
        // MopacResult::AtomicOrbitalPopulations() accessor this block
        // cross-checked the emitted AO NPYs against. Value correctness vs
        // the full-precision/probe oracle + NaN semantics now live in
        // tests/validate_mopac_npy.py; here we verify directly from the
        // emitted arrays the invariant 10 owns — [s,p,d] shell totals equal
        // the shell-sums of the 9 diagonal populations, NaN-aware. (cf. the
        // test_mopac_result.cpp merge resolution.)
        ASSERT_EQ(ao_npys[frame].shape,
                  (std::vector<std::size_t>{conf.AtomCount(), 9u}));
        ASSERT_EQ(ao_total_npys[frame].shape,
                  (std::vector<std::size_t>{conf.AtomCount(), 3u}));
        for (std::size_t row = 0; row < conf.AtomCount(); ++row) {
            const double* pop = &ao_npys[frame].values[row * 9u];
            const std::array<double, 3> expected_totals{{
                pop[0],
                pop[1] + pop[2] + pop[3],
                pop[4] + pop[5] + pop[6] + pop[7] + pop[8],
            }};
            for (std::size_t shell = 0; shell < 3; ++shell) {
                const double emitted =
                    ao_total_npys[frame].values[row * 3u + shell];
                if (std::isnan(expected_totals[shell]))
                    EXPECT_TRUE(std::isnan(emitted));
                else
                    EXPECT_DOUBLE_EQ(emitted, expected_totals[shell]);
            }
        }

        ASSERT_EQ(sidechain_bo_npys[frame].shape,
                  (std::vector<std::size_t>{conf.AtomCount(), 9u}));
        const std::size_t sidechain_category = static_cast<std::size_t>(
            McConnellSourceCategory::SidechainCO);
        ASSERT_EQ(category_bo_npys[frame][sidechain_category].shape,
                  sidechain_bo_npys[frame].shape);
        EXPECT_EQ(category_bo_npys[frame][sidechain_category].values,
                  sidechain_bo_npys[frame].values)
            << "typed sidechain owner must serialize the same production "
               "BO tensor as mc_sidechain_co_bo.npy";

        ASSERT_EQ(aromatic_fixed_npys[frame].shape,
                  (std::vector<std::size_t>{conf.AtomCount(), 9u}));
        ASSERT_EQ(category_bo_npys[frame][aromatic_category].shape,
                  aromatic_fixed_npys[frame].shape);
        ASSERT_EQ(mopac_aromatic_sum_npys[frame].shape,
                  (std::vector<std::size_t>{conf.AtomCount()}));
        ASSERT_EQ(mopac_aromatic_total_npys[frame].shape,
                  aromatic_fixed_npys[frame].shape);
        ASSERT_EQ(mopac_shielding_npys[frame].shape,
                  aromatic_fixed_npys[frame].shape);
        for (size_t atom_index = 0;
             atom_index < conf.AtomCount(); ++atom_index) {
            std::array<double, 9> expected_fixed{};
            std::array<double, 9> expected_bo{};
            SphericalTensor::Decompose(
                aromatic_oracles[frame].fixed[atom_index])
                .PackFull9(expected_fixed.data());
            SphericalTensor::Decompose(
                aromatic_oracles[frame].bond_order[atom_index])
                .PackFull9(expected_bo.data());
            for (size_t component = 0; component < 9; ++component) {
                EXPECT_NEAR(
                    aromatic_fixed_npys[frame]
                        .values[atom_index * 9 + component],
                    expected_fixed[component], 1e-12);
                EXPECT_NEAR(
                    category_bo_npys[frame][aromatic_category]
                        .values[atom_index * 9 + component],
                    expected_bo[component], 1e-11);
                EXPECT_NEAR(
                    mopac_aromatic_total_npys[frame]
                        .values[atom_index * 9 + component],
                    expected_bo[component], 1e-11);
            }
            EXPECT_NEAR(
                mopac_aromatic_sum_npys[frame].values[atom_index],
                expected_bo[0], 1e-11);

            std::array<double, 9> expected_shielding{};
            conf.AtomAt(atom_index).mopac_mc_shielding_contribution
                .PackFull9(expected_shielding.data());
            for (size_t component = 0; component < 9; ++component) {
                EXPECT_DOUBLE_EQ(
                    mopac_shielding_npys[frame]
                        .values[atom_index * 9 + component],
                    expected_shielding[component]);
            }
        }
    }

    EXPECT_GT(max_aromatic_expected_fixed, 1e-8);
    EXPECT_GT(max_aromatic_expected_bo, 1e-8);
    EXPECT_LT(max_aromatic_fixed_oracle_error, 1e-12);
    EXPECT_LT(max_aromatic_bo_oracle_error, 1e-11);
    EXPECT_LT(max_aromatic_fixed_aggregate_error, 1e-12);
    EXPECT_LT(max_aromatic_bo_aggregate_error, 1e-11);

    double max_bond_order_error = 0.0;
    double max_tensor_error = 0.0;
    double max_t0_error = 0.0;
    double max_t1_error = 0.0;
    double max_t2_error = 0.0;
    double max_t1_signal = 0.0;
    bool any_nonzero_sidechain_bo = false;
    std::array<bool, kMcConnellSourceCategoryCount>
        category_has_nonzero_signal{};
    for (std::size_t moved_frame = 1; moved_frame < 3; ++moved_frame) {
        const auto& transform = transforms[moved_frame - 1u];
        for (std::size_t bond = 0; bond < protein->BondCount(); ++bond) {
            const double error = std::abs(
                topology_bond_orders[moved_frame][bond] -
                topology_bond_orders[0][bond]);
            max_bond_order_error = std::max(max_bond_order_error, error);
            EXPECT_LE(error, kBondOrderAbsTolerance)
                << "MOPAC Wiberg order bond=" << bond
                << " frame=" << moved_frame;
        }

        for (std::size_t category = 0;
             category < kMcConnellSourceCategoryCount; ++category) {
            const auto& source = category_bo_npys[0][category];
            const auto& moved = category_bo_npys[moved_frame][category];
            ASSERT_EQ(source.shape,
                      (std::vector<std::size_t>{protein->AtomCount(), 9u}));
            ASSERT_EQ(moved.shape, source.shape);
            for (std::size_t atom = 0; atom < protein->AtomCount(); ++atom) {
                const SphericalTensor source_tensor =
                    UnpackFull9(&source.values[atom * 9u]);
                const SphericalTensor moved_tensor =
                    UnpackFull9(&moved.values[atom * 9u]);
                for (std::size_t component = 0; component < 9;
                     ++component) {
                    category_has_nonzero_signal[category] =
                        category_has_nonzero_signal[category] ||
                        std::abs(source.values[atom * 9u + component]) >
                            1.0e-12;
                }
                const Mat3 expected_matrix =
                    EvenRank2(transform, source_tensor.Reconstruct());
                const double matrix_error =
                    (moved_tensor.Reconstruct() - expected_matrix).norm();
                max_tensor_error = std::max(max_tensor_error, matrix_error);
                EXPECT_TRUE(NearMatrix(
                    moved_tensor.Reconstruct(), expected_matrix,
                    kTensorAbsTolerance, kTensorRelTolerance))
                    << bo_filenames[category] << " atom=" << atom
                    << " moved_frame=" << moved_frame;

                max_t0_error = std::max(
                    max_t0_error,
                    std::abs(moved_tensor.T0 - source_tensor.T0));
                EXPECT_TRUE(Near(moved_tensor.T0, source_tensor.T0,
                                 kTensorAbsTolerance,
                                 kTensorRelTolerance));

                const Vec3 expected_t1 = Axial(
                    transform, T1Vector(source_tensor));
                max_t1_error = std::max(
                    max_t1_error,
                    (T1Vector(moved_tensor) - expected_t1).norm());
                max_t1_signal = std::max(
                    max_t1_signal, T1Vector(source_tensor).norm());
                EXPECT_TRUE(NearVector(
                    T1Vector(moved_tensor), expected_t1,
                    kTensorAbsTolerance, kTensorRelTolerance));

                const SphericalTensor expected_t2 =
                    RotateNativeT2(transform, source_tensor);
                for (std::size_t component = 0; component < 5;
                     ++component) {
                    max_t2_error = std::max(
                        max_t2_error,
                        std::abs(moved_tensor.T2[component] -
                                 expected_t2.T2[component]));
                    EXPECT_TRUE(Near(
                        moved_tensor.T2[component],
                        expected_t2.T2[component],
                        kTensorAbsTolerance, kTensorRelTolerance))
                        << bo_filenames[category] << " atom=" << atom
                        << " native_T2=" << component
                        << " moved_frame=" << moved_frame;
                }

                if (category == static_cast<std::size_t>(
                                    McConnellSourceCategory::SidechainCO)) {
                    for (std::size_t component = 0; component < 9;
                         ++component) {
                        any_nonzero_sidechain_bo =
                            any_nonzero_sidechain_bo ||
                            std::abs(source.values[atom * 9u + component]) >
                                1.0e-12;
                    }
                }
            }
        }

        // Sidechain source-frame parity comes from the owning cross-product
        // body: origin is a position, x/y are polar, and z is axial.
        const auto& source_frames = sidechain_frame_npys[0];
        const auto& moved_frames = sidechain_frame_npys[moved_frame];
        ASSERT_EQ(source_frames.shape, moved_frames.shape);
        ASSERT_EQ(source_frames.shape.size(), 2u);
        ASSERT_EQ(source_frames.shape[1], 12u);
        for (std::size_t row = 0; row < source_frames.shape[0]; ++row) {
            const double* s = &source_frames.values[row * 12u];
            const double* m = &moved_frames.values[row * 12u];
            EXPECT_TRUE(NearVector(Vec3(m[0], m[1], m[2]),
                                   Position(transform,
                                            Vec3(s[0], s[1], s[2])),
                                   kGeometryAbsTolerance, 0.0));
            EXPECT_TRUE(NearVector(Vec3(m[3], m[4], m[5]),
                                   Polar(transform,
                                         Vec3(s[3], s[4], s[5])),
                                   kGeometryAbsTolerance, 0.0));
            EXPECT_TRUE(NearVector(Vec3(m[6], m[7], m[8]),
                                   Polar(transform,
                                         Vec3(s[6], s[7], s[8])),
                                   kGeometryAbsTolerance, 0.0));
            EXPECT_TRUE(NearVector(Vec3(m[9], m[10], m[11]),
                                   Axial(transform,
                                         Vec3(s[9], s[10], s[11])),
                                   kGeometryAbsTolerance, 0.0));
        }
        ASSERT_EQ(sidechain_quality_npys[moved_frame].shape,
                  sidechain_quality_npys[0].shape);
        EXPECT_LE(MaxFiniteDifference(
                      sidechain_quality_npys[moved_frame].values,
                      sidechain_quality_npys[0].values),
                  kGeometryAbsTolerance);
    }
    EXPECT_TRUE(any_nonzero_sidechain_bo)
        << "sidechain_co_bo_T2.npy must force a finite nonzero external BO";
    EXPECT_TRUE(category_has_nonzero_signal[static_cast<size_t>(
        McConnellSourceCategory::Aromatic)])
        << "real MOPAC fixture must force a finite nonzero aromatic BO channel";
    EXPECT_GT(max_t1_signal, 1.0e-10)
        << "fixture must exercise the axial T1 branch of nonsymmetric DQ";

    // The molecule is neutral (net_charge=0), so the molecular dipole is a
    // polar vector and is translation invariant.  A charged calculation
    // would instead obey p' = Qp + q_total*t and is outside this assertion.
    ASSERT_GT(molecular_dipoles[0].norm(), 5.0e-2)
        << "mopac_global.npy dipole parser did not produce a finite signal";
    double max_dipole_error = 0.0;
    double max_heat_error = 0.0;
    for (std::size_t moved_frame = 1; moved_frame < 3; ++moved_frame) {
        const Vec3 expected = Polar(
            transforms[moved_frame - 1u], molecular_dipoles[0]);
        max_dipole_error = std::max(
            max_dipole_error,
            (molecular_dipoles[moved_frame] - expected).norm());
        EXPECT_TRUE(NearVector(
            molecular_dipoles[moved_frame], expected,
            kDipoleAbsToleranceDebye, kDipoleRelTolerance));
        max_heat_error = std::max(
            max_heat_error, std::abs(heats[moved_frame] - heats[0]));
        EXPECT_NEAR(heats[moved_frame], heats[0],
                    kHeatAbsToleranceKcalMol);
    }

    // Printed AO populations are diagonal lab-axis occupations, not vector
    // components.  Their p/d shell traces are closed scalars; the nine raw
    // cells are not O(3)-closed because the omitted off-diagonal density
    // coherences are required to rotate the diagonal.  The external rerun
    // must therefore show orientation dependence without assigning a false
    // vector/tensor law to mopac_atomic_orbital_populations.npy.
    double max_ao_raw_change = 0.0;
    double max_ao_shell_error = 0.0;
    double min_false_polar_error = std::numeric_limits<double>::infinity();
    for (std::size_t moved_frame = 1; moved_frame < 3; ++moved_frame) {
        ASSERT_EQ(ao_npys[moved_frame].shape, ao_npys[0].shape);
        ASSERT_EQ(ao_total_npys[moved_frame].shape,
                  ao_total_npys[0].shape);
        max_ao_raw_change = std::max(
            max_ao_raw_change,
            MaxFiniteDifference(ao_npys[moved_frame].values,
                                ao_npys[0].values));
        max_ao_shell_error = std::max(
            max_ao_shell_error,
            MaxFiniteDifference(ao_total_npys[moved_frame].values,
                                ao_total_npys[0].values));

        const std::size_t rows = ao_npys[0].shape[0];
        for (std::size_t row = 0; row < rows; ++row) {
            const Vec3 source_p(
                ao_npys[0].values[row * 9u + 1u],
                ao_npys[0].values[row * 9u + 2u],
                ao_npys[0].values[row * 9u + 3u]);
            const Vec3 moved_p(
                ao_npys[moved_frame].values[row * 9u + 1u],
                ao_npys[moved_frame].values[row * 9u + 2u],
                ao_npys[moved_frame].values[row * 9u + 3u]);
            if (!source_p.allFinite() || !moved_p.allFinite()) continue;
            min_false_polar_error = std::min(
                min_false_polar_error,
                (moved_p - Polar(transforms[moved_frame - 1u],
                                 source_p)).norm());
        }
    }
    EXPECT_GT(max_ao_raw_change, 1.0e-4)
        << "fixture must demonstrate lab-axis AO orientation dependence";
    EXPECT_LE(max_ao_shell_error, kAOShellAbsToleranceElectron)
        << "only [s_total,p_total,d_total] is the closed companion";
    EXPECT_GT(min_false_polar_error, 1.0e-2)
        << "px/py/pz occupations must not be mislabeled as a polar vector";

    // Feed the three independently rerun owner stacks into the real sparse
    // trajectory source and validate the exact serialized H5 payload.
    auto trajectory_protein =
        TrajectoryProtein::CreateForTesting(std::move(protein));
    ASSERT_NE(trajectory_protein, nullptr);
    auto time_series =
        MopacMcConnellShieldingTimeSeriesTrajectoryResult::Create(
            *trajectory_protein);
    ASSERT_NE(time_series, nullptr);
    Trajectory dummy("", "", "");
    for (std::size_t frame = 0; frame < 3; ++frame) {
        time_series->Compute(
            trajectory_protein->ProteinRef().ConformationAt(frame),
            *trajectory_protein, dummy, 700u + frame, 0.25 * frame);
    }
    time_series->Finalize(*trajectory_protein, dummy);
    const fs::path h5_path = output_root / "mopac_mc_directional.h5";
    {
        HighFive::File file(h5_path.string(), HighFive::File::Truncate);
        time_series->WriteH5Group(*trajectory_protein, file);
    }

    const std::string h5_group =
        "/trajectory/mopac_mc_shielding_time_series/";
    std::vector<std::size_t> h5_dimensions;
    const std::vector<double> h5_xyz = ReadH5Flat<double>(
        h5_path, h5_group + "xyz", &h5_dimensions);
    ASSERT_EQ(h5_dimensions,
              (std::vector<std::size_t>{
                  trajectory_protein->AtomCount(), 3u, 9u}));
    EXPECT_EQ(ReadH5Flat<std::uint8_t>(
                  h5_path, h5_group + "source_attached_per_frame"),
              (std::vector<std::uint8_t>{1u, 1u, 1u}));
    EXPECT_EQ(ReadH5Flat<std::uint64_t>(
                  h5_path, h5_group + "frame_indices"),
              (std::vector<std::uint64_t>{700u, 701u, 702u}));
    EXPECT_EQ(ReadH5Flat<double>(h5_path, h5_group + "frame_times"),
              (std::vector<double>{0.0, 0.25, 0.5}));

    double max_h5_exact_error = 0.0;
    double max_h5_covariance_error = 0.0;
    double max_h5_aromatic_residual_error = 0.0;
    const size_t aromatic_category = static_cast<size_t>(
        McConnellSourceCategory::Aromatic);
    for (std::size_t atom = 0;
         atom < trajectory_protein->AtomCount(); ++atom) {
        std::array<double, 9> source_pack{};
        trajectory_protein->ProteinRef().ConformationAt(0)
            .AtomAt(atom).mopac_mc_shielding_contribution
            .PackFull9(source_pack.data());
        const SphericalTensor source_tensor =
            UnpackFull9(source_pack.data());
        for (std::size_t frame = 0; frame < 3; ++frame) {
            std::array<double, 9> expected_pack{};
            trajectory_protein->ProteinRef().ConformationAt(frame)
                .AtomAt(atom).mopac_mc_shielding_contribution
                .PackFull9(expected_pack.data());
            const std::size_t base = (atom * 3u + frame) * 9u;
            for (std::size_t component = 0; component < 9; ++component) {
                max_h5_exact_error = std::max(
                    max_h5_exact_error,
                    std::abs(h5_xyz[base + component] -
                             expected_pack[component]));
                EXPECT_DOUBLE_EQ(h5_xyz[base + component],
                                 expected_pack[component]);
            }
            const SphericalTensor emitted = UnpackFull9(&h5_xyz[base]);
            const ConformationAtom& frame_atom =
                trajectory_protein->ProteinRef().ConformationAt(frame)
                    .AtomAt(atom);
            const Mat3 aromatic_residual = emitted.Reconstruct() -
                SumMcConnellCategoriesExcept(
                    frame_atom, aromatic_category,
                    McConnellChannel::BondOrder);
            max_h5_aromatic_residual_error = std::max(
                max_h5_aromatic_residual_error,
                MaxAbs(aromatic_residual -
                       aromatic_oracles[frame].bond_order[atom]));
            if (frame > 0) {
                const Mat3 expected = EvenRank2(
                    transforms[frame - 1u],
                    source_tensor.Reconstruct());
                max_h5_covariance_error = std::max(
                    max_h5_covariance_error,
                    (emitted.Reconstruct() - expected).norm());
                EXPECT_TRUE(NearMatrix(
                    emitted.Reconstruct(), expected,
                    kTensorAbsTolerance, kTensorRelTolerance));
            }
        }
    }
    EXPECT_LT(max_h5_aromatic_residual_error, 1e-11);

    std::string nonzero_categories;
    std::string zero_categories;
    for (std::size_t category = 0;
         category < kMcConnellSourceCategoryCount; ++category) {
        std::string& destination = category_has_nonzero_signal[category]
            ? nonzero_categories : zero_categories;
        if (!destination.empty()) destination += ',';
        destination += McConnellSourceCategoryStem(
            static_cast<McConnellSourceCategory>(category));
    }

    std::cout
        << "MOPAC external O(3) seed=0x" << std::hex << kSeed << std::dec
        << " max_BO_abs=" << max_bond_order_error
        << " tol_BO_abs=" << kBondOrderAbsTolerance
        << " max_tensor_Frobenius=" << max_tensor_error
        << " tensor_abs=" << kTensorAbsTolerance
        << " tensor_rel=" << kTensorRelTolerance
        << " max_T0_abs=" << max_t0_error
        << " max_T1_L2=" << max_t1_error
        << " max_native_T2_abs=" << max_t2_error
        << " aromatic_bonds=" << aromatic_oracles[0].aromatic_bond_count
        << " retained_aromatic_BO_bonds="
        << aromatic_oracles[0].retained_bo_bond_count
        << " max_aromatic_fixed=" << max_aromatic_expected_fixed
        << " max_aromatic_BO=" << max_aromatic_expected_bo
        << " max_aromatic_fixed_oracle_error="
        << max_aromatic_fixed_oracle_error
        << " max_aromatic_BO_oracle_error="
        << max_aromatic_bo_oracle_error
        << " max_aromatic_fixed_aggregate_error="
        << max_aromatic_fixed_aggregate_error
        << " max_aromatic_BO_aggregate_error="
        << max_aromatic_bo_aggregate_error
        << " max_dipole_L2_D=" << max_dipole_error
        << " dipole_abs_D=" << kDipoleAbsToleranceDebye
        << " dipole_rel=" << kDipoleRelTolerance
        << " max_heat_abs_kcal_mol=" << max_heat_error
        << " heat_abs=" << kHeatAbsToleranceKcalMol
        << " max_AO_raw_change_e=" << max_ao_raw_change
        << " max_AO_shell_abs_e=" << max_ao_shell_error
        << " AO_shell_abs=" << kAOShellAbsToleranceElectron
        << " min_false_polar_L2_e=" << min_false_polar_error
        << " max_H5_exact_abs=" << max_h5_exact_error
        << " max_H5_cov_Frobenius=" << max_h5_covariance_error
        << " max_H5_aromatic_residual="
        << max_h5_aromatic_residual_error
        << " nonzero_BO_categories=" << nonzero_categories
        << " zero_BO_categories=" << zero_categories
        << std::endl;

    RemoveTempTree(output_root);
}


TEST(McConnellImplementationChecks, PeptideCORhombicFallsBackToAxial) {
    {
        auto protein = BuildSyntheticPeptideCOProtein();
        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        const size_t nc_bond = FirstBondBetween(*protein, 0, 1);
        ASSERT_NE(nc_bond, SIZE_MAX);
        ASSERT_NE(protein->BondAt(nc_bond).category, BondCategory::PeptideCO);

        const auto axial = McConnellResult::ComputePairKernel(
            conf.PositionAt(3),
            conf.bond_midpoints[nc_bond],
            conf.bond_directions[nc_bond]);
        const auto rhombic = McConnellResult::ComputePeptideCORhombicPairKernel(
            conf, nc_bond, conf.PositionAt(3));

        EXPECT_LT(MaxAbs(rhombic.source_shape - axial.source_shape), 1e-12);
        EXPECT_LT(MaxAbs(rhombic.response - axial.response), 1e-12);
    }

    {
        auto protein = BuildSyntheticPeptideCOProtein(true);
        auto& conf = protein->Conformation();
        conf.AttachResult(GeometryResult::Compute(conf));
        const size_t co_bond = FirstBondWithCategory(*protein, BondCategory::PeptideCO);
        ASSERT_NE(co_bond, SIZE_MAX);

        const auto axial = McConnellResult::ComputePairKernel(
            conf.PositionAt(3),
            conf.bond_midpoints[co_bond],
            conf.bond_directions[co_bond]);
        const auto rhombic = McConnellResult::ComputePeptideCORhombicPairKernel(
            conf, co_bond, conf.PositionAt(3));

        EXPECT_LT(MaxAbs(rhombic.source_shape - axial.source_shape), 1e-12);
        EXPECT_LT(MaxAbs(rhombic.response - axial.response), 1e-12);

        conf.AttachResult(SpatialIndexResult::Compute(conf));
        conf.AttachResult(McConnellResult::Compute(conf));
        EXPECT_LT(MaxAbs(conf.AtomAt(3)
                             .mcconnell_peptide_co_rhombic.Reconstruct()),
                  1e-12);
    }
}


TEST_F(McConnellProteinTest,
       AromaticCategoryMatchesIndependentBondKernelSum) {
    auto& conf = protein->Conformation();
    const IndependentAromaticOracle oracle =
        ComputeIndependentAromaticOracle(conf);
    ASSERT_GT(oracle.aromatic_bond_count, 0u);
    ASSERT_TRUE(conf.AttachResult(McConnellResult::Compute(conf)));

    const size_t cat = static_cast<size_t>(
        McConnellSourceCategory::Aromatic);
    EXPECT_STREQ(McConnellSourceCategoryStem(
                     McConnellSourceCategory::Aromatic),
                 "aromatic");

    double max_expected = 0.0;
    double max_category_error = 0.0;
    double max_scalar_error = 0.0;
    double max_total_error = 0.0;
    double max_aggregate_residual_error = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const ConformationAtom& atom = conf.AtomAt(ai);
        const Mat3& expected = oracle.fixed[ai];
        max_expected = std::max(max_expected, MaxAbs(expected));

        const Mat3 category = atom.mcconnell_source_tensors[cat]
            [static_cast<size_t>(McConnellChannel::Fixed)].Reconstruct();
        max_category_error = std::max(
            max_category_error, MaxAbs(category - expected));

        const SphericalTensor expected_spherical =
            SphericalTensor::Decompose(expected);
        max_scalar_error = std::max(
            max_scalar_error,
            std::abs(atom.mcconnell_aromatic_sum - expected_spherical.T0));
        max_scalar_error = std::max(
            max_scalar_error,
            std::abs(conf.Result<McConnellResult>().CategoryScalarSum(
                         ai, BondCategory::Aromatic) -
                     expected_spherical.T0));
        max_total_error = std::max(
            max_total_error,
            MaxAbs(atom.T2_aromatic_total.Reconstruct() - expected));

        const Mat3 aromatic_residual =
            atom.mc_shielding_contribution.Reconstruct() -
            SumMcConnellCategoriesExcept(
                atom, cat, McConnellChannel::Fixed);
        max_aggregate_residual_error = std::max(
            max_aggregate_residual_error,
            MaxAbs(aromatic_residual - expected));
    }

    const Protein& topology = conf.ProteinRef();
    const size_t first_aromatic =
        FirstBondWithCategory(topology, BondCategory::Aromatic);
    ASSERT_NE(first_aromatic, SIZE_MAX);
    const Bond& forcing_bond = topology.BondAt(first_aromatic);
    const Vec3 forcing_a = conf.PositionAt(forcing_bond.atom_index_a);
    const Vec3 forcing_b = conf.PositionAt(forcing_bond.atom_index_b);
    const Vec3 forcing_u = (forcing_b - forcing_a).normalized();
    const Vec3 perpendicular = forcing_u.cross(
        std::abs(forcing_u.x()) < 0.9 ? Vec3::UnitX() : Vec3::UnitY())
        .normalized();
    const Vec3 sample_point =
        0.5 * (forcing_a + forcing_b) + 3.0 * perpendicular;

    const bool include_xh_sources =
        CalculatorConfig::Get("mcconnell_include_xh_sources") != 0.0;
    const double cutoff =
        CalculatorConfig::Get("mcconnell_bond_anisotropy_cutoff");
    const double singularity_guard =
        CalculatorConfig::Get("singularity_guard_distance");
    const double near_field_ratio =
        CalculatorConfig::Get("near_field_exclusion_ratio");
    Mat3 expected_aromatic_sample = Mat3::Zero();
    Mat3 unchanged_sample = Mat3::Zero();
    for (size_t bi = 0; bi < topology.BondCount(); ++bi) {
        const Bond& bond = topology.BondAt(bi);
        if (!include_xh_sources &&
            mcconnell_result_detail::IsXHBond(topology, bond)) {
            continue;
        }

        const Vec3 endpoint_a = conf.PositionAt(bond.atom_index_a);
        const Vec3 endpoint_b = conf.PositionAt(bond.atom_index_b);
        const Vec3 bond_vector = endpoint_b - endpoint_a;
        const double bond_length = bond_vector.norm();
        if (bond_length < 1e-15) continue;
        const Vec3 midpoint = 0.5 * (endpoint_a + endpoint_b);
        const Vec3 displacement = sample_point - midpoint;
        const double distance = displacement.norm();
        if (distance < singularity_guard || distance > cutoff ||
            distance < near_field_ratio * bond_length) {
            continue;
        }

        if (bond.category == BondCategory::Aromatic) {
            const Vec3 u = bond_vector / bond_length;
            const Vec3 n = displacement / distance;
            const Mat3 source_shape =
                u * u.transpose() - Mat3::Identity() / 3.0;
            const Mat3 dipolar =
                (3.0 * n * n.transpose() - Mat3::Identity()) /
                (distance * distance * distance);
            expected_aromatic_sample += dipolar * source_shape;
        } else {
            const McConnellPairKernel kernel =
                bond.category == BondCategory::PeptideCO
                ? McConnellResult::ComputePeptideCORhombicPairKernel(
                      conf, bi, sample_point)
                : McConnellResult::ComputePairKernel(
                      sample_point, midpoint, bond_vector / bond_length);
            unchanged_sample += kernel.response;
        }
    }
    const Mat3 sampled_aromatic_residual =
        conf.Result<McConnellResult>().SampleKernelAt(sample_point)
            .Reconstruct() -
        unchanged_sample;
    const double sample_error = MaxAbs(
        sampled_aromatic_residual - expected_aromatic_sample);

    EXPECT_GT(max_expected, 1e-8)
        << "fixture must force a real aromatic bond response";
    EXPECT_GT(MaxAbs(expected_aromatic_sample), 1e-8)
        << "sample point must force a real aromatic bond response";
    EXPECT_LT(max_category_error, 1e-12);
    EXPECT_LT(max_scalar_error, 1e-12);
    EXPECT_LT(max_total_error, 1e-12);
    EXPECT_LT(max_aggregate_residual_error, 1e-12);
    EXPECT_LT(sample_error, 1e-12);

    std::cout << "McConnell aromatic fixed oracle: bonds="
              << oracle.aromatic_bond_count
              << " max|expected|=" << max_expected
              << " max_category_error=" << max_category_error
              << " max_scalar_error=" << max_scalar_error
              << " max_total_error=" << max_total_error
              << " max_aggregate_residual_error="
              << max_aggregate_residual_error
              << " sample_error=" << sample_error << "\n";
}


TEST_F(McConnellProteinTest, MissingMopacZerosOnlyBoChannel) {
    auto& conf = protein->Conformation();
    ASSERT_FALSE(conf.HasResult<MopacResult>());
    conf.AttachResult(McConnellResult::Compute(conf));

    double max_fixed_t2 = 0.0;
    double max_bo_abs = 0.0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        for (size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
            max_fixed_t2 = std::max(
                max_fixed_t2,
                ca.mcconnell_source_tensors[c][0].T2Magnitude());
            const auto& bo = ca.mcconnell_source_tensors[c][1];
            max_bo_abs = std::max(max_bo_abs, std::abs(bo.T0));
            for (double v : bo.T1) max_bo_abs = std::max(max_bo_abs, std::abs(v));
            for (double v : bo.T2) max_bo_abs = std::max(max_bo_abs, std::abs(v));
        }
    }
    EXPECT_GT(max_fixed_t2, 0.0);
    EXPECT_EQ(max_bo_abs, 0.0);
}


TEST_F(McConnellProteinTest, NearFieldAuditCountsAreReported) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int accepted = 0;
    int rejected = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        accepted += conf.AtomAt(ai).mcconnell_near_field_accepted_lt3A;
        rejected += conf.AtomAt(ai).mcconnell_near_field_rejected_lt3A;
    }
    EXPECT_GE(accepted, 0);
    EXPECT_GE(rejected, 0);
    std::cout << "  McConnell near-field audit (<3 A): accepted="
              << accepted << " rejected=" << rejected << "\n";
}


TEST_F(McConnellProteinTest, WriteFeaturesEmitsTwentyEightArraysAndManifest) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    const fs::path out_dir = fs::temp_directory_path() /
        ("mcconnell_features_" + std::to_string(::getpid()));
    fs::create_directories(out_dir);
    const auto& mc = conf.Result<McConnellResult>();
    EXPECT_EQ(mc.WriteFeatures(conf, out_dir.string()), 28);

    for (size_t c = 0; c < kMcConnellSourceCategoryCount; ++c) {
        const auto cat = static_cast<McConnellSourceCategory>(c);
        for (size_t h = 0; h < kMcConnellChannelCount; ++h) {
            const auto channel = static_cast<McConnellChannel>(h);
            const fs::path p = out_dir / (
                std::string("mc_") + McConnellSourceCategoryStem(cat) +
                "_" + McConnellChannelStem(channel) + ".npy");
            EXPECT_TRUE(fs::exists(p)) << p;
        }
    }
    EXPECT_TRUE(fs::exists(out_dir / "mc_peptide_co_rhombic.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_dir.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_midpoint.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_T2.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_cn_T2.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_co_bond_index.npy"));
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearest_cn_bond_index.npy"));

    const fs::path manifest = out_dir / "extraction_manifest.json";
    ASSERT_TRUE(fs::exists(manifest));
    const std::string text = ReadText(manifest);
    EXPECT_NE(text.find("\"source_model\": \"unit susceptibility shape; axial scale learned; peptide C=O fixed/BO responses use the full pinned rhombic shape\""),
              std::string::npos);
    EXPECT_NE(text.find("\"bo_source\": \"MOPAC Wiberg bond order\""),
              std::string::npos);
    EXPECT_NE(text.find("\"aromatic_sources_included\": true"),
              std::string::npos);
    EXPECT_NE(text.find(
        "\"aromatic_source_model\": \"aromatic bond D(r)Qhat responses accumulate independently of BS/HM ring-current results\""),
        std::string::npos);
    EXPECT_EQ(text.find("aromatic_zeroed"), std::string::npos);
    EXPECT_NE(text.find(
        "\"tensor_basis\": \"project_native_full9_spherical_tensor_v1\""),
        std::string::npos);
    EXPECT_NE(text.find(
        "\"tensor_component_order\": \"T0,T1_x,T1_y,T1_z,T2_m-2,T2_m-1,T2_m0,T2_m+1,T2_m+2\""),
        std::string::npos);
    EXPECT_NE(text.find(
        "\"tensor_frame\": \"conformation_cartesian_xyz\""),
        std::string::npos);
    EXPECT_NE(text.find(
        "\"e3nn_export\": \"explicit project-basis to e3nn conversion required before use\""),
        std::string::npos);
    EXPECT_EQ(text.find("0e,1e_x,1e_y,1e_z,2e_m-2..+2"),
              std::string::npos)
        << "project-basis full9 arrays must never be advertised as e3nn irreps";
    EXPECT_NE(text.find("\"units\": \"Angstrom^-3\""),
              std::string::npos);
    EXPECT_NE(text.find("\"rhombic_status\": \"peptide_co_pinned_present\""),
              std::string::npos);
    EXPECT_NE(text.find("\"rhombic_array\": \"mc_peptide_co_rhombic.npy\""),
              std::string::npos);
    EXPECT_NE(text.find("\"source\": \"Hooper & Kaiser 1965 Table III, EF-corrected acetamide A, Abraham-anchored sign\""),
              std::string::npos);
    EXPECT_TRUE(fs::exists(out_dir / "mc_nearfield_counts.npy"));

    RemoveMcConnellOutputs(out_dir);
}


TEST_F(McConnellProteinTest, NearestBondTrackingWorks) {
    auto& conf = protein->Conformation();
    conf.AttachResult(McConnellResult::Compute(conf));

    int has_co = 0;
    int has_cn = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& atom = conf.AtomAt(ai);
        if (atom.nearest_CO_dist < NO_DATA_SENTINEL) {
            has_co++;
            EXPECT_GT(atom.nearest_CO_dist, 0.0);
            EXPECT_LT(atom.nearest_CO_dist, MCCONNELL_CUTOFF_A);
            ASSERT_LT(atom.nearest_CO_bond_index, protein->BondCount());
            EXPECT_EQ(mcconnell_result_detail::SourceCategory(
                          *protein,
                          protein->BondAt(atom.nearest_CO_bond_index)),
                      McConnellSourceCategory::PeptideCO);
        } else {
            EXPECT_EQ(atom.nearest_CO_bond_index, SIZE_MAX);
        }
        if (atom.nearest_CN_dist < NO_DATA_SENTINEL) {
            has_cn++;
            EXPECT_GT(atom.nearest_CN_dist, 0.0);
            EXPECT_LT(atom.nearest_CN_dist, MCCONNELL_CUTOFF_A);
            ASSERT_LT(atom.nearest_CN_bond_index, protein->BondCount());
            EXPECT_EQ(mcconnell_result_detail::SourceCategory(
                          *protein,
                          protein->BondAt(atom.nearest_CN_bond_index)),
                      McConnellSourceCategory::PeptideCN);
        } else {
            EXPECT_EQ(atom.nearest_CN_bond_index, SIZE_MAX);
        }
    }

    EXPECT_GT(has_co, static_cast<int>(conf.AtomCount() / 2))
        << "Most atoms should have a nearest CO bond within range";
    EXPECT_GT(has_cn, static_cast<int>(conf.AtomCount() / 2))
        << "Most atoms should have a nearest CN bond within range";

    std::cout << "  Atoms with nearest CO: " << has_co
              << " / " << conf.AtomCount()
              << "; nearest CN: " << has_cn
              << " / " << conf.AtomCount() << "\n";
}


// ============================================================================
// ORCA protein test (protonated, more atoms)
// ============================================================================

TEST(McConnellOrcaTest, RunOnProtonatedProtein) {
    OrcaRunFiles files;
    files.pdb_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.pdb";
    files.xyz_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.xyz";
    files.prmtop_path = std::string(nmr::test::TestEnvironment::OrcaDir()) + "A0A7C5FAR6_WT.prmtop";

    if (!fs::exists(files.xyz_path) || !fs::exists(files.prmtop_path))
        GTEST_SKIP() << "ORCA test data not found";

    auto load = BuildFromOrca(files);
    ASSERT_TRUE(load.Ok()) << load.error;

    auto& conf = load.protein->Conformation();
    conf.AttachResult(GeometryResult::Compute(conf));
    conf.AttachResult(SpatialIndexResult::Compute(conf));
    auto mc = McConnellResult::Compute(conf);
    ASSERT_NE(mc, nullptr);
    conf.AttachResult(std::move(mc));

    // Summary: bond neighbour counts, T0/T2 ranges
    double min_t0 = 1e30, max_t0 = -1e30;
    double max_t2 = 0;
    int total_bn = 0;
    for (size_t ai = 0; ai < conf.AtomCount(); ++ai) {
        const auto& ca = conf.AtomAt(ai);
        total_bn += static_cast<int>(ca.bond_neighbours.size());
        double t0 = ca.mc_shielding_contribution.T0;
        min_t0 = std::min(min_t0, t0);
        max_t0 = std::max(max_t0, t0);
        max_t2 = std::max(max_t2, ca.mc_shielding_contribution.T2Magnitude());
    }

    std::cout << "  ORCA protein McConnell summary:\n"
              << "    atoms=" << conf.AtomCount()
              << " bonds=" << load.protein->BondCount()
              << " total_bond_neighbours=" << total_bn << "\n"
              << "    T0 range: [" << min_t0 << ", " << max_t0 << "]\n"
              << "    max |T2| = " << max_t2 << "\n";

    EXPECT_GT(total_bn, 0);
    EXPECT_GT(max_t0 - min_t0, 0.001) << "Should have a range of T0 values";
    EXPECT_GT(max_t2, 0.001) << "T2 should be non-zero";
}
