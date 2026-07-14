//
// test_dihedral_bin_transition: discipline + integration for
// DihedralBinTransitionTrajectoryResult. AV companion to DihedralTS;
// stats-only (no per-frame bin labels — those are derivable from
// DihedralTS phi/psi/chi via the same binning function copied here).
//

#include "AminoAcidType.h"
#include "CalculatorConfig.h"
#include "ConformationAtom.h"
#include "DirectionalTestHelpers.h"
#include "DihedralBinTransitionTrajectoryResult.h"
#include "GeometryResult.h"
#include "OperationLog.h"
#include "PdbFileReader.h"
#include "Protein.h"
#include "ProteinConformation.h"
#include "Residue.h"
#include "RunConfiguration.h"
#include "Session.h"
#include "SpatialIndexResult.h"
#include "TestEnvironment.h"
#include "Trajectory.h"
#include "TrajectoryProtein.h"
#include "Types.h"

#include <gtest/gtest.h>
#include <highfive/H5DataSet.hpp>
#include <highfive/H5File.hpp>
#include <highfive/H5Group.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <functional>
#include <iostream>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <string>
#include <vector>

namespace fs = std::filesystem;

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

namespace {

constexpr const char* kFixtureProtein = "1P9J_5801";

std::string TrrPathFor(const std::string& p) {
    return fs::path(p).replace_extension(".trr").string();
}
std::string ProductionDirFor(const std::string& p) {
    return fs::path(p).parent_path().string();
}
bool FixtureAvailable(const nmr::test::AmberTrajectoryFixture& fix) {
    return !fix.tpr_path.empty() && fs::exists(fix.tpr_path)
        && fs::exists(TrrPathFor(fix.tpr_path)) && fs::exists(fix.edr_path);
}
nmr::RunConfiguration BuildConfig(unsigned stride) {
    nmr::RunConfiguration config;
    auto& opts = config.MutablePerFrameRunOptions();
    opts.skip_mopac = true; opts.skip_coulomb = true; opts.skip_apbs = true;
    opts.skip_dssp = true;
    config.RequireConformationResult(typeid(nmr::GeometryResult));
    config.RequireConformationResult(typeid(nmr::SpatialIndexResult));
    config.AddTrajectoryResultFactory(
        [](const nmr::TrajectoryProtein& tp_in)
        -> std::unique_ptr<nmr::TrajectoryResult> {
        return nmr::DihedralBinTransitionTrajectoryResult::Create(tp_in);
    });
    config.SetStride(stride);
    return config;
}

template <typename T>
std::vector<T> ReadDihedralBinFlat(const std::string& path,
                                   const std::string& dataset) {
    HighFive::File file(path, HighFive::File::ReadOnly);
    auto data_set = file.getDataSet(dataset);
    const auto dims = data_set.getSpace().getDimensions();
    const std::size_t count = std::accumulate(
        dims.begin(), dims.end(), std::size_t{1},
        std::multiplies<std::size_t>());
    std::vector<T> values(count);
    if (!values.empty()) data_set.read(values.data());
    return values;
}

std::string FreshDihedralBinPath(const std::string& stem) {
    const std::string path = nmr::test::TestEnvironment::TempPath(stem);
    (void)std::remove(path.c_str());
    return path;
}

std::uint8_t ReflectedChiBin(std::uint8_t source) {
    using R = nmr::DihedralBinTransitionTrajectoryResult;
    if (source == R::kChiBinGplus) return R::kChiBinGminus;
    if (source == R::kChiBinGminus) return R::kChiBinGplus;
    return source;  // trans and unassigned are fixed by chi -> -chi.
}

std::vector<nmr::Vec3> ChiEndpointFrame(double chi_radians) {
    // For these coordinates the production Dihedral() evaluates to
    // atan2(-y, 0.5).  Setting y=0 therefore exercises the calculator's
    // exact chi==0 endpoint; y=-0.5*tan(epsilon) gives +epsilon.
    const double y = chi_radians == 0.0
        ? 0.0
        : -0.5 * std::tan(chi_radians);
    return {
        nmr::Vec3(1.0, 0.0, 0.0),
        nmr::Vec3(0.0, 0.0, 0.0),
        nmr::Vec3(0.0, 0.0, 1.0),
        nmr::Vec3(0.5, y, 1.0),
    };
}

std::unique_ptr<nmr::Protein> BuildChiEndpointProtein(
        const std::vector<nmr::Vec3>& canonical_positions) {
    auto protein = std::make_unique<nmr::Protein>();
    nmr::Residue residue;
    residue.type = nmr::AminoAcid::Unknown;
    residue.chain_id = "A";
    residue.sequence_number = 1;
    const std::size_t residue_index = protein->AddResidue(std::move(residue));

    std::array<std::size_t, 4> chi_atoms{};
    for (std::size_t i = 0; i < chi_atoms.size(); ++i) {
        auto atom = nmr::Atom::Create(nmr::Element::C);
        atom->residue_index = residue_index;
        chi_atoms[i] = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(residue_index).atom_indices.push_back(
            chi_atoms[i]);
        protein->MutableResidueAt(residue_index).chi[0].a[i] = chi_atoms[i];
    }

    protein->FinalizeConstruction(canonical_positions);
    protein->AddConformation(canonical_positions,
                             "chi endpoint transition source");
    return protein;
}

std::vector<nmr::Vec3> BackboneForcingFrame(double phi_degrees,
                                             double psi_degrees) {
    // Five production positions in the order
    //   C(i-1), N(i), CA(i), C(i), N(i+1).
    // With N=(0,0,0), CA=(1,0,0), C=(1,1,0), the calculator's exact
    // atan2 construction evaluates to the requested phi/psi values.
    const double radians_per_degree = std::acos(-1.0) / 180.0;
    const double phi = phi_degrees * radians_per_degree;
    const double psi = psi_degrees * radians_per_degree;
    return {
        nmr::Vec3(0.0, std::cos(phi), std::sin(phi)),
        nmr::Vec3(0.0, 0.0, 0.0),
        nmr::Vec3(1.0, 0.0, 0.0),
        nmr::Vec3(1.0, 1.0, 0.0),
        nmr::Vec3(1.0 - std::cos(psi), 1.0, -std::sin(psi)),
    };
}

std::unique_ptr<nmr::Protein> BuildBackboneForcingProtein(
        const std::vector<nmr::Vec3>& canonical_positions) {
    auto protein = std::make_unique<nmr::Protein>();
    for (int sequence_number = 1; sequence_number <= 3; ++sequence_number) {
        nmr::Residue residue;
        residue.type = nmr::AminoAcid::Unknown;
        residue.chain_id = "A";
        residue.sequence_number = sequence_number;
        protein->AddResidue(std::move(residue));
    }

    auto add_atom = [&](std::size_t residue_index, nmr::Element element) {
        auto atom = nmr::Atom::Create(element);
        atom->residue_index = residue_index;
        const std::size_t atom_index = protein->AddAtom(std::move(atom));
        protein->MutableResidueAt(residue_index).atom_indices.push_back(
            atom_index);
        return atom_index;
    };

    protein->MutableResidueAt(0).C = add_atom(0, nmr::Element::C);
    protein->MutableResidueAt(1).N = add_atom(1, nmr::Element::N);
    protein->MutableResidueAt(1).CA = add_atom(1, nmr::Element::C);
    protein->MutableResidueAt(1).C = add_atom(1, nmr::Element::C);
    protein->MutableResidueAt(2).N = add_atom(2, nmr::Element::N);

    protein->FinalizeConstruction(canonical_positions);
    protein->AddConformation(canonical_positions,
                             "backbone transition forcing source");
    return protein;
}

}  // namespace


TEST(DihedralBinTransition,
     DirectionalCategoricalRerunO3SerializedH5) {
    using nmr::test::directional::Positions;
    using nmr::test::directional::SeededTransform;
    using R = nmr::DihedralBinTransitionTrajectoryResult;

    nmr::test::TestEnvironment::LoadCalculatorConfig();
    const std::string pdb = nmr::test::TestEnvironment::UbqProtonated();
    if (pdb.empty() || !fs::exists(pdb))
        GTEST_SKIP() << "1UBQ protonated fixture unavailable";

    auto source_loaded = nmr::BuildFromProtonatedPdb(pdb);
    ASSERT_TRUE(source_loaded.Ok()) << source_loaded.error;
    auto source_tp = nmr::TrajectoryProtein::CreateForTesting(
        std::move(source_loaded.protein));
    ASSERT_NE(source_tp, nullptr);
    const auto source_positions =
        source_tp->CanonicalConformation().Positions();
    auto source_conf = source_tp->TickConformation(source_positions);
    auto source_result = R::Create(*source_tp);
    ASSERT_NE(source_result, nullptr);
    nmr::Trajectory dummy("", "", "");
    source_result->Compute(*source_conf, *source_tp, dummy, 97, 6.0);
    source_result->Finalize(*source_tp, dummy);
    const std::string source_path =
        FreshDihedralBinPath("dihedral_bin_directional_source.h5");
    {
        HighFive::File file(source_path, HighFive::File::Truncate);
        source_result->WriteH5Group(*source_tp, file);
    }

    const std::string group = "/trajectory/dihedral_bin_transition/";
    {
        HighFive::File file(source_path, HighFive::File::ReadOnly);
        for (const char* name : {
                 "backbone_transition_count", "backbone_dominant_region",
                 "backbone_bin_occupancy", "chi_transition_count",
                 "chi_dominant_rotamer", "chi_rotamer_occupancy"}) {
            auto dataset = file.getDataSet(group + name);
            std::string frame;
            std::string parity;
            std::string law;
            dataset.getAttribute("coordinate_frame").read(frame);
            dataset.getAttribute("parity").read(parity);
            dataset.getAttribute("transformation").read(law);
            EXPECT_EQ(frame, "intrinsic_signed_dihedral_bins") << name;
            EXPECT_EQ(parity, "mixed") << name;
            EXPECT_NE(law.find("no global improper-transform map"),
                      std::string::npos) << name;
        }
        for (const char* name : {"n_frames_observed",
                                 "chi_n_frames_observed"}) {
            auto dataset = file.getDataSet(group + name);
            std::string parity;
            std::string law;
            dataset.getAttribute("parity").read(parity);
            dataset.getAttribute("transformation").read(law);
            EXPECT_EQ(parity, "even") << name;
            EXPECT_NE(law.find("reflection-invariant"),
                      std::string::npos) << name;
        }
    }
    const std::size_t residue_count = source_tp->ProteinRef().ResidueCount();
    const auto source_backbone_count = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "backbone_transition_count");
    const auto source_backbone_dominant = ReadDihedralBinFlat<std::uint8_t>(
        source_path, group + "backbone_dominant_region");
    const auto source_backbone_observed = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "n_frames_observed");
    const auto source_backbone_occupancy = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "backbone_bin_occupancy");
    const auto source_chi_count = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "chi_transition_count");
    const auto source_chi_dominant = ReadDihedralBinFlat<std::uint8_t>(
        source_path, group + "chi_dominant_rotamer");
    const auto source_chi_observed = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "chi_n_frames_observed");
    const auto source_chi_occupancy = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "chi_rotamer_occupancy");
    const auto source_frame_indices = ReadDihedralBinFlat<std::size_t>(
        source_path, group + "frame_indices");
    const auto source_frame_times = ReadDihedralBinFlat<double>(
        source_path, group + "frame_times");
    ASSERT_EQ(source_backbone_dominant.size(), residue_count);
    ASSERT_EQ(source_backbone_occupancy.size(),
              residue_count * R::kBackboneBinCount);
    ASSERT_EQ(source_chi_dominant.size(), residue_count * R::kChiCount);
    ASSERT_EQ(source_chi_occupancy.size(),
              residue_count * R::kChiCount * R::kChiBinCount);
    ASSERT_EQ(source_frame_indices, (std::vector<std::size_t>{97u}));
    ASSERT_EQ(source_frame_times, (std::vector<double>{6.0}));

    constexpr std::uint64_t kTransformSeed = 0x444948454452414CULL;
    std::size_t improper_backbone_changes = 0;
    std::array<std::set<std::uint8_t>, R::kBackboneBinCount>
        reflected_destinations_by_source_bin;

    for (const bool improper : {false, true}) {
        auto moved_loaded = nmr::BuildFromProtonatedPdb(pdb);
        ASSERT_TRUE(moved_loaded.Ok()) << moved_loaded.error;
        auto moved_tp = nmr::TrajectoryProtein::CreateForTesting(
            std::move(moved_loaded.protein));
        ASSERT_NE(moved_tp, nullptr);
        const auto transform = SeededTransform(kTransformSeed, improper);
        auto moved_conf = moved_tp->TickConformation(
            Positions(transform, source_positions));
        auto moved_result = R::Create(*moved_tp);
        ASSERT_NE(moved_result, nullptr);
        moved_result->Compute(*moved_conf, *moved_tp, dummy, 97, 6.0);
        moved_result->Finalize(*moved_tp, dummy);
        const std::string moved_path = FreshDihedralBinPath(
            improper ? "dihedral_bin_directional_improper.h5" :
                       "dihedral_bin_directional_proper.h5");
        {
            HighFive::File file(moved_path, HighFive::File::Truncate);
            moved_result->WriteH5Group(*moved_tp, file);
        }

        const auto moved_backbone_count =
            ReadDihedralBinFlat<std::uint32_t>(
                moved_path, group + "backbone_transition_count");
        const auto moved_backbone_dominant =
            ReadDihedralBinFlat<std::uint8_t>(
                moved_path, group + "backbone_dominant_region");
        const auto moved_backbone_observed =
            ReadDihedralBinFlat<std::uint32_t>(
                moved_path, group + "n_frames_observed");
        const auto moved_backbone_occupancy =
            ReadDihedralBinFlat<std::uint32_t>(
                moved_path, group + "backbone_bin_occupancy");
        const auto moved_chi_count = ReadDihedralBinFlat<std::uint32_t>(
            moved_path, group + "chi_transition_count");
        const auto moved_chi_dominant = ReadDihedralBinFlat<std::uint8_t>(
            moved_path, group + "chi_dominant_rotamer");
        const auto moved_chi_observed = ReadDihedralBinFlat<std::uint32_t>(
            moved_path, group + "chi_n_frames_observed");
        const auto moved_chi_occupancy = ReadDihedralBinFlat<std::uint32_t>(
            moved_path, group + "chi_rotamer_occupancy");

        if (!improper) {
            EXPECT_EQ(moved_backbone_count, source_backbone_count);
            EXPECT_EQ(moved_backbone_dominant,
                      source_backbone_dominant);
            EXPECT_EQ(moved_backbone_observed, source_backbone_observed);
            EXPECT_EQ(moved_backbone_occupancy,
                      source_backbone_occupancy);
            EXPECT_EQ(moved_chi_count, source_chi_count);
            EXPECT_EQ(moved_chi_dominant, source_chi_dominant);
            EXPECT_EQ(moved_chi_observed, source_chi_observed);
            EXPECT_EQ(moved_chi_occupancy, source_chi_occupancy);
        } else {
            // Signed phi/psi both negate under reflection. The literature
            // Rama regions are not sign-symmetric and a coarse source bin
            // does not determine one destination bin. Validate the fresh
            // serialized owner result without asserting a false invariant.
            EXPECT_EQ(moved_backbone_observed, source_backbone_observed);
            for (std::size_t ri = 0; ri < residue_count; ++ri) {
                const auto source_bin = source_backbone_dominant[ri];
                const auto moved_bin = moved_backbone_dominant[ri];
                ASSERT_LT(source_bin, R::kBackboneBinCount);
                ASSERT_LT(moved_bin, R::kBackboneBinCount);
                reflected_destinations_by_source_bin[source_bin].insert(
                    moved_bin);
                if (source_bin != moved_bin)
                    ++improper_backbone_changes;
                for (std::size_t bin = 0;
                     bin < R::kBackboneBinCount; ++bin) {
                    const auto value = moved_backbone_occupancy[
                        ri * R::kBackboneBinCount + bin];
                    EXPECT_EQ(value, bin == moved_bin ? 1u : 0u);
                }
            }

            // Away from the exact-zero endpoint, each chi angle negates and
            // the rotamer bins obey g+ <-> g-, trans fixed.  This one-frame
            // run cannot make a transition-count claim; the two-frame
            // endpoint forcing test below covers that non-closed case.
            EXPECT_EQ(moved_chi_observed, source_chi_observed);
            for (std::size_t row = 0;
                 row < residue_count * R::kChiCount; ++row) {
                EXPECT_EQ(moved_chi_dominant[row],
                          ReflectedChiBin(source_chi_dominant[row]));
                for (std::size_t source_bin = 0;
                     source_bin < R::kChiBinCount; ++source_bin) {
                    const std::size_t moved_bin =
                        ReflectedChiBin(static_cast<std::uint8_t>(
                            source_bin));
                    EXPECT_EQ(
                        moved_chi_occupancy[
                            row * R::kChiBinCount + moved_bin],
                        source_chi_occupancy[
                            row * R::kChiBinCount + source_bin]);
                }
            }
        }

        EXPECT_EQ(ReadDihedralBinFlat<std::size_t>(
                      moved_path, group + "frame_indices"),
                  source_frame_indices);
        EXPECT_EQ(ReadDihedralBinFlat<double>(
                      moved_path, group + "frame_times"),
                  source_frame_times);
        EXPECT_EQ(ReadDihedralBinFlat<std::uint8_t>(
                      moved_path, group + "source_attached_per_frame"),
                  ReadDihedralBinFlat<std::uint8_t>(
                      source_path, group + "source_attached_per_frame"));
        EXPECT_EQ(std::remove(moved_path.c_str()), 0);
    }

    EXPECT_GT(improper_backbone_changes, 0u);
    bool source_bin_has_multiple_reflected_destinations = false;
    for (const auto& destinations : reflected_destinations_by_source_bin) {
        if (destinations.size() > 1u)
            source_bin_has_multiple_reflected_destinations = true;
    }
    EXPECT_TRUE(source_bin_has_multiple_reflected_destinations)
        << "A single categorical permutation would be a false Rama law";
    EXPECT_EQ(std::remove(source_path.c_str()), 0);
}


TEST(DihedralBinTransition,
     BackboneTransitionCountImproperMergeRerunSerializedH5) {
    using nmr::test::directional::OrthogonalTransform;
    using nmr::test::directional::Positions;
    using R = nmr::DihedralBinTransitionTrajectoryResult;

    nmr::test::TestEnvironment::LoadCalculatorConfig();

    // Source: beta -> other.  Reflection negates both signed dihedrals,
    // mapping (+120,+170) and (+120,+120) to "other" and therefore
    // merging the source transition.  All angles are well inside their
    // bins, so this exercises the categorical physics rather than a
    // floating-point boundary.
    const std::array<std::vector<nmr::Vec3>, 2> source_frames = {
        BackboneForcingFrame(-120.0, -170.0),
        BackboneForcingFrame(-120.0, -120.0),
    };
    auto tp = nmr::TrajectoryProtein::CreateForTesting(
        BuildBackboneForcingProtein(source_frames[0]));
    ASSERT_NE(tp, nullptr);
    ASSERT_EQ(tp->ProteinRef().BackbonePredecessor(1),
              std::optional<std::size_t>(0));
    ASSERT_EQ(tp->ProteinRef().BackboneSuccessor(1),
              std::optional<std::size_t>(2));

    // Exact signed permutations preserve the requested dihedrals under
    // the proper transform and negate them under the improper transform,
    // while the nonzero translation checks coordinate-origin independence.
    OrthogonalTransform proper;
    proper.Q << 0.0, -1.0, 0.0,
                1.0,  0.0, 0.0,
                0.0,  0.0, 1.0;
    proper.translation = nmr::Vec3(4.0, -2.0, 7.0);
    OrthogonalTransform improper = proper;
    improper.Q.col(0) *= -1.0;
    ASSERT_DOUBLE_EQ(proper.Determinant(), 1.0);
    ASSERT_DOUBLE_EQ(improper.Determinant(), -1.0);

    const std::string group = "/trajectory/dihedral_bin_transition/";
    nmr::Trajectory dummy("", "", "");
    auto run_and_write = [&](const std::array<std::vector<nmr::Vec3>, 2>& frames,
                             const std::string& path) {
        auto result = R::Create(*tp);
        EXPECT_NE(result, nullptr);
        if (!result) return;
        for (std::size_t frame = 0; frame < frames.size(); ++frame) {
            auto conf = tp->TickConformation(frames[frame]);
            result->Compute(*conf, *tp, dummy, 211u + frame,
                            0.25 * static_cast<double>(frame));
        }
        result->Finalize(*tp, dummy);
        HighFive::File file(path, HighFive::File::Truncate);
        result->WriteH5Group(*tp, file);
    };

    const std::string source_path = FreshDihedralBinPath(
        "dihedral_bin_backbone_forcing_source.h5");
    const std::string proper_path = FreshDihedralBinPath(
        "dihedral_bin_backbone_forcing_proper.h5");
    const std::string improper_path = FreshDihedralBinPath(
        "dihedral_bin_backbone_forcing_improper.h5");
    run_and_write(source_frames, source_path);

    std::array<std::vector<nmr::Vec3>, 2> proper_frames;
    std::array<std::vector<nmr::Vec3>, 2> improper_frames;
    for (std::size_t frame = 0; frame < source_frames.size(); ++frame) {
        proper_frames[frame] = Positions(proper, source_frames[frame]);
        improper_frames[frame] = Positions(improper, source_frames[frame]);
    }
    run_and_write(proper_frames, proper_path);
    run_and_write(improper_frames, improper_path);

    const auto source_count = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "backbone_transition_count");
    const auto proper_count = ReadDihedralBinFlat<std::uint32_t>(
        proper_path, group + "backbone_transition_count");
    const auto improper_count = ReadDihedralBinFlat<std::uint32_t>(
        improper_path, group + "backbone_transition_count");
    ASSERT_EQ(source_count.size(), 3u);
    ASSERT_EQ(proper_count.size(), source_count.size());
    ASSERT_EQ(improper_count.size(), source_count.size());
    EXPECT_EQ(source_count[1], 1u);
    EXPECT_EQ(proper_count[1], 1u);
    EXPECT_EQ(improper_count[1], 0u);
    EXPECT_NE(improper_count[1], source_count[1]);

    const auto source_occupancy = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "backbone_bin_occupancy");
    const auto proper_occupancy = ReadDihedralBinFlat<std::uint32_t>(
        proper_path, group + "backbone_bin_occupancy");
    const auto improper_occupancy = ReadDihedralBinFlat<std::uint32_t>(
        improper_path, group + "backbone_bin_occupancy");
    ASSERT_EQ(source_occupancy.size(), 3u * R::kBackboneBinCount);
    EXPECT_EQ(proper_occupancy, source_occupancy);
    const std::size_t central_row = R::kBackboneBinCount;
    EXPECT_EQ(source_occupancy[central_row + R::kBinBeta], 1u);
    EXPECT_EQ(source_occupancy[central_row + R::kBinOther], 1u);
    EXPECT_EQ(improper_occupancy[central_row + R::kBinBeta], 0u);
    EXPECT_EQ(improper_occupancy[central_row + R::kBinOther], 2u);

    const auto source_observed = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "n_frames_observed");
    const auto proper_observed = ReadDihedralBinFlat<std::uint32_t>(
        proper_path, group + "n_frames_observed");
    const auto improper_observed = ReadDihedralBinFlat<std::uint32_t>(
        improper_path, group + "n_frames_observed");
    ASSERT_EQ(source_observed.size(), 3u);
    EXPECT_EQ(proper_observed, source_observed);
    EXPECT_EQ(improper_observed, source_observed);
    EXPECT_EQ(source_observed[1], 2u);

    EXPECT_EQ(std::remove(source_path.c_str()), 0);
    EXPECT_EQ(std::remove(proper_path.c_str()), 0);
    EXPECT_EQ(std::remove(improper_path.c_str()), 0);
}


TEST(DihedralBinTransition,
     ChiTransitionCountExactZeroEndpointRerunSerializedH5) {
    using nmr::test::directional::OrthogonalTransform;
    using nmr::test::directional::Positions;
    using R = nmr::DihedralBinTransitionTrajectoryResult;

    nmr::test::TestEnvironment::LoadCalculatorConfig();

    constexpr double kEpsilonRadians = 1.0e-4;
    const std::array<std::vector<nmr::Vec3>, 2> source_frames = {
        ChiEndpointFrame(0.0),
        ChiEndpointFrame(kEpsilonRadians),
    };
    auto tp = nmr::TrajectoryProtein::CreateForTesting(
        BuildChiEndpointProtein(source_frames[0]));
    ASSERT_NE(tp, nullptr);

    // Exact signed-permutation transforms keep the chi==0 forcing endpoint
    // at floating-point zero while still exercising rotation, reflection,
    // and translation of the calculator input.
    OrthogonalTransform proper;
    proper.Q << 0.0, -1.0, 0.0,
                1.0,  0.0, 0.0,
                0.0,  0.0, 1.0;
    proper.translation = nmr::Vec3(2.0, -3.0, 5.0);
    OrthogonalTransform improper = proper;
    improper.Q.col(0) *= -1.0;
    ASSERT_DOUBLE_EQ(proper.Determinant(), 1.0);
    ASSERT_DOUBLE_EQ(improper.Determinant(), -1.0);

    const std::string group = "/trajectory/dihedral_bin_transition/";
    nmr::Trajectory dummy("", "", "");
    auto run_and_write = [&](const std::array<std::vector<nmr::Vec3>, 2>& frames,
                             const std::string& path) {
        auto result = R::Create(*tp);
        EXPECT_NE(result, nullptr);
        if (!result) return;
        for (std::size_t frame = 0; frame < frames.size(); ++frame) {
            auto conf = tp->TickConformation(frames[frame]);
            result->Compute(*conf, *tp, dummy, 101u + frame,
                            0.5 * static_cast<double>(frame));
        }
        result->Finalize(*tp, dummy);
        HighFive::File file(path, HighFive::File::Truncate);
        result->WriteH5Group(*tp, file);
    };

    const std::string source_path =
        FreshDihedralBinPath("dihedral_bin_chi_endpoint_source.h5");
    const std::string proper_path =
        FreshDihedralBinPath("dihedral_bin_chi_endpoint_proper.h5");
    const std::string improper_path =
        FreshDihedralBinPath("dihedral_bin_chi_endpoint_improper.h5");
    run_and_write(source_frames, source_path);

    std::array<std::vector<nmr::Vec3>, 2> proper_frames;
    std::array<std::vector<nmr::Vec3>, 2> improper_frames;
    for (std::size_t frame = 0; frame < source_frames.size(); ++frame) {
        proper_frames[frame] = Positions(proper, source_frames[frame]);
        improper_frames[frame] = Positions(improper, source_frames[frame]);
    }
    run_and_write(proper_frames, proper_path);
    run_and_write(improper_frames, improper_path);

    const auto source_count = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "chi_transition_count");
    const auto proper_count = ReadDihedralBinFlat<std::uint32_t>(
        proper_path, group + "chi_transition_count");
    const auto improper_count = ReadDihedralBinFlat<std::uint32_t>(
        improper_path, group + "chi_transition_count");
    ASSERT_EQ(source_count.size(), R::kChiCount);
    ASSERT_EQ(proper_count.size(), source_count.size());
    ASSERT_EQ(improper_count.size(), source_count.size());

    // Source/proper: 0 (g-) -> +epsilon (g+) is one transition.
    // Improper: -0 (g-) -> -epsilon (g-) is no transition.  Thus the
    // serialized count is non-closed at the calculator's documented exact
    // endpoint and cannot be globally classified as reflection-invariant.
    EXPECT_EQ(source_count[0], 1u);
    EXPECT_EQ(proper_count[0], 1u);
    EXPECT_EQ(improper_count[0], 0u);
    EXPECT_NE(improper_count[0], source_count[0]);
    for (std::size_t chi = 1; chi < R::kChiCount; ++chi) {
        EXPECT_EQ(source_count[chi], 0u);
        EXPECT_EQ(proper_count[chi], 0u);
        EXPECT_EQ(improper_count[chi], 0u);
    }

    const auto source_occupancy = ReadDihedralBinFlat<std::uint32_t>(
        source_path, group + "chi_rotamer_occupancy");
    const auto proper_occupancy = ReadDihedralBinFlat<std::uint32_t>(
        proper_path, group + "chi_rotamer_occupancy");
    const auto improper_occupancy = ReadDihedralBinFlat<std::uint32_t>(
        improper_path, group + "chi_rotamer_occupancy");
    ASSERT_EQ(source_occupancy.size(), R::kChiCount * R::kChiBinCount);
    EXPECT_EQ(proper_occupancy, source_occupancy);
    EXPECT_EQ(source_occupancy[R::kChiBinGplus], 1u);
    EXPECT_EQ(source_occupancy[R::kChiBinGminus], 1u);
    EXPECT_EQ(improper_occupancy[R::kChiBinGplus], 0u);
    EXPECT_EQ(improper_occupancy[R::kChiBinGminus], 2u);

    EXPECT_EQ(std::remove(source_path.c_str()), 0);
    EXPECT_EQ(std::remove(proper_path.c_str()), 0);
    EXPECT_EQ(std::remove(improper_path.c_str()), 0);
}


TEST(DihedralBinTransition, Frame0Semantics) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::DihedralBinTransitionTrajectoryResult>();
    EXPECT_EQ(tr.NumFrames(), 1u);
}


TEST(DihedralBinTransition, FinalizeIdempotency) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    auto& tr = tp.Result<nmr::DihedralBinTransitionTrajectoryResult>();
    const std::size_t T = tr.NumFrames();
    tr.Finalize(tp, traj);
    EXPECT_EQ(tr.NumFrames(), T);
}


TEST(DihedralBinTransition, H5RoundTrip) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(99999);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::DihedralBinTransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();

    const std::string h5_path = (fs::temp_directory_path() /
        ("dihedral_bin_trans_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    ASSERT_TRUE(reopen.exist("/trajectory/dihedral_bin_transition"));
    auto grp = reopen.getGroup("/trajectory/dihedral_bin_transition");

    // Per-residue stats datasets.
    for (const std::string& name : {"backbone_transition_count",
                                     "backbone_dominant_region",
                                     "n_frames_observed"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 1u);
        EXPECT_EQ(dims[0], R);
    }
    {
        ASSERT_TRUE(grp.exist("backbone_bin_occupancy"));
        const auto dims =
            grp.getDataSet("backbone_bin_occupancy").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 6u);
    }
    for (const std::string& name : {"chi_transition_count",
                                     "chi_dominant_rotamer",
                                     "chi_n_frames_observed"}) {
        ASSERT_TRUE(grp.exist(name)) << name;
        const auto dims = grp.getDataSet(name).getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 2u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 4u);
    }
    {
        ASSERT_TRUE(grp.exist("chi_rotamer_occupancy"));
        const auto dims =
            grp.getDataSet("chi_rotamer_occupancy").getSpace().getDimensions();
        ASSERT_EQ(dims.size(), 3u);
        EXPECT_EQ(dims[0], R);
        EXPECT_EQ(dims[1], 4u);
        EXPECT_EQ(dims[2], 3u);
    }
    EXPECT_TRUE(grp.exist("frame_indices"));
    EXPECT_TRUE(grp.exist("frame_times"));
    EXPECT_TRUE(grp.exist("source_attached_per_frame"));

    // Convention attrs
    std::string legend, gate, source_policy;
    grp.getAttribute("backbone_bin_legend").read(legend);
    grp.getAttribute("transition_gate").read(gate);
    grp.getAttribute("source_attached_policy").read(source_policy);
    EXPECT_NE(legend.find("alphaR"), std::string::npos);
    EXPECT_NE(gate.find("Both prev and curr"), std::string::npos);
    EXPECT_NE(source_policy.find("always_attached"), std::string::npos);

    fs::remove(h5_path);
}


TEST(DihedralBinTransition, Integration1P9J) {
    nmr::test::TestEnvironment::LoadCalculatorConfig();
    nmr::test::TestEnvironment::Load();
    auto fix = nmr::test::TestEnvironment::FleetAmberTrajectory(kFixtureProtein);
    if (!FixtureAvailable(fix)) GTEST_SKIP() << "fixture not on disk";

    auto config = BuildConfig(300);
    nmr::TrajectoryProtein tp;
    ASSERT_TRUE(tp.BuildFromTrajectory(ProductionDirFor(fix.tpr_path)))
        << tp.Error();
    nmr::Trajectory traj(TrrPathFor(fix.tpr_path), fix.tpr_path, fix.edr_path);
    nmr::Session session;
    ASSERT_EQ(traj.Run(tp, config, session), nmr::kOk);

    const auto& tr = tp.Result<nmr::DihedralBinTransitionTrajectoryResult>();
    const std::size_t R = tp.ProteinRef().ResidueCount();
    const std::size_t T = tr.NumFrames();
    EXPECT_GE(T, 2u);

    const std::string h5_path = (fs::temp_directory_path() /
        ("dihedral_bin_trans_int_" + std::to_string(::getpid()) + ".h5")).string();
    { HighFive::File file(h5_path, HighFive::File::Truncate);
      tr.WriteH5Group(tp, file); }
    HighFive::File reopen(h5_path, HighFive::File::ReadOnly);
    auto grp = reopen.getGroup("/trajectory/dihedral_bin_transition");

    std::vector<std::uint32_t> bb_trans_count, n_obs;
    grp.getDataSet("backbone_transition_count").read(bb_trans_count);
    grp.getDataSet("n_frames_observed").read(n_obs);
    ASSERT_EQ(bb_trans_count.size(), R);
    ASSERT_EQ(n_obs.size(), R);

    std::vector<std::vector<std::uint32_t>> bb_occ;
    grp.getDataSet("backbone_bin_occupancy").read(bb_occ);
    ASSERT_EQ(bb_occ.size(), R);

    // Invariants:
    // (1) Sum of bin occupancy across non-unassigned bins equals
    //     n_frames_observed for each residue.
    // (2) Sum of ALL bin occupancy (including unassigned bin 0) equals
    //     T — every frame accounted for (codex-review-2026-05-19 fix).
    // (3) transition_count <= n_frames_observed (transitions can occur
    //     at most once per consecutive observed-frame pair).
    // (4) For non-N-term, non-C-term residues on 1P9J (linear single-
    //     chain monotonic), n_frames_observed should equal T.
    std::size_t residues_fully_observed = 0;
    std::size_t total_transitions = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        std::uint32_t sum_non_unassigned = 0;
        for (std::size_t b = 1; b < 6; ++b) sum_non_unassigned += bb_occ[ri][b];
        EXPECT_EQ(sum_non_unassigned, n_obs[ri])
            << "non-unassigned occupancy sum != n_frames_observed at ri=" << ri;
        std::uint32_t sum_all = bb_occ[ri][0] + sum_non_unassigned;
        EXPECT_EQ(sum_all, T)
            << "total occupancy sum (including unassigned bin 0) != T at ri=" << ri;
        EXPECT_LE(bb_trans_count[ri], n_obs[ri])
            << "transition_count > n_frames_observed at ri=" << ri;
        if (n_obs[ri] == T) ++residues_fully_observed;
        total_transitions += bb_trans_count[ri];
    }
    EXPECT_GE(residues_fully_observed, R - 5u)
        << "Expected most residues to have phi+psi observed every frame";

    std::cout << "DihedralBinTransition: " << R << " residues, " << T
              << " frames; residues_fully_observed=" << residues_fully_observed
              << "; total_backbone_transitions=" << total_transitions << "\n";

    // Chi rotamer invariants: occupancy sum equals
    // chi_n_frames_observed per (residue, chi).
    std::vector<std::vector<std::vector<std::uint32_t>>> chi_occ;
    grp.getDataSet("chi_rotamer_occupancy").read(chi_occ);
    ASSERT_EQ(chi_occ.size(), R);
    std::vector<std::vector<std::uint32_t>> chi_n_obs;
    grp.getDataSet("chi_n_frames_observed").read(chi_n_obs);
    ASSERT_EQ(chi_n_obs.size(), R);

    std::size_t chi_observations = 0;
    for (std::size_t ri = 0; ri < R; ++ri) {
        for (std::size_t k = 0; k < 4; ++k) {
            const auto sum = chi_occ[ri][k][0] + chi_occ[ri][k][1]
                + chi_occ[ri][k][2];
            EXPECT_EQ(sum, chi_n_obs[ri][k])
                << "chi occupancy sum mismatch at ri=" << ri << " k=" << k;
            if (chi_n_obs[ri][k] > 0) ++chi_observations;
        }
    }
    EXPECT_GT(chi_observations, 0u);
    std::cout << "  chi observations (residue × chi pairs with >0 obs): "
              << chi_observations << "\n";

    // dominant_region == argmax of occupancy on residues with >0 obs.
    std::vector<std::uint8_t> dom_region;
    grp.getDataSet("backbone_dominant_region").read(dom_region);
    ASSERT_EQ(dom_region.size(), R);
    for (std::size_t ri = 0; ri < R; ++ri) {
        if (n_obs[ri] == 0) {
            EXPECT_EQ(dom_region[ri], 0u);  // unassigned
            continue;
        }
        // Verify argmax: dom_region[ri] is one of the bins with max count.
        std::uint32_t max_count = 0;
        for (std::size_t b = 0; b < 6; ++b)
            max_count = std::max(max_count, bb_occ[ri][b]);
        EXPECT_EQ(bb_occ[ri][dom_region[ri]], max_count)
            << "dominant_region not argmax at ri=" << ri;
    }

    fs::remove(h5_path);
}
