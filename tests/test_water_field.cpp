//
// test_water_field: first test of explicit solvent calculators.
//
// Reads a full-system XTC (protein + water + ions), splits via
// FullSystemReader, creates a conformation from protein positions,
// runs WaterFieldResult and HydrationShellResult.
//

#include "FullSystemReader.h"
#include "BuildResult.h"
#include "WaterFieldResult.h"
#include "HydrationShellResult.h"
#include "GromacsEnergyResult.h"
#include "OperationRunner.h"
#include "Trajectory.h"
#include "GromacsEnsembleLoader.h"
#include "PdbFileReader.h"
#include "ChargeSource.h"
#include "RuntimeEnvironment.h"
#include "OperationLog.h"
#include "xtc_reader.h"

#include <gtest/gtest.h>
#include <iostream>

#ifndef NMR_TEST_DATA_DIR
#error "NMR_TEST_DATA_DIR must be defined"
#endif

static const std::string FULLSYS_DIR =
    std::string(NMR_TEST_DATA_DIR) + "/fleet_test_fullsys/1ZR7_6721";


TEST(WaterField, ReadTopologyAndExtractFrame) {
    nmr::OperationLog::SetChannelMask(0xFFFFFFFF);

    // 1. Read full-system topology from TPR
    nmr::FullSystemReader reader;
    std::string tpr = FULLSYS_DIR + "/run_params/prod_fullsys.tpr";
    ASSERT_TRUE(reader.ReadTopology(tpr)) << reader.error();

    auto& topo = reader.Topology();
    std::cout << "\n=== Topology ==="
              << "\n  total atoms:    " << topo.total_atoms
              << "\n  protein atoms:  " << topo.protein_count
              << "\n  water molecules:" << topo.water_count
              << "\n  ions:           " << topo.ion_count
              << "\n  water O charge: " << topo.water_O_charge
              << "\n  water H charge: " << topo.water_H_charge
              << std::endl;

    ASSERT_GT(topo.protein_count, 0u);
    ASSERT_GT(topo.water_count, 0u);
    EXPECT_EQ(topo.total_atoms, 11079u);

    // 2. Read one XTC frame (full system)
    std::string xtc = FULLSYS_DIR + "/walker_0/md.xtc";
    auto frames = read_all_xtc_frames(xtc);
    ASSERT_GT(frames.size(), 0u);
    EXPECT_EQ(static_cast<size_t>(frames[0].natoms), topo.total_atoms);

    std::cout << "  XTC frames: " << frames.size()
              << "\n  first frame t=" << frames[0].time << "ps"
              << std::endl;

    // 3. Extract protein positions + solvent from frame 0
    std::vector<nmr::Vec3> protein_pos;
    nmr::SolventEnvironment solvent;
    ASSERT_TRUE(reader.ExtractFrame(frames[0].x, protein_pos, solvent));

    EXPECT_EQ(protein_pos.size(), topo.protein_count);
    EXPECT_EQ(solvent.WaterCount(), topo.water_count);

    std::cout << "\n=== Frame 0 ==="
              << "\n  protein positions: " << protein_pos.size()
              << "\n  water molecules:   " << solvent.WaterCount()
              << "\n  ions:              " << solvent.IonCount()
              << std::endl;

    // 4. Build the protein from the same TPR slice, then seed it with
    // the extracted frame positions.
    auto build = reader.BuildProtein("1ZR7_6721");
    ASSERT_TRUE(build.Ok()) << build.error;

    build.protein->FinalizeConstruction(protein_pos);
    auto& conf = build.protein->AddMDFrame(
        std::move(protein_pos), 0, frames[0].time, 1.0, 0.0, 0.0);

    // EDR energy is loaded by Trajectory and handed to OperationRunner
    // per frame in the production trajectory path. Mirror that here.
    nmr::Trajectory traj(xtc, tpr, FULLSYS_DIR + "/walker_0/md.edr");
    ASSERT_TRUE(traj.HasEdr());
    const nmr::GromacsEnergy* frame_energy = traj.EnergyAtTime(frames[0].time);
    ASSERT_NE(frame_energy, nullptr);

    // 5. Run standard calculators (no MOPAC, no Coulomb, with APBS)
    nmr::RunOptions opts;
    opts.skip_mopac = true;
    opts.skip_coulomb = true;
    opts.charge_source = build.charges.get();
    opts.solvent = &solvent;
    opts.frame_energy = frame_energy;
    if (reader.HasBondedParams())
        opts.bonded_params = &reader.BondedParams();

    nmr::RunResult rr = nmr::OperationRunner::Run(conf, opts);
    ASSERT_TRUE(rr.Ok()) << rr.error;

    std::cout << "\n=== Results ==="
              << "\n  attached (" << rr.attached.size() << "):";
    for (const auto& name : rr.attached)
        std::cout << "\n    " << name;
    std::cout << std::endl;

    // 6. Verify water calculators ran
    bool has_water = false;
    bool has_hydration = false;
    bool has_energy = false;
    for (const auto& name : rr.attached) {
        if (name == "WaterFieldResult") has_water = true;
        if (name == "HydrationShellResult") has_hydration = true;
        if (name == "GromacsEnergyResult") has_energy = true;
    }
    EXPECT_TRUE(has_water) << "WaterFieldResult did not attach";
    EXPECT_TRUE(has_hydration) << "HydrationShellResult did not attach";
    EXPECT_TRUE(has_energy) << "GromacsEnergyResult did not attach";

    // 7. Spot-check: some atoms should have non-zero water fields
    int atoms_with_water = 0;
    for (size_t i = 0; i < conf.AtomCount(); ++i) {
        if (conf.AtomAt(i).water_n_first > 0)
            ++atoms_with_water;
    }
    std::cout << "\n  atoms with first-shell water: "
              << atoms_with_water << " / " << conf.AtomCount()
              << std::endl;
    EXPECT_GT(atoms_with_water, 0) << "No atoms have first-shell water";

    // 8. Print a few atoms' water data
    std::cout << "\n=== Sample atoms ===" << std::endl;
    for (size_t i = 0; i < std::min<size_t>(10, conf.AtomCount()); ++i) {
        const auto& a = conf.AtomAt(i);
        std::cout << "  atom " << i
                  << "  n_first=" << a.water_n_first
                  << "  n_second=" << a.water_n_second
                  << "  |E_water|=" << a.water_efield.norm()
                  << "  half_shell=" << a.half_shell_asymmetry
                  << "  dipole_cos=" << a.mean_water_dipole_cos
                  << std::endl;
    }
}
