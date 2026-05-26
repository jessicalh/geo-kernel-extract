#include "ReduceProtonation.h"
#include "OperationLog.h"

// reducelib exposes mutable globals; set them before each invocation.
#include "reduce.h"
#include "CTab.h"
#include "AtomPositions.h"
#include "DotSph.h"

#include <sstream>
#include <filesystem>

namespace fs = std::filesystem;

namespace nmr {

// CMake supplies the default reduce_wwPDB_het_dict.txt path; REDUCE_HET_DICT
// can override it at runtime.
#ifndef HET_DICTIONARY
#error "HET_DICTIONARY must be defined by CMake (set REDUCE_SRC in CMakeLists.txt)"
#endif


// Set reduce globals to BUILD mode (equivalent to `reduce -BUILD`).
// This must be called before every invocation because reduce uses
// mutable global state.
static void SetReduceBuildMode() {
    BuildHisHydrogens      = true;   // add HIS sidechain NH
    SaveOHetcHydrogens     = true;   // keep OH/SH hydrogens
    AddOtherHydrogens      = true;   // add non-water hydrogens
    AddWaterHydrogens      = false;  // skip water (we strip waters)

    RotExistingOH          = true;   // rotate existing OH
    DemandFlipAllHNQs      = true;   // flip NQH groups
    StopBeforeOptimizing   = false;  // do full optimization
    OKtoAdjust             = true;   // allow adjustments

    // Existing ATOM hydrogens are stripped before reduce rebuilds them.
    RemoveATOMHydrogens    = true;
    RemoveOtherHydrogens   = false;

    Verbose                = false;
    KeepConnections        = true;
    StandardizeRHBondLengths = true;
    ProcessConnHydOnHets   = false;
    UseXplorNames          = false;
    UseOldNames            = false;
    BackBoneModel          = false;
    NeutralTermini         = false;
    DemandRotNH3           = false;
    DemandRotExisting      = false;
    DoOnlyAltA             = true;
    ShowCliqueTicks        = false;
    ShowOrientScore        = false;
    StringInput            = false;
    ShowCharges            = false;
    UseNuclearDistances    = false;
    RenameFlip             = false;
    GenerateFinalFlip      = false;

    NBondCutoff            = 4;
    ExhaustiveLimit        = 600;
    ProbeRadius            = 0.0;
    VdwDotDensity          = 16.0;
    OccupancyCutoff        = 0.01;
    WaterBcutoff           = 40.0;
    WaterOCCcutoff         = 0.66;
    PenaltyMagnitude       = 1.0;
    MinRegHBgap            = 0.6;
    MinChargedHBgap        = 0.8;
    BadBumpGapCut          = 0.4;
    NonMetalBumpBias       = 0.125;
    MetalBumpBias          = 0.865;
    GapWidth               = 0.3;
    MinNTermResNo          = 1;
    MaxAromRingDih          = 10;

    Tally = SummaryStats();
}


std::string ProtonateWithReduce(const std::string& pdb_content) {
    OperationLog::Scope scope("ProtonateWithReduce", "");

    if (pdb_content.empty()) {
        OperationLog::Error("ProtonateWithReduce", "empty PDB content");
        return {};
    }

    const char* env_dict = std::getenv("REDUCE_HET_DICT");
    std::string het_dict_path = env_dict ? env_dict : HET_DICTIONARY;
    if (!fs::exists(het_dict_path)) {
        OperationLog::Error("ProtonateWithReduce",
            "het dictionary not found: " + het_dict_path +
            " (set REDUCE_HET_DICT env var or REDUCE_SRC in CMake)");
        return {};
    }

    SetReduceBuildMode();

    CTab hetdatabase(het_dict_path);

    auto models = inputModels(pdb_content);
    if (models.empty()) {
        OperationLog::Error("ProtonateWithReduce", "no models parsed from PDB");
        return {};
    }

    auto& m = models[0];

    if (RemoveATOMHydrogens || RemoveOtherHydrogens) {
        dropHydrogens(m, RemoveATOMHydrogens, RemoveOtherHydrogens);
    }

    UseSEGIDasChain = checkSEGIDs(m);

    DotSphManager dotBucket(VdwDotDensity);
    AtomPositions xyz(2000, DoOnlyAltA, UseXplorNames, UseOldNames,
        BackBoneModel, NBondCutoff, MinRegHBgap, MinChargedHBgap,
        BadBumpGapCut, dotBucket, ProbeRadius,
        PenaltyMagnitude, OccupancyCutoff,
        Verbose, ShowOrientScore, ShowCliqueTicks, std::cerr);

    auto infoPtr = m.begin();
    scanAndGroupRecords(m, xyz, infoPtr);

    std::vector<std::string> adjNotes;
    Tally._num_adj = 0;
    reduceList(hetdatabase, m, xyz, adjNotes);

    int ret = optimize(xyz, adjNotes);
    if (OKtoAdjust && xyz.numChanges() > 0) {
        xyz.describeChanges(m, infoPtr, adjNotes);
    }

    std::string result = outputRecords_all_string(models);

    OperationLog::Info(LogCharges, "ProtonateWithReduce",
        "H found=" + std::to_string(Tally._H_found) +
        " added=" + std::to_string(Tally._H_added) +
        " removed=" + std::to_string(Tally._H_removed) +
        " adjusted=" + std::to_string(Tally._num_adj) +
        " optimize_ret=" + std::to_string(ret));

    return result;
}

}  // namespace nmr
