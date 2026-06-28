#pragma once

#include "ConformationGeometry.h"
#include "CsaShape.h"
#include "DftShielding.h"
#include "QtResultBlocks.h"

#include "../io/CalcsetManifest.h"

#include <QString>

#include <cstddef>
#include <optional>
#include <vector>

namespace h5reader::model {

class Conformation;
class QtProtein;

struct RingNullCollarOptions {
    std::optional<std::size_t> atom;
    std::optional<std::size_t> ring;
    std::optional<int> startFrame;
    std::optional<int> endFrame;
    double surfaceToleranceA = 1e-9;
    bool includeSaturatedRings = false;
    bool includeSignalStamps = true;
    int stampRadiusDft = 2;
};

struct RingNullOrcaSnapshot {
    int frameIndex = -1;
    double timePs = 0.0;
    RingNullMeasurement null;
    DftAtomShielding shielding;
    CsaShape totalShape;
    CsaShape diaShape;
    CsaShape paraShape;
};

struct RingNullAtomIdentity {
    QString atomLabelAmber;
    QString atomLabelIupac;
    QString atomLabelBmrb;
    int residueIndex = -1;
    int residueNumber = 0;
    QString residueLabelAmber;
    QString residueLabelIupac;
    QString residueLabelBmrb;
};

struct RingNullRingIdentity {
    QString typeName;
    int typeIndex = -1;
    QString kind;
    int parentResidueIndex = -1;
    int parentResidueNumber = 0;
    QString parentResidueLabelAmber;
    QString parentResidueLabelIupac;
    QString parentResidueLabelBmrb;
    int fusedPartnerRingId = -1;
    std::vector<int> atomIndices;
};

struct RingNullEventFrame {
    double fromSignedNullMarginA = 0.0;
    double toSignedNullMarginA = 0.0;
    double signedNullMarginStepA = 0.0;
    double zeroFraction = 0.0;
    double zeroTimePs = 0.0;
    Vec3 zeroAtomPosition = Vec3::Zero();
};

struct RingNullMotion {
    Vec3 worldVectorA = Vec3::Zero();
    double distanceA = 0.0;
    double timeStepPs = 0.0;
    double radialChangeA = 0.0;
    double absAxialChangeA = 0.0;
    double distanceChangeA = 0.0;
    double angleChangeDeg = 0.0;
    double angularFactorChange = 0.0;
};

struct RingNullMopacSignals {
    bool present = false;

    bool chargePresent = false;
    double charge = 0.0;
    bool coreScalarsPresent = false;
    MopacScalars coreScalars;

    bool coulombEPresent = false;
    Vec3 coulombE = Vec3::Zero();
    bool coulombScalarsPresent = false;
    CoulombScalars coulombScalars;
    bool coulombShieldingPresent = false;
    SphericalTensor coulombShielding;
    bool coulombEfgBackbonePresent = false;
    QtEfg coulombEfgBackbone;
    bool coulombEfgAromaticPresent = false;
    QtEfg coulombEfgAromatic;

    bool mcShieldingPresent = false;
    SphericalTensor mcShielding;
    bool mcCategoryT2Present = false;
    PerBondCategoryT2 mcCategoryT2;
    bool mcScalarsPresent = false;
    McConnellScalars mcScalars;
};

struct RingNullSignalStamp {
    int frameIndex = -1;
    double timePs = 0.0;
    double timeOffsetFromZeroPs = 0.0;
    double dftOrdinalOffsetFromZero = 0.0;
    RingNullMeasurement null;

    bool orcaPresent = false;
    DftAtomShielding orca;
    CsaShape orcaTotalShape;
    CsaShape orcaDiaShape;
    CsaShape orcaParaShape;

    bool snapshotPresent = false;
    RingNullMopacSignals mopac;
};

struct RingNullCollarEntry {
    std::size_t atom = 0;
    std::size_t ring = 0;
    RingNullAtomIdentity atomIdentity;
    RingNullRingIdentity ringIdentity;
    RingNullEventFrame eventFrame;
    RingNullMotion motion;
    RingNullOrcaSnapshot from;
    RingNullOrcaSnapshot to;
    std::vector<RingNullSignalStamp> signalStamps;
};

struct RingNullCollarSummary {
    bool complete = false;
    int dftFramesDeclared = 0;
    int dftFramesLoaded = 0;
    int dftFramesSkipped = 0;
    int dftPairsScanned = 0;
    int atomsScanned = 0;
    int ringsScanned = 0;
    int entryCount = 0;
    int signalStampCount = 0;
};

// A run-level operational object: hold a collar definition, scan the DFT frame
// sequence once, and retain the whole per-atom ORCA shielding records for the
// DFT pairs that cross the ring-null surface.
class RingNullCollar {
public:
    explicit RingNullCollar(RingNullCollarOptions options = {});

    bool collect(const QtProtein& protein,
                 Conformation& conformation,
                 const std::vector<h5reader::io::DftFrame>& dftFrames,
                 QString* error = nullptr);

    const RingNullCollarOptions& options() const { return options_; }
    const RingNullCollarSummary& summary() const { return summary_; }
    const std::vector<RingNullCollarEntry>& entries() const { return entries_; }

private:
    RingNullCollarOptions options_;
    RingNullCollarSummary summary_;
    std::vector<RingNullCollarEntry> entries_;
};

}  // namespace h5reader::model
