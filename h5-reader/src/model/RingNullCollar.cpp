#include "RingNullCollar.h"

#include "Conformation.h"
#include "QtConformationSnapshot.h"
#include "QtMopacCoreGroup.h"
#include "QtMopacCoulombGroup.h"
#include "QtMopacMcConnellGroup.h"
#include "QtProtein.h"

#include "../io/DftShieldingLoader.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <utility>

namespace h5reader::model {
namespace {

struct LoadedDftFrame {
    int frameIndex = -1;
    double timePs = 0.0;
    std::shared_ptr<const DftShieldingFrame> shielding;
};

std::vector<std::size_t> atomScanList(const QtProtein& protein,
                                      const RingNullCollarOptions& options) {
    std::vector<std::size_t> out;
    if (options.atom) {
        out.push_back(*options.atom);
        return out;
    }
    out.reserve(protein.atomCount());
    for (std::size_t i = 0; i < protein.atomCount(); ++i)
        out.push_back(i);
    return out;
}

std::vector<std::size_t> ringScanList(const QtProtein& protein,
                                      const RingNullCollarOptions& options) {
    std::vector<std::size_t> out;
    if (options.ring) {
        out.push_back(*options.ring);
        return out;
    }
    out.reserve(protein.ringCount());
    for (std::size_t i = 0; i < protein.ringCount(); ++i) {
        if (!options.includeSaturatedRings && !protein.ring(i).IsAromatic())
            continue;
        out.push_back(i);
    }
    return out;
}

std::optional<LoadedDftFrame> loadDftFrame(const h5reader::io::DftFrame& declared,
                                           const QtProtein& protein) {
    QString metaError;
    if (!declared.LoadMeta(&metaError))
        return std::nullopt;

    std::shared_ptr<const DftShieldingFrame> frame =
        h5reader::io::DftShieldingLoader::LoadAndValidate(
            declared.meta_json_abspath,
            &protein);
    if (!frame || !frame->valid || frame->atoms.size() < protein.atomCount())
        return std::nullopt;

    return LoadedDftFrame{
        declared.frame_index,
        declared.framePs(),
        std::move(frame),
    };
}

QString ringKindName(RingKind kind) {
    switch (kind) {
    case RingKind::Aromatic:
        return QStringLiteral("aromatic");
    case RingKind::Saturated:
        return QStringLiteral("saturated");
    }
    return QStringLiteral("unknown");
}

RingNullAtomIdentity makeAtomIdentity(const QtProtein& protein, std::size_t atom) {
    RingNullAtomIdentity out;
    out.atomLabelAmber = protein.atomLabel(atom, NamingConvention::Amber);
    out.atomLabelIupac = protein.atomLabel(atom, NamingConvention::Iupac);
    out.atomLabelBmrb = protein.atomLabel(atom, NamingConvention::Bmrb);

    const QtAtom& qtAtom = protein.atom(atom);
    out.residueIndex = qtAtom.residueIndex;
    if (qtAtom.residueIndex >= 0 &&
        static_cast<std::size_t>(qtAtom.residueIndex) < protein.residueCount()) {
        const std::size_t residue = static_cast<std::size_t>(qtAtom.residueIndex);
        out.residueNumber = protein.residue(residue).address.residueNumber;
        out.residueLabelAmber =
            protein.residueLabel(residue, NamingConvention::Amber, NamingSource::Verbatim);
        out.residueLabelIupac =
            protein.residueLabel(residue, NamingConvention::Iupac, NamingSource::Verbatim);
        out.residueLabelBmrb =
            protein.residueLabel(residue, NamingConvention::Bmrb, NamingSource::Verbatim);
    }
    return out;
}

RingNullRingIdentity makeRingIdentity(const QtProtein& protein, std::size_t ring) {
    RingNullRingIdentity out;
    const QtRing& qtRing = protein.ring(ring);
    out.typeName = QString::fromLatin1(qtRing.TypeName());
    out.typeIndex = qtRing.TypeIndexAsInt();
    out.kind = ringKindName(qtRing.ringKind);
    out.parentResidueIndex = qtRing.parentResidueIndex;
    out.parentResidueNumber = qtRing.parentResidueNumber;
    out.fusedPartnerRingId = qtRing.fusedPartnerRingId;
    out.atomIndices.reserve(qtRing.atomIndices.size());
    for (int32_t atomIndex : qtRing.atomIndices)
        out.atomIndices.push_back(static_cast<int>(atomIndex));

    if (qtRing.parentResidueIndex >= 0 &&
        static_cast<std::size_t>(qtRing.parentResidueIndex) < protein.residueCount()) {
        const std::size_t residue = static_cast<std::size_t>(qtRing.parentResidueIndex);
        out.parentResidueLabelAmber =
            protein.residueLabel(residue, NamingConvention::Amber, NamingSource::Verbatim);
        out.parentResidueLabelIupac =
            protein.residueLabel(residue, NamingConvention::Iupac, NamingSource::Verbatim);
        out.parentResidueLabelBmrb =
            protein.residueLabel(residue, NamingConvention::Bmrb, NamingSource::Verbatim);
    }
    return out;
}

double crossingFraction(const RingNullMeasurement& from,
                        const RingNullMeasurement& to,
                        double toleranceA) {
    const double tol = std::max(0.0, toleranceA);
    if (std::abs(from.nullMarginA) <= tol)
        return 0.0;
    if (std::abs(to.nullMarginA) <= tol)
        return 1.0;
    const double denom = std::abs(from.nullMarginA) + std::abs(to.nullMarginA);
    if (denom > 1e-12)
        return std::abs(from.nullMarginA) / denom;
    return 0.0;
}

RingNullEventFrame makeEventFrame(const RingNullOrcaSnapshot& from,
                                  const RingNullOrcaSnapshot& to,
                                  double toleranceA) {
    RingNullEventFrame out;
    out.fromSignedNullMarginA = from.null.nullMarginA;
    out.toSignedNullMarginA = to.null.nullMarginA;
    out.signedNullMarginStepA = to.null.nullMarginA - from.null.nullMarginA;
    out.zeroFraction = crossingFraction(from.null, to.null, toleranceA);
    out.zeroTimePs = from.timePs + out.zeroFraction * (to.timePs - from.timePs);
    out.zeroAtomPosition =
        from.null.atomPosition + out.zeroFraction * (to.null.atomPosition - from.null.atomPosition);
    return out;
}

RingNullMotion makeMotion(const RingNullOrcaSnapshot& from,
                          const RingNullOrcaSnapshot& to) {
    RingNullMotion out;
    out.worldVectorA = to.null.atomPosition - from.null.atomPosition;
    out.distanceA = out.worldVectorA.norm();
    out.timeStepPs = to.timePs - from.timePs;
    out.radialChangeA = to.null.radialA - from.null.radialA;
    out.absAxialChangeA = to.null.absAxialA - from.null.absAxialA;
    out.distanceChangeA = to.null.distanceA - from.null.distanceA;
    out.angleChangeDeg = to.null.angleDeg - from.null.angleDeg;
    out.angularFactorChange = to.null.angularFactor - from.null.angularFactor;
    return out;
}

RingNullMopacSignals makeMopacSignals(const QtConformationSnapshot& snapshot,
                                      std::size_t atom) {
    RingNullMopacSignals out;

    const QtMopacCoreGroup core(snapshot);
    if (const std::optional<double> charge = core.charge(atom)) {
        out.chargePresent = true;
        out.charge = *charge;
    }
    if (const std::optional<MopacScalars> scalars = core.scalars(atom)) {
        out.coreScalarsPresent = true;
        out.coreScalars = *scalars;
    }

    const QtMopacCoulombGroup coulomb(snapshot);
    if (const std::optional<Vec3> e = coulomb.E(atom)) {
        out.coulombEPresent = true;
        out.coulombE = *e;
    }
    if (const std::optional<CoulombScalars> scalars = coulomb.scalars(atom)) {
        out.coulombScalarsPresent = true;
        out.coulombScalars = *scalars;
    }
    if (const std::optional<SphericalTensor> shielding = coulomb.shielding(atom)) {
        out.coulombShieldingPresent = true;
        out.coulombShielding = *shielding;
    }
    if (const std::optional<QtEfg> efg = coulomb.efgBackbone(atom)) {
        out.coulombEfgBackbonePresent = true;
        out.coulombEfgBackbone = *efg;
    }
    if (const std::optional<QtEfg> efg = coulomb.efgAromatic(atom)) {
        out.coulombEfgAromaticPresent = true;
        out.coulombEfgAromatic = *efg;
    }

    const QtMopacMcConnellGroup mc(snapshot);
    if (const std::optional<SphericalTensor> shielding = mc.shielding(atom)) {
        out.mcShieldingPresent = true;
        out.mcShielding = *shielding;
    }
    if (const std::optional<PerBondCategoryT2> category = mc.categoryT2(atom)) {
        out.mcCategoryT2Present = true;
        out.mcCategoryT2 = *category;
    }
    if (const std::optional<McConnellScalars> scalars = mc.scalars(atom)) {
        out.mcScalarsPresent = true;
        out.mcScalars = *scalars;
    }

    out.present = out.chargePresent || out.coreScalarsPresent ||
                  out.coulombEPresent || out.coulombScalarsPresent ||
                  out.coulombShieldingPresent ||
                  out.coulombEfgBackbonePresent || out.coulombEfgAromaticPresent ||
                  out.mcShieldingPresent || out.mcCategoryT2Present ||
                  out.mcScalarsPresent;
    return out;
}

RingNullSignalStamp makeSignalStamp(const QtProtein& protein,
                                    Conformation& conformation,
                                    const LoadedDftFrame& dft,
                                    std::size_t dftOrdinal,
                                    double zeroDftOrdinal,
                                    double zeroTimePs,
                                    std::size_t atom,
                                    std::size_t ring,
                                    double toleranceA) {
    RingNullSignalStamp out;
    out.frameIndex = dft.frameIndex;
    out.timePs = dft.timePs;
    out.timeOffsetFromZeroPs = dft.timePs - zeroTimePs;
    out.dftOrdinalOffsetFromZero = static_cast<double>(dftOrdinal) - zeroDftOrdinal;

    const std::size_t frame = static_cast<std::size_t>(dft.frameIndex);
    out.null = MeasureRingNull(RingGeometryAt(conformation, ring, frame),
                               conformation.atomPosition(frame, atom),
                               toleranceA);

    if (dft.shielding && atom < dft.shielding->atoms.size()) {
        out.orcaPresent = true;
        out.orca = dft.shielding->atoms[atom];
        out.orcaTotalShape = ComputeCsaShape(out.orca.total_raw);
        out.orcaDiaShape = ComputeCsaShape(out.orca.dia_raw);
        out.orcaParaShape = ComputeCsaShape(out.orca.para_raw);
    }

    conformation.requestSnapshot(frame);
    const std::shared_ptr<const QtConformationSnapshot> snapshot = conformation.snapshot(frame);
    if (snapshot && snapshot->protein() == &protein) {
        out.snapshotPresent = true;
        out.mopac = makeMopacSignals(*snapshot, atom);
    }
    return out;
}

std::vector<RingNullSignalStamp> makeSignalStamps(
    const QtProtein& protein,
    Conformation& conformation,
    const std::vector<const h5reader::io::DftFrame*>& candidates,
    std::size_t fromOrdinal,
    const LoadedDftFrame& fromDft,
    std::size_t toOrdinal,
    const LoadedDftFrame& toDft,
    const RingNullEventFrame& eventFrame,
    std::size_t atom,
    std::size_t ring,
    int radius,
    double toleranceA) {
    std::vector<RingNullSignalStamp> out;
    if (radius < 0)
        return out;

    const std::size_t first =
        fromOrdinal > static_cast<std::size_t>(radius)
            ? fromOrdinal - static_cast<std::size_t>(radius)
            : 0;
    const std::size_t last = std::min(candidates.size() - 1,
                                      toOrdinal + static_cast<std::size_t>(radius));
    out.reserve(last - first + 1);

    const double zeroDftOrdinal = static_cast<double>(fromOrdinal) + eventFrame.zeroFraction;
    for (std::size_t ordinal = first; ordinal <= last; ++ordinal) {
        std::optional<LoadedDftFrame> loaded;
        if (ordinal == fromOrdinal) {
            loaded = fromDft;
        } else if (ordinal == toOrdinal) {
            loaded = toDft;
        } else {
            loaded = loadDftFrame(*candidates[ordinal], protein);
        }
        if (!loaded)
            continue;
        out.push_back(makeSignalStamp(protein, conformation, *loaded, ordinal,
                                      zeroDftOrdinal, eventFrame.zeroTimePs,
                                      atom, ring, toleranceA));
    }
    return out;
}

RingNullOrcaSnapshot makeSnapshot(int frameIndex,
                                  double timePs,
                                  const RingNullMeasurement& null,
                                  const DftAtomShielding& shielding) {
    RingNullOrcaSnapshot out;
    out.frameIndex = frameIndex;
    out.timePs = timePs;
    out.null = null;
    out.shielding = shielding;
    out.totalShape = ComputeCsaShape(shielding.total_raw);
    out.diaShape = ComputeCsaShape(shielding.dia_raw);
    out.paraShape = ComputeCsaShape(shielding.para_raw);
    return out;
}

}  // namespace

RingNullCollar::RingNullCollar(RingNullCollarOptions options)
    : options_(std::move(options)) {}

bool RingNullCollar::collect(const QtProtein& protein,
                             Conformation& conformation,
                             const std::vector<h5reader::io::DftFrame>& dftFrames,
                             QString* error) {
    summary_ = {};
    entries_.clear();
    summary_.dftFramesDeclared = static_cast<int>(dftFrames.size());

    if (protein.atomCount() == 0 || protein.ringCount() == 0) {
        if (error)
            *error = QStringLiteral("protein has no atoms or rings");
        return false;
    }
    if (conformation.frameCount() == 0) {
        if (error)
            *error = QStringLiteral("conformation has no frames");
        return false;
    }
    if (options_.atom && *options_.atom >= protein.atomCount()) {
        if (error)
            *error = QStringLiteral("atom out of range");
        return false;
    }
    if (options_.ring && *options_.ring >= protein.ringCount()) {
        if (error)
            *error = QStringLiteral("ring out of range");
        return false;
    }
    if (!std::isfinite(options_.surfaceToleranceA) || options_.surfaceToleranceA < 0.0) {
        if (error)
            *error = QStringLiteral("surfaceToleranceA must be finite and >= 0");
        return false;
    }
    if (options_.stampRadiusDft < 0) {
        if (error)
            *error = QStringLiteral("stampRadiusDft must be >= 0");
        return false;
    }

    const int frameCount = static_cast<int>(conformation.frameCount());
    const int start = options_.startFrame.value_or(0);
    const int end = options_.endFrame.value_or(frameCount - 1);
    if (start < 0 || end < 0 || start >= frameCount || end >= frameCount || end < start) {
        if (error)
            *error = QStringLiteral("frame range out of range");
        return false;
    }

    const std::vector<std::size_t> atoms = atomScanList(protein, options_);
    const std::vector<std::size_t> rings = ringScanList(protein, options_);
    summary_.atomsScanned = static_cast<int>(atoms.size());
    summary_.ringsScanned = static_cast<int>(rings.size());
    if (atoms.empty() || rings.empty()) {
        summary_.complete = true;
        return true;
    }

    std::vector<const h5reader::io::DftFrame*> candidates;
    candidates.reserve(dftFrames.size());
    for (const h5reader::io::DftFrame& declared : dftFrames) {
        if (declared.frame_index < start || declared.frame_index > end)
            continue;
        if (declared.frame_index < 0 || declared.frame_index >= frameCount) {
            ++summary_.dftFramesSkipped;
            continue;
        }
        candidates.push_back(&declared);
    }
    std::sort(candidates.begin(), candidates.end(),
              [](const h5reader::io::DftFrame* a, const h5reader::io::DftFrame* b) {
        return a->frame_index < b->frame_index;
    });

    std::optional<int> lastCandidateFrame;
    std::optional<LoadedDftFrame> prevDft;
    std::optional<std::size_t> prevOrdinal;
    for (std::size_t ordinal = 0; ordinal < candidates.size(); ++ordinal) {
        const h5reader::io::DftFrame* declared = candidates[ordinal];
        if (lastCandidateFrame && *lastCandidateFrame == declared->frame_index)
            continue;
        lastCandidateFrame = declared->frame_index;

        std::optional<LoadedDftFrame> currentDft = loadDftFrame(*declared, protein);
        if (!currentDft) {
            ++summary_.dftFramesSkipped;
            continue;
        }
        ++summary_.dftFramesLoaded;
        if (!prevDft) {
            prevDft = std::move(currentDft);
            prevOrdinal = ordinal;
            continue;
        }

        const LoadedDftFrame& fromDft = *prevDft;
        const LoadedDftFrame& toDft = *currentDft;
        const std::size_t fromOrdinal = prevOrdinal.value_or(ordinal - 1);
        ++summary_.dftPairsScanned;

        const std::size_t fromFrame = static_cast<std::size_t>(fromDft.frameIndex);
        const std::size_t toFrame = static_cast<std::size_t>(toDft.frameIndex);

        for (std::size_t ring : rings) {
            const RingGeometry fromRing = RingGeometryAt(conformation, ring, fromFrame);
            const RingGeometry toRing = RingGeometryAt(conformation, ring, toFrame);

            for (std::size_t atom : atoms) {
                const RingNullMeasurement fromNull =
                    MeasureRingNull(fromRing, conformation.atomPosition(fromFrame, atom),
                                    options_.surfaceToleranceA);
                const RingNullMeasurement toNull =
                    MeasureRingNull(toRing, conformation.atomPosition(toFrame, atom),
                                    options_.surfaceToleranceA);
                if (!RingNullCrosses(fromNull, toNull, options_.surfaceToleranceA))
                    continue;

                RingNullCollarEntry entry;
                entry.atom = atom;
                entry.ring = ring;
                entry.atomIdentity = makeAtomIdentity(protein, atom);
                entry.ringIdentity = makeRingIdentity(protein, ring);
                entry.from = makeSnapshot(fromDft.frameIndex, fromDft.timePs, fromNull,
                                          fromDft.shielding->atoms[atom]);
                entry.to = makeSnapshot(toDft.frameIndex, toDft.timePs, toNull,
                                        toDft.shielding->atoms[atom]);
                entry.eventFrame =
                    makeEventFrame(entry.from, entry.to, options_.surfaceToleranceA);
                entry.motion = makeMotion(entry.from, entry.to);
                if (options_.includeSignalStamps) {
                    entry.signalStamps = makeSignalStamps(
                        protein, conformation, candidates, fromOrdinal, fromDft,
                        ordinal, toDft, entry.eventFrame, atom, ring,
                        options_.stampRadiusDft, options_.surfaceToleranceA);
                    summary_.signalStampCount += static_cast<int>(entry.signalStamps.size());
                }
                entries_.push_back(std::move(entry));
            }
        }

        prevDft = std::move(currentDft);
        prevOrdinal = ordinal;
    }

    summary_.entryCount = static_cast<int>(entries_.size());
    summary_.complete = true;
    return true;
}

}  // namespace h5reader::model
