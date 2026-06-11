#include "ResidentIndexes.h"

#include "AnalysisBody.h"
#include "Catalog.h"
#include "RunData.h"

#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace h5reader::rediscover {
namespace {

QString chemicalCategoryForAtom(const model::QtAtom& atom) {
    if (atom.IsBackbone()) return QStringLiteral("backbone");
    if (atom.IsAromaticRingHydrogen() || atom.aromatic) return QStringLiteral("aromatic");
    if (atom.IsPolarH()) return QStringLiteral("polar_h");
    if (atom.residueIndex >= 0) return QStringLiteral("sidechain");
    return QStringLiteral("unknown");
}

std::size_t indexFrameCount(const RunData& run) {
    if (run.frameMap.frameCount() > 0) return run.frameMap.frameCount();
    if (run.conformation) return run.conformation->frameCount();
    return 1;
}

std::optional<std::array<double, 8>> readSs8(const Body& body,
                                             std::size_t atom,
                                             std::size_t frame) {
    std::array<double, 8> row = {};
    for (int c = 0; c < 8; ++c) {
        QString reason;
        const std::optional<double> v =
            body.catalog.value(body, io::FieldKind::DSSPSs8, atom, frame, c, &reason);
        if (!v) return std::nullopt;
        row[static_cast<std::size_t>(c)] = *v;
    }
    return row;
}

std::optional<double> readFieldValue(const Body& body,
                                     io::FieldKind kind,
                                     std::size_t nativeRow,
                                     std::size_t frame,
                                     int component = 0) {
    QString reason;
    return body.catalog.value(body, kind, nativeRow, frame, component, &reason);
}

double dihedralRadians(const model::Vec3& p0,
                       const model::Vec3& p1,
                       const model::Vec3& p2,
                       const model::Vec3& p3) {
    model::Vec3 b0 = -(p1 - p0);
    model::Vec3 b1 = p2 - p1;
    model::Vec3 b2 = p3 - p2;
    const double b1Norm = b1.norm();
    if (!(b1Norm > 1.0e-12)) return std::numeric_limits<double>::quiet_NaN();
    b1 /= b1Norm;

    const model::Vec3 v = b0 - b0.dot(b1) * b1;
    const model::Vec3 w = b2 - b2.dot(b1) * b1;
    const double vNorm = v.norm();
    const double wNorm = w.norm();
    if (!(vNorm > 1.0e-12) || !(wNorm > 1.0e-12))
        return std::numeric_limits<double>::quiet_NaN();
    const double x = v.dot(w);
    const double y = b1.cross(v).dot(w);
    return WrapRadians(std::atan2(y, x));
}

BackboneAngleSignPolicy resolveBackboneSignPolicy(const Body& body) {
    if (!body.run.protein || !body.run.conformation) return BackboneAngleSignPolicy::Unresolved;
    const model::QtProtein& protein = *body.run.protein;
    const model::Conformation& conformation = *body.run.conformation;
    const std::size_t frames = std::min<std::size_t>(indexFrameCount(body.run), 8);

    double straightError = 0.0;
    double negatedError = 0.0;
    std::size_t samples = 0;
    for (std::size_t frame = 0; frame < frames; ++frame) {
        for (std::size_t ri = 0; ri < protein.residueCount(); ++ri) {
            const model::QtResidue& residue = protein.residue(ri);
            if (!residue.HasCA()) continue;
            const std::size_t atomRow = static_cast<std::size_t>(residue.CA);

            if (residue.prevResidueIndex >= 0
                && static_cast<std::size_t>(residue.prevResidueIndex) < protein.residueCount()
                && residue.HasN() && residue.HasC()) {
                const model::QtResidue& prev =
                    protein.residue(static_cast<std::size_t>(residue.prevResidueIndex));
                if (prev.HasC()) {
                    const std::optional<double> npy =
                        readFieldValue(body, io::FieldKind::DSSPBackbone, atomRow, frame, 0);
                    if (npy && std::abs(*npy) > 1.0e-12) {
                        const double computed = dihedralRadians(
                            conformation.atomPosition(frame, static_cast<std::size_t>(prev.C)),
                            conformation.atomPosition(frame, static_cast<std::size_t>(residue.N)),
                            conformation.atomPosition(frame, static_cast<std::size_t>(residue.CA)),
                            conformation.atomPosition(frame, static_cast<std::size_t>(residue.C)));
                        if (std::isfinite(computed)) {
                            straightError += AngularDistance(computed, *npy);
                            negatedError += AngularDistance(computed, -*npy);
                            ++samples;
                        }
                    }
                }
            }

            if (residue.nextResidueIndex >= 0
                && static_cast<std::size_t>(residue.nextResidueIndex) < protein.residueCount()
                && residue.HasN() && residue.HasCA() && residue.HasC()) {
                const model::QtResidue& next =
                    protein.residue(static_cast<std::size_t>(residue.nextResidueIndex));
                if (next.HasN()) {
                    const std::optional<double> npy =
                        readFieldValue(body, io::FieldKind::DSSPBackbone, atomRow, frame, 1);
                    if (npy && std::abs(*npy) > 1.0e-12) {
                        const double computed = dihedralRadians(
                            conformation.atomPosition(frame, static_cast<std::size_t>(residue.N)),
                            conformation.atomPosition(frame, static_cast<std::size_t>(residue.CA)),
                            conformation.atomPosition(frame, static_cast<std::size_t>(residue.C)),
                            conformation.atomPosition(frame, static_cast<std::size_t>(next.N)));
                        if (std::isfinite(computed)) {
                            straightError += AngularDistance(computed, *npy);
                            negatedError += AngularDistance(computed, -*npy);
                            ++samples;
                        }
                    }
                }
            }
        }
    }
    return ChooseBackboneSignPolicy(straightError, negatedError, samples);
}

DihedralState dihedralState(double radians) {
    if (!std::isfinite(radians)) return {};
    return DihedralState{WrapRadians(radians), FixedDihedralBin(radians), true};
}

void buildNameIndexes(const RunData& run, ResidentIndexes* idx) {
    if (!run.protein) return;
    const model::QtProtein& protein = *run.protein;
    idx->iupacNames.reset(protein.atomCount());
    idx->chemicalCategories.reset(protein.atomCount());
    for (std::size_t atom = 0; atom < protein.atomCount(); ++atom) {
        idx->iupacNames.add(atom, protein.atomLabel(atom, model::NamingConvention::Iupac));
        idx->chemicalCategories.add(atom, chemicalCategoryForAtom(protein.atom(atom)));
    }
}

void buildSecondaryStructureIndex(const Body& body, ResidentIndexes* idx) {
    if (!body.run.protein) return;
    const std::size_t atoms = body.run.protein->atomCount();
    const std::size_t frames = indexFrameCount(body.run);
    idx->secondaryStructure.reset(atoms, frames);
    for (std::size_t frame = 0; frame < frames; ++frame) {
        for (std::size_t atom = 0; atom < atoms; ++atom) {
            const std::optional<std::array<double, 8>> row = readSs8(body, atom, frame);
            idx->secondaryStructure.set(atom, frame, row ? DecodeSs8(*row) : SecondaryStructureState{});
        }
    }
}

void buildDihedralIndex(const Body& body, ResidentIndexes* idx) {
    if (!body.run.protein) return;
    const std::size_t atoms = body.run.protein->atomCount();
    const std::size_t frames = indexFrameCount(body.run);
    idx->dihedrals.reset(atoms, frames);
    const BackboneAngleSignPolicy signPolicy = resolveBackboneSignPolicy(body);
    idx->dihedrals.setBackboneSignPolicy(signPolicy);

    for (std::size_t frame = 0; frame < frames; ++frame) {
        for (std::size_t atom = 0; atom < atoms; ++atom) {
            const std::optional<double> phiRaw =
                readFieldValue(body, io::FieldKind::DSSPBackbone, atom, frame, 0);
            const std::optional<double> psiRaw =
                readFieldValue(body, io::FieldKind::DSSPBackbone, atom, frame, 1);
            if (phiRaw) {
                idx->dihedrals.set(DihedralKind::Phi, atom, frame,
                                   dihedralState(ApplyBackboneSignPolicy(*phiRaw, signPolicy)));
            }
            if (psiRaw) {
                idx->dihedrals.set(DihedralKind::Psi, atom, frame,
                                   dihedralState(ApplyBackboneSignPolicy(*psiRaw, signPolicy)));
            }

            const int32_t residueIndex = body.run.protein->atom(atom).residueIndex;
            if (residueIndex >= 0) {
                const std::optional<double> omega =
                    readFieldValue(body, io::FieldKind::OmegaActual,
                                   static_cast<std::size_t>(residueIndex), frame);
                if (omega)
                    idx->dihedrals.set(DihedralKind::Omega, atom, frame, dihedralState(*omega));
            }

            for (int k = 0; k < 4; ++k) {
                const int base = k * 3;
                const std::optional<double> c =
                    readFieldValue(body, io::FieldKind::DSSPChi, atom, frame, base);
                const std::optional<double> s =
                    readFieldValue(body, io::FieldKind::DSSPChi, atom, frame, base + 1);
                const std::optional<double> exists =
                    readFieldValue(body, io::FieldKind::DSSPChi, atom, frame, base + 2);
                if (c && s && exists && *exists != 0.0) {
                    idx->dihedrals.set(static_cast<DihedralKind>(static_cast<int>(DihedralKind::Chi1) + k),
                                       atom, frame, dihedralState(std::atan2(*s, *c)));
                }
            }
        }
    }
}

}  // namespace

ResidentIndexes BuildResidentIndexes(const RunData& run) {
    ResidentIndexes idx;
    if (run.protein) idx.typedAtoms = TypedAtomIndex(*run.protein);
    buildNameIndexes(run, &idx);
    Catalog catalog(run);
    Body body{run, idx, catalog};
    buildSecondaryStructureIndex(body, &idx);
    buildDihedralIndex(body, &idx);
    idx.spatial = SpatialIndexSet(run);
    idx.ringGeometry = RingGeometryCache(run);
    idx.temporal = TemporalIndex(run.conformation ? run.conformation->frameCount() : 0);
    return idx;
}

}  // namespace h5reader::rediscover
