#include "NearbySignalModel.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../io/QtTrajectoryH5.h"
#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtBondVectorBuffers.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"
#include "../model/TrajectoryConformation.h"

#include <QBrush>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <utility>
#include <variant>

namespace h5reader::app {

namespace {
double distance(const model::Vec3& a, const model::Vec3& b) {
    const double dx = a.x() - b.x();
    const double dy = a.y() - b.y();
    const double dz = a.z() - b.z();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

QString kindLabel(NearbySignalModel::CandidateKind kind) {
    switch (kind) {
    case NearbySignalModel::CandidateKind::Atom:
        return QStringLiteral("atom");
    case NearbySignalModel::CandidateKind::Residue:
        return QStringLiteral("residue");
    case NearbySignalModel::CandidateKind::Bond:
        return QStringLiteral("bond");
    case NearbySignalModel::CandidateKind::BondVector:
        return QStringLiteral("bond vector");
    case NearbySignalModel::CandidateKind::Ring:
        return QStringLiteral("ring");
    case NearbySignalModel::CandidateKind::AromaticRing:
        return QStringLiteral("aromatic ring");
    case NearbySignalModel::CandidateKind::SaturatedRing:
        return QStringLiteral("saturated ring");
    case NearbySignalModel::CandidateKind::RingMembership:
        return QStringLiteral("ring member");
    }
    return QStringLiteral("candidate");
}

int kindSortRank(NearbySignalModel::CandidateKind kind) {
    switch (kind) {
    case NearbySignalModel::CandidateKind::Residue:
        return 0;
    // Bond vectors sit alongside atoms (both are residue
    // substructures). Same rank — distance + label tiebreak.
    case NearbySignalModel::CandidateKind::BondVector:
    case NearbySignalModel::CandidateKind::Atom:
        return 1;
    case NearbySignalModel::CandidateKind::Bond:
        return 2;
    case NearbySignalModel::CandidateKind::Ring:
    case NearbySignalModel::CandidateKind::AromaticRing:
    case NearbySignalModel::CandidateKind::SaturatedRing:
        return 3;
    case NearbySignalModel::CandidateKind::RingMembership:
        return 4;
    }
    return 9;
}

QString bondVectorKindSuffix(std::uint8_t kind) {
    switch (kind) {
    case 1: return QStringLiteral(" N-H");
    case 2: return QString::fromUtf8(" Cα-Hα");
    case 3: return QStringLiteral(" C=O");
    default: return QStringLiteral(" vector");
    }
}
}  // namespace

NearbySignalModel::NearbySignalModel(QObject* parent)
    : QAbstractTableModel(parent)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("NearbySignalModel"));
}

int NearbySignalModel::rowCount(const QModelIndex& parent) const {
    if (parent.isValid())
        return 0;
    return static_cast<int>(candidates_.size());
}

int NearbySignalModel::columnCount(const QModelIndex& parent) const {
    if (parent.isValid())
        return 0;
    return 3;
}

QVariant NearbySignalModel::data(const QModelIndex& index, int role) const {
    if (!index.isValid())
        return {};
    const int row = index.row();
    if (row < 0 || static_cast<std::size_t>(row) >= candidates_.size())
        return {};
    const Candidate& c = candidates_[static_cast<std::size_t>(row)];

    if (role == Qt::DisplayRole) {
        switch (index.column()) {
        case 0:
            return kindLabel(c.kind);
        case 1:
            return c.label;
        case 2:
            return QStringLiteral("%1 Å").arg(c.distanceAngstrom, 0, 'f', 2);
        default:
            return {};
        }
    }
    if (role == Qt::TextAlignmentRole && index.column() == 2)
        return QVariant::fromValue(Qt::AlignRight | Qt::AlignVCenter);
    if (role == Qt::ToolTipRole)
        return QStringLiteral("%1 within %2 Å at frame %3")
            .arg(c.label)
            .arg(radiusAngstrom_, 0, 'f', 1)
            .arg(frame_ + 1);
    if (role == KindRole)
        return kindLabel(c.kind);
    if (role == AtomIndexRole) {
        if (const auto* atom = std::get_if<model::AtomAnchor>(&c.anchor))
            return QVariant::fromValue(static_cast<qulonglong>(atom->atom));
        if (const auto* membership = std::get_if<model::RingMembershipAnchor>(&c.anchor)) {
            if (protein_ && membership->membership < protein_->ringMembershipCount()) {
                const auto& row = protein_->topology().ringMembershipAt(membership->membership);
                if (row.atomIndex >= 0)
                    return QVariant::fromValue(static_cast<qulonglong>(row.atomIndex));
            }
        }
        return {};
    }
    if (role == ResidueIndexRole) {
        if (const auto* residue = std::get_if<model::ResidueAnchor>(&c.anchor))
            return QVariant::fromValue(static_cast<qulonglong>(residue->residue));
        return c.residueContext ? QVariant::fromValue(static_cast<qulonglong>(*c.residueContext)) : QVariant();
    }
    if (role == BondIndexRole) {
        if (const auto* bond = std::get_if<model::BondAnchor>(&c.anchor))
            return QVariant::fromValue(static_cast<qulonglong>(bond->bond));
        return {};
    }
    if (role == RingIndexRole) {
        if (c.absoluteRing)
            return QVariant::fromValue(static_cast<qulonglong>(*c.absoluteRing));
        if (const auto* ring = std::get_if<model::RingAnchor>(&c.anchor))
            return QVariant::fromValue(static_cast<qulonglong>(ring->ring));
        return {};
    }
    if (role == RingMembershipIndexRole) {
        if (const auto* membership = std::get_if<model::RingMembershipAnchor>(&c.anchor))
            return QVariant::fromValue(static_cast<qulonglong>(membership->membership));
        return {};
    }
    if (role == DistanceRole)
        return c.distanceAngstrom;
    if (role == AnchorKindRole) {
        switch (c.kind) {
        case CandidateKind::Atom:
            return static_cast<int>(model::SignalAnchorKind::Atom);
        case CandidateKind::Residue:
            return static_cast<int>(model::SignalAnchorKind::Residue);
        case CandidateKind::Bond:
            return static_cast<int>(model::SignalAnchorKind::Bond);
        case CandidateKind::BondVector:
            return static_cast<int>(model::SignalAnchorKind::BondVector);
        case CandidateKind::Ring:
            return static_cast<int>(model::SignalAnchorKind::Ring);
        case CandidateKind::AromaticRing:
            return static_cast<int>(model::SignalAnchorKind::AromaticRing);
        case CandidateKind::SaturatedRing:
            return static_cast<int>(model::SignalAnchorKind::SaturatedRing);
        case CandidateKind::RingMembership:
            return static_cast<int>(model::SignalAnchorKind::RingMembership);
        }
    }
    return {};
}

QVariant NearbySignalModel::headerData(int section, Qt::Orientation orientation, int role) const {
    if (orientation != Qt::Horizontal || role != Qt::DisplayRole)
        return {};
    switch (section) {
    case 0: return QStringLiteral("Type");
    case 1: return QStringLiteral("Candidate");
    case 2: return QStringLiteral("Distance");
    default: return {};
    }
}

QHash<int, QByteArray> NearbySignalModel::roleNames() const {
    QHash<int, QByteArray> roles = QAbstractTableModel::roleNames();
    roles[KindRole] = "kind";
    roles[AtomIndexRole] = "atomIndex";
    roles[ResidueIndexRole] = "residueIndex";
    roles[BondIndexRole] = "bondIndex";
    roles[RingIndexRole] = "ringIndex";
    roles[RingMembershipIndexRole] = "ringMembershipIndex";
    roles[DistanceRole] = "distance";
    roles[AnchorKindRole] = "anchorKind";
    return roles;
}

void NearbySignalModel::setContext(const model::QtProtein* protein, model::Conformation* conformation) {
    ASSERT_THREAD(this);
    protein_ = protein;
    conformation_ = conformation;
    rebuild();
}

void NearbySignalModel::setRadiusAngstrom(double radius) {
    ASSERT_THREAD(this);
    const double clamped = std::clamp(radius, 1.0, 30.0);
    if (std::abs(clamped - radiusAngstrom_) < 1e-9)
        return;
    radiusAngstrom_ = clamped;
    rebuild();
}

void NearbySignalModel::setAnchor(std::size_t atom, int frame) {
    ASSERT_THREAD(this);
    anchorAtom_ = atom;
    frame_ = std::max(0, frame);
    rebuild();
}

void NearbySignalModel::clear() {
    ASSERT_THREAD(this);
    anchorAtom_.reset();
    beginResetModel();
    candidates_.clear();
    endResetModel();
}

const NearbySignalModel::Candidate* NearbySignalModel::candidateAt(const QModelIndex& index) const {
    if (!index.isValid())
        return nullptr;
    const int row = index.row();
    if (row < 0 || static_cast<std::size_t>(row) >= candidates_.size())
        return nullptr;
    return &candidates_[static_cast<std::size_t>(row)];
}

QString NearbySignalModel::residueLabel(std::size_t residueIdx) const {
    if (!protein_ || residueIdx >= protein_->residueCount())
        return QStringLiteral("residue %1").arg(residueIdx + 1);
    const auto& residue = protein_->residue(residueIdx);
    const QString chain = residue.address.chainId.isEmpty()
                              ? QString()
                              : QStringLiteral("%1:").arg(residue.address.chainId);
    const QString name = protein_->residueLabel(residueIdx,
                                                model::NamingConvention::Iupac,
                                                model::NamingSource::Verbatim);
    return QStringLiteral("%1%2%3").arg(chain, name).arg(residue.address.residueNumber);
}

QString NearbySignalModel::atomLabel(std::size_t atomIdx) const {
    if (!protein_ || atomIdx >= protein_->atomCount())
        return QStringLiteral("#%1").arg(atomIdx);
    const auto& atom = protein_->atom(atomIdx);
    const QString atomName = protein_->atomNames(atomIdx).iupac;
    if (atom.residueIndex >= 0
        && static_cast<std::size_t>(atom.residueIndex) < protein_->residueCount()) {
        return QStringLiteral("%1:%2").arg(
            residueLabel(static_cast<std::size_t>(atom.residueIndex)),
            atomName);
    }
    return QStringLiteral("#%1:%2").arg(atomIdx).arg(atomName);
}

QString NearbySignalModel::bondLabel(std::size_t bondIdx) const {
    if (!protein_ || bondIdx >= protein_->bondCount())
        return QStringLiteral("bond %1").arg(bondIdx + 1);
    const auto& bond = protein_->bond(bondIdx);
    const QString a = bond.atomIndexA >= 0
                          ? atomLabel(static_cast<std::size_t>(bond.atomIndexA))
                          : QStringLiteral("?");
    const QString b = bond.atomIndexB >= 0
                          ? atomLabel(static_cast<std::size_t>(bond.atomIndexB))
                          : QStringLiteral("?");
    return QStringLiteral("%1 - %2").arg(a, b);
}

QString NearbySignalModel::ringLabel(std::size_t ringIdx) const {
    if (!protein_ || ringIdx >= protein_->ringCount())
        return QStringLiteral("ring %1").arg(ringIdx + 1);
    const auto& ring = protein_->ring(ringIdx);
    QString owner = QStringLiteral("ring %1").arg(ringIdx + 1);
    if (ring.parentResidueIndex >= 0
        && static_cast<std::size_t>(ring.parentResidueIndex) < protein_->residueCount()) {
        owner = residueLabel(static_cast<std::size_t>(ring.parentResidueIndex));
    }
    return QStringLiteral("%1 %2").arg(owner, QString::fromLatin1(ring.TypeName()));
}

QString NearbySignalModel::ringMembershipLabel(std::size_t membershipIdx) const {
    if (!protein_ || membershipIdx >= protein_->ringMembershipCount())
        return QStringLiteral("ring member %1").arg(membershipIdx + 1);
    const auto& membership = protein_->topology().ringMembershipAt(membershipIdx);
    const QString atom = membership.atomIndex >= 0
                             ? atomLabel(static_cast<std::size_t>(membership.atomIndex))
                             : QStringLiteral("?");
    const QString ring = membership.ringId >= 0
                             ? ringLabel(static_cast<std::size_t>(membership.ringId))
                             : QStringLiteral("ring ?");
    return QStringLiteral("%1 in %2").arg(atom, ring);
}

void NearbySignalModel::rebuild() {
    ASSERT_THREAD(this);
    std::vector<Candidate> next;
    if (protein_ && conformation_ && anchorAtom_ && *anchorAtom_ < protein_->atomCount()
        && frame_ >= 0 && static_cast<std::size_t>(frame_) < conformation_->frameCount()) {
        const model::Vec3 focus = conformation_->atomPosition(static_cast<std::size_t>(frame_), *anchorAtom_);
        std::vector<double> atomDistances(protein_->atomCount(), std::numeric_limits<double>::infinity());
        std::vector<double> residueDistances(protein_->residueCount(), std::numeric_limits<double>::infinity());

        for (std::size_t atomIdx = 0; atomIdx < protein_->atomCount(); ++atomIdx) {
            const double d = distance(focus, conformation_->atomPosition(static_cast<std::size_t>(frame_), atomIdx));
            if (!std::isfinite(d))
                continue;
            atomDistances[atomIdx] = d;
            if (d > radiusAngstrom_)
                continue;
            const auto& atom = protein_->atom(atomIdx);
            if (atom.residueIndex >= 0
                && static_cast<std::size_t>(atom.residueIndex) < residueDistances.size()) {
                double& best = residueDistances[static_cast<std::size_t>(atom.residueIndex)];
                best = std::min(best, d);
            }
            next.push_back(Candidate{
                CandidateKind::Atom,
                model::AtomAnchor{atomIdx},
                atom.residueIndex >= 0 ? std::optional<std::size_t>(static_cast<std::size_t>(atom.residueIndex))
                                        : std::nullopt,
                std::nullopt,
                d,
                atomLabel(atomIdx),
            });
        }

        for (std::size_t residueIdx = 0; residueIdx < residueDistances.size(); ++residueIdx) {
            const double d = residueDistances[residueIdx];
            if (!std::isfinite(d) || d > radiusAngstrom_)
                continue;
            next.push_back(Candidate{
                CandidateKind::Residue,
                model::ResidueAnchor{residueIdx},
                residueIdx,
                std::nullopt,
                d,
                residueLabel(residueIdx),
            });
        }

        for (std::size_t bondIdx = 0; bondIdx < protein_->bondCount(); ++bondIdx) {
            const auto& bond = protein_->bond(bondIdx);
            double d = std::numeric_limits<double>::infinity();
            if (bond.atomIndexA >= 0 && static_cast<std::size_t>(bond.atomIndexA) < atomDistances.size())
                d = std::min(d, atomDistances[static_cast<std::size_t>(bond.atomIndexA)]);
            if (bond.atomIndexB >= 0 && static_cast<std::size_t>(bond.atomIndexB) < atomDistances.size())
                d = std::min(d, atomDistances[static_cast<std::size_t>(bond.atomIndexB)]);
            if (!std::isfinite(d) || d > radiusAngstrom_)
                continue;
            next.push_back(Candidate{
                CandidateKind::Bond,
                model::BondAnchor{bondIdx},
                std::nullopt,
                std::nullopt,
                d,
                bondLabel(bondIdx),
            });
        }

        for (std::size_t ringIdx = 0; ringIdx < protein_->ringCount(); ++ringIdx) {
            const auto& ring = protein_->ring(ringIdx);
            double d = std::numeric_limits<double>::infinity();
            for (int32_t atom : ring.atomIndices) {
                if (atom >= 0 && static_cast<std::size_t>(atom) < atomDistances.size())
                    d = std::min(d, atomDistances[static_cast<std::size_t>(atom)]);
            }
            if (!std::isfinite(d) || d > radiusAngstrom_)
                continue;
            const bool aromatic = ring.IsAromatic();
            model::SignalAnchor ringAnchor = model::RingAnchor{ringIdx};
            if (ring.nativeAxisIndex >= 0) {
                const auto nativeIndex = static_cast<std::size_t>(ring.nativeAxisIndex);
                ringAnchor = aromatic ? model::SignalAnchor{model::AromaticRingAnchor{nativeIndex}}
                                      : model::SignalAnchor{model::SaturatedRingAnchor{nativeIndex}};
            }
            next.push_back(Candidate{
                aromatic ? CandidateKind::AromaticRing : CandidateKind::SaturatedRing,
                std::move(ringAnchor),
                ring.parentResidueIndex >= 0
                    ? std::optional<std::size_t>(static_cast<std::size_t>(ring.parentResidueIndex))
                    : std::nullopt,
                ringIdx,
                d,
                ringLabel(ringIdx),
            });
        }

        for (std::size_t membershipIdx = 0; membershipIdx < protein_->ringMembershipCount(); ++membershipIdx) {
            const auto& membership = protein_->topology().ringMembershipAt(membershipIdx);
            if (membership.atomIndex < 0
                || static_cast<std::size_t>(membership.atomIndex) >= atomDistances.size())
                continue;
            const double d = atomDistances[static_cast<std::size_t>(membership.atomIndex)];
            if (!std::isfinite(d) || d > radiusAngstrom_)
                continue;
            next.push_back(Candidate{
                CandidateKind::RingMembership,
                model::RingMembershipAnchor{membershipIdx},
                std::nullopt,
                membership.ringId >= 0 ? std::optional<std::size_t>(static_cast<std::size_t>(membership.ringId))
                                       : std::nullopt,
                d,
                ringMembershipLabel(membershipIdx),
            });
        }

        // Bond-vector candidates from the iRED + Reorient identity tables.
        // Each table has its own row set (iRED is N-H only; Reorient is
        // N-H + Cα-Hα + C=O); the picker emits one candidate per
        // (residue, kind), distance = min over the row's tail/head atoms,
        // and dedups across the two tables since BondVectorAnchor is
        // semantic (residue, kind), not a per-table row index. Without
        // this enumeration the dialog has no way to surface a specific
        // bond-vector pick — and the controller-side wildcard would only
        // ever return row 0.
        const model::TrajectoryConformation* trajectory = conformation_->asTrajectory();
        const auto* h5 = trajectory ? trajectory->h5() : nullptr;
        if (h5) {
            std::set<std::pair<std::size_t, std::uint8_t>> seenVectors;
            auto walkTable = [&](const model::QtBondVectorTable* table) {
                if (!table || table->n_vectors == 0)
                    return;
                for (std::size_t i = 0; i < table->n_vectors; ++i) {
                    if (i >= table->residue_index.size()
                        || i >= table->kind.size()
                        || i >= table->tail_atom.size()
                        || i >= table->head_atom.size())
                        break;
                    if (table->residue_index[i] < 0)
                        continue;
                    const auto residueIdx = static_cast<std::size_t>(table->residue_index[i]);
                    const auto kind = table->kind[i];
                    const auto key = std::make_pair(residueIdx, kind);
                    if (seenVectors.count(key))
                        continue;

                    double d = std::numeric_limits<double>::infinity();
                    if (table->tail_atom[i] >= 0
                        && static_cast<std::size_t>(table->tail_atom[i]) < atomDistances.size())
                        d = std::min(d, atomDistances[static_cast<std::size_t>(table->tail_atom[i])]);
                    if (table->head_atom[i] >= 0
                        && static_cast<std::size_t>(table->head_atom[i]) < atomDistances.size())
                        d = std::min(d, atomDistances[static_cast<std::size_t>(table->head_atom[i])]);
                    if (!std::isfinite(d) || d > radiusAngstrom_)
                        continue;

                    seenVectors.insert(key);
                    const QString label = residueLabel(residueIdx) + bondVectorKindSuffix(kind);
                    next.push_back(Candidate{
                        CandidateKind::BondVector,
                        model::BondVectorAnchor{residueIdx, kind},
                        residueIdx,
                        std::nullopt,
                        d,
                        label,
                    });
                }
            };
            const model::QtIRedOrderParameters* ired = h5->iredOrderParameters();
            const model::QtReorientationalDynamics* rd = h5->reorientationalDynamics();
            walkTable(ired ? &ired->identity : nullptr);
            walkTable(rd ? &rd->identity : nullptr);
        }

        std::sort(next.begin(), next.end(), [](const Candidate& a, const Candidate& b) {
            if (std::abs(a.distanceAngstrom - b.distanceAngstrom) > 1e-9)
                return a.distanceAngstrom < b.distanceAngstrom;
            if (a.kind != b.kind)
                return kindSortRank(a.kind) < kindSortRank(b.kind);
            return a.label < b.label;
        });
    }

    beginResetModel();
    candidates_ = std::move(next);
    endResetModel();
}

}  // namespace h5reader::app
