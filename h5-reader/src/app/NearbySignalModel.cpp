#include "NearbySignalModel.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"
#include "../model/Conformation.h"
#include "../model/QtAtom.h"
#include "../model/QtProtein.h"
#include "../model/QtResidue.h"

#include <QBrush>

#include <algorithm>
#include <cmath>
#include <limits>

namespace h5reader::app {

namespace {
double distance(const model::Vec3& a, const model::Vec3& b) {
    const double dx = a.x() - b.x();
    const double dy = a.y() - b.y();
    const double dz = a.z() - b.z();
    return std::sqrt(dx * dx + dy * dy + dz * dz);
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
            return c.kind == CandidateKind::Atom ? QStringLiteral("atom") : QStringLiteral("residue");
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
        return c.kind == CandidateKind::Atom ? QStringLiteral("atom") : QStringLiteral("residue");
    if (role == AtomIndexRole)
        return c.atom ? QVariant::fromValue(static_cast<qulonglong>(*c.atom)) : QVariant();
    if (role == ResidueIndexRole)
        return c.residue ? QVariant::fromValue(static_cast<qulonglong>(*c.residue)) : QVariant();
    if (role == DistanceRole)
        return c.distanceAngstrom;
    if (role == AnchorKindRole) {
        return c.kind == CandidateKind::Atom
            ? static_cast<int>(model::SignalAnchorKind::Atom)
            : static_cast<int>(model::SignalAnchorKind::Residue);
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

void NearbySignalModel::rebuild() {
    ASSERT_THREAD(this);
    std::vector<Candidate> next;
    if (protein_ && conformation_ && anchorAtom_ && *anchorAtom_ < protein_->atomCount()
        && frame_ >= 0 && static_cast<std::size_t>(frame_) < conformation_->frameCount()) {
        const model::Vec3 focus = conformation_->atomPosition(static_cast<std::size_t>(frame_), *anchorAtom_);
        std::vector<double> residueDistances(protein_->residueCount(), std::numeric_limits<double>::infinity());

        for (std::size_t atomIdx = 0; atomIdx < protein_->atomCount(); ++atomIdx) {
            const double d = distance(focus, conformation_->atomPosition(static_cast<std::size_t>(frame_), atomIdx));
            if (!std::isfinite(d) || d > radiusAngstrom_)
                continue;
            const auto& atom = protein_->atom(atomIdx);
            if (atom.residueIndex >= 0
                && static_cast<std::size_t>(atom.residueIndex) < residueDistances.size()) {
                double& best = residueDistances[static_cast<std::size_t>(atom.residueIndex)];
                best = std::min(best, d);
            }
            next.push_back(Candidate{
                CandidateKind::Atom,
                atomIdx,
                atom.residueIndex >= 0 ? std::optional<std::size_t>(static_cast<std::size_t>(atom.residueIndex))
                                        : std::nullopt,
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
                std::nullopt,
                residueIdx,
                d,
                residueLabel(residueIdx),
            });
        }

        std::sort(next.begin(), next.end(), [](const Candidate& a, const Candidate& b) {
            if (std::abs(a.distanceAngstrom - b.distanceAngstrom) > 1e-9)
                return a.distanceAngstrom < b.distanceAngstrom;
            if (a.kind != b.kind)
                return a.kind == CandidateKind::Residue;
            return a.label < b.label;
        });
    }

    beginResetModel();
    candidates_ = std::move(next);
    endResetModel();
}

}  // namespace h5reader::app
