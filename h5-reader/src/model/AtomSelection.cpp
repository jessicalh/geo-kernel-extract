// AtomSelection — gesture interpretation, model mutation, and the role
// surface. The comments record WHY each branch fires which signal, and which
// begin/end pair guards each structural change (get those wrong and a bound
// view dereferences a stale row).

#include "AtomSelection.h"

#include "QtProtein.h"

#include "../diagnostics/ObjectCensus.h"
#include "../diagnostics/ThreadGuard.h"

#include <QLoggingCategory>

#include <algorithm>

namespace h5reader::model {

namespace {
Q_LOGGING_CATEGORY(cSel, "h5reader.selection")
}

const char* NameForGeometryKind(GeometryKind k) {
    switch (k) {
    case GeometryKind::None:     return "—";
    case GeometryKind::Distance: return "Distance";
    case GeometryKind::Angle:    return "Angle";
    case GeometryKind::Dihedral: return "Dihedral";
    }
    return "—";
}

AtomSelection::AtomSelection(const QtProtein* protein, QObject* parent)
    : QAbstractListModel(parent),
      protein_(protein)
{
    CENSUS_REGISTER(this);
    setObjectName(QStringLiteral("AtomSelection"));
}

// ----- QAbstractListModel ---------------------------------------------------

int AtomSelection::rowCount(const QModelIndex& parent) const {
    // Flat list: children of a valid index do not exist.
    if (parent.isValid())
        return 0;
    return static_cast<int>(atoms_.size());
}

QVariant AtomSelection::data(const QModelIndex& index, int role) const {
    if (!index.isValid())
        return {};
    const int row = index.row();
    if (row < 0 || static_cast<std::size_t>(row) >= atoms_.size())
        return {};
    const std::size_t atomIdx = atoms_[static_cast<std::size_t>(row)];

    switch (role) {
    case Qt::DisplayRole:
        return rowLabel(static_cast<std::size_t>(row));
    case Qt::DecorationRole:
    case SlotColorRole:
        return QVariant::fromValue(SlotColor(static_cast<std::size_t>(row)));
    case Qt::ToolTipRole:
        return QStringLiteral("slot %1 · atom #%2%3")
            .arg(row)
            .arg(atomIdx)
            .arg(focus_.has_value() && focus_.value() == atomIdx
                     ? QStringLiteral(" · focus")
                     : QString());
    case AtomIndexRole:
        return static_cast<int>(atomIdx);
    case SlotRole:
        return row;
    case IsFocusRole:
        return focus_.has_value() && focus_.value() == atomIdx;
    default:
        return {};
    }
}

QHash<int, QByteArray> AtomSelection::roleNames() const {
    QHash<int, QByteArray> r = QAbstractListModel::roleNames();
    r[AtomIndexRole] = "atomIndex";
    r[SlotRole]      = "slot";
    r[SlotColorRole] = "slotColor";
    r[IsFocusRole]   = "isFocus";
    return r;
}

// ----- Typed surface --------------------------------------------------------

bool AtomSelection::contains(std::size_t atomIdx) const {
    return std::find(atoms_.begin(), atoms_.end(), atomIdx) != atoms_.end();
}

int AtomSelection::slotOf(std::size_t atomIdx) const {
    for (std::size_t i = 0; i < atoms_.size(); ++i) {
        if (atoms_[i] == atomIdx)
            return static_cast<int>(i);
    }
    return -1;
}

GeometryKind AtomSelection::geometryKind() const {
    switch (atoms_.size()) {
    case 2:  return GeometryKind::Distance;
    case 3:  return GeometryKind::Angle;
    case 4:  return GeometryKind::Dihedral;
    default: return GeometryKind::None;  // 0 or 1 atoms define no measurement
    }
}

std::array<double, 3> AtomSelection::SlotColorRgb(std::size_t slot) {
    switch (slot) {
    case 0:  return {0.902, 0.624, 0.000};  // orange         #E69F00
    case 1:  return {0.337, 0.706, 0.914};  // sky blue       #56B4E9
    case 2:  return {0.000, 0.620, 0.451};  // bluish green   #009E73
    case 3:  return {0.800, 0.475, 0.655};  // reddish purple #CC79A7
    default: return {0.600, 0.600, 0.600};  // grey — out of palette
    }
}

QColor AtomSelection::SlotColor(std::size_t slot) {
    const std::array<double, 3> c = SlotColorRgb(slot);
    return QColor::fromRgbF(c[0], c[1], c[2]);
}

// ----- Mutation -------------------------------------------------------------

void AtomSelection::applyPick(std::size_t atomIdx, Qt::KeyboardModifiers mods) {
    ASSERT_THREAD(this);

    // The protein spine is the authority on the atom count; a stale or
    // out-of-range index is dropped loudly rather than corrupting the set.
    if (protein_ && atomIdx >= protein_->atomCount()) {
        qCWarning(cSel).noquote()
            << "pick out of range:" << atomIdx << ">= atomCount" << protein_->atomCount();
        return;
    }

    const bool addMode = mods.testFlag(Qt::ShiftModifier);

    if (!addMode) {
        // Plain pick — replace the whole set with this one atom. Keeps
        // today's "double-click inspects an atom" muscle memory: a bare click
        // always lands on exactly one focused atom. The membership change is
        // wholesale, so a model reset is the correct (and simplest) signal.
        beginResetModel();
        atoms_ = {atomIdx};
        focus_ = atomIdx;
        endResetModel();
        qCInfo(cSel).noquote() << "focus" << atomIdx << "| replace | n=1";
        emit changed();
        emit focusChanged(atomIdx);
        return;
    }

    // Shift+pick — toggle membership in the ordered set.
    const int slot = slotOf(atomIdx);
    if (slot >= 0) {
        // Already selected → remove it (Shift-click again to deselect).
        beginRemoveRows(QModelIndex(), slot, slot);
        atoms_.erase(atoms_.begin() + slot);
        endRemoveRows();

        if (atoms_.empty()) {
            focus_.reset();
            qCInfo(cSel).noquote() << "removed" << atomIdx << "| n=0 | cleared";
            emit changed();
            emit cleared();
            return;
        }
        // Focus moves to the most recently added survivor (back of the set).
        focus_ = atoms_.back();
        emitFocusRolesChanged();
        qCInfo(cSel).noquote()
            << "removed" << atomIdx << "| n=" << atoms_.size() << "| focus" << *focus_;
        emit changed();
        emit focusChanged(*focus_);
        return;
    }

    if (full()) {
        // The geometric cap is reached; don't silently evict a member.
        qCInfo(cSel).noquote()
            << "selection full (" << kMaxAtoms
            << "); plain-click to reset, or Shift-click a member to remove";
        return;
    }

    // Append at the next slot; the new atom becomes the focus.
    const int newRow = static_cast<int>(atoms_.size());
    beginInsertRows(QModelIndex(), newRow, newRow);
    atoms_.push_back(atomIdx);
    focus_ = atomIdx;
    endInsertRows();
    emitFocusRolesChanged();  // the previous focus row lost its marker
    qCInfo(cSel).noquote()
        << "added" << atomIdx << "| slot" << newRow << "| n=" << atoms_.size();
    emit changed();
    emit focusChanged(atomIdx);
}

void AtomSelection::clear() {
    ASSERT_THREAD(this);
    if (atoms_.empty())
        return;
    beginResetModel();
    atoms_.clear();
    focus_.reset();
    endResetModel();
    qCInfo(cSel).noquote() << "cleared";
    emit changed();
    emit cleared();
}

void AtomSelection::bulkSet(const std::vector<std::size_t>& atoms) {
    ASSERT_THREAD(this);

    // Validate against the protein's atomCount and the geometric cap; drop
    // out-of-range entries loudly rather than silently corrupting the set.
    // The REST handler validates first; this is the belt+suspenders pass for
    // any future caller.
    std::vector<std::size_t> kept;
    kept.reserve(std::min(atoms.size(), kMaxAtoms));
    for (std::size_t atomIdx : atoms) {
        if (kept.size() >= kMaxAtoms) {
            qCWarning(cSel).noquote()
                << "bulkSet truncated to" << kMaxAtoms << "atoms; ignored" << atomIdx;
            break;
        }
        if (protein_ && atomIdx >= protein_->atomCount()) {
            qCWarning(cSel).noquote()
                << "bulkSet out of range:" << atomIdx
                << ">= atomCount" << protein_->atomCount() << "— dropped";
            continue;
        }
        kept.push_back(atomIdx);
    }

    // Reset model: the whole membership is being replaced wholesale. Per the
    // QAbstractItemModel contract a beginResetModel/endResetModel pair is the
    // right signal for full-content replacement and is the same shape
    // applyPick() uses for the plain-pick branch (line 149).
    beginResetModel();
    atoms_ = std::move(kept);
    if (atoms_.empty())
        focus_.reset();
    else
        focus_ = atoms_.back();
    endResetModel();

    qCInfo(cSel).noquote()
        << "bulkSet | n=" << atoms_.size()
        << "| focus=" << (focus_.has_value() ? QString::number(*focus_) : QStringLiteral("none"));

    emit changed();
    if (atoms_.empty())
        emit cleared();
    else
        emit focusChanged(*focus_);
}

// ----- Helpers --------------------------------------------------------------

QString AtomSelection::rowLabel(std::size_t slot) const {
    if (!protein_ || slot >= atoms_.size())
        return QStringLiteral("—");
    const std::size_t a = atoms_[slot];
    if (a >= protein_->atomCount())
        return QStringLiteral("#%1").arg(a);

    const QtAtom& atom     = protein_->atom(a);
    const QString atomName = protein_->atomNames(a).amber;

    if (atom.residueIndex >= 0
        && static_cast<std::size_t>(atom.residueIndex) < protein_->residueCount()) {
        const QtResidue& res = protein_->residue(static_cast<std::size_t>(atom.residueIndex));
        const QString res3   = protein_->residueLabel(
            static_cast<std::size_t>(atom.residueIndex),
            NamingConvention::Amber, NamingSource::Verbatim);
        return QStringLiteral("%1%2:%3").arg(res3).arg(res.address.residueNumber).arg(atomName);
    }
    return QStringLiteral("#%1:%2").arg(a).arg(atomName);
}

void AtomSelection::emitFocusRolesChanged() {
    // The focus marker is a per-row boolean; when focus moves WITHOUT a
    // structural change, tell views the IsFocusRole of every row may have
    // flipped. ≤4 rows, so a blanket dataChanged is trivial and correct.
    if (atoms_.empty())
        return;
    const QModelIndex top = index(0, 0);
    const QModelIndex bot = index(static_cast<int>(atoms_.size()) - 1, 0);
    emit dataChanged(top, bot, {IsFocusRole, Qt::ToolTipRole});
}

}  // namespace h5reader::model
