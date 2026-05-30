// AtomSelection — the live UI atom selection, modelled as a proper Qt list
// model: an ordered set of up to four atoms plus a "focus" member.
//
// Design (decided 2026-05-26, memory
// project_h5reader_killer_app_multiatom_compare_20260526 — ONE unified
// selection, focus-tracked; user chose "full Qt model now" over a value type):
// the group IS the QAbstractListModel — one source of truth, not a manager
// plus a wrapper model to keep in sync. The picker stays dumb and reports a
// pick + the keyboard modifiers at click time; THIS object interprets the
// gesture, mutates the model with the proper begin/end row signals, and fans
// typed changes out to the consumers:
//   * focus  → the single-atom consumers (inspector dock, time-series dock)
//              that show one atom's full state. focus follows the most
//              recently touched member.
//   * the ordered set → the measurement overlay (≤4 colour-coded spheres,
//              and later distance/angle/dihedral) and the comparison plots.
//   * the model rows → any QAbstractItemView (the SelectionDock panel today;
//              the comparison table tomorrow).
//
// Gesture policy:
//   * plain pick  → replace the whole set with {idx} (focus = idx).
//   * Shift+pick  → toggle idx in the ordered set: append if absent and there
//                   is room (max 4), REMOVE if already present.
// Slot order = insertion order, so the overlay colours spheres 0..3 stably and
// the geometry helpers read the set as an ordered tuple (a-b distance,
// a-b-c angle, a-b-c-d dihedral).
//
// Two surfaces, deliberately: the role-based data()/roleNames() for views, and
// a plain typed C++ API (atoms()/focus()/geometryKind()/SlotColor) for the
// renderer + geometry, which want raw indices, not QModelIndexes.
//
// NOT QtSelectionBag: that is the H5 *event* bag (DftPoseCoordinator /
// RmsdSpikeSelection / ChiRotamerSelection curation records). This is the
// interactive UI selection. "Selection" is overloaded in this codebase by
// history; the two never meet.

#pragma once

#include <QAbstractListModel>
#include <QColor>
#include <QString>

#include <array>
#include <cstddef>
#include <optional>
#include <vector>

namespace h5reader::model {

class QtProtein;

// What a k-tuple of selected atoms MEASURES — purely a function of the
// selection's cardinality. The VALUE (which needs positions) is computed by
// model::Measure() in ConformationGeometry; this is just the group answering
// "what am I" from its size.
enum class GeometryKind { None, Distance, Angle, Dihedral };
const char* NameForGeometryKind(GeometryKind k);

class AtomSelection final : public QAbstractListModel {
    Q_OBJECT

public:
    // The killer app measures up to a dihedral (four atoms). The cap is a
    // physics statement, not a UI nicety: distance(2) / angle(3) /
    // dihedral(4) are the geometric observables a selection can define.
    static constexpr std::size_t kMaxAtoms = 4;

    // Custom roles for views. DisplayRole = "TRP42:CA"; DecorationRole =
    // the slot colour swatch (matches the 3-D sphere). The typed roles feed
    // a delegate / the comparison plots that want the raw values.
    enum Roles {
        AtomIndexRole = Qt::UserRole + 1,  // int    — global atom index
        SlotRole,                          // int    — 0-based slot
        SlotColorRole,                     // QColor — Okabe-Ito slot colour
        IsFocusRole,                       // bool   — is this row the focus
    };

    // protein supplies the atom-count bound + the row labels; the selection
    // never outlives it (both owned by ReaderMainWindow for the window's life).
    explicit AtomSelection(const QtProtein* protein, QObject* parent = nullptr);

    // ----- QAbstractListModel -----
    int rowCount(const QModelIndex& parent = QModelIndex()) const override;
    QVariant data(const QModelIndex& index, int role = Qt::DisplayRole) const override;
    QHash<int, QByteArray> roleNames() const override;

    // ----- Typed C++ surface (renderer + geometry; no QModelIndex tax) -----
    const std::vector<std::size_t>& atoms() const { return atoms_; }  // slot order
    std::size_t count() const { return atoms_.size(); }
    bool        empty() const { return atoms_.empty(); }
    bool        full()  const { return atoms_.size() >= kMaxAtoms; }
    bool        contains(std::size_t atomIdx) const;
    int         slotOf(std::size_t atomIdx) const;  // 0-based, or -1 if absent

    bool        hasFocus() const { return focus_.has_value(); }
    std::size_t focus()    const { return focus_.value_or(0); }

    GeometryKind geometryKind() const;

    // Okabe-Ito colour-universal palette (Okabe & Ito 2008,
    // "Color Universal Design") — the ONE source for slot colours, read by
    // BOTH the model's DecorationRole and the MeasurementOverlay spheres, so a
    // panel swatch always matches its sphere. RGB components in 0..1.
    static std::array<double, 3> SlotColorRgb(std::size_t slot);
    static QColor                SlotColor(std::size_t slot);

public slots:
    // Interpret one pick. mods is the keyboard state captured at click time.
    // Out-of-range idx, and an append when the set is full, are no-ops
    // (logged at info — a wrongly-empty selection should be a question the
    // log answers, not a silent state).
    void applyPick(std::size_t atomIdx, Qt::KeyboardModifiers mods);

    // Bulk replace the entire selection with `atoms`, in order, focus =
    // atoms.back() if non-empty. Validates each index against the protein's
    // atomCount (out-of-range are dropped with a warning; the REST handler
    // does primary validation, this is belt+suspenders). One `changed()` and
    // one `focusChanged()` (or `cleared()`) emit at the end, NOT per atom —
    // consumers that rebuild from the full set (MeasurementOverlay,
    // comparison plots) only need one update.
    void bulkSet(const std::vector<std::size_t>& atoms);

    // Empty the selection (a future Clear action / Esc).
    void clear();

signals:
    // Set membership changed (atom appended / removed / replaced). The
    // measurement overlay + comparison plots rebuild from atoms(). Emitted
    // alongside the QAbstractItemModel structural signals, as the single
    // "something changed" hook for non-view consumers.
    void changed();

    // The focus member changed. The inspector + time-series dock retarget.
    void focusChanged(std::size_t atomIdx);

    // The selection became empty. Focus consumers clear their views.
    void cleared();

private:
    QString rowLabel(std::size_t slot) const;  // "TRP42:CA"
    void    emitFocusRolesChanged();           // dataChanged(IsFocusRole) for all rows

    const QtProtein*           protein_ = nullptr;
    std::vector<std::size_t>   atoms_;   // insertion order == slot order
    std::optional<std::size_t> focus_;
};

}  // namespace h5reader::model
