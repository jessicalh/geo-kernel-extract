// NewmanDock — a small dockable Newman projection of the focused residue's
// backbone phi/psi torsion, redrawn live as frames play. Mirrors
// QtAtomInspectorDock's reveal-on-focus + frameChanged update pattern. The
// geometry comes from NewmanProjection (pure compute); this is the view.
//
// Source is ASCII-only on purpose (no Greek phi/psi or degree glyphs) until the
// build opts into /utf-8 for reader-authored files; labels read "phi"/"psi".

#pragma once

#include "NewmanProjection.h"

#include <QDockWidget>
#include <QString>
#include <QWidget>

#include <cstddef>

class QLabel;
class QToolButton;

namespace h5reader::model {
class QtProtein;
class Conformation;
}

namespace h5reader::app {

// The drawing surface: paints one NewmanProjection (front hub + spokes from the
// centre, back atom as a circle with spokes from the rim, reference bonds
// accented, torsion read out). Holds the drawn torsion plus its companion for a
// numeric line.
class NewmanView final : public QWidget {
    Q_OBJECT
public:
    explicit NewmanView(QWidget* parent = nullptr);
    void setData(const NewmanProjection& drawn, const NewmanProjection& other);
    void setEmpty(const QString& hint);

protected:
    void paintEvent(QPaintEvent*) override;
    QSize minimumSizeHint() const override { return QSize(210, 230); }

private:
    NewmanProjection drawn_;
    NewmanProjection other_;
    QString          hint_ = QStringLiteral("Select an atom to see its backbone torsion.");
    bool             has_  = false;
};

class NewmanDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit NewmanDock(QWidget* parent = nullptr);

    // Bind the typed model (base Conformation is fine: torsions are
    // rigid-transform invariant). Pass (nullptr, nullptr) to unbind on unload.
    void setContext(const model::QtProtein* protein, const model::Conformation* conf);

public slots:
    void setFocusAtom(std::size_t atomIdx);  // maps atom -> residue, recomputes
    void setFrame(int frame);
    void clear();

private:
    void recompute();

    const model::QtProtein*    protein_ = nullptr;
    const model::Conformation* conf_    = nullptr;
    int        residue_ = -1;            // -1 == none
    int        frame_   = 0;
    NewmanKind kind_    = NewmanKind::Phi;

    NewmanView*  view_    = nullptr;
    QToolButton* toggle_  = nullptr;
    QLabel*      caption_ = nullptr;
};

}  // namespace h5reader::app
