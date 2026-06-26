// MeasurementsDock -- a small dockable readout of the current AtomSelection's
// geometry measurement (distance / angle / dihedral) for the focused frame,
// redrawn live as frames play. The value comes from model::Measure (pure compute
// over Conformation positions); this is just the view.
//
// Realises the "text info lives in a window, not on the molecule" rule (see the
// scene-text-belongs-in-windows memory): the MeasurementOverlay draws only the
// spheres + connecting lines in the 3-D scene; the measured VALUE shows here.
//
// Source stays ASCII; the degree + Angstrom glyphs come via QChar code points.

#pragma once

#include <QDockWidget>

#include <cstddef>
#include <vector>

class QLabel;

namespace h5reader::model {
class QtProtein;
class Conformation;
}

namespace h5reader::app {

class MeasurementsDock final : public QDockWidget {
    Q_OBJECT
public:
    explicit MeasurementsDock(QWidget* parent = nullptr);

    // Bind the typed model (base Conformation is fine: geometry is a pure
    // function of positions). Pass (nullptr, nullptr) to unbind on unload.
    void setContext(const model::QtProtein* protein, const model::Conformation* conf);

public slots:
    // The ordered selection tuple (slot order): 2 atoms -> distance, 3 -> angle
    // (vertex = middle), 4 -> dihedral; any other count -> the empty hint.
    void setAtoms(const std::vector<std::size_t>& atoms);
    void setFrame(int frame);
    void clear();

private:
    void recompute();

    const model::QtProtein*    protein_ = nullptr;
    const model::Conformation* conf_    = nullptr;
    std::vector<std::size_t>   atoms_;
    int                        frame_ = 0;

    QLabel* kindLabel_  = nullptr;  // "Distance"/"Angle"/"Dihedral" or the hint
    QLabel* valueLabel_ = nullptr;  // the measured value, large
    QLabel* atomsLabel_ = nullptr;  // the ordered atoms, one per line
};

}  // namespace h5reader::app
