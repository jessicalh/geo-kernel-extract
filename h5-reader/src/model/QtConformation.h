// QtConformation — the trajectory. Mirrors `nmr::ProteinConformation`
// in role (per-protein computed/projected state) plus
// `nmr::TrajectoryProtein` in scope (T frames + per-TR buffers).
//
// Owns the QtTrajectoryH5 (which holds all eager-loaded per-TR
// buffers). Back-pointer to QtProtein for typed topology lookups.
//
// QtFrame is constructed cheaply on demand — see frame(t). All per-
// frame data lives in QtTrajectoryH5's per-TR buffers; QtFrame is a
// thin (conformation, frame_index) view that delegates accessor calls
// to the right buffer.

#pragma once

#include "../io/QtTrajectoryH5.h"
#include "QtFrame.h"

#include <cstddef>
#include <memory>

namespace h5reader::model {

class QtProtein;

class QtConformation {
public:
    QtConformation(const QtProtein* protein, std::unique_ptr<h5reader::io::QtTrajectoryH5> h5);

    ~QtConformation() = default;
    QtConformation(const QtConformation&) = delete;
    QtConformation& operator=(const QtConformation&) = delete;
    QtConformation(QtConformation&&) = delete;
    QtConformation& operator=(QtConformation&&) = delete;

    // ----- Back-references -----
    const QtProtein* protein() const { return protein_; }
    const h5reader::io::QtTrajectoryH5* h5() const { return h5_.get(); }

    // ----- Trajectory parameters -----
    std::size_t frameCount() const { return h5_->frameCount(); }
    double startTimePicoseconds() const;
    double endTimePicoseconds() const;

    // ----- Frame access -----
    QtFrame frame(std::size_t t) const { return QtFrame(this, t); }

    // ----- Convenience (mirrors old QtConformation surface) -----
    std::size_t ringCount() const;

private:
    const QtProtein* protein_;
    std::unique_ptr<h5reader::io::QtTrajectoryH5> h5_;
};

}  // namespace h5reader::model
