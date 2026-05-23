// QtSelectionBag implementation — mangled-name dispatch table + simple
// filtering queries.
//
// The mangled C++ type_index names in /trajectory/selections/ subgroups
// are itanium-ABI mangled (g++ / clang on Linux/macOS). The reader
// hardcodes the known mappings here — adding a new selection kind
// means appending one entry to kSelectionMangledNames + one enum value.
// Cross-platform without platform-specific demangling code (design
// §11.F).

#include "QtSelectionBag.h"

namespace h5reader::model {

namespace {

struct MangledMapping {
    const char* mangled;
    QtSelectionKind kind;
};

// The three known selection kinds from the 1P9J fixture. Agent 3's
// inspection reported them as:
//   N3nmr34DftPoseCoordinatorTrajectoryResultE
//   N3nmr34RmsdSpikeSelectionTrajectoryResultE
//   N3nmr35ChiRotamerSelectionTrajectoryResultE
//
// The number after "N3nmr" is the length of the class name minus 1
// (itanium mangling convention for the class-name part of a nested
// name). New kinds: add a row here + a QtSelectionKind enum value.
constexpr MangledMapping kSelectionMangledNames[] = {
    {"N3nmr34DftPoseCoordinatorTrajectoryResultE", QtSelectionKind::DftPoseCoordinator},
    {"N3nmr34RmsdSpikeSelectionTrajectoryResultE", QtSelectionKind::RmsdSpikeSelection},
    {"N3nmr35ChiRotamerSelectionTrajectoryResultE", QtSelectionKind::ChiRotamerSelection},
};

}  // namespace


const char* NameForSelectionKind(QtSelectionKind k) {
    switch (k) {
    case QtSelectionKind::DftPoseCoordinator:
        return "DftPoseCoordinator";
    case QtSelectionKind::RmsdSpikeSelection:
        return "RmsdSpikeSelection";
    case QtSelectionKind::ChiRotamerSelection:
        return "ChiRotamerSelection";
    case QtSelectionKind::Unknown:
        return "Unknown";
    }
    return "Unknown";
}


QtSelectionKind ParseSelectionGroupName(const QString& mangled_name) {
    const QByteArray needle = mangled_name.toLatin1();
    for (const auto& row : kSelectionMangledNames) {
        if (needle == row.mangled)
            return row.kind;
    }
    return QtSelectionKind::Unknown;
}


std::vector<std::size_t> QtSelectionBag::indicesByKind(QtSelectionKind k) const {
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < events_.size(); ++i) {
        if (events_[i].kind == k)
            out.push_back(i);
    }
    return out;
}


std::vector<std::size_t> QtSelectionBag::indicesInTimeRange(double t_lo_ps, double t_hi_ps) const {
    std::vector<std::size_t> out;
    for (std::size_t i = 0; i < events_.size(); ++i) {
        const double t = events_[i].time_ps;
        if (t >= t_lo_ps && t <= t_hi_ps)
            out.push_back(i);
    }
    return out;
}


std::size_t QtSelectionBag::countByKind(QtSelectionKind k) const {
    std::size_t n = 0;
    for (const auto& ev : events_) {
        if (ev.kind == k)
            ++n;
    }
    return n;
}


}  // namespace h5reader::model
