// QtSelectionBag — typed event bag for /trajectory/selections/.
//
// The H5 carries one subgroup per selection kind, keyed by the
// mangled C++ std::type_index name (e.g.
// "N3nmr34DftPoseCoordinatorTrajectoryResultE"). The reader translates
// mangled names → typed QtSelectionKind enum at load via a static
// lookup table rather than platform-specific runtime demangling.
//
// Per-record fields per subgroup: frame_idx (uint64), time_ps
// (float64), reason (object/string), metadata_json (object/string).
//
// The reason string is free-form diagnostic text used only for display and is
// never dispatched on.
//
// metadata_json is a per-kind structured blob; we PARSE it at load
// into typed std::variant<QtDftPoseMeta, QtRmsdSpikeMeta,
// QtChiRotamerMeta, std::monostate>. Failed optional metadata parsing leaves
// std::monostate and logs a warning.
//
// Filtering queries (indicesByKind, indicesInTimeRange) are linear scans over
// the loaded event list.

#pragma once

#include <QString>
#include <cstddef>
#include <cstdint>
#include <variant>
#include <vector>

namespace h5reader::model {

// ──────────────────────────────────────────────────────────────────
// Typed selection kinds. Mirror the mangled C++ type_index names in
// the H5; mapping table lives in the .cpp.
// ──────────────────────────────────────────────────────────────────

enum class QtSelectionKind : int8_t {
    Unknown = -1,
    DftPoseCoordinator = 0,
    RmsdSpikeSelection = 1,
    ChiRotamerSelection = 2,
};

const char* NameForSelectionKind(QtSelectionKind k);

// Static mangled-name → typed-kind lookup. Returns QtSelectionKind::Unknown
// when the mangled name is not in the lookup table; loader logs Warn.
QtSelectionKind ParseSelectionGroupName(const QString& mangled_name);


// ──────────────────────────────────────────────────────────────────
// Typed metadata variants per selection kind. v1: minimal fields the
// metadata_json is known to carry; loader parses best-effort.
// ──────────────────────────────────────────────────────────────────

struct QtDftPoseMeta {
    // Optional metadata_json fields; missing keys retain these defaults.
    double score = 0.0;
    QString method;
};

struct QtRmsdSpikeMeta {
    double rmsd_angstroms = 0.0;
    double threshold = 0.0;
};

struct QtChiRotamerMeta {
    int32_t residue_index = -1;
    int8_t chi_axis = -1;      // 0..3
    int8_t from_rotamer = -1;  // bin index before transition
    int8_t to_rotamer = -1;    // bin index after
};


// ──────────────────────────────────────────────────────────────────
// QtSelectionEvent — one record from a selections subgroup.
// ──────────────────────────────────────────────────────────────────

struct QtSelectionEvent {
    QtSelectionKind kind = QtSelectionKind::Unknown;
    uint64_t frame_idx = 0;
    double time_ps = 0.0;

    QString reason;             // display only
    QString metadata_json_raw;  // unparsed JSON; refinement extracts typed

    std::variant<std::monostate, QtDftPoseMeta, QtRmsdSpikeMeta, QtChiRotamerMeta> meta;
};


// ──────────────────────────────────────────────────────────────────
// QtSelectionBag — flat vector of events with typed filtering queries.
// ──────────────────────────────────────────────────────────────────

class QtSelectionBag {
public:
    void push(QtSelectionEvent ev) { events_.push_back(std::move(ev)); }

    std::size_t count() const { return events_.size(); }
    const QtSelectionEvent& at(std::size_t i) const { return events_[i]; }
    const std::vector<QtSelectionEvent>& events() const { return events_; }

    // O(N) linear scans; sized for thousands of events at most.
    std::vector<std::size_t> indicesByKind(QtSelectionKind k) const;
    std::vector<std::size_t> indicesInTimeRange(double t_lo_ps, double t_hi_ps) const;
    std::size_t countByKind(QtSelectionKind k) const;

private:
    std::vector<QtSelectionEvent> events_;
};


}  // namespace h5reader::model
