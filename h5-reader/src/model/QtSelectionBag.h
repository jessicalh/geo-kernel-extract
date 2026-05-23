// QtSelectionBag — typed event bag for /trajectory/selections/.
//
// The H5 carries one subgroup per selection kind, keyed by the
// mangled C++ std::type_index name (e.g.
// "N3nmr34DftPoseCoordinatorTrajectoryResultE"). The reader translates
// mangled names → typed QtSelectionKind enum at load via a static
// lookup table (design §11.F — static table, not runtime demangling,
// for cross-platform).
//
// Per-record fields per subgroup: frame_idx (uint64), time_ps
// (float64), reason (object/string), metadata_json (object/string).
//
// The reason string is genuine free-form diagnostic text (per Agent 3:
// "chi_transition_phi_side_chain_4_TRP" etc.) — display only, never
// dispatched on. We store it as QString.
//
// metadata_json is a per-kind structured blob; we PARSE it at load
// into typed std::variant<QtDftPoseMeta, QtRmsdSpikeMeta,
// QtChiRotamerMeta, std::monostate>. v1: minimal typed parsing (best-
// effort); failures keep std::monostate but log Warn. Refinement in a
// later session as the metadata shape stabilises.
//
// Filtering queries (indicesByKind, indicesInTimeRange) are O(N)
// linear scans — adequate for the 8-6031 event range in 1P9J fixture.

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
    // Fields populated from metadata_json best-effort; left as defaults
    // when JSON parse misses a key. The "score" / "method" pattern from
    // Agent 3's example is the current expectation.
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
