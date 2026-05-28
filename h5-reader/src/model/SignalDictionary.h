// SignalDictionary -- typed vocabulary for chartable trajectory signals.
//
// A SignalKey identifies DATA, not a widget. A strip then decides how to show
// that signal: time-domain, FFT, histogram, etc. This keeps dashboard state as
// pinned bindings to data signals rather than as incidental UI selection.

#pragma once

#include <QString>

#include <cstdint>
#include <optional>
#include <vector>

namespace h5reader::model {

enum class SignalAnchorKind : uint8_t {
    None = 0,     // system-wide signal
    Atom,         // one atom: DFT sigma, atom-axis TR scalar/component
    Residue,      // one residue: phi/psi/omega/chi, DSSP, residue stats
    AtomTuple,    // ordered atoms: distance/angle/dihedral
};

enum class SignalDomainKind : uint8_t {
    FrameTime = 0,   // x-axis is trajectory frame/time
    Frequency,       // x-axis is frequency, e.g. power spectrum
};

enum class SignalValueKind : uint8_t {
    Scalar = 0,
    Angle,
    Distance,
    Power,
    Category,
};

enum class SignalFamily : uint8_t {
    Geometry = 0,
    DftShielding,
    ResidueDihedral,
    AtomTimeSeries,
    ResidueTimeSeries,
    SystemTimeSeries,
};

struct SignalKey {
    SignalFamily family = SignalFamily::Geometry;
    QString name;       // stable dictionary id inside the family, e.g. "chi1"
    QString component;  // optional component, e.g. "total.T0" / "isotropic"

    friend bool operator==(const SignalKey& a, const SignalKey& b) {
        return a.family == b.family && a.name == b.name && a.component == b.component;
    }
};

struct SignalSpec {
    SignalKey key;
    QString label;
    QString unit;
    SignalAnchorKind anchorKind = SignalAnchorKind::None;
    SignalValueKind valueKind = SignalValueKind::Scalar;
    bool finiteScalar = true;
};

struct SignalBinding {
    SignalKey key;
    SignalAnchorKind anchorKind = SignalAnchorKind::None;
    std::optional<std::size_t> atom;
    std::optional<std::size_t> residue;
    std::vector<std::size_t> atomTuple;
    bool followsFocus = false;  // false == pinned dashboard strip
};

std::vector<SignalSpec> CoreSignalDictionary();
std::optional<SignalSpec> FindSignalSpec(const SignalKey& key);

}  // namespace h5reader::model
