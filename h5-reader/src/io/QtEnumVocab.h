// QtEnumVocab — typed projection of the manifest's `enum_vocab` block.
//
// The manifest declares the int8 ordinal → name mapping for six enums
// used in the sidecar NPYs: terminal_state, bond_order, bond_category,
// ring_kind, ring_type_index, residue_type. The reader has typed
// counterparts (TerminalState, BondOrder, BondCategory, RingKind,
// RingTypeIndex, AminoAcid) with hardcoded ordinals.
//
// QtEnumVocab loads the manifest entries and validates them against
// the reader's typed enums via per-enum NameForX() helpers. A drift
// (manifest says ordinal 5 = "HisImidazole" but the reader's
// RingTypeIndex(5) NameForRingType()'s "HieImidazole") logs an
// ErrorBus warning at load time — the silent mis-decode failure mode
// from the dead-format era is structurally impossible here.
//
// QtEnumVocab also provides DisplayName(kind, ordinal) for inspector
// tooltips that want the manifest's sanctioned name rather than the
// reader's hardcoded literal.
//
// No-strings discipline: the QString values are display-projection
// surfaces only. Chemistry dispatch goes through the reader's typed
// enums; the vocab block exists for drift detection + display.

#pragma once

#include <QString>
#include <array>
#include <cstdint>
#include <unordered_map>

namespace h5reader::io {

class QtEnumVocab {
public:
    enum class Kind : int8_t {
        TerminalState = 0,
        BondOrder = 1,
        BondCategory = 2,
        RingKind = 3,
        RingTypeIndex = 4,
        AminoAcid = 5,
    };
    static constexpr int kKindCount = 6;

    // Load + parse the manifest's enum_vocab block. The manifest_path
    // is passed for diagnostics only; the QJsonObject is the parsed
    // value of `root["enum_vocab"]`.
    static QtEnumVocab LoadFromManifest(const QString& manifest_path);

    bool ok() const { return ok_; }
    const QString& error() const { return error_; }

    // Validate the loaded ordinals against the reader's typed enum
    // NameForX() helpers. Returns true if all manifest entries agree
    // with the reader. Each disagreement is logged via ErrorBus
    // (Severity::Warning) with the specific (kind, ordinal, expected,
    // got) tuple. Called once at load by QtTopologySidecar.
    bool ValidateAgainstReaderEnums() const;

    // Display lookup for inspector tooltips. Returns the manifest's
    // canonical name for the given ordinal, or empty QString if the
    // ordinal is not present in this vocab.
    QString DisplayName(Kind kind, int ordinal) const;

private:
    bool ok_ = false;
    QString error_;
    // entries_[kind_int][ordinal] = name. Outer array is fixed-size
    // by Kind enum; inner is sparse map (most enums have <16 entries).
    std::array<std::unordered_map<int, QString>, kKindCount> entries_;

    friend class QtEnumVocabLoader;  // file-local helper in .cpp
};

}  // namespace h5reader::io
