// QtManifest — typed projection of extraction_manifest.json.
//
// The manifest is emitted by `src/TopologySidecar.cpp::WriteManifest`
// alongside the 5 sidecar NPYs. It carries:
//   - schema_version + extractor metadata
//   - protein_id (the canonical identifier)
//   - topology block (has_atom_semantic / has_ff_atom_types / has_ff_mass)
//   - enum_vocab (6 enums × {ordinal → name} maps, the load-time
//     contract checked by QtEnumVocab)
//   - axis_sizes (atom, residue, bond, aromatic_ring, saturated_ring,
//     ring, ring_membership) — cross-check against actual NPY row counts
//   - axis_alignment (informational; documents the cross-NPY row
//     correspondence contract)
//
// Strings live here only as projection metadata — schema version,
// extractor name, axis-alignment narrative strings. Zero dispatch on
// these strings.
//
// Cross-platform via QJsonDocument (Qt-native, no POSIX).

#pragma once

#include <QString>
#include <cstddef>

namespace h5reader::io {

// Topology presence block — `topology` object in the manifest.
struct QtManifestTopology {
    QString source;  // e.g. "amber-ff14SB+cifpp"
    bool hasAtomSemantic = false;
    bool hasFfAtomTypes = false;
    bool hasFfMass = false;
};

// Axis sizes — `axis_sizes` object. Used for cross-checking NPY shapes
// at load time. A mismatch is a structural error.
struct QtManifestAxisSizes {
    std::size_t atom = 0;
    std::size_t residue = 0;
    std::size_t bond = 0;
    std::size_t aromaticRing = 0;
    std::size_t saturatedRing = 0;
    std::size_t ring = 0;  // == aromaticRing + saturatedRing
    std::size_t ringMembership = 0;
};

struct QtManifest {
    bool ok = false;
    QString error;  // empty on success

    QString schemaVersion;     // e.g. "1.0"
    QString extractor;         // "nmr_extract"
    QString extractorVersion;  // "0.2.0"
    QString generatedAtUtc;    // ISO 8601 timestamp
    QString proteinId;         // canonical protein identifier

    QtManifestTopology topology;
    QtManifestAxisSizes axisSizes;

    // axis_alignment is informational; loaded as a single
    // human-readable summary block, not parsed into typed structures.
    QString axisAlignmentSummary;

    // Load + parse. On failure: returns a QtManifest with ok=false,
    // error populated, and logs via ErrorBus.
    static QtManifest Load(const QString& path);
};

}  // namespace h5reader::io
