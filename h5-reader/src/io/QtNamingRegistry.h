// QtNamingRegistry — the string/enum boundary for the reader.
//
// Decodes three-letter amino-acid codes into AminoAcid enum values.
// Translates atom names between force-field conventions (CHARMM "HN"
// vs PDB "H", "HSD"/"HSE"/"HSP" vs "HID"/"HIE"/"HIP" residue labels).
// Everything else in the reader is typed; this file and the loader
// are the only places strings turn into enum values or are parsed for
// physics meaning.
//
// Do not call the translator from rendering, overlay, or inspector
// code — use typed properties (AminoAcid, AtomRole, ring subclass
// identity) instead. String output for UI is derived from typed
// properties via TypeName()/NameFor*() helpers in Types.h.

#pragma once

#include "../model/Types.h"

#include <QString>
#include <string>

namespace h5reader::io {

// ---------------------------------------------------------------------------
// Amino-acid code → enum.
//
// Accepts the 20 standard codes, plus force-field protonation variants:
//   HIS:   HID HIE HIP (AMBER), HSD HSE HSP (CHARMM)
//   ASP:   ASH
//   GLU:   GLH
//   CYS:   CYX CYM
//   LYS:   LYN
//   ARG:   ARN
//   TYR:   TYM
//   MET:   MSE (selenomethionine → treated as MET)
//
// Unknown codes → AminoAcid::Unknown, caller should ErrorBus::Report.
// ---------------------------------------------------------------------------
model::AminoAcid ThreeLetterCodeToAminoAcid(const std::string& code);

// Canonical three-letter code for display purposes. Does NOT reproduce
// protonation-variant suffixes (e.g. HIS, not HID) — use
// ProtonationVariantName() for that.
const char* ThreeLetterCodeForAminoAcid(model::AminoAcid aa);

// ---------------------------------------------------------------------------
// Protonation variant helpers.
// ---------------------------------------------------------------------------
const char* ProtonationVariantName(model::ProtonationVariant v);

// ---------------------------------------------------------------------------
// Display-only: translate a CHARMM-flavoured atom name to canonical PDB
// naming. The library's NamingRegistry is the authority for the mapping;
// we carry a subset covering what appears in analysis H5 files from
// CHARMM36m trajectories.
//
// Used by inspector widgets that want to show "H" rather than "HN".
// Never used for dispatch.
// ---------------------------------------------------------------------------
QString CharmmAtomNameToCanonical(const QString& charmmName,
                                   const QString& residueName);

}  // namespace h5reader::io
