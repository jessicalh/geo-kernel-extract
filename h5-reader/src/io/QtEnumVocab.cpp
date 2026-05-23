// QtEnumVocab implementation — manifest enum_vocab parser + reader-
// side enum drift validator.
//
// The reader-side NameForX() helpers in this file produce the
// canonical names the manifest is expected to declare. If a future
// extractor change renames an enum value, ValidateAgainstReaderEnums()
// logs the mismatch at load.

#include "QtEnumVocab.h"

#include "../diagnostics/ErrorBus.h"
#include "../model/QtSemanticEnums.h"  // TerminalState, etc.
#include "../model/Types.h"            // BondOrder, BondCategory, RingKind, RingTypeIndex, AminoAcid

#include <QFile>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonValue>

namespace h5reader::io {

namespace {

// ── Reader-side canonical names per enum (must match manifest) ──
//
// These are the strings the manifest's enum_vocab declares for each
// ordinal. The reader hardcodes them here as the validation target;
// drift between the manifest and these literals trips an ErrorBus
// warning at load.

const char* NameForReaderTerminalState(int ord) {
    switch (static_cast<h5reader::model::TerminalState>(ord)) {
    case h5reader::model::TerminalState::Internal:
        return "Internal";
    case h5reader::model::TerminalState::NTerminus:
        return "NTerminus";
    case h5reader::model::TerminalState::CTerminus:
        return "CTerminus";
    case h5reader::model::TerminalState::NAndCTerminus:
        return "NAndCTerminus";
    case h5reader::model::TerminalState::Unknown:
        return "Unknown";
    }
    return nullptr;
}

const char* NameForReaderBondOrder(int ord) {
    switch (static_cast<h5reader::model::BondOrder>(ord)) {
    case h5reader::model::BondOrder::Single:
        return "Single";
    case h5reader::model::BondOrder::Double:
        return "Double";
    case h5reader::model::BondOrder::Triple:
        return "Triple";
    case h5reader::model::BondOrder::Aromatic:
        return "Aromatic";
    case h5reader::model::BondOrder::Peptide:
        return "Peptide";
    case h5reader::model::BondOrder::Unknown:
        return "Unknown";
    }
    return nullptr;
}

const char* NameForReaderBondCategory(int ord) {
    switch (static_cast<h5reader::model::BondCategory>(ord)) {
    case h5reader::model::BondCategory::PeptideCO:
        return "PeptideCO";
    case h5reader::model::BondCategory::PeptideCN:
        return "PeptideCN";
    case h5reader::model::BondCategory::BackboneOther:
        return "BackboneOther";
    case h5reader::model::BondCategory::SidechainCO:
        return "SidechainCO";
    case h5reader::model::BondCategory::Aromatic:
        return "Aromatic";
    case h5reader::model::BondCategory::Disulfide:
        return "Disulfide";
    case h5reader::model::BondCategory::SidechainOther:
        return "SidechainOther";
    case h5reader::model::BondCategory::Unknown:
        return "Unknown";
    }
    return nullptr;
}

const char* NameForReaderRingKind(int ord) {
    switch (static_cast<h5reader::model::RingKind>(ord)) {
    case h5reader::model::RingKind::Aromatic:
        return "aromatic";
    case h5reader::model::RingKind::Saturated:
        return "saturated";
    }
    return nullptr;
}

const char* NameForReaderRingTypeIndex(int ord) {
    switch (static_cast<h5reader::model::RingTypeIndex>(ord)) {
    case h5reader::model::RingTypeIndex::PheBenzene:
        return "PheBenzene";
    case h5reader::model::RingTypeIndex::TyrPhenol:
        return "TyrPhenol";
    case h5reader::model::RingTypeIndex::TrpBenzene:
        return "TrpBenzene";
    case h5reader::model::RingTypeIndex::TrpPyrrole:
        return "TrpPyrrole";
    case h5reader::model::RingTypeIndex::TrpPerimeter:
        return "TrpPerimeter";
    case h5reader::model::RingTypeIndex::HisImidazole:
        return "HisImidazole";
    case h5reader::model::RingTypeIndex::HidImidazole:
        return "HidImidazole";
    case h5reader::model::RingTypeIndex::HieImidazole:
        return "HieImidazole";
    case h5reader::model::RingTypeIndex::ProPyrrolidine:
        return "ProPyrrolidine";
    }
    return nullptr;
}

const char* NameForReaderAminoAcid(int ord) {
    switch (static_cast<h5reader::model::AminoAcid>(ord)) {
    case h5reader::model::AminoAcid::ALA:
        return "ALA";
    case h5reader::model::AminoAcid::ARG:
        return "ARG";
    case h5reader::model::AminoAcid::ASN:
        return "ASN";
    case h5reader::model::AminoAcid::ASP:
        return "ASP";
    case h5reader::model::AminoAcid::CYS:
        return "CYS";
    case h5reader::model::AminoAcid::GLN:
        return "GLN";
    case h5reader::model::AminoAcid::GLU:
        return "GLU";
    case h5reader::model::AminoAcid::GLY:
        return "GLY";
    case h5reader::model::AminoAcid::HIS:
        return "HIS";
    case h5reader::model::AminoAcid::ILE:
        return "ILE";
    case h5reader::model::AminoAcid::LEU:
        return "LEU";
    case h5reader::model::AminoAcid::LYS:
        return "LYS";
    case h5reader::model::AminoAcid::MET:
        return "MET";
    case h5reader::model::AminoAcid::PHE:
        return "PHE";
    case h5reader::model::AminoAcid::PRO:
        return "PRO";
    case h5reader::model::AminoAcid::SER:
        return "SER";
    case h5reader::model::AminoAcid::THR:
        return "THR";
    case h5reader::model::AminoAcid::TRP:
        return "TRP";
    case h5reader::model::AminoAcid::TYR:
        return "TYR";
    case h5reader::model::AminoAcid::VAL:
        return "VAL";
    case h5reader::model::AminoAcid::Unknown:
        return "Unknown";
    }
    return nullptr;
}

// Dispatch table for the per-kind reader-side NameForX. Indexed by
// QtEnumVocab::Kind ordinal.
using NameFn = const char* (*)(int);
constexpr NameFn kReaderNameFn[QtEnumVocab::kKindCount] = {
    &NameForReaderTerminalState,  // 0
    &NameForReaderBondOrder,      // 1
    &NameForReaderBondCategory,   // 2
    &NameForReaderRingKind,       // 3
    &NameForReaderRingTypeIndex,  // 4
    &NameForReaderAminoAcid,      // 5
};

constexpr const char* kKindManifestKey[QtEnumVocab::kKindCount] = {
    "terminal_state",
    "bond_order",
    "bond_category",
    "ring_kind",
    "ring_type_index",
    "residue_type",
};

constexpr const char* kKindDisplayName[QtEnumVocab::kKindCount] = {
    "TerminalState",
    "BondOrder",
    "BondCategory",
    "RingKind",
    "RingTypeIndex",
    "AminoAcid",
};

}  // namespace


QtEnumVocab QtEnumVocab::LoadFromManifest(const QString& manifest_path) {
    QtEnumVocab v;

    QFile f(manifest_path);
    if (!f.open(QIODevice::ReadOnly)) {
        v.error_ = QStringLiteral("QtEnumVocab: could not open %1").arg(manifest_path);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtEnumVocab"),
                                                v.error_,
                                                manifest_path);
        return v;
    }
    const QByteArray bytes = f.readAll();
    f.close();

    QJsonParseError perr{};
    const QJsonDocument doc = QJsonDocument::fromJson(bytes, &perr);
    if (perr.error != QJsonParseError::NoError || !doc.isObject()) {
        v.error_ = QStringLiteral("QtEnumVocab: JSON parse failed at offset %1: %2").arg(perr.offset).arg(perr.errorString());
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtEnumVocab"),
                                                v.error_,
                                                manifest_path);
        return v;
    }
    const QJsonObject root = doc.object();
    if (!root.contains(QStringLiteral("enum_vocab"))) {
        // Vocab block is optional in principle; reader-side enums still
        // work without it (we trust the int8 ordinals blindly), but log
        // Warn so this is visible.
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                QStringLiteral("QtEnumVocab"),
                                                QStringLiteral("manifest lacks enum_vocab block; reader-side "
                                                               "ordinals will be trusted blindly"),
                                                manifest_path);
        v.ok_ = true;
        return v;
    }
    const QJsonObject vocab = root.value(QStringLiteral("enum_vocab")).toObject();

    for (int ki = 0; ki < kKindCount; ++ki) {
        const QString key = QString::fromLatin1(kKindManifestKey[ki]);
        if (!vocab.contains(key))
            continue;  // missing kind is non-fatal
        const QJsonObject kindObj = vocab.value(key).toObject();
        for (auto it = kindObj.constBegin(); it != kindObj.constEnd(); ++it) {
            bool ok = false;
            const int ord = it.key().toInt(&ok);
            if (!ok || !it.value().isString())
                continue;
            v.entries_[ki][ord] = it.value().toString();
        }
    }
    v.ok_ = true;
    return v;
}


bool QtEnumVocab::ValidateAgainstReaderEnums() const {
    bool all_ok = true;
    for (int ki = 0; ki < kKindCount; ++ki) {
        const auto& table = entries_[ki];
        const NameFn fn = kReaderNameFn[ki];
        for (const auto& [ord, name] : table) {
            const char* expected = fn(ord);
            if (!expected) {
                h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                        QStringLiteral("QtEnumVocab"),
                                                        QStringLiteral("manifest %1[%2]=\"%3\" but reader's "
                                                                       "typed enum has no value at this ordinal")
                                                            .arg(QString::fromLatin1(kKindDisplayName[ki]))
                                                            .arg(ord)
                                                            .arg(name),
                                                        QString());
                all_ok = false;
                continue;
            }
            if (QString::fromLatin1(expected) != name) {
                h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                        QStringLiteral("QtEnumVocab"),
                                                        QStringLiteral("manifest %1[%2]=\"%3\" but reader expects \"%4\"")
                                                            .arg(QString::fromLatin1(kKindDisplayName[ki]))
                                                            .arg(ord)
                                                            .arg(name)
                                                            .arg(QString::fromLatin1(expected)),
                                                        QString());
                all_ok = false;
            }
        }
    }
    return all_ok;
}


QString QtEnumVocab::DisplayName(Kind kind, int ordinal) const {
    const int ki = static_cast<int>(kind);
    if (ki < 0 || ki >= kKindCount)
        return QString();
    const auto& table = entries_[ki];
    auto it = table.find(ordinal);
    return (it == table.end()) ? QString() : it->second;
}

}  // namespace h5reader::io
