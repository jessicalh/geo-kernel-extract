// QtTopologySidecar implementation.
//
// Sequence:
//   1. Load extraction_manifest.json → QtManifest + QtEnumVocab
//   2. QtEnumVocab.ValidateAgainstReaderEnums() (drift detection)
//   3. Load atoms_category_info.npy → decode rows into QtAtom + QtAtomNames
//   4. Load residues.npy → decode rows into QtResidue (no atom membership yet)
//   5. Load bonds.npy → decode rows into QtBond
//   6. Load rings.npy → CreateQtRing per row, populate identity fields
//   7. Load ring_membership.npy → decode rows into QtRingMembership
//   8. Bucket atom indices into QtResidue.atomIndices via residue_index walk
//   9. Build per-residue backbone cache from atoms' typed BackboneRole + Locant
//  10. Populate each QtRing's atomIndices from ring_membership join
//      (canonical-walk order via ring_atom_order)
//  11. Cross-check axis sizes against manifest.axisSizes
//
// Failure modes log via ErrorBus. Hard fails: missing manifest,
// missing atoms_category_info.npy, missing residues.npy. Soft fails:
// missing bonds.npy / rings.npy / ring_membership.npy (the reader can
// function with reduced topology — just no bond / ring visualizations).

#include "QtTopologySidecar.h"

#include "QtNpyReader.h"
#include "QtNpyRecords.h"

#include "../diagnostics/ErrorBus.h"
#include "../model/QtAminoAcidType.h"

#include <QDir>
#include <QString>
#include <algorithm>
#include <cstring>
#include <string_view>

namespace h5reader::io {

namespace {

using namespace h5reader::model;

// ── Fixed-width string extractor ────────────────────────────────────
//
// Sidecar S<N> columns are zero-padded byte arrays of fixed width;
// strlen-equivalent stops at first NUL or column edge.
QString PackedStringToQString(const char* data, std::size_t maxLen) {
    std::size_t len = 0;
    while (len < maxLen && data[len] != '\0')
        ++len;
    return QString::fromUtf8(data, static_cast<int>(len));
}

std::string_view PackedStringView(const char* data, std::size_t maxLen) {
    std::size_t len = 0;
    while (len < maxLen && data[len] != '\0')
        ++len;
    return std::string_view(data, len);
}

// ── QtFfAtomType string → enum ─────────────────────────────────────
//
// AMBER ff14SB atom-type vocabulary. The S4 column carries values like
// "CT", "HC", "N", "N3", "OW", etc. Maps to typed enum at the loader
// boundary so runtime stays string-free per design §11.B.
//
// Unknown ff types fall to QtFfAtomType::Unknown; the loader logs a
// Warn. Adding a new ff14SB / ff19SB type means adding an enum value
// here + an enum case to the switch.
QtFfAtomType ParseFfAtomType(std::string_view sv) {
    // Carbons
    if (sv == "C")
        return QtFfAtomType::C;
    if (sv == "CA")
        return QtFfAtomType::CA;
    if (sv == "CB")
        return QtFfAtomType::CB;
    if (sv == "CC")
        return QtFfAtomType::CC;
    if (sv == "CD")
        return QtFfAtomType::CD;
    if (sv == "CK")
        return QtFfAtomType::CK;
    if (sv == "CM")
        return QtFfAtomType::CM;
    if (sv == "CN")
        return QtFfAtomType::CN;
    if (sv == "CO")
        return QtFfAtomType::CO;
    if (sv == "CP")
        return QtFfAtomType::CP;
    if (sv == "CQ")
        return QtFfAtomType::CQ;
    if (sv == "CR")
        return QtFfAtomType::CR;
    if (sv == "CT")
        return QtFfAtomType::CT;
    if (sv == "CV")
        return QtFfAtomType::CV;
    if (sv == "CW")
        return QtFfAtomType::CW;
    if (sv == "CX")
        return QtFfAtomType::CX;
    if (sv == "CY")
        return QtFfAtomType::CY;
    if (sv == "CZ")
        return QtFfAtomType::CZ;
    if (sv == "C*")
        return QtFfAtomType::Cstar;
    // Nitrogens
    if (sv == "N")
        return QtFfAtomType::N;
    if (sv == "N2")
        return QtFfAtomType::N2;
    if (sv == "N3")
        return QtFfAtomType::N3;
    if (sv == "NA")
        return QtFfAtomType::NA;
    if (sv == "NB")
        return QtFfAtomType::NB;
    if (sv == "NC")
        return QtFfAtomType::NC;
    if (sv == "NP")
        return QtFfAtomType::NP;
    if (sv == "NT")
        return QtFfAtomType::NT;
    if (sv == "NY")
        return QtFfAtomType::NY;
    if (sv == "N*")
        return QtFfAtomType::Nstar;
    // Hydrogens
    if (sv == "H")
        return QtFfAtomType::H;
    if (sv == "H1")
        return QtFfAtomType::H1;
    if (sv == "H2")
        return QtFfAtomType::H2;
    if (sv == "H3")
        return QtFfAtomType::H3;
    if (sv == "H4")
        return QtFfAtomType::H4;
    if (sv == "H5")
        return QtFfAtomType::H5;
    if (sv == "HA")
        return QtFfAtomType::HA;
    if (sv == "HC")
        return QtFfAtomType::HC;
    if (sv == "HO")
        return QtFfAtomType::HO;
    if (sv == "HP")
        return QtFfAtomType::HP;
    if (sv == "HS")
        return QtFfAtomType::HS;
    if (sv == "HW")
        return QtFfAtomType::HW;
    if (sv == "HZ")
        return QtFfAtomType::HZ;
    // Oxygens
    if (sv == "O")
        return QtFfAtomType::O;
    if (sv == "O2")
        return QtFfAtomType::O2;
    if (sv == "OH")
        return QtFfAtomType::OH;
    if (sv == "OS")
        return QtFfAtomType::OS;
    if (sv == "OW")
        return QtFfAtomType::OW;
    if (sv == "OP")
        return QtFfAtomType::OP;
    if (sv == "OD")
        return QtFfAtomType::OD;
    // Sulfurs
    if (sv == "S")
        return QtFfAtomType::S;
    if (sv == "SH")
        return QtFfAtomType::SH;
    return QtFfAtomType::Unknown;
}

// ── Decoders: NPY row → typed Qt model ─────────────────────────────

QtAtom DecodeAtom(const QtNpyAtomCategoryRow& row) {
    QtAtom a;
    a.atomIndex = row.atom_index;
    a.residueIndex = row.residue_index;
    a.parentAtomIndex = row.parent_atom_index;

    a.element = ElementFromAtomicNumber(static_cast<int>(row.element));
    a.locant = static_cast<Locant>(row.locant);
    a.branch = {row.branch_outer, row.branch_inner};
    a.diIndex = static_cast<DiastereotopicIndex>(row.di_index);
    a.backboneRole = static_cast<BackboneRole>(row.backbone_role);

    a.prochiral = static_cast<ProchiralStereo>(row.prochiral);
    a.planarGroup = static_cast<PlanarGroupKind>(row.planar_group);
    a.planarStereo = static_cast<PlanarStereo>(row.planar_stereo);
    a.polarH = static_cast<PolarHKind>(row.polar_h_kind);

    a.ringPositionPrimary = static_cast<RingPositionLabel>(row.ring_position_primary);
    a.ringPositionSecondary = static_cast<RingPositionLabel>(row.ring_position_secondary);

    a.pseudoatomKind = static_cast<PseudoatomKind>(row.pseudoatom_kind);
    a.inSuperGroup = (row.in_super_group != 0);
    a.aromatic = (row.aromatic != 0);
    a.formalCharge = row.formal_charge;
    a.isExchangeable = (row.is_exchangeable != 0);
    a.equivalenceClass = row.equivalence_class;

    a.ffAtomType = ParseFfAtomType(PackedStringView(row.ff_atom_type_string, 4));

    return a;
}

// Note: atomIndex parameter for diagnostics on ParseFfAtomType warns is
// passed by the caller, not stored on QtAtom (it's already on QtAtom.atomIndex).

QtAtomNames DecodeAtomNames(const QtNpyAtomCategoryRow& row) {
    QtAtomNames n;
    n.amber = PackedStringToQString(row.amber_atom_name, 8);
    n.iupac = PackedStringToQString(row.iupac_atom_name, 8);
    n.bmrb = PackedStringToQString(row.bmrb_atom_name, 8);
    n.iupacProvenance = static_cast<NamingProvenance>(row.iupac_naming_provenance);
    n.bmrbProvenance = static_cast<NamingProvenance>(row.bmrb_naming_provenance);
    return n;
}

QtResidue DecodeResidue(const QtNpyResidueRow& row) {
    QtResidue r;
    r.residueIndex = row.residue_index;

    r.address.chainId = PackedStringToQString(row.chain_id, 2);
    r.address.residueNumber = row.residue_number;
    r.address.insertionCode = PackedStringToQString(row.insertion_code, 1);

    r.aminoAcid = static_cast<AminoAcid>(row.residue_type);
    r.protonationVariantIndex = row.protonation_variant_index;
    r.terminalState = static_cast<TerminalState>(row.terminal_state);

    r.prevResidueIndex = row.prev_residue_index;
    r.nextResidueIndex = row.next_residue_index;
    r.prevResidueType = static_cast<AminoAcid>(row.prev_residue_type);
    r.nextResidueType = static_cast<AminoAcid>(row.next_residue_type);

    r.atomCount = row.atom_count;
    r.isProline = (row.is_proline != 0);
    r.isAromatic = (row.is_aromatic != 0);
    r.isTitratable = (row.is_titratable != 0);
    r.hasAmideH = (row.has_amide_h != 0);
    r.isXProContext = (row.is_xpro_context != 0);

    return r;
}

QtBond DecodeBond(const QtNpyBondRow& row) {
    QtBond b;
    b.bondIndex = row.bond_index;
    b.atomIndexA = row.atom_index_a;
    b.atomIndexB = row.atom_index_b;
    b.order = static_cast<BondOrder>(row.bond_order);
    b.category = static_cast<BondCategory>(row.bond_category);
    b.isRotatable = (row.is_rotatable != 0);
    b.isAromatic = (row.is_aromatic != 0);
    b.isPeptide = (row.is_peptide != 0);
    b.isBackbone = (row.is_backbone != 0);
    return b;
}

std::unique_ptr<QtRing> DecodeRing(const QtNpyRingRow& row, bool& warned_unknown_type) {
    const auto typeIdx = static_cast<RingTypeIndex>(row.ring_type_index);
    auto ring = CreateQtRing(typeIdx);
    if (!ring) {
        if (!warned_unknown_type) {
            h5reader::diagnostics::ErrorBus::Report(
                h5reader::diagnostics::Severity::Warning,
                QStringLiteral("QtTopologySidecar"),
                QStringLiteral("rings.npy: unknown ring_type_index=%1 at ring_id=%2; ring skipped")
                    .arg(static_cast<int>(row.ring_type_index))
                    .arg(row.ring_id),
                QString());
            warned_unknown_type = true;
        }
        return nullptr;
    }
    ring->ringId = row.ring_id;
    ring->nativeAxisIndex = row.native_axis_index;
    ring->ringKind = static_cast<RingKind>(row.ring_kind);
    ring->parentResidueIndex = row.parent_residue_index;
    ring->parentResidueNumber = row.parent_residue_number;
    ring->fusedPartnerRingId = row.fused_partner_ring_id;
    return ring;
}

QtRingMembership DecodeRingMembership(const QtNpyRingMembershipRow& row) {
    QtRingMembership m;
    m.ringId = row.ring_id;
    m.atomIndex = row.atom_index;
    m.ringAtomOrder = row.ring_atom_order;
    m.isVertex = (row.is_vertex != 0);
    m.isSubstituent = (row.is_substituent != 0);
    return m;
}

// ── Backbone cache builder ─────────────────────────────────────────
//
// For each residue, walk its atoms and assign backbone slots from the
// typed BackboneRole + Locant (no string scan).
//
// Special case: GLY's HA2/HA3 carry BackboneRole::None + Locant::Alpha
// (Markley convention). We pick HA == DiastereotopicIndex::Position2
// for GLY by convention so the reader has a deterministic HA slot.
void PopulateBackboneCache(QtResidue& res, const std::vector<QtAtom>& atoms) {
    const bool is_gly = (res.aminoAcid == AminoAcid::GLY);

    for (int32_t a_idx : res.atomIndices) {
        if (a_idx < 0 || static_cast<std::size_t>(a_idx) >= atoms.size())
            continue;
        const QtAtom& a = atoms[static_cast<std::size_t>(a_idx)];

        switch (a.backboneRole) {
        case BackboneRole::Nitrogen:
            res.N = a_idx;
            break;
        case BackboneRole::AlphaCarbon:
            res.CA = a_idx;
            break;
        case BackboneRole::CarbonylCarbon:
            res.C = a_idx;
            break;
        case BackboneRole::CarbonylOxygen:
            res.O = a_idx;
            break;
        case BackboneRole::AmideHydrogen:
            res.H = a_idx;
            break;
        case BackboneRole::AlphaHydrogen:
            res.HA = a_idx;
            break;
        case BackboneRole::None:
            break;
        }

        // GLY HA: BackboneRole::None + Locant::Alpha + DI::Position2.
        if (is_gly && a.element == Element::H && a.locant == Locant::Alpha && a.backboneRole == BackboneRole::None
            && a.diIndex == DiastereotopicIndex::Position2 && res.HA == QtResidue::NONE) {
            res.HA = a_idx;
        }

        // CB: substrate-typed (Element::C, Locant::Beta, BackboneRole::None).
        // No "first non-aromatic sidechain C" heuristic — the substrate
        // says which atom is at beta-position outright.
        if (a.element == Element::C && a.locant == Locant::Beta && a.backboneRole == BackboneRole::None
            && res.CB == QtResidue::NONE) {
            res.CB = a_idx;
        }
    }
}

}  // namespace


QtTopologySidecar::LoadResult QtTopologySidecar::Load(const QString& sidecar_dir,
                                                       const QString& manifest_path) {
    LoadResult r;
    QDir d(sidecar_dir);
    if (!d.exists()) {
        r.error = QStringLiteral("QtTopologySidecar: sidecar directory does not exist: %1").arg(sidecar_dir);
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Error,
                                                QStringLiteral("QtTopologySidecar"),
                                                r.error,
                                                sidecar_dir);
        return r;
    }

    // The `.LGS` (CalcsetManifest) supplies the extraction_manifest.json
    // path explicitly; the legacy sidecar/extraction_manifest.json
    // convention is used when no override is passed.
    const QString manifestPath = manifest_path.isEmpty()
        ? d.filePath(QStringLiteral("extraction_manifest.json"))
        : manifest_path;
    const QString atomsPath = d.filePath(QStringLiteral("atoms_category_info.npy"));
    const QString residuesPath = d.filePath(QStringLiteral("residues.npy"));
    const QString bondsPath = d.filePath(QStringLiteral("bonds.npy"));
    const QString ringsPath = d.filePath(QStringLiteral("rings.npy"));
    const QString membershipPath = d.filePath(QStringLiteral("ring_membership.npy"));

    // 1. Manifest + EnumVocab
    r.manifest = QtManifest::Load(manifestPath);
    if (!r.manifest.ok) {
        r.error = r.manifest.error;
        return r;
    }
    r.vocab = QtEnumVocab::LoadFromManifest(manifestPath);
    if (r.vocab.ok()) {
        if (!r.vocab.ValidateAgainstReaderEnums()) {
            // ValidateAgainstReaderEnums already logs per-mismatch; bump
            // the warning count.
            ++r.warningCount;
        }
    }

    // 2. Atoms (required)
    std::vector<QtNpyAtomCategoryRow> atomRows;
    {
        auto res = QtNpyReader::ReadStructured<QtNpyAtomCategoryRow>(atomsPath, "atom_index", atomRows);
        if (!res.ok) {
            r.error = res.error;
            return r;
        }
    }
    if (r.manifest.axisSizes.atom != 0 && r.manifest.axisSizes.atom != atomRows.size()) {
        h5reader::diagnostics::ErrorBus::Report(
            h5reader::diagnostics::Severity::Warning,
            QStringLiteral("QtTopologySidecar"),
            QStringLiteral("atoms_category_info row count %1 disagrees with manifest.axisSizes.atom %2")
                .arg(atomRows.size())
                .arg(r.manifest.axisSizes.atom),
            atomsPath);
        ++r.warningCount;
    }
    r.atoms.reserve(atomRows.size());
    r.atomNames.reserve(atomRows.size());
    for (const auto& row : atomRows) {
        r.atoms.push_back(DecodeAtom(row));
        r.atomNames.push_back(DecodeAtomNames(row));
    }

    // 3. Residues (required)
    std::vector<QtNpyResidueRow> residueRows;
    {
        auto res = QtNpyReader::ReadStructured<QtNpyResidueRow>(residuesPath, "residue_index", residueRows);
        if (!res.ok) {
            r.error = res.error;
            return r;
        }
    }
    if (r.manifest.axisSizes.residue != 0 && r.manifest.axisSizes.residue != residueRows.size()) {
        h5reader::diagnostics::ErrorBus::Report(
            h5reader::diagnostics::Severity::Warning,
            QStringLiteral("QtTopologySidecar"),
            QStringLiteral("residues row count %1 disagrees with manifest.axisSizes.residue %2")
                .arg(residueRows.size())
                .arg(r.manifest.axisSizes.residue),
            residuesPath);
        ++r.warningCount;
    }
    r.residues.reserve(residueRows.size());
    for (const auto& row : residueRows) {
        r.residues.push_back(DecodeResidue(row));
    }

    // 3b. Verbatim residue-name projection (#7). amber/iupac/one_letter come
    //     from residues.npy; bmrb_residue_3letter lives only per-atom in
    //     atoms_category_info, so take it from any atom of the residue (it is
    //     a residue property replicated on each of its atoms). Labels are
    //     projections, never identity — see QtResidueNames.h.
    std::vector<QString> residueBmrb(residueRows.size());
    for (const auto& arow : atomRows) {
        if (arow.residue_index >= 0
            && static_cast<std::size_t>(arow.residue_index) < residueBmrb.size()
            && residueBmrb[static_cast<std::size_t>(arow.residue_index)].isEmpty()) {
            residueBmrb[static_cast<std::size_t>(arow.residue_index)] =
                PackedStringToQString(arow.bmrb_residue_3letter, 4);
        }
    }
    r.residueNames.reserve(residueRows.size());
    for (std::size_t i = 0; i < residueRows.size(); ++i) {
        const auto& row = residueRows[i];
        h5reader::model::QtResidueNames n;
        n.amber = PackedStringToQString(row.amber_residue_3letter, 4);
        n.iupac = PackedStringToQString(row.iupac_residue_3letter, 4);
        n.bmrb = residueBmrb[i].isEmpty() ? n.amber : residueBmrb[i];  // amber fallback
        n.oneLetter = (row.one_letter[0] != '\0') ? row.one_letter[0] : '?';
        r.residueNames.push_back(std::move(n));
    }

    // 4. Bonds (soft fail — reader can function without)
    std::vector<QtNpyBondRow> bondRows;
    if (QFile::exists(bondsPath)) {
        auto res = QtNpyReader::ReadStructured<QtNpyBondRow>(bondsPath, "bond_index", bondRows);
        if (res.ok) {
            r.bonds.reserve(bondRows.size());
            for (const auto& row : bondRows)
                r.bonds.push_back(DecodeBond(row));
        } else {
            ++r.warningCount;
        }
    } else {
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                QStringLiteral("QtTopologySidecar"),
                                                QStringLiteral("bonds.npy missing; bond overlays will be empty"),
                                                bondsPath);
        ++r.warningCount;
    }

    // 5. Rings + ring_membership (soft fail)
    std::vector<QtNpyRingRow> ringRows;
    if (QFile::exists(ringsPath)) {
        auto res = QtNpyReader::ReadStructured<QtNpyRingRow>(ringsPath, "ring_id", ringRows);
        if (res.ok) {
            bool warned_unknown = false;
            r.rings.reserve(ringRows.size());
            for (const auto& row : ringRows) {
                auto ring = DecodeRing(row, warned_unknown);
                if (ring) {
                    if (ring->ringKind == RingKind::Aromatic)
                        ++r.aromaticRingCount;
                    else
                        ++r.saturatedRingCount;
                }
                r.rings.push_back(std::move(ring));
            }
        } else {
            ++r.warningCount;
        }
    } else {
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                QStringLiteral("QtTopologySidecar"),
                                                QStringLiteral("rings.npy missing; ring overlays will be empty"),
                                                ringsPath);
        ++r.warningCount;
    }

    std::vector<QtNpyRingMembershipRow> membershipRows;
    if (QFile::exists(membershipPath)) {
        auto res = QtNpyReader::ReadStructured<QtNpyRingMembershipRow>(membershipPath, "ring_id", membershipRows);
        if (res.ok) {
            r.ringMemberships.reserve(membershipRows.size());
            for (const auto& row : membershipRows)
                r.ringMemberships.push_back(DecodeRingMembership(row));
        } else {
            ++r.warningCount;
        }
    } else if (!r.rings.empty()) {
        h5reader::diagnostics::ErrorBus::Report(h5reader::diagnostics::Severity::Warning,
                                                QStringLiteral("QtTopologySidecar"),
                                                QStringLiteral("ring_membership.npy missing despite rings.npy present; "
                                                               "ring vertices unrecoverable"),
                                                membershipPath);
        ++r.warningCount;
    }

    // 8. Bucket atom indices into QtResidue.atomIndices
    for (const auto& a : r.atoms) {
        if (a.residueIndex >= 0 && static_cast<std::size_t>(a.residueIndex) < r.residues.size()) {
            r.residues[static_cast<std::size_t>(a.residueIndex)].atomIndices.push_back(a.atomIndex);
        }
    }

    // 9. Populate backbone cache per residue from typed BackboneRole + Locant
    for (auto& res : r.residues) {
        PopulateBackboneCache(res, r.atoms);
    }

    // 10. Populate each ring's atomIndices from ring_membership join
    //     (canonical-walk order via ring_atom_order).
    if (!r.rings.empty() && !r.ringMemberships.empty()) {
        // Bucket memberships by ring_id, sorted by ring_atom_order.
        std::vector<std::vector<int32_t>> byRing(r.rings.size());
        for (std::size_t mi = 0; mi < r.ringMemberships.size(); ++mi) {
            const auto& m = r.ringMemberships[mi];
            if (m.ringId >= 0 && static_cast<std::size_t>(m.ringId) < r.rings.size())
                byRing[static_cast<std::size_t>(m.ringId)].push_back(static_cast<int32_t>(mi));
        }
        for (std::size_t ri = 0; ri < r.rings.size(); ++ri) {
            if (!r.rings[ri])
                continue;
            auto& bucket = byRing[ri];
            std::sort(bucket.begin(), bucket.end(), [&r](int32_t a, int32_t b) {
                return r.ringMemberships[static_cast<std::size_t>(a)].ringAtomOrder
                       < r.ringMemberships[static_cast<std::size_t>(b)].ringAtomOrder;
            });
            r.rings[ri]->atomIndices.clear();
            r.rings[ri]->atomIndices.reserve(bucket.size());
            for (int32_t mi : bucket) {
                r.rings[ri]->atomIndices.push_back(r.ringMemberships[static_cast<std::size_t>(mi)].atomIndex);
            }
        }
    }

    r.ok = true;
    return r;
}

}  // namespace h5reader::io
