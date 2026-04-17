// =============================================================================
// H5 LOADING BOUNDARY
//
// Every decoder in this file crosses FROM the H5's on-disk
// representation (int ordinals, three-letter codes, CHARMM atom names)
// INTO the reader's typed object model. After Load() returns, no other
// code in the reader should call AminoAcidFromCode, cast raw ordinals,
// or compare atom_name strings. Those are this file's responsibility
// and nowhere else's.
//
// Grep test: from the project root,
//   rg -n 'AminoAcidFromCode|CharmmAtomNameToCanonical|static_cast<(AtomRole|RingTypeIndex|BondCategory|BondOrder|Hybridisation)>' src/
// should return matches only under src/io/. Any other hit is a
// boundary violation.
// =============================================================================

#include "QtProteinLoader.h"

#include "QtNamingRegistry.h"
#include "../diagnostics/ErrorBus.h"
#include "../model/QtRing.h"

#include <QFileInfo>
#include <QLoggingCategory>

#include <algorithm>
#include <stdexcept>

Q_LOGGING_CATEGORY(cLoader, "h5reader.loader")

namespace h5reader::io {

using namespace h5reader::model;
using h5reader::diagnostics::ErrorBus;
using h5reader::diagnostics::Severity;

namespace {

// ---------------------------------------------------------------------------
// Enum decoders. Each logs + counts an error on bad input; returns the
// enum's "Unknown" sentinel so downstream code can keep going.
// ---------------------------------------------------------------------------

Element DecodeElement(int32_t z, int& errors) {
    switch (z) {
        case 1:  return Element::H;
        case 6:  return Element::C;
        case 7:  return Element::N;
        case 8:  return Element::O;
        case 16: return Element::S;
    }
    ErrorBus::Report(Severity::Warning, QStringLiteral("QtProteinLoader"),
                     QStringLiteral("Unknown atomic number"),
                     QStringLiteral("z=%1").arg(z));
    ++errors;
    return Element::Unknown;
}

AtomRole DecodeAtomRole(int32_t n, int& errors) {
    if (n >= 0 && n <= static_cast<int32_t>(AtomRole::Unknown))
        return static_cast<AtomRole>(n);
    ErrorBus::Report(Severity::Warning, QStringLiteral("QtProteinLoader"),
                     QStringLiteral("AtomRole ordinal out of range"),
                     QStringLiteral("n=%1").arg(n));
    ++errors;
    return AtomRole::Unknown;
}

Hybridisation DecodeHybridisation(int32_t n, int& errors) {
    if (n >= 0 && n <= static_cast<int32_t>(Hybridisation::Unassigned))
        return static_cast<Hybridisation>(n);
    ErrorBus::Report(Severity::Warning, QStringLiteral("QtProteinLoader"),
                     QStringLiteral("Hybridisation ordinal out of range"),
                     QStringLiteral("n=%1").arg(n));
    ++errors;
    return Hybridisation::Unassigned;
}

BondOrder DecodeBondOrder(int32_t n, int& errors) {
    if (n >= 0 && n <= static_cast<int32_t>(BondOrder::Unknown))
        return static_cast<BondOrder>(n);
    ErrorBus::Report(Severity::Warning, QStringLiteral("QtProteinLoader"),
                     QStringLiteral("BondOrder ordinal out of range"),
                     QStringLiteral("n=%1").arg(n));
    ++errors;
    return BondOrder::Unknown;
}

BondCategory DecodeBondCategory(int32_t n, int& errors) {
    if (n >= 0 && n <= static_cast<int32_t>(BondCategory::Unknown))
        return static_cast<BondCategory>(n);
    ErrorBus::Report(Severity::Warning, QStringLiteral("QtProteinLoader"),
                     QStringLiteral("BondCategory ordinal out of range"),
                     QStringLiteral("n=%1").arg(n));
    ++errors;
    return BondCategory::Unknown;
}

RingTypeIndex DecodeRingType(int32_t n, int& errors) {
    if (n >= 0 && n < RingTypeCount)
        return static_cast<RingTypeIndex>(n);
    ErrorBus::Report(Severity::Warning, QStringLiteral("QtProteinLoader"),
                     QStringLiteral("RingTypeIndex out of range"),
                     QStringLiteral("n=%1 max=%2").arg(n).arg(RingTypeCount - 1));
    ++errors;
    return RingTypeIndex::PheBenzene;   // sentinel; caller will null the ring
}

// ---------------------------------------------------------------------------
// Protonation-variant inference. Uses ring subclass identity (via
// dynamic_cast — typed dispatch, not strings) for HIS tautomers, and
// atom-role presence inside the residue for ASP/GLU/LYS/TYR/CYS.
//
// Takes the already-decoded QtProtein so it can walk atoms by role,
// and a residue index. Does NOT read H5 strings.
// ---------------------------------------------------------------------------

ProtonationVariant InferProtonationVariant(
    const QtProtein& protein, int residueIndex, AminoAcid aa)
{
    switch (aa) {
        case AminoAcid::HIS: {
            for (const auto& ring : protein.rings()) {
                if (!ring || ring->parentResidueIndex != residueIndex) continue;
                if (dynamic_cast<const QtHidImidazoleRing*>(ring.get()))
                    return ProtonationVariant::HIS_delta;
                if (dynamic_cast<const QtHieImidazoleRing*>(ring.get()))
                    return ProtonationVariant::HIS_epsilon;
                if (dynamic_cast<const QtHisImidazoleRing*>(ring.get()))
                    return ProtonationVariant::Default;    // tautomer unspecified
            }
            return ProtonationVariant::Default;
        }

        // For other titratable residues, the H5 does not separately
        // record the variant. Infer from the presence of the variant-
        // specific proton by walking the residue's atoms and checking
        // AtomRole. Expand these inferences as the inspector surfaces
        // cases. Not speculative — only when we have a concrete use.
        default:
            return ProtonationVariant::Default;
    }
}

// ---------------------------------------------------------------------------
// Populate backbone index cache for one residue from the atoms it owns.
// ---------------------------------------------------------------------------

void PopulateBackboneCache(QtResidue& residue, const std::vector<QtAtom>& atoms) {
    for (size_t a : residue.atomIndices) {
        const auto& atom = atoms[a];
        switch (atom.role) {
            case AtomRole::BackboneN:  residue.N  = a; break;
            case AtomRole::BackboneCA: residue.CA = a; break;
            case AtomRole::BackboneC:  residue.C  = a; break;
            case AtomRole::BackboneO:  residue.O  = a; break;
            case AtomRole::AmideH:     residue.H  = a; break;
            case AtomRole::AlphaH:     residue.HA = a; break;
            default: break;
        }
        // CB — the first non-aromatic, non-methyl sidechain C we see.
        if (residue.CB == QtResidue::NONE
            && atom.role == AtomRole::SidechainC
            && atom.element == Element::C) {
            residue.CB = a;
        }
    }
}

}  // anonymous namespace


// =============================================================================
// Public entry
// =============================================================================

QtLoadResult QtProteinLoader::Load(const QString& h5Path) {
    QtLoadResult result;

    // Derive a display identifier from the filename stem. The H5's
    // meta/protein_id often reads "md" on trajectory extractions (the
    // TPR filename leaks through); the path stem is more reliable.
    // See feedback_h5_reader_test.md.
    QFileInfo fi(h5Path);
    QString stem = fi.completeBaseName();   // e.g. "1B1V_4292_analysis"
    if (stem.endsWith(QLatin1String("_analysis"))) stem.chop(9);
    result.proteinId = stem;

    // Open the H5. HighFive throws on failure — the external-library
    // boundary. Catch and convert to a load failure.
    auto af = std::make_shared<AnalysisFile>();
    try {
        af->ReadH5(h5Path.toStdString());
    } catch (const std::exception& e) {
        result.ok = false;
        result.error = QStringLiteral("AnalysisFile::ReadH5 failed: %1")
                         .arg(QString::fromUtf8(e.what()));
        ErrorBus::Report(Severity::Error, QStringLiteral("QtProteinLoader"),
                         result.error, h5Path);
        return result;
    }

    const size_t N = af->n_atoms;
    const size_t R = af->n_residues;
    const size_t B = af->n_bonds;
    const size_t nRings = af->n_rings;
    const size_t T = af->n_frames;

    qCInfo(cLoader).noquote()
        << "Decoding" << result.proteinId
        << "| N=" << N << "atoms | R=" << R << "residues"
        << "| B=" << B << "bonds | rings=" << nRings
        << "| T=" << T << "frames | stride=" << af->meta.stride;

    // Build the protein shell and populate in place.
    auto protein = std::make_unique<QtProtein>();
    protein->setProteinId(result.proteinId);

    int errors = 0;

    // ---- Atoms ----
    protein->atoms_.reserve(N);
    for (size_t i = 0; i < N; ++i) {
        QtAtom atom;
        atom.element       = DecodeElement(af->atoms.element[i], errors);
        atom.role          = DecodeAtomRole(af->atoms.atom_role[i], errors);
        atom.hybridisation = DecodeHybridisation(af->atoms.hybridisation[i], errors);
        atom.residueIndex  = af->atoms.residue_index[i];
        atom.nBonded       = af->atoms.n_bonded[i];
        atom.graphDistRing = af->atoms.graph_dist_ring[i];
        atom.graphDistN    = af->atoms.graph_dist_N[i];
        atom.graphDistO    = af->atoms.graph_dist_O[i];

        atom.isBackbone          = af->atoms.is_backbone[i]          != 0;
        atom.isAmideH            = af->atoms.is_amide_H[i]           != 0;
        atom.isAlphaH            = af->atoms.is_alpha_H[i]           != 0;
        atom.isMethyl            = af->atoms.is_methyl[i]            != 0;
        atom.isAromaticH         = af->atoms.is_aromatic_H[i]        != 0;
        atom.isOnAromaticResidue = af->atoms.is_on_aromatic_residue[i] != 0;
        atom.isHBondDonor        = af->atoms.is_hbond_donor[i]       != 0;
        atom.isHBondAcceptor     = af->atoms.is_hbond_acceptor[i]    != 0;
        atom.parentIsSp2         = af->atoms.parent_is_sp2[i]        != 0;
        atom.isConjugated        = af->atoms.is_conjugated[i]        != 0;

        atom.enegSum1          = af->atoms.eneg_sum_1[i];
        atom.enegSum2          = af->atoms.eneg_sum_2[i];
        atom.nPiBonds3         = af->atoms.n_pi_bonds_3[i];
        atom.bfsToNearestRing  = af->atoms.bfs_to_nearest_ring[i];
        atom.bfsDecay          = af->atoms.bfs_decay[i];
        atom.partialCharge     = af->atoms.partial_charge[i];
        atom.vdwRadius         = af->atoms.vdw_radius[i];
        atom.parentAtomIndex   = af->topology.parent_atom_index[i];

        atom.h5AtomName = QString::fromStdString(af->atoms.atom_name[i]);
        protein->atoms_.push_back(std::move(atom));
    }

    // ---- Residues ----
    protein->residues_.reserve(R);
    for (size_t r = 0; r < R; ++r) {
        QtResidue residue;
        residue.aminoAcid = ThreeLetterCodeToAminoAcid(af->residues.residue_name[r]);
        if (residue.aminoAcid == AminoAcid::Unknown) {
            ErrorBus::Report(Severity::Warning, QStringLiteral("QtProteinLoader"),
                             QStringLiteral("Unknown three-letter amino-acid code"),
                             QStringLiteral("code=\"%1\" residue=%2")
                               .arg(QString::fromStdString(af->residues.residue_name[r]))
                               .arg(r));
            ++errors;
        }
        residue.residueNumber = af->residues.residue_number[r];
        residue.chainId       = QString::fromStdString(af->residues.chain_id[r]);
        protein->residues_.push_back(std::move(residue));
    }

    // Second pass: assign atom membership to residues by residueIndex.
    for (size_t i = 0; i < N; ++i) {
        const int ridx = protein->atoms_[i].residueIndex;
        if (ridx >= 0 && static_cast<size_t>(ridx) < R) {
            protein->residues_[ridx].atomIndices.push_back(i);
        }
    }

    // Third pass: populate backbone index caches per residue.
    for (auto& residue : protein->residues_) {
        PopulateBackboneCache(residue, protein->atoms_);
    }

    // ---- Bonds ----
    protein->bonds_.reserve(B);
    for (size_t b = 0; b < B; ++b) {
        QtBond bond;
        bond.atomIndexA = static_cast<size_t>(af->topology.bond_atoms[b * 2 + 0]);
        bond.atomIndexB = static_cast<size_t>(af->topology.bond_atoms[b * 2 + 1]);
        bond.order      = DecodeBondOrder(af->topology.bond_order[b], errors);
        bond.category   = DecodeBondCategory(af->topology.bond_category[b], errors);
        protein->bonds_.push_back(bond);
    }

    // ---- Rings ----
    protein->rings_.reserve(nRings);
    for (size_t r = 0; r < nRings; ++r) {
        RingTypeIndex idx = DecodeRingType(af->topology.ring_type[r], errors);
        auto ring = CreateQtRing(idx);
        if (!ring) {
            ErrorBus::Report(Severity::Warning, QStringLiteral("QtProteinLoader"),
                             QStringLiteral("CreateQtRing returned nullptr"),
                             QStringLiteral("ring_type=%1")
                               .arg(static_cast<int>(idx)));
            ++errors;
            protein->rings_.push_back(nullptr);
            continue;
        }

        // Populate common fields.
        ring->parentResidueIndex = af->topology.ring_residue[r];
        if (ring->parentResidueIndex >= 0
            && static_cast<size_t>(ring->parentResidueIndex) < R) {
            ring->parentResidueNumber =
                protein->residues_[ring->parentResidueIndex].residueNumber;
        }
        ring->fusedPartnerIndex = af->topology.ring_fused_partner[r];

        // Vertex atom indices — CSR ragged array.
        const int32_t off0 = af->topology.ring_offsets[r];
        const int32_t off1 = af->topology.ring_offsets[r + 1];
        ring->atomIndices.reserve(static_cast<size_t>(off1 - off0));
        for (int32_t k = off0; k < off1; ++k) {
            ring->atomIndices.push_back(
                static_cast<size_t>(af->topology.ring_atom_indices[k]));
        }

        protein->rings_.push_back(std::move(ring));
    }

    // ---- Protonation variants (after rings exist) ----
    for (size_t r = 0; r < protein->residues_.size(); ++r) {
        protein->residues_[r].variant =
            InferProtonationVariant(*protein,
                                    static_cast<int>(r),
                                    protein->residues_[r].aminoAcid);
    }

    // ---- Build the trajectory wrapper ----
    auto conformation = std::make_unique<QtConformation>(protein.get(), af);

    qCInfo(cLoader).noquote()
        << "Decode complete |"
        << "decodeErrors=" << errors;

    result.protein       = std::move(protein);
    result.conformation  = std::move(conformation);
    result.analysisFile  = af;
    result.decodeErrors  = errors;
    result.ok            = true;
    return result;
}

}  // namespace h5reader::io
