// QtProtein implementation — ring-type counting helper.

#include "QtProtein.h"

#include "QtResidueNames.h"  // derived residue-name helpers (residueLabel)

#include <QString>

namespace h5reader::model {

int QtProtein::ringCountByType(RingTypeIndex t) const {
    if (!topology_)
        return 0;
    int n = 0;
    for (std::size_t i = 0; i < topology_->ringCount(); ++i) {
        const QtRing& r = topology_->ringAt(i);
        if (r.TypeIndex() == t)
            ++n;
    }
    return n;
}

// ---- Name projections (selectable; see Types.h NamingConvention/Source).
//      Labels are projections, never identity — the reader is not
//      label-driven; these only pick which display / ML-join label to show.

QString QtProtein::atomLabel(std::size_t i, NamingConvention conv) const {
    // Atom names are verbatim-only (the library derived them via atom_nom.tbl;
    // the reader reads the result). Provenance is on QtAtomNames if needed.
    const QtAtomNames& n = atomNames_[i];
    switch (conv) {
    case NamingConvention::Amber: return n.amber;
    case NamingConvention::Iupac: return n.iupac;
    case NamingConvention::Bmrb:  return n.bmrb;
    }
    return n.amber;
}

QString QtProtein::residueLabel(std::size_t i, NamingConvention conv, NamingSource src) const {
    // Verbatim: the sidecar projection as CategoryInfoProjection wrote it —
    // the label to match BMRB records as-deposited.
    if (src == NamingSource::Verbatim && i < residueNames_.size()) {
        const QtResidueNames& n = residueNames_[i];
        switch (conv) {
        case NamingConvention::Amber: return n.amber;
        case NamingConvention::Iupac: return n.iupac;
        case NamingConvention::Bmrb:  return n.bmrb;
        }
    }
    // Derived: recomputed from the typed AminoAcid (+ variant). Also the
    // fallback when the verbatim projection wasn't loaded (degraded sidecar).
    const QtResidue& res = residues_[i];
    switch (conv) {
    case NamingConvention::Amber:
        return QString::fromLatin1(AmberResidue3LetterFor(res.aminoAcid, res.protonationVariantIndex));
    case NamingConvention::Iupac:
        return QString::fromLatin1(IupacResidue3LetterFor(res.aminoAcid));
    case NamingConvention::Bmrb:
        // No reader-side BMRB derivation; the canonical IUPAC 3-letter is the
        // variant-blind equivalent (HID/HIE/HIP -> HIS). Documented choice.
        return QString::fromLatin1(IupacResidue3LetterFor(res.aminoAcid));
    }
    return QString();
}

}  // namespace h5reader::model
