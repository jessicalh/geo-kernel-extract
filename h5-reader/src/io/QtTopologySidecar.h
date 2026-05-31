// QtTopologySidecar — orchestrates the 5 sidecar NPY files +
// extraction_manifest.json into the typed Qt model objects.
//
// Reads:
//   - extraction_manifest.json → QtManifest + QtEnumVocab
//   - atoms_category_info.npy → vector<QtAtom> + vector<QtAtomNames>
//   - residues.npy → vector<QtResidue> (atom-membership + backbone cache
//     populated from atom-side BackboneRole/Locant)
//   - bonds.npy → vector<QtBond>
//   - rings.npy → vector<unique_ptr<QtRing>> (via CreateQtRing factory)
//   - ring_membership.npy → vector<QtRingMembership>
//
// Cross-validates axis sizes (manifest.axisSizes.atom == atoms NPY row
// count, etc.) and logs Warn on mismatch with fallback to actual row
// counts.
//
// Produces a LoadResult struct that QtProteinLoader hands to QtProtein
// for ownership transfer. The sidecar reader itself is a one-shot
// static method, no instance state held past the call.
//
// Cross-platform: QFile + QJsonDocument + QtNpyReader. No POSIX.

#pragma once

#include "QtEnumVocab.h"
#include "QtManifest.h"

#include "../model/QtAtom.h"
#include "../model/QtAtomNames.h"
#include "../model/QtResidueNames.h"
#include "../model/QtBond.h"
#include "../model/QtResidue.h"
#include "../model/QtRing.h"
#include "../model/QtRingMembership.h"

#include <QString>
#include <cstddef>
#include <memory>
#include <vector>

namespace h5reader::io {

class QtTopologySidecar {
public:
    struct LoadResult {
        bool ok = false;
        QString error;         // empty on success
        int warningCount = 0;  // count of non-fatal warnings logged

        QtManifest manifest;
        QtEnumVocab vocab;

        std::vector<h5reader::model::QtAtom> atoms;
        std::vector<h5reader::model::QtAtomNames> atomNames;
        std::vector<h5reader::model::QtResidueNames> residueNames;  // verbatim projection
        std::vector<h5reader::model::QtResidue> residues;
        std::vector<h5reader::model::QtBond> bonds;
        std::vector<std::unique_ptr<h5reader::model::QtRing>> rings;
        std::vector<h5reader::model::QtRingMembership> ringMemberships;

        std::size_t aromaticRingCount = 0;
        std::size_t saturatedRingCount = 0;
    };

    // sidecar_dir is the directory containing the 5 NPYs; typically
    // dirname(trajectory.h5). manifest_path is the explicit path to
    // extraction_manifest.json (the `.LGS` carries it as
    // trajectory.extraction_manifest; the legacy
    // sidecar_dir/extraction_manifest.json convention is implied when
    // omitted). Files are looked up by their canonical names (no
    // globbing, no discovery — per feedback_no_file_discovery).
    static LoadResult Load(const QString& sidecar_dir,
                           const QString& manifest_path = QString());
};

}  // namespace h5reader::io
