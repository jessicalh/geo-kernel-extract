// QtProteinLoader — the H5 → typed-model boundary.
//
// This is the ONE file where H5 ints become enums and H5 strings become
// AminoAcid / QtRing subclasses / ProtonationVariant. After Load()
// returns, no code in the reader touches AnalysisFile's enum ordinals
// or name strings for physics dispatch.
//
// Load() reads the AnalysisFile once, decodes every enum at the
// boundary, builds the QtProtein (identity/topology) and
// QtConformation (trajectory) together, and returns them. Decoder
// failures (unknown ordinal, unknown AA code) are reported through
// ErrorBus and the loader's BuildResult.ok flag; the caller decides
// whether to continue with a partially-decoded model or abort.
//
// NEVER add string-dispatch logic outside this file. If a future
// feature needs a typed property that isn't in the H5, add it to
// QtAtom/QtResidue as a decoded enum and populate it here.

#pragma once

#include "../model/QtProtein.h"
#include "../model/QtConformation.h"
#include "analysis_file.h"

#include <QString>
#include <memory>
#include <string>

namespace h5reader::io {

struct QtLoadResult {
    std::unique_ptr<model::QtProtein>           protein;
    std::unique_ptr<model::QtConformation>      conformation;
    std::shared_ptr<const AnalysisFile>         analysisFile;
    QString                                     proteinId;     // derived from H5 path stem
    bool                                        ok = false;
    int                                         decodeErrors = 0;
    QString                                     error;         // empty if ok
};

class QtProteinLoader {
public:
    // Open the H5 file, decode the typed model. On HighFive exceptions
    // or I/O failures, returns a result with ok=false and a populated
    // error string. Individual decode errors (bad enum ordinal, unknown
    // AA code) are logged via ErrorBus::Report and counted in
    // decodeErrors — load continues so the inventory can still report.
    static QtLoadResult Load(const QString& h5Path);
};

}  // namespace h5reader::io
