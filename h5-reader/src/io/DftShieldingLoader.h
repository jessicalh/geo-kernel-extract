#pragma once

#include "../model/DftShielding.h"

#include <QString>

#include <cstddef>
#include <memory>

namespace h5reader::model {
class QtProtein;
}

namespace h5reader::io {

class DftShieldingLoader {
public:
    // Parse + validate one DFT job. meta_json_abspath is the path to the
    // job's `<job_id>_meta.json` — the calcset's `.LGS` carries it as
    // `dft.frames[].meta_json`, so there is no dir-name parsing.
    // Reads files.out_primary, parses the .out, validates against the
    // protein topology. Returns null and logs at the seam on any failure
    // (cannot open meta.json, no files.out_primary, parser hole,
    // atom-count mismatch, dia+para != total).
    static std::shared_ptr<const h5reader::model::DftShieldingFrame>
    LoadAndValidate(const QString& meta_json_abspath,
                    const h5reader::model::QtProtein* protein);
};

}  // namespace h5reader::io
