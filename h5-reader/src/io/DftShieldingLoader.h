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
    // Parse + validate one DFT job. Reads meta.json -> files.out_primary,
    // parses the .out, validates against the protein topology. Returns null
    // and logs at the seam on any failure (missing meta.json, parser hole,
    // atom-count mismatch, dia+para != total).
    //
    // jobsDir is the dft/jobs root; originalIndex is the original XTC frame
    // index used in the job-dir naming (_fNNNNNN_t<ps>).
    static std::shared_ptr<const h5reader::model::DftShieldingFrame>
    LoadAndValidate(std::size_t originalIndex,
                    const QString& jobsDir,
                    const h5reader::model::QtProtein* protein);
};

}  // namespace h5reader::io
