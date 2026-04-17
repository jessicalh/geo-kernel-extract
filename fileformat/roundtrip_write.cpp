// Roundtrip all H5 files: read -> write to <basename>_rt.h5
// Does NOT delete the output — leaves them for NPY cross-validation.

#include "analysis_file.h"
#include <cstdio>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

int main(int argc, char** argv) {
    std::string dir = argc > 1 ? argv[1] : "test";
    std::string out_dir = argc > 2 ? argv[2] : (dir + "/rt");

    fs::create_directories(out_dir);

    for (const auto& entry : fs::directory_iterator(dir)) {
        if (entry.path().extension() != ".h5") continue;
        std::string stem = entry.path().stem().string();
        if (stem.find("_rt") != std::string::npos) continue;  // skip rt files

        std::string src = entry.path().string();
        std::string dst = out_dir + "/" + stem + ".h5";

        printf("%s -> %s ... ", stem.c_str(), dst.c_str());
        fflush(stdout);

        AnalysisFile af;
        af.ReadH5(src);
        af.WriteH5(dst);

        printf("OK (%zu frames, %zu atoms)\n", af.n_frames, af.n_atoms);
    }
    return 0;
}
