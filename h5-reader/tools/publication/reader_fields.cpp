#include "io/FrameFieldPolicy.h"
#include <iostream>

int main() {
    for (const auto& field : h5reader::io::kFieldCatalog)
        if (h5reader::io::ShouldLoadFrameField(field.kind))
            std::cout << field.stem << ".npy\n";
    return std::cout.good() ? 0 : 1;
}
