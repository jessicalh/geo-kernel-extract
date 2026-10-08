#pragma once

#include <array>

namespace h5reader::app {

using TensorAxisColours = std::array<std::array<double, 3>, 3>;

// Principal-value order, shared by the scene arrows and inspector swatches.
inline constexpr TensorAxisColours kDefaultTensorColours{{
    {0.96, 0.66, 0.16},
    {0.18, 0.74, 0.74},
    {0.74, 0.36, 0.86},
}};
inline constexpr TensorAxisColours kShieldingTensorColours{{
    {0.96, 0.72, 0.18}, // gold
    {1.00, 0.40, 0.28}, // coral
    {0.93, 0.36, 0.66}, // pink
}};
inline constexpr TensorAxisColours kOrientationTensorColours{{
    {0.14, 0.78, 0.64}, // green-teal
    {0.25, 0.68, 0.96}, // sky blue
    {0.57, 0.50, 1.00}, // lavender
}};

} // namespace h5reader::app
