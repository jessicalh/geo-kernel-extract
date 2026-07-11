#pragma once
//
// Raw local-frame geometric pi-quadrupole tensor companion.
//
// PiQuadrupoleResult preserves the historical axial A-term scalar.  This
// result emits the full symmetric/traceless geometry tensor in a deterministic
// ring-local frame without changing that scalar channel.
//

#include "ConformationResult.h"
#include "Ring.h"
#include "Types.h"

#include <typeindex>
#include <vector>

namespace nmr {

class ProteinConformation;

namespace pi_quadrupole_local_tensor_detail {

struct LocalFrame {
    Vec3 x_axis = Vec3::Zero();
    Vec3 y_axis = Vec3::Zero();
    Vec3 z_axis = Vec3::Zero();
    bool valid = false;
};

struct TensorEvaluation {
    Mat3 tensor = Mat3::Zero();
    bool valid = false;
};

LocalFrame BuildLocalFrame(const RingGeometry& geom);

// Exact production formula for the local geometric EFG tensor (A^-5).
TensorEvaluation ComputeLocalTensor(const Vec3& d_local);

}  // namespace pi_quadrupole_local_tensor_detail


class PiQuadrupoleLocalTensorResult : public ConformationResult {
public:
    std::string Name() const override {
        return "PiQuadrupoleLocalTensorResult";
    }

    std::vector<std::type_index> Dependencies() const override;

    static std::unique_ptr<PiQuadrupoleLocalTensorResult> Compute(
        ProteinConformation& conf);

    int WriteFeatures(const ProteinConformation& conf,
                      const std::string& output_dir) const override;
};

}  // namespace nmr
