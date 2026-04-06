#pragma once
//
// DemoResult: a trivial ConformationResult that exercises the complete pattern.
//
// Computes per-atom distance to nearest ring center.
// Stores a double and a Vec3 on ConformationAtom.
// Declares dependency on GeometryResult.
// Also decomposes a test Mat3 via SphericalTensor to prove the pipeline.
//

#include "ConformationResult.h"
#include "GeometryResult.h"
#include "ProteinConformation.h"

namespace nmr {

class DemoResult : public ConformationResult {
public:
    std::string Name() const override { return "DemoResult"; }

    std::vector<std::type_index> Dependencies() const override {
        return { std::type_index(typeid(GeometryResult)) };
    }

    // Factory
    static std::unique_ptr<DemoResult> Compute(ProteinConformation& conf);

    // Physics query methods
    double NearestRingDistance(size_t atom_index) const;
    Vec3 NearestRingDirection(size_t atom_index) const;

    // SphericalTensor demonstration
    const SphericalTensor& TestDecomposition() const { return test_decomposition_; }
    const Mat3& TestReconstructed() const { return test_reconstructed_; }

private:
    const ProteinConformation* conf_ = nullptr;
    SphericalTensor test_decomposition_;
    Mat3 test_reconstructed_ = Mat3::Zero();
};

}  // namespace nmr
