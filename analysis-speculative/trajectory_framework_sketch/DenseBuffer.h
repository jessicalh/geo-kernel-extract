// DenseBuffer.h
//
// Dense buffer template owned by TrajectoryProtein at trajectory scope.
// Per WIP_OBJECT_MODEL.md §4 ("Dense buffers on TrajectoryProtein") and
// Appendix C ("Dense buffer layouts for non-scalar payloads").
//
// Atom-major layout: one contiguous region per atom, `stride_per_atom_`
// elements long. Used by TrajectoryResults whose output is per-atom x N
// (where N is a frame count, lag count, or similar), so per-atom time
// series stays contiguous.
//
// Template parameter T is a native C++ type (double, Vec3, Mat3,
// SphericalTensor) per Appendix C. H5 emission flattens the trailing
// dimensions; that is the writer's concern, not this class's.
//
// This class exists to keep per-atom dense data out of
// std::vector<...> per-atom fields (which would fragment). It does not
// know anything about frames, lags, or frequencies; the meaning of
// "stride" is documented by the owning TrajectoryResult's WriteH5Group.

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_DENSEBUFFER_H
#define TRAJECTORY_FRAMEWORK_SKETCH_DENSEBUFFER_H

#include <cstddef>
#include <vector>

// Non-template base so TrajectoryProtein can hold a polymorphic map
// of buffers keyed by type_index. Per WIP §3 TrajectoryProtein structure:
//     std::unordered_map<std::type_index,
//           std::unique_ptr<DenseBufferBase>> dense_buffers_;
class DenseBufferBase {
public:
    virtual ~DenseBufferBase() = default;
};

template <typename T>
class DenseBuffer : public DenseBufferBase {
public:
    DenseBuffer(std::size_t atom_count, std::size_t stride_per_atom)
        : atom_count_(atom_count),
          stride_per_atom_(stride_per_atom),
          storage_(atom_count * stride_per_atom) {}

    // Per Appendix C: typed element access by (atom, offset-within-stride).
    T& At(std::size_t atom_idx, std::size_t offset) {
        return storage_[atom_idx * stride_per_atom_ + offset];
    }
    const T& At(std::size_t atom_idx, std::size_t offset) const {
        return storage_[atom_idx * stride_per_atom_ + offset];
    }

    // Raw data + element count for H5 emission (cast to const double* by
    // the writer if T is a contiguous-doubles struct, per Appendix C).
    const T* RawData() const { return storage_.data(); }
    std::size_t TotalElementCount() const { return storage_.size(); }

    std::size_t AtomCount() const { return atom_count_; }
    std::size_t StridePerAtom() const { return stride_per_atom_; }

private:
    std::size_t atom_count_;
    std::size_t stride_per_atom_;
    std::vector<T> storage_;  // atom-major
};

#endif // TRAJECTORY_FRAMEWORK_SKETCH_DENSEBUFFER_H
