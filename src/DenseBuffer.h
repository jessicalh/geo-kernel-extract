#pragma once
//
// DenseBuffer<T>: contiguous per-atom storage for TrajectoryResult
// outputs whose indexing is (atom, stride) with stride meaning frame,
// lag, frequency, or similar. Owned by TrajectoryProtein after a
// TrajectoryResult transfers ownership at Finalize.
//
// Layout: atom-major. Storage is a single flat vector<T> of length
// (atom_count × stride_per_atom). Per-atom slices are contiguous in
// memory, which is cache-friendly for per-atom autocorrelation and
// time-series access patterns.
//
// Typed with the native C++ type (T = double, Vec3, Mat3,
// SphericalTensor). These are memory-contiguous structs whose
// sizeof(T) corresponds to the natural serialisation size; no
// intermediate flattening to double arrays is needed. H5 emission
// flattens to trailing dimensions with metadata attributes describing
// irrep layout (see spec/WIP_OBJECT_MODEL.md Appendix C).
//
// For the pattern, see spec/WIP_OBJECT_MODEL.md §4 (dense buffers on
// TrajectoryProtein) + Appendix C (layout conventions).
//

#include <cstddef>
#include <stdexcept>
#include <vector>

namespace nmr {

// Non-template base so TrajectoryProtein can own buffers of mixed T
// via unique_ptr<DenseBufferBase>.
class DenseBufferBase {
public:
    virtual ~DenseBufferBase() = default;
    virtual size_t AtomCount() const = 0;
    virtual size_t StridePerAtom() const = 0;
    virtual size_t ElementSizeBytes() const = 0;
};


template <typename T>
class DenseBuffer : public DenseBufferBase {
public:
    DenseBuffer(size_t atom_count, size_t stride_per_atom)
        : atom_count_(atom_count),
          stride_per_atom_(stride_per_atom),
          storage_(atom_count * stride_per_atom) {}

    // Typed element access: atom index + offset within the stride.
    T& At(size_t atom_idx, size_t offset) {
        return storage_[atom_idx * stride_per_atom_ + offset];
    }
    const T& At(size_t atom_idx, size_t offset) const {
        return storage_[atom_idx * stride_per_atom_ + offset];
    }

    // Per-atom slice pointer + length (C++17; std::span would simplify
    // but is C++20 per PATTERNS.md).
    T* AtomSlicePtr(size_t atom_idx) {
        return storage_.data() + atom_idx * stride_per_atom_;
    }
    const T* AtomSlicePtr(size_t atom_idx) const {
        return storage_.data() + atom_idx * stride_per_atom_;
    }

    // Raw storage access for H5 emission. Layout is atom-major,
    // contiguous per-atom.
    T* RawData() { return storage_.data(); }
    const T* RawData() const { return storage_.data(); }
    size_t TotalElementCount() const { return storage_.size(); }

    size_t AtomCount() const override { return atom_count_; }
    size_t StridePerAtom() const override { return stride_per_atom_; }
    size_t ElementSizeBytes() const override { return sizeof(T); }

private:
    size_t atom_count_;
    size_t stride_per_atom_;
    std::vector<T> storage_;
};

}  // namespace nmr
