// DftShielding — QM ground-truth shielding for one atom, ingested from a DFT
// NMR run (ORCA .out now; a LeanSCF tensors file — /shared/dft-ex1 — later).
//
// Kept WHOLE, not fluffed (user, 2026-05-27): the TOTAL tensor AND its
// DIAMAGNETIC and PARAMAGNETIC parts, each decomposed to a SphericalTensor
// (T0/T1/T2). This is downstream of the library (nmr-extract emits the PDB ->
// ORCA runs on it -> .out), so ingesting it is reading a QM artifact, not
// re-deriving kernels — and the reader keeps more than ORCA's summary or the
// dft-ex1 JSON bother to (which drop the dia/para split and the antisymmetric
// T1). It is its own result category, source-agnostic behind this type.

#pragma once

#include "Types.h"  // SphericalTensor, Element

#include <vector>

namespace h5reader::model {

struct DftAtomShielding {
    SphericalTensor total;  // T0 (iso, ppm) / T1 (antisym) / T2 (sym-traceless)
    SphericalTensor dia;    // diamagnetic part
    SphericalTensor para;   // paramagnetic part  (total == dia + para, checked)
    Element         element = Element::Unknown;  // cross-check vs topology order
};

// One frame's DFT shielding: one entry per atom, indexed by ORCA nucleus index
// (== emitted-PDB atom order == topology order). `valid` is false if the
// shielding section was absent or empty.
struct DftShieldingFrame {
    std::vector<DftAtomShielding> atoms;
    bool                          valid = false;
};

}  // namespace h5reader::model
