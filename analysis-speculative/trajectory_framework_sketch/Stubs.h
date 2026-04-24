// Stubs.h
//
// Forward declarations and minimal stub types for classes the spec
// treats as pre-existing: Protein, ProteinConformation, ConformationAtom,
// ConformationResult, RunOptions, OperationRunner, AIMNet2Model,
// GromacsEnergy, SolventEnvironment.
//
// These are NOT reimplementations. They exist only so the sketch code
// compiles against a consistent set of names. Real bodies live elsewhere
// (src/Protein.{h,cpp}, src/ProteinConformation.{h,cpp}, etc.).
//
// Per sandbox scope: do not model anything after the existing src/. The
// shapes here are literal from OBJECT_MODEL.md and WIP_OBJECT_MODEL.md
// references only (e.g. "Protein has an AtomAt(i)", "ProteinConformation
// has AtomAt(i) returning ConformationAtom").
//
// HighFive::File is also forward-declared. Real code #includes highfive.

#ifndef TRAJECTORY_FRAMEWORK_SKETCH_STUBS_H
#define TRAJECTORY_FRAMEWORK_SKETCH_STUBS_H

#include <cstddef>
#include <string>
#include <vector>

// External libraries — forward declarations only.
namespace HighFive { class File; }

// Eigen types — real library provides these.
namespace Eigen { template <typename, int, int, int, int, int> class Matrix; }

// Stand-ins for typed Eigen aliases used in the spec.
struct Vec3 { double x = 0.0, y = 0.0, z = 0.0; };
struct Mat3 { double m[3][3] = {{0}}; };

// Stand-in for SphericalTensor (OBJECT_MODEL.md §"SphericalTensor"):
// T0 scalar, T1[3] antisymmetric, T2[5] traceless symmetric.
struct SphericalTensor {
    double T0 = 0.0;
    double T1[3] = {0.0, 0.0, 0.0};
    double T2[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
};

// Minimal stand-ins for library classes. Real definitions are in
// the main src/ tree; this sandbox does not recreate their behaviour.
// Only the methods invoked by the framework sketch are declared here,
// with the exact signatures the spec cites.
//
// Each of these stubs exists because the framework sketch calls one
// or two methods on it; their real bodies live elsewhere.

// Protein: identity + topology. Framework sketch calls AtomCount()
// (from §3: "TrajectoryProtein::AtomCount() const { ... == protein().
// AtomCount() }"). Real class per OBJECT_MODEL.md §"Protein".
class Protein {
public:
    std::size_t AtomCount() const { return 0; }
};

// ProteinConformation: per-frame geometric instance of a Protein.
// Framework sketch calls AtomAt(i) in TrajectoryResult::Compute
// reading ConformationAtom fields (§4 worked example). Real class
// per OBJECT_MODEL.md §"ProteinConformation".
class ConformationAtom {
public:
    // bs_shielding_contribution is the field BiotSavartResult writes
    // per OBJECT_MODEL.md Atom field tables; BsWelfordTrajectoryResult
    // reads its .T0 component per §4 example.
    SphericalTensor bs_shielding_contribution;
};

class ProteinConformation {
public:
    const ConformationAtom& AtomAt(std::size_t) const {
        static ConformationAtom sentinel;
        return sentinel;
    }
};

// ConformationResult: only referenced by Dependencies() returning
// type_index values (e.g. typeid(BiotSavartResult)). The sketch
// forward-declares the specific result types it cites.
class ConformationResult;
class BiotSavartResult;  // cited by BsWelfordTrajectoryResult Dependencies

// Operation plumbing. The framework references these only by pointer /
// type_index; full bodies are elsewhere.
class AIMNet2Model;
class GromacsEnergy;
class SolventEnvironment;
class OperationRunner;

// A minimal RunOptions shape. Per WIP §6, RunConfiguration stores a
// RunOptions value and Trajectory::Run hands it to the per-frame work.
// Fields here are the skip flags and aimnet2 pointer the spec names
// explicitly in §6's factory bodies; anything else belongs to the
// existing src/RunOptions.h which this sandbox does not redefine.
struct RunOptions {
    bool skip_dssp = false;
    bool skip_mopac = false;
    bool skip_apbs = false;
    bool skip_coulomb = false;
    AIMNet2Model* aimnet2_model = nullptr;
};

// Forward declaration for the per-frame reader. The spec treats this
// as process-scope infrastructure (§9: "stays Gromacs-named per format").
// Spec §5 shows Trajectory::Run driving it via Open / Next /
// LastConformation / LastTimePs / LastXtcIndex. Full body elsewhere.
class GromacsFrameHandler;

#endif // TRAJECTORY_FRAMEWORK_SKETCH_STUBS_H
