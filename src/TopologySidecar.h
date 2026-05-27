#pragma once
//
// TopologySidecar -- output-side per-protein topology projection.
// Emits five sibling files alongside the existing NPY output:
//
//   residues.npy              structured NPY, one row per Residue
//   bonds.npy                 structured NPY, one row per Bond
//   rings.npy                 structured NPY, one row per Ring
//                             (aromatic rings first, then saturated)
//   ring_membership.npy       structured NPY, one row per (ring, atom)
//   extraction_manifest.json  topology declarations + axis sizes
//
// Primary reads (no model mutation):
//   - protein residues and atoms                  (residue rows + axis sizes)
//   - protein.LegacyAmber().BondList()            (Bond struct)
//   - protein.LegacyAmber().Rings()               (RingTopology surface)
//   - protein.LegacyAmber().HasAtomSemantic()     (manifest flags)
//
// Same architectural shape as CategoryInfoProjection and FramePdbEmitter:
// singleton static-method, fixed shape, no virtuals, no
// ConformationResult / TrajectoryResult lifecycle. Holds no per-frame
// state. Called at output boundaries alongside the
// CategoryInfoProjection call.
//
// The four NPY outputs are invariant per-protein -- topology does not
// change between conformations or trajectory frames. Calling
// WriteFeatures multiple times per Protein on the same output_dir
// overwrites the files; the manifest carries a fresh generated_at_utc
// timestamp.
//
// Architectural rule (memory feedback_naming_input_output_asymmetry +
// the OBJECT_MODEL "Protein is identity and topology only" wall):
// this projection exists at the output boundary. Strings are an
// EXPORTED projection surface (chain_id, atom-type names); the
// internal substrate stays typed.
//

#include <cstddef>
#include <string>

namespace nmr {

class Protein;

class TopologySidecar {
public:
    // Emit residues.npy, bonds.npy, rings.npy, ring_membership.npy, and
    // extraction_manifest.json to ``output_dir``. Returns the number of
    // files written (5 on full success). Idempotent; safe to call
    // multiple times.
    //
    // ``protein_id`` is recorded in the manifest. Empty string is
    // acceptable; callers typically pass fs::path(output_dir).filename().
    //
    // Always emits the five files when writes succeed. The manifest's
    // populated flags let consumers distinguish absent substrate fields
    // from omitted topology output.
    static int WriteFeatures(const Protein& protein,
                              const std::string& output_dir,
                              const std::string& protein_id = "");

    TopologySidecar() = delete;
    TopologySidecar(const TopologySidecar&) = delete;
    TopologySidecar& operator=(const TopologySidecar&) = delete;
};

}  // namespace nmr
