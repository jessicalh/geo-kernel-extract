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
// Reads ONLY (no model mutation):
//   - protein.LegacyAmber().BondList()            (Bond struct)
//   - protein.LegacyAmber().Rings()               (RingTopology surface)
//   - protein.LegacyAmber().HasAtomSemantic()     (gates substrate fields)
//   - protein.AtomCount(), ResidueCount()         (manifest axis sizes)
//
// Same architectural shape as CategoryInfoProjection and FramePdbEmitter:
// singleton static-method, fixed shape, no virtuals, no
// ConformationResult / TrajectoryResult lifecycle. Holds no per-frame
// state. Called once per Protein from each entry point alongside the
// CategoryInfoProjection call.
//
// All five output files are INVARIANT per-protein -- topology does not
// change between conformations or trajectory frames. Calling
// WriteFeatures multiple times per Protein on the same output_dir is
// idempotent; the second write overwrites with bit-identical content.
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
    // Inert when the protein has no LegacyAmber substrate
    // (HasAtomSemantic() false AND zero bonds AND zero rings) -- mirrors
    // CategoryInfoProjection's stub-fixture behavior. Manifest is still
    // emitted with the populated flags so a consumer can distinguish
    // "no topology" from "topology omitted".
    static int WriteFeatures(const Protein& protein,
                              const std::string& output_dir,
                              const std::string& protein_id = "");

    TopologySidecar() = delete;
    TopologySidecar(const TopologySidecar&) = delete;
    TopologySidecar& operator=(const TopologySidecar&) = delete;
};

}  // namespace nmr
