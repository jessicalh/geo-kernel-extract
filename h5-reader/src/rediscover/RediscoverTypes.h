// RediscoverTypes — the typed primitives the substrate emits: the DFT target,
// one source slot, and the per-(atom, frame) neighbourhood record. Plain data
// classes (no QObject), per DESIGN.md.
//
// The record is row-kind-agnostic carrier data; the CSV sink (RecordSink)
// turns it into the two emitted row kinds — un-summed per-source rows and one
// aggregated per-(atom, frame) row — per the lead's decision.

#pragma once

#include "../model/Types.h"  // Vec3, Mat3, SphericalTensor, enums

#include <QString>

#include <cstdint>
#include <vector>

namespace h5reader::rediscover {

using model::Mat3;
using model::SphericalTensor;
using model::Vec3;

// ── The DFT target ────────────────────────────────────────────────────────
// Raw 3x3 (total/dia/para) kept verbatim from ORCA + the library-basis
// decomposition of the raw total tensor (DecomposeLibrary) + σ_iso (= T0,
// rotation-invariant, the safe scalar). dia/para are diagnostics.
//
// The raw total tensor is ALSO carried in the per-atom local frame
// (total_local) so a tensor-valued / equivariant fitter has the target in the
// same rotation-stable frame as the source vectors. The T2-component caveat
// (cross-frame validity) is recorded as a flag, not silently assumed.
struct DftTarget {
    bool present = false;  // false ⇒ no DFT job for this (atom, frame)

    Mat3 total_raw = Mat3::Zero();
    Mat3 dia_raw   = Mat3::Zero();
    Mat3 para_raw  = Mat3::Zero();

    // Library-basis decomposition of total_raw (and dia/para for diagnostics).
    SphericalTensor total_decomp;  // T0 == σ_iso
    SphericalTensor dia_decomp;
    SphericalTensor para_decomp;

    // total_raw rotated into the atom's local frame (only meaningful when the
    // local frame is valid). For an equivariant fit the source vectors and the
    // target tensor are then both in {x,y,z}.
    Mat3 total_local = Mat3::Zero();
    bool local_frame_valid = false;
};

// ── One source in the neighbourhood ───────────────────────────────────────
// A source is an aromatic ring, anisotropic bond, or charge site. The fields
// are a superset; the per-extraction schema names which columns it fills.
enum class SourceKind : int { Ring = 0, Bond = 1, Charge = 2 };

struct SourceSlot {
    SourceKind kind = SourceKind::Ring;

    // Geometry, all in the target atom's LOCAL frame unless noted.
    Vec3   disp_local = Vec3::Zero();  // displacement target→source-reference (Å)
    double r          = 0.0;           // |displacement| (Å)
    double cos_theta  = 0.0;           // cosθ about the source axis (ring: z/r)
    double dipolar    = 0.0;           // (3cos²θ − 1) / r³  (this source's term)

    // Ring-only geometry, read straight from the H5 ring-neighbourhood TS.
    double ring_z              = 0.0;  // signed axial distance (Å)
    double ring_rho            = 0.0;  // in-plane radial distance (Å)
    double ring_in_plane_angle = 0.0;  // radians
    // The SOURCE ring's unit normal in the target's local frame — the dipole
    // axis the ring-current l=2 tensor is oriented by (the scalar cosθ alone
    // can't reconstruct the tensor). Sign is irrelevant (the l=2 form is even
    // under n → −n). Ring sources only; zero for bonds.
    Vec3   source_normal_local = Vec3::Zero();

    // Ring identity / physics (frame-0 membership snapshot; QtRing virtuals).
    int32_t ring_index       = -1;     // absolute rings.npy row (source provenance)
    bool   is_self_or_bonded = false;  // source ring is the H's own ring (or fused
                                       // with it) — the producer excludes these
                                       // from its kernel; we emit both sums so the
                                       // fitter decides.
    int    ring_type_index   = -1;     // RingTypeIndex ordinal
    double ring_intensity    = 0.0;    // LiteratureIntensity (nA/T)
    int    ring_nitrogen     = 0;      // NitrogenCount
    double ring_jb_offset    = 0.0;    // JohnsonBoveyLobeOffset (Å)
    int    ring_aromaticity  = -1;     // RingAromaticity ordinal
    int    ring_size         = 0;      // RingSizeValue
    bool   ring_fused        = false;  // IsFused
    bool   ring_jb_kernel_present = false;  // source-level Johnson-Bovey kernel emitted by C++
    SphericalTensor ring_jb_unit_kernel;  // unit-current source kernel in target local frame (ppm_T_per_nA)
    SphericalTensor ring_jb_kernel;    // fixed-literature source kernel in target local frame (ppm)

    // Bond-only identity + orientation (McConnell). bond_axis_local is the UNIT
    // bond axis (B-A) expressed in the target atom's local frame. A/B is the
    // index-oriented sidecar endpoint order, not a chemical C->O/C->N direction.
    int    bond_category   = -1;       // BondCategory ordinal
    int    bond_order      = -1;       // BondOrder ordinal
    int    bond_elem_a     = -1;       // Element ordinal of endpoint A
    int    bond_elem_b     = -1;       // Element ordinal of endpoint B
    int32_t bond_index     = -1;       // row in bonds.npy (provenance)
    int32_t bond_atom_a    = -1;       // endpoint A atom index (provenance)
    int32_t bond_atom_b    = -1;       // endpoint B atom index (provenance)
    Vec3   bond_axis_local = Vec3::Zero();  // unit (B−A) in the local frame
    bool   mc_source_is_self_or_bonded = false;  // producer McConnell exclusion:
                                                // endpoint self-source OR inside
                                                // bond-length near-field radius.
    bool   bond_mc_lit_kernel_present = false;
    SphericalTensor bond_mc_lit_kernel;  // Δχ-scaled McConnell PCS kernel, local frame (ppm)

    // Charge-site identity + weight (charge_dipole / charge_quadrupole). The
    // displacement fields above are still the target→source vector in the
    // target's local frame; source_q_e is the selected charge source's value.
    QString source_charge_source;
    double  source_q_e = 0.0;
    int32_t source_atom_index = -1;
    int32_t source_residue_index = -1;
    int32_t source_residue_number = 0;
    int     source_amino_acid = -1;
    int     source_element = -1;
    QString source_atom_name;
};

// ── The per-(atom, frame) neighbourhood record ────────────────────────────
struct NeighborhoodRecord {
    // Identity / frame.
    int32_t atom_index     = -1;
    int32_t residue_index  = -1;
    int32_t residue_number = 0;       // PDB numbering, display
    int     amino_acid     = -1;      // AminoAcid ordinal
    int     element        = -1;      // Element ordinal
    QString atom_name;                // verbatim projection (display only)
    QString stratum;                  // "ring_current" | "mcconnell"
    int32_t h5_row         = -1;      // H5 frame row
    int32_t original_index = -1;      // original (trr) frame index
    double  time_ps        = 0.0;

    // Local frame recorded with the record (the basis the vectors live in).
    Vec3 frame_z = Vec3::Zero();
    Vec3 frame_x = Vec3::Zero();
    Vec3 frame_y = Vec3::Zero();
    int  frame_variant = 0;           // FrameVariant ordinal
    bool frame_valid   = false;
    // The chemistry-typed atom whose centroid→atom direction fixed the frame x
    // azimuth (aromatic-H: CG/CD2 per substrate_conventions). -1 for HN / none.
    int32_t frame_anchor_atom_index = -1;

    // Sources (un-summed).
    std::vector<SourceSlot> sources;

    // The producer's bare per-atom kernel cross-check (this stratum's kernel:
    // bs for ring current, mc for McConnell), read from the H5 at this frame.
    SphericalTensor bare_kernel;
    bool bare_kernel_present = false;

    // The DFT target (same for both row kinds).
    DftTarget target;
};

}  // namespace h5reader::rediscover
