#include "MutationDeltaResult.h"
#include "Protein.h"
#include "LegacyAmberTopology.h"
#include "OrcaShieldingResult.h"
#include "ApbsFieldResult.h"
#include "MopacResult.h"
#include "DsspResult.h"
#include "MolecularGraphResult.h"
#include "GeometryResult.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"
#include "SemanticEnums.h"
#include "generated/LegacyAmberSemanticTables.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <vector>

namespace nmr {

const Mat3 MutationDeltaResult::zero_mat3_ = Mat3::Zero();
const SphericalTensor MutationDeltaResult::zero_spherical_ = {};
const MatchedAtomData MutationDeltaResult::empty_match_ = {};


std::vector<std::type_index> MutationDeltaResult::Dependencies() const {
    return { std::type_index(typeid(OrcaShieldingResult)) };
}


// ============================================================================
// Atom matching: substrate-typed-identity, no spatial dependence
//
// Replaces the pre-2026-05-08 KD-tree+greedy-NN matcher (`MatchCloud`,
// `MatchTree`, `MATCH_TOLERANCE`, second-chance-pass, evicted-claims
// bookkeeping) with binding by `(residue_index, AtomMechanicalIdentity)`
// against the LegacyAmber substrate. The replaced matcher silently
// rebound across mutation sites, variant differences, and rotamer
// flips — fine for mechanical-swap PDB mutants where coordinates are
// pinned, indefensible as soon as coordinates drift. See the design
// pass at spec/plan/category-info-projection-implementation-plan-2026-05-08.md
// §8 (MutationDeltaResult matchup rewrite).
//
// Within a residue, atoms with collision-equal identity tuples (methyl
// HB1/HB2/HB3 etc.) bind by consume-in-order: each mut atom is taken
// at most once, in residue-atom-index order. The choice within an
// equivalence class is arbitrary, deterministic, and chemically
// immaterial — the substrate has stamped these atoms as equivalent.
// No atom-name string crosses the boundary.
//
// Spatial nearest-neighbor distance is computed AFTER binding as a
// reality check (drift on bound, sanity-check on rejected). It does
// not enter the binding decision.
// ============================================================================


// ============================================================================
// Ring proximity: cylindrical coordinates in ring frame
//
// From old project EquivariantExtractor: z = d.n, rho = |d - z*n|,
// theta = atan2(rho, z). McConnell factor = (3cos^2 theta - 1)/r^3.
// Exponential decay = exp(-r/4.0) with tau = 4A characteristic length.
// ============================================================================

static RingProximity ComputeRingProximity(
        const Vec3& atom_pos,
        size_t ring_index,
        RingTypeIndex ring_type,
        const RingGeometry& geom) {

    RingProximity rp;
    rp.ring_index = ring_index;
    rp.ring_type = ring_type;

    Vec3 d = atom_pos - geom.center;
    rp.distance = d.norm();

    if (rp.distance < MIN_DISTANCE) {
        rp.distance = MIN_DISTANCE;
        return rp;
    }

    // Cylindrical decomposition in ring frame
    rp.z = d.dot(geom.normal);
    Vec3 in_plane = d - rp.z * geom.normal;
    rp.rho = in_plane.norm();
    rp.theta = std::atan2(rp.rho, rp.z);

    // McConnell factor: (3cos^2 theta - 1) / r^3
    double cos_theta = rp.z / rp.distance;
    double r3 = rp.distance * rp.distance * rp.distance;
    rp.mcconnell_factor = (3.0 * cos_theta * cos_theta - 1.0) / r3;

    // Exponential decay with tau = 4A (ring current characteristic length)
    rp.exp_decay = std::exp(-rp.distance / EXP_DECAY_LENGTH);

    return rp;
}


// ============================================================================
// Detect mutation sites and identify which rings were removed
// ============================================================================

static std::vector<MutationSite> DetectMutationSites(
        const Protein& wt_protein,
        const Protein& mut_protein) {

    std::vector<MutationSite> sites;
    size_t n = std::min(wt_protein.ResidueCount(), mut_protein.ResidueCount());

    for (size_t ri = 0; ri < n; ++ri) {
        AminoAcid wt_type = wt_protein.ResidueAt(ri).type;
        AminoAcid mut_type = mut_protein.ResidueAt(ri).type;
        if (wt_type != mut_type) {
            MutationSite site;
            site.residue_index = ri;
            site.wt_type = wt_type;
            site.mut_type = mut_type;

            // Find WT rings at this residue
            for (size_t rgi = 0; rgi < wt_protein.RingCount(); ++rgi) {
                if (wt_protein.RingAt(rgi).parent_residue_index == ri) {
                    site.wt_ring_indices.push_back(rgi);
                }
            }
            sites.push_back(std::move(site));
        }
    }
    return sites;
}


// ============================================================================
// Build summary statistics
// ============================================================================

static DeltaSummary BuildSummary(
        const std::vector<MatchedAtomData>& matched) {

    DeltaSummary summary;

    // --- By element ---
    std::map<Element, std::vector<const MatchedAtomData*>> by_elem;
    for (const auto& m : matched)
        by_elem[m.element].push_back(&m);

    for (auto& [elem, atoms] : by_elem) {
        DeltaSummary::ElementBin bin;
        bin.element = elem;
        bin.count = static_cast<int>(atoms.size());
        double sum_t0 = 0, sum_abs_t0 = 0, max_abs = 0, sum_t2 = 0;
        for (const auto* a : atoms) {
            double t0 = a->delta_shielding_spherical.T0;
            sum_t0 += t0;
            sum_abs_t0 += std::abs(t0);
            max_abs = std::max(max_abs, std::abs(t0));
            sum_t2 += a->delta_shielding_spherical.T2Magnitude();
        }
        bin.mean_delta_t0 = sum_t0 / bin.count;
        bin.mean_abs_delta_t0 = sum_abs_t0 / bin.count;
        bin.max_abs_delta_t0 = max_abs;
        bin.mean_t2_magnitude = sum_t2 / bin.count;
        summary.by_element.push_back(bin);
    }

    // --- By distance to nearest removed ring (1A bins, 0-15A) ---
    for (int bin_start = 0; bin_start < 15; ++bin_start) {
        DeltaSummary::DistanceBin bin;
        bin.bin_start = static_cast<double>(bin_start);
        bin.bin_end = static_cast<double>(bin_start + 1);
        double sum_abs_t0 = 0, sum_t2 = 0;
        int count = 0;
        for (const auto& m : matched) {
            double d = m.nearest_removed_ring_dist;
            if (d >= bin.bin_start && d < bin.bin_end) {
                sum_abs_t0 += std::abs(m.delta_shielding_spherical.T0);
                sum_t2 += m.delta_shielding_spherical.T2Magnitude();
                count++;
            }
        }
        bin.count = count;
        if (count > 0) {
            bin.mean_abs_delta_t0 = sum_abs_t0 / count;
            bin.mean_t2_magnitude = sum_t2 / count;
        }
        summary.by_distance.push_back(bin);
    }

    // --- Backbone vs sidechain ---
    double bb_sum = 0, sc_sum = 0;
    for (const auto& m : matched) {
        double abs_t0 = std::abs(m.delta_shielding_spherical.T0);
        if (m.is_backbone) {
            summary.backbone_count++;
            bb_sum += abs_t0;
        } else {
            summary.sidechain_count++;
            sc_sum += abs_t0;
        }
    }
    if (summary.backbone_count > 0)
        summary.backbone_mean_abs_t0 = bb_sum / summary.backbone_count;
    if (summary.sidechain_count > 0)
        summary.sidechain_mean_abs_t0 = sc_sum / summary.sidechain_count;

    return summary;
}


// ============================================================================
// MutationDeltaResult::Compute
// ============================================================================

std::unique_ptr<MutationDeltaResult> MutationDeltaResult::Compute(
        ProteinConformation& wt_conf,
        const ProteinConformation& mut_conf) {

    OperationLog::Scope scope("MutationDeltaResult::Compute",
        "wt_atoms=" + std::to_string(wt_conf.AtomCount()) +
        " mut_atoms=" + std::to_string(mut_conf.AtomCount()));

    const Protein& wt_protein = wt_conf.ProteinRef();
    const Protein& mut_protein = mut_conf.ProteinRef();

    // ---- Precondition checks ----

    if (!wt_conf.HasResult<OrcaShieldingResult>()) {
        OperationLog::Error("MutationDeltaResult::Compute",
            "WT conformation missing OrcaShieldingResult");
        return nullptr;
    }
    if (!mut_conf.HasResult<OrcaShieldingResult>()) {
        OperationLog::Error("MutationDeltaResult::Compute",
            "mutant conformation missing OrcaShieldingResult");
        return nullptr;
    }
    if (wt_protein.ResidueCount() != mut_protein.ResidueCount()) {
        OperationLog::Error("MutationDeltaResult::Compute",
            "residue count mismatch: WT=" +
            std::to_string(wt_protein.ResidueCount()) +
            " mut=" + std::to_string(mut_protein.ResidueCount()));
        return nullptr;
    }

    // ---- Detect mutation sites and removed rings ----

    auto mutation_sites = DetectMutationSites(wt_protein, mut_protein);

    // Collect all removed ring indices
    std::vector<size_t> removed_ring_indices;
    for (const auto& site : mutation_sites) {
        for (size_t ri : site.wt_ring_indices)
            removed_ring_indices.push_back(ri);
        OperationLog::Info(LogAtomMapping, "MutationDeltaResult",
            "mutation at residue " + std::to_string(site.residue_index) +
            ": " + ThreeLetterCodeForAminoAcid(site.wt_type) +
            " -> " + ThreeLetterCodeForAminoAcid(site.mut_type) +
            " (" + std::to_string(site.wt_ring_indices.size()) + " rings removed)");
    }

    // ---- Check optional data sources ----

    bool has_apbs = wt_conf.HasResult<ApbsFieldResult>() &&
                    mut_conf.HasResult<ApbsFieldResult>();
    bool has_mopac = wt_conf.HasResult<MopacResult>() &&
                     mut_conf.HasResult<MopacResult>();
    bool has_dssp = wt_conf.HasResult<DsspResult>() &&
                    mut_conf.HasResult<DsspResult>();
    bool has_graph = wt_conf.HasResult<MolecularGraphResult>() &&
                     mut_conf.HasResult<MolecularGraphResult>();
    bool has_geom = wt_conf.HasResult<GeometryResult>();

    // ---- Substrate-required check ----

    if (!wt_protein.LegacyAmber().HasAtomSemantic() ||
        !mut_protein.LegacyAmber().HasAtomSemantic()) {
        OperationLog::Error("MutationDeltaResult::Compute",
            "typed-identity matching requires LegacyAmber substrate on "
            "both proteins (HasAtomSemantic() returned false). The "
            "spatial matcher this Result used to fall back on has been "
            "retired; no fallback path exists. Verify the upstream "
            "loader populated AtomSemantic via ComposeAtomSemantic.");
        return nullptr;
    }

    // ---- Per-residue lookup over mut atoms ----

    const size_t wt_count = wt_conf.AtomCount();
    const size_t mut_count = mut_conf.AtomCount();

    std::vector<bool> residue_is_mutation(wt_protein.ResidueCount(), false);
    for (const auto& site : mutation_sites)
        residue_is_mutation[site.residue_index] = true;

    // Per non-mutation residue: ordered list of (mut_atom_index, identity).
    // Atoms get consumed in residue-atom-index order as WT atoms claim
    // matches — equivalent-H sets bind position-by-position.
    std::map<size_t, std::vector<std::pair<size_t, AtomMechanicalIdentity>>>
        mut_residue_pool;
    for (size_t mi = 0; mi < mut_count; ++mi) {
        const size_t ri = mut_protein.AtomAt(mi).residue_index;
        if (ri >= residue_is_mutation.size() || residue_is_mutation[ri])
            continue;
        const auto& sem = mut_protein.LegacyAmber().SemanticAt(mi);
        mut_residue_pool[ri].push_back({mi,
            AtomMechanicalIdentity{sem.element, sem.locant, sem.branch,
                                    sem.di_index, sem.backbone_role}});
    }

    // ---- Atom matching: typed identity, residue-local, consume-in-order ----

    std::vector<size_t> wt_to_mut(wt_count, SIZE_MAX);

    // Diagnostic counters. Mutation-skipped is counted at the atom
    // level (not residue) so the totals add up to wt_count.
    size_t mutation_skipped = 0;
    size_t variant_unmatched = 0;
    size_t no_id_match = 0;
    std::set<size_t> variant_residues_seen;

    for (size_t wi = 0; wi < wt_count; ++wi) {
        const Atom& wt_a = wt_protein.AtomAt(wi);
        const size_t ri = wt_a.residue_index;

        if (residue_is_mutation[ri]) {
            ++mutation_skipped;
            continue;
        }

        // Variant difference at this residue is per-atom — common atoms
        // still bind, only variant-specific atoms (HID HD1 vs. HIE HE2)
        // surface as no-match. Track residue-level so the diagnostic
        // log can name them.
        const bool variant_diff =
            (wt_protein.ResidueAt(ri).protonation_variant_index !=
             mut_protein.ResidueAt(ri).protonation_variant_index);
        if (variant_diff) variant_residues_seen.insert(ri);

        const auto& wt_sem = wt_protein.LegacyAmber().SemanticAt(wi);
        const AtomMechanicalIdentity wt_id{
            wt_sem.element, wt_sem.locant, wt_sem.branch,
            wt_sem.di_index, wt_sem.backbone_role};

        // Consume the first matching mut atom in this residue's pool.
        // Equivalent-H sets land in residue-atom-index order — this is
        // chemically arbitrary but deterministic, and the substrate has
        // already declared these atoms equivalent.
        auto& pool = mut_residue_pool[ri];
        bool bound = false;
        for (auto it = pool.begin(); it != pool.end(); ++it) {
            if (it->second == wt_id) {
                wt_to_mut[wi] = it->first;
                pool.erase(it);
                bound = true;
                break;
            }
        }
        if (!bound) {
            if (variant_diff) ++variant_unmatched;
            else              ++no_id_match;
        }
    }

    // ---- Spatial-NN sanity check (post-binding) ----
    //
    // For bound atoms: drift = ||pos_wt - pos_mut||. Mechanical-swap
    // mutants put non-mutation atoms at near-identical coordinates; large
    // drift surfaces rotamer flips or coordinate drift, which the typed
    // matcher handles correctly while this metric makes them visible.
    //
    // For rejected atoms: same-element nearest mut atom + distance. A
    // nearby same-element neighbour confirms the mechanical-mutation
    // pipeline's structural integrity (rejection is chemistry-only). The
    // distance does NOT enter the binding decision; it's the methodology
    // section's reality check.

    auto SameElementSpatialNn = [&](size_t wi) -> double {
        const Vec3 wp = wt_conf.PositionAt(wi);
        const Element we = wt_protein.AtomAt(wi).element;
        double best = std::numeric_limits<double>::infinity();
        for (size_t mi = 0; mi < mut_count; ++mi) {
            if (mut_protein.AtomAt(mi).element != we) continue;
            const double d = (wp - mut_conf.PositionAt(mi)).norm();
            if (d < best) best = d;
        }
        return best;
    };

    std::vector<double> bound_drifts;
    bound_drifts.reserve(wt_count);
    for (size_t wi = 0; wi < wt_count; ++wi) {
        if (wt_to_mut[wi] == SIZE_MAX) continue;
        bound_drifts.push_back(
            (wt_conf.PositionAt(wi)
             - mut_conf.PositionAt(wt_to_mut[wi])).norm());
    }
    double drift_med = 0.0, drift_max = 0.0;
    if (!bound_drifts.empty()) {
        std::vector<double> sorted = bound_drifts;
        std::sort(sorted.begin(), sorted.end());
        drift_med = sorted[sorted.size() / 2];
        drift_max = sorted.back();
    }

    size_t reject_within_0_5 = 0, reject_within_2_0 = 0, reject_far = 0;
    for (size_t wi = 0; wi < wt_count; ++wi) {
        if (wt_to_mut[wi] != SIZE_MAX) continue;
        const double d = SameElementSpatialNn(wi);
        if      (d < 0.5) ++reject_within_0_5;
        else if (d < 2.0) ++reject_within_2_0;
        else              ++reject_far;
    }

    // ---- DSSP residue-level lookup helpers ----

    const DsspResult* wt_dssp = has_dssp ? &wt_conf.Result<DsspResult>() : nullptr;
    const DsspResult* mut_dssp = has_dssp ? &mut_conf.Result<DsspResult>() : nullptr;

    // ---- Compute deltas for matched pairs ----

    auto result = std::make_unique<MutationDeltaResult>();
    result->wt_conf_ = &wt_conf;
    result->mutation_sites_ = std::move(mutation_sites);
    result->has_apbs_delta_ = has_apbs;
    result->has_mopac_delta_ = has_mopac;
    result->has_dssp_delta_ = has_dssp;
    result->has_graph_delta_ = has_graph;
    result->wt_to_matched_.assign(wt_count, SIZE_MAX);

    size_t matched = 0;
    size_t unmatched = 0;

    for (size_t wi = 0; wi < wt_count; ++wi) {
        size_t mi = wt_to_mut[wi];
        if (mi == SIZE_MAX) { unmatched++; continue; }

        MatchedAtomData data;
        data.wt_index = wi;
        data.mut_index = mi;
        // Drift between identity-bound atoms — diagnostic only (the
        // binding criterion was substrate identity, not proximity).
        data.match_distance = (wt_conf.PositionAt(wi)
                               - mut_conf.PositionAt(mi)).norm();

        // WT atom identity (typed, not strings)
        data.element = wt_protein.AtomAt(wi).element;
        data.residue_index = wt_protein.AtomAt(wi).residue_index;
        const auto& wt_ca = wt_conf.AtomAt(wi);
        data.role = wt_ca.role;
        data.is_backbone = wt_ca.is_backbone;

        const auto& mut_ca = mut_conf.AtomAt(mi);

        // DFT shielding total delta: WT - mutant.
        data.delta_shielding = wt_ca.orca_shielding_total - mut_ca.orca_shielding_total;
        data.delta_shielding_spherical = SphericalTensor::Decompose(data.delta_shielding);

        // DFT shielding component decomposition (dia + para = total).
        // Sides are direct copies; deltas are wt - mut.
        data.wt_shielding_diamagnetic            = wt_ca.orca_shielding_diamagnetic;
        data.wt_shielding_diamagnetic_spherical  = wt_ca.orca_shielding_diamagnetic_spherical;
        data.wt_shielding_paramagnetic           = wt_ca.orca_shielding_paramagnetic;
        data.wt_shielding_paramagnetic_spherical = wt_ca.orca_shielding_paramagnetic_spherical;

        data.mut_shielding_diamagnetic            = mut_ca.orca_shielding_diamagnetic;
        data.mut_shielding_diamagnetic_spherical  = mut_ca.orca_shielding_diamagnetic_spherical;
        data.mut_shielding_paramagnetic           = mut_ca.orca_shielding_paramagnetic;
        data.mut_shielding_paramagnetic_spherical = mut_ca.orca_shielding_paramagnetic_spherical;

        data.delta_shielding_diamagnetic =
            wt_ca.orca_shielding_diamagnetic - mut_ca.orca_shielding_diamagnetic;
        data.delta_shielding_diamagnetic_spherical =
            SphericalTensor::Decompose(data.delta_shielding_diamagnetic);
        data.delta_shielding_paramagnetic =
            wt_ca.orca_shielding_paramagnetic - mut_ca.orca_shielding_paramagnetic;
        data.delta_shielding_paramagnetic_spherical =
            SphericalTensor::Decompose(data.delta_shielding_paramagnetic);

        // ff14SB charge delta
        data.delta_partial_charge = wt_ca.partial_charge - mut_ca.partial_charge;

        // APBS delta
        if (has_apbs) {
            data.delta_efield = wt_ca.apbs_efield - mut_ca.apbs_efield;
            data.delta_efg = wt_ca.apbs_efg - mut_ca.apbs_efg;
            data.delta_efg_spherical = SphericalTensor::Decompose(data.delta_efg);
            data.has_apbs_delta = true;
        }

        // MOPAC charge delta
        if (has_mopac) {
            data.delta_mopac_charge = wt_ca.mopac_charge - mut_ca.mopac_charge;
            data.has_mopac_delta = true;
        }

        // DSSP delta (per residue, projected to atom)
        if (has_dssp) {
            size_t wt_ri = data.residue_index;
            size_t mut_ri = mut_protein.AtomAt(mi).residue_index;
            data.delta_phi = wt_dssp->Phi(wt_ri) - mut_dssp->Phi(mut_ri);
            data.delta_psi = wt_dssp->Psi(wt_ri) - mut_dssp->Psi(mut_ri);
            data.delta_sasa = wt_dssp->SASA(wt_ri) - mut_dssp->SASA(mut_ri);
            data.has_dssp_delta = true;
        }

        // Graph delta
        if (has_graph) {
            data.delta_graph_dist_ring = wt_ca.graph_dist_ring - mut_ca.graph_dist_ring;
            data.delta_bfs_decay = wt_ca.bfs_decay - mut_ca.bfs_decay;
            data.has_graph_delta = true;
        }

        // Ring proximity to removed rings
        if (has_geom && !removed_ring_indices.empty()) {
            Vec3 atom_pos = wt_conf.PositionAt(wi);
            data.nearest_removed_ring_dist = 99.0;

            for (size_t rri : removed_ring_indices) {
                const Ring& ring = wt_protein.RingAt(rri);
                const RingGeometry& geom = wt_conf.ring_geometries[rri];

                RingProximity rp = ComputeRingProximity(
                    atom_pos, rri, ring.type_index, geom);
                data.removed_ring_proximity.push_back(rp);

                if (rp.distance < data.nearest_removed_ring_dist) {
                    data.nearest_removed_ring_dist = rp.distance;
                    data.nearest_removed_ring = data.removed_ring_proximity.size() - 1;
                }
            }
        }

        result->wt_to_matched_[wi] = result->matched_atoms_.size();
        result->matched_atoms_.push_back(std::move(data));
        matched++;
    }

    // ---- Build summary statistics ----

    result->summary_ = BuildSummary(result->matched_atoms_);

    // ---- Log ----

    OperationLog::Info(LogAtomMapping, "MutationDeltaResult::Compute",
        "matched=" + std::to_string(matched) +
        " mutation_skipped=" + std::to_string(mutation_skipped) +
        " variant_residues=" + std::to_string(variant_residues_seen.size()) +
        " variant_unmatched=" + std::to_string(variant_unmatched) +
        " no_id_match=" + std::to_string(no_id_match) +
        " mutation_sites=" + std::to_string(result->mutation_sites_.size()) +
        " removed_rings=" + std::to_string(removed_ring_indices.size()) +
        " has_apbs=" + std::to_string(has_apbs) +
        " has_dssp=" + std::to_string(has_dssp) +
        " has_graph=" + std::to_string(has_graph));

    // Drift on bound atoms (rotamer-flip / coordinate-drift indicator).
    OperationLog::Info(LogAtomMapping, "MutationDeltaResult::Compute",
        "drift_med=" + std::to_string(drift_med) +
        " drift_max=" + std::to_string(drift_max) +
        " bound_count=" + std::to_string(bound_drifts.size()));

    // Spatial-NN sanity check on rejections — confirms mechanical-swap
    // structural integrity without the binding criterion depending on
    // it. See methodology paragraph in the spec/plan/category-info-
    // projection-implementation-plan-2026-05-08.md slice for the
    // intended thesis prose.
    OperationLog::Info(LogAtomMapping, "MutationDeltaResult::Compute",
        "rejection_spatial_nn"
        " within_0.5A=" + std::to_string(reject_within_0_5) +
        " within_2A=" + std::to_string(reject_within_2_0) +
        " far=" + std::to_string(reject_far));

    // Log summary by element
    for (const auto& bin : result->summary_.by_element) {
        OperationLog::Info(LogAtomMapping, "MutationDeltaResult",
            "element=" + SymbolForElement(bin.element) +
            " n=" + std::to_string(bin.count) +
            " mean_|dT0|=" + std::to_string(bin.mean_abs_delta_t0) +
            " max_|dT0|=" + std::to_string(bin.max_abs_delta_t0) +
            " mean_|T2|=" + std::to_string(bin.mean_t2_magnitude));
    }

    // Log distance decay (the 1/r^3 diagnostic)
    for (const auto& bin : result->summary_.by_distance) {
        if (bin.count > 0) {
            OperationLog::Info(LogAtomMapping, "MutationDeltaResult",
                "dist=[" + std::to_string((int)bin.bin_start) + "-" +
                std::to_string((int)bin.bin_end) + "A] n=" +
                std::to_string(bin.count) +
                " mean_|dT0|=" + std::to_string(bin.mean_abs_delta_t0) +
                " mean_|T2|=" + std::to_string(bin.mean_t2_magnitude));
        }
    }

    return result;
}


// ============================================================================
// Query methods
// ============================================================================

size_t MutationDeltaResult::UnmatchedWtAtomCount() const {
    size_t n = 0;
    for (size_t idx : wt_to_matched_) if (idx == SIZE_MAX) n++;
    return n;
}

bool MutationDeltaResult::HasMatch(size_t i) const {
    return i < wt_to_matched_.size() && wt_to_matched_[i] != SIZE_MAX;
}

const MatchedAtomData& MutationDeltaResult::MatchedDataAt(size_t i) const {
    if (i >= wt_to_matched_.size() || wt_to_matched_[i] == SIZE_MAX)
        return empty_match_;
    return matched_atoms_[wt_to_matched_[i]];
}

const Mat3& MutationDeltaResult::DeltaShieldingAt(size_t i) const {
    if (!HasMatch(i)) return zero_mat3_;
    return matched_atoms_[wt_to_matched_[i]].delta_shielding;
}

const SphericalTensor& MutationDeltaResult::DeltaShieldingSphericalAt(size_t i) const {
    if (!HasMatch(i)) return zero_spherical_;
    return matched_atoms_[wt_to_matched_[i]].delta_shielding_spherical;
}

double MutationDeltaResult::DeltaT0At(size_t i) const {
    return DeltaShieldingSphericalAt(i).T0;
}

Vec3 MutationDeltaResult::DeltaEFieldAt(size_t i) const {
    if (!HasMatch(i)) return Vec3::Zero();
    return matched_atoms_[wt_to_matched_[i]].delta_efield;
}

Mat3 MutationDeltaResult::DeltaEFGAt(size_t i) const {
    if (!HasMatch(i)) return Mat3::Zero();
    return matched_atoms_[wt_to_matched_[i]].delta_efg;
}

double MutationDeltaResult::DeltaPartialChargeAt(size_t i) const {
    if (!HasMatch(i)) return 0.0;
    return matched_atoms_[wt_to_matched_[i]].delta_partial_charge;
}

double MutationDeltaResult::DeltaMopacChargeAt(size_t i) const {
    if (!HasMatch(i)) return 0.0;
    return matched_atoms_[wt_to_matched_[i]].delta_mopac_charge;
}

double MutationDeltaResult::NearestRemovedRingDistance(size_t i) const {
    if (!HasMatch(i)) return 99.0;
    return matched_atoms_[wt_to_matched_[i]].nearest_removed_ring_dist;
}

size_t MutationDeltaResult::MutantAtomFor(size_t i) const {
    if (!HasMatch(i)) return SIZE_MAX;
    return matched_atoms_[wt_to_matched_[i]].mut_index;
}

double MutationDeltaResult::MatchDistanceAt(size_t i) const {
    if (!HasMatch(i)) return 0.0;
    return matched_atoms_[wt_to_matched_[i]].match_distance;
}

// ============================================================================
// WriteFeatures: export delta arrays indexed by WT atom index.
//
// Unmatched atoms get zeros. Arrays are (N, K) where N = WT atom count.
// ============================================================================

static void PackST(const SphericalTensor& st, double* out) {
    out[0] = st.T0;
    for (int i = 0; i < 3; ++i) out[1+i] = st.T1[i];
    for (int i = 0; i < 5; ++i) out[4+i] = st.T2[i];
}

int MutationDeltaResult::WriteFeatures(const ProteinConformation& conf,
                                        const std::string& output_dir) const {
    const size_t N = conf.AtomCount();
    int written = 0;

    // delta_shielding: (N, 9) — DFT shielding total delta as SphericalTensor.
    {
        std::vector<double> data(N * 9, 0.0);
        for (size_t i = 0; i < N; ++i)
            if (HasMatch(i))
                PackST(matched_atoms_[wt_to_matched_[i]].delta_shielding_spherical, &data[i*9]);
        NpyWriter::WriteFloat64(output_dir + "/delta_shielding.npy", data.data(), N, 9);
        written++;
    }

    // DFT shielding component decomposition: WT side, mut side, and
    // delta for both diamagnetic and paramagnetic. All (N, 9) packed
    // SphericalTensors. Total = dia + para; the existing delta_shielding
    // and these two component deltas satisfy delta_shielding ==
    // delta_shielding_diamagnetic + delta_shielding_paramagnetic
    // analytically.
    auto WriteShieldingComponent = [&](const std::string& filename,
                                        SphericalTensor MatchedAtomData::* field) {
        std::vector<double> data(N * 9, 0.0);
        for (size_t i = 0; i < N; ++i)
            if (HasMatch(i))
                PackST(matched_atoms_[wt_to_matched_[i]].*field, &data[i*9]);
        NpyWriter::WriteFloat64(output_dir + "/" + filename, data.data(), N, 9);
    };

    WriteShieldingComponent("wt_shielding_diamagnetic.npy",
        &MatchedAtomData::wt_shielding_diamagnetic_spherical);
    WriteShieldingComponent("wt_shielding_paramagnetic.npy",
        &MatchedAtomData::wt_shielding_paramagnetic_spherical);
    WriteShieldingComponent("mut_shielding_diamagnetic.npy",
        &MatchedAtomData::mut_shielding_diamagnetic_spherical);
    WriteShieldingComponent("mut_shielding_paramagnetic.npy",
        &MatchedAtomData::mut_shielding_paramagnetic_spherical);
    WriteShieldingComponent("delta_shielding_diamagnetic.npy",
        &MatchedAtomData::delta_shielding_diamagnetic_spherical);
    WriteShieldingComponent("delta_shielding_paramagnetic.npy",
        &MatchedAtomData::delta_shielding_paramagnetic_spherical);
    written += 6;

    // delta_scalars: (N, 6) — [matched, delta_T0, nearest_ring_dist, delta_charge, delta_mopac_charge, match_dist]
    {
        std::vector<double> data(N * 6, 0.0);
        for (size_t i = 0; i < N; ++i) {
            if (HasMatch(i)) {
                const auto& m = matched_atoms_[wt_to_matched_[i]];
                data[i*6 + 0] = 1.0;
                data[i*6 + 1] = m.delta_shielding_spherical.T0;
                data[i*6 + 2] = m.nearest_removed_ring_dist;
                data[i*6 + 3] = m.delta_partial_charge;
                data[i*6 + 4] = m.delta_mopac_charge;
                data[i*6 + 5] = m.match_distance;
            }
        }
        NpyWriter::WriteFloat64(output_dir + "/delta_scalars.npy", data.data(), N, 6);
        written++;
    }

    // delta_apbs: (N, 12) — [delta_E(3), delta_EFG_spherical(9)]
    if (has_apbs_delta_) {
        std::vector<double> data(N * 12, 0.0);
        for (size_t i = 0; i < N; ++i) {
            if (HasMatch(i)) {
                const auto& m = matched_atoms_[wt_to_matched_[i]];
                data[i*12 + 0] = m.delta_efield.x();
                data[i*12 + 1] = m.delta_efield.y();
                data[i*12 + 2] = m.delta_efield.z();
                PackST(m.delta_efg_spherical, &data[i*12 + 3]);
            }
        }
        NpyWriter::WriteFloat64(output_dir + "/delta_apbs.npy", data.data(), N, 12);
        written++;
    }

    // delta_ring_proximity: (N, R, 6) flattened to (N, R*6) where R = removed ring count
    // Per removed ring: [distance, z, rho, theta, mcconnell_factor, exp_decay]
    if (!mutation_sites_.empty()) {
        size_t total_removed = 0;
        for (const auto& site : mutation_sites_)
            total_removed += site.wt_ring_indices.size();

        if (total_removed > 0) {
            const size_t cols = total_removed * 6;
            std::vector<double> data(N * cols, 0.0);
            for (size_t i = 0; i < N; ++i) {
                if (HasMatch(i)) {
                    const auto& m = matched_atoms_[wt_to_matched_[i]];
                    for (size_t r = 0; r < m.removed_ring_proximity.size() && r < total_removed; ++r) {
                        const auto& rp = m.removed_ring_proximity[r];
                        size_t base = i * cols + r * 6;
                        data[base + 0] = rp.distance;
                        data[base + 1] = rp.z;
                        data[base + 2] = rp.rho;
                        data[base + 3] = rp.theta;
                        data[base + 4] = rp.mcconnell_factor;
                        data[base + 5] = rp.exp_decay;
                    }
                }
            }
            NpyWriter::WriteFloat64(output_dir + "/delta_ring_proximity.npy", data.data(), N, cols);
            written++;
        }
    }

    return written;
}

}  // namespace nmr
