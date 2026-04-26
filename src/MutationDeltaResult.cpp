#include "MutationDeltaResult.h"
#include "Protein.h"
#include "AtomReference.h"
#include "OrcaShieldingResult.h"
#include "ApbsFieldResult.h"
#include "MopacResult.h"
#include "DsspResult.h"
#include "MolecularGraphResult.h"
#include "GeometryResult.h"
#include "NpyWriter.h"
#include "OperationLog.h"
#include "PhysicalConstants.h"

#include <cmath>
#include <algorithm>
#include <map>

namespace nmr {

const Mat3 MutationDeltaResult::zero_mat3_ = Mat3::Zero();
const SphericalTensor MutationDeltaResult::zero_spherical_ = {};
const MatchedAtomData MutationDeltaResult::empty_match_ = {};


std::vector<std::type_index> MutationDeltaResult::Dependencies() const {
    return { std::type_index(typeid(OrcaShieldingResult)) };
}


// ============================================================================
// Atom matching: typed cross-protein key (AtomLocator).
//
// Pre-2026-04-26 used a KD-tree on positions with a 0.5A tolerance + element
// filter. Replaced with AtomLocator (chain, residue_position, atom_name) —
// a typed cross-protein identity. Backbone atoms at mutation sites match
// across the residue type change (PHE5 N → ALA5 N) because their locators
// are identical; side-chain atoms unique to one side (PHE's CG, CD1, ... vs
// ALA's HB1) naturally fail to match. No geometric tolerance, no fall-back
// to position-based matching: identity is symbolic and the proteins are
// expected to have come from the same tleap setup.
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
            double t0 = a->delta_shielding_total_spherical.T0;
            sum_t0 += t0;
            sum_abs_t0 += std::abs(t0);
            max_abs = std::max(max_abs, std::abs(t0));
            sum_t2 += a->delta_shielding_total_spherical.T2Magnitude();
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
                sum_abs_t0 += std::abs(m.delta_shielding_total_spherical.T0);
                sum_t2 += m.delta_shielding_total_spherical.T2Magnitude();
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
        double abs_t0 = std::abs(m.delta_shielding_total_spherical.T0);
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

    // ---- Atom matching: typed AtomLocator (chain, position, atom_name) ----
    //
    // Build a map over the mutant protein keyed by AtomLocator. Each WT
    // atom looks itself up by its own AtomLocator. A match is exact at the
    // typed level — no geometric tolerance, no nearest-neighbor search.
    // Side-chain atoms unique to one side (PHE's CG vs ALA's HB1) fail
    // the lookup. Element is cross-checked as a sanity guard against
    // naming weirdness; mismatch logs a warning and the pair is skipped.

    const size_t wt_count = wt_conf.AtomCount();
    const size_t mut_count = mut_conf.AtomCount();

    auto mut_locator_map = BuildAtomLocatorMap(mut_protein);
    std::vector<size_t> wt_to_mut(wt_count, SIZE_MAX);
    int element_mismatches = 0;

    for (size_t wi = 0; wi < wt_count; ++wi) {
        AtomLocator loc = MakeAtomLocator(wt_protein, wi);
        auto it = mut_locator_map.find(loc);
        if (it == mut_locator_map.end()) continue;       // unmatched (mutation-site sidechain etc.)
        size_t mi = it->second;

        const Element wt_elem = wt_protein.AtomAt(wi).element;
        const Element mut_elem = mut_protein.AtomAt(mi).element;
        if (wt_elem != mut_elem) {
            element_mismatches++;
            const Residue& wt_res = wt_protein.ResidueAt(loc.residue_index);
            const char* res_code = wt_res.AminoAcidInfo().three_letter_code;
            OperationLog::Warn("MutationDeltaResult::Compute",
                "element mismatch at " + std::string(res_code) +
                std::to_string(wt_res.sequence_number) + ":" +
                loc.atom_name.AsString() + " — WT=" +
                SymbolForElement(wt_elem) + " mut=" +
                SymbolForElement(mut_elem) + "; skipping");
            continue;
        }

        wt_to_mut[wi] = mi;
    }

    if (element_mismatches > 0) {
        OperationLog::Warn("MutationDeltaResult::Compute",
            std::to_string(element_mismatches) + " atom locators had element disagreement");
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
        // match_distance is now Euclidean position delta — typed match
        // happened at the locator level; the position delta is reported
        // for diagnostic purposes (e.g., backbone stability across the
        // mutation, MD-frame relaxation).
        data.match_distance = (wt_conf.PositionAt(wi) - mut_conf.PositionAt(mi)).norm();

        // WT atom identity (typed, not strings)
        data.element = wt_protein.AtomAt(wi).element;
        data.residue_index = wt_protein.AtomAt(wi).residue_index;
        const auto& wt_ca = wt_conf.AtomAt(wi);
        data.role = wt_ca.role;
        data.is_backbone = wt_ca.is_backbone;

        const auto& mut_ca = mut_conf.AtomAt(mi);

        // DFT shielding — six full tensors per channel, copied from
        // ConformationAtom (OrcaShieldingResult populated them at protein
        // load). Three deltas computed by subtraction; the spherical
        // decomposition is recomputed from the delta Mat3 to preserve the
        // T0+T1+T2 structure.
        data.wt_shielding_total                  = wt_ca.orca_shielding_total;
        data.wt_shielding_total_spherical        = wt_ca.orca_shielding_total_spherical;
        data.wt_shielding_diamagnetic            = wt_ca.orca_shielding_diamagnetic;
        data.wt_shielding_diamagnetic_spherical  = wt_ca.orca_shielding_diamagnetic_spherical;
        data.wt_shielding_paramagnetic           = wt_ca.orca_shielding_paramagnetic;
        data.wt_shielding_paramagnetic_spherical = wt_ca.orca_shielding_paramagnetic_spherical;

        data.mut_shielding_total                  = mut_ca.orca_shielding_total;
        data.mut_shielding_total_spherical        = mut_ca.orca_shielding_total_spherical;
        data.mut_shielding_diamagnetic            = mut_ca.orca_shielding_diamagnetic;
        data.mut_shielding_diamagnetic_spherical  = mut_ca.orca_shielding_diamagnetic_spherical;
        data.mut_shielding_paramagnetic           = mut_ca.orca_shielding_paramagnetic;
        data.mut_shielding_paramagnetic_spherical = mut_ca.orca_shielding_paramagnetic_spherical;

        data.delta_shielding_total = data.wt_shielding_total - data.mut_shielding_total;
        data.delta_shielding_total_spherical =
            SphericalTensor::Decompose(data.delta_shielding_total);
        data.delta_shielding_diamagnetic =
            data.wt_shielding_diamagnetic - data.mut_shielding_diamagnetic;
        data.delta_shielding_diamagnetic_spherical =
            SphericalTensor::Decompose(data.delta_shielding_diamagnetic);
        data.delta_shielding_paramagnetic =
            data.wt_shielding_paramagnetic - data.mut_shielding_paramagnetic;
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
        " unmatched_wt=" + std::to_string(unmatched) +
        " mutation_sites=" + std::to_string(result->mutation_sites_.size()) +
        " removed_rings=" + std::to_string(removed_ring_indices.size()) +
        " has_apbs=" + std::to_string(has_apbs) +
        " has_dssp=" + std::to_string(has_dssp) +
        " has_graph=" + std::to_string(has_graph));

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
    return matched_atoms_[wt_to_matched_[i]].delta_shielding_total;
}

const SphericalTensor& MutationDeltaResult::DeltaShieldingSphericalAt(size_t i) const {
    if (!HasMatch(i)) return zero_spherical_;
    return matched_atoms_[wt_to_matched_[i]].delta_shielding_total_spherical;
}

double MutationDeltaResult::DeltaT0At(size_t i) const {
    return DeltaShieldingSphericalAt(i).T0;
}

const Mat3& MutationDeltaResult::DeltaShieldingDiamagneticAt(size_t i) const {
    if (!HasMatch(i)) return zero_mat3_;
    return matched_atoms_[wt_to_matched_[i]].delta_shielding_diamagnetic;
}

const SphericalTensor& MutationDeltaResult::DeltaShieldingDiamagneticSphericalAt(size_t i) const {
    if (!HasMatch(i)) return zero_spherical_;
    return matched_atoms_[wt_to_matched_[i]].delta_shielding_diamagnetic_spherical;
}

const Mat3& MutationDeltaResult::DeltaShieldingParamagneticAt(size_t i) const {
    if (!HasMatch(i)) return zero_mat3_;
    return matched_atoms_[wt_to_matched_[i]].delta_shielding_paramagnetic;
}

const SphericalTensor& MutationDeltaResult::DeltaShieldingParamagneticSphericalAt(size_t i) const {
    if (!HasMatch(i)) return zero_spherical_;
    return matched_atoms_[wt_to_matched_[i]].delta_shielding_paramagnetic_spherical;
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

    // Helper: pack a SphericalTensor field selected per matched atom into
    // a (N, 9) array, zero rows for unmatched atoms, write NPY.
    auto write_st_array = [&](const std::string& filename,
                               const SphericalTensor MatchedAtomData::*field) {
        std::vector<double> data(N * 9, 0.0);
        for (size_t i = 0; i < N; ++i) {
            if (!HasMatch(i)) continue;
            PackST(matched_atoms_[wt_to_matched_[i]].*field, &data[i*9]);
        }
        if (NpyWriter::WriteFloat64(output_dir + "/" + filename,
                                     data.data(), N, 9))
            written++;
    };

    // delta_shielding.npy preserves its filename (downstream consumers).
    // Content unchanged: total-channel delta as SphericalTensor.
    write_st_array("delta_shielding.npy",
                   &MatchedAtomData::delta_shielding_total_spherical);

    // 2026-04-26: per-channel deltas — paramagnetic isolates the
    // aromatic-ring and heavy-atom contributions, diamagnetic the
    // local-density baseline. Both required for residual analysis.
    write_st_array("delta_shielding_diamagnetic.npy",
                   &MatchedAtomData::delta_shielding_diamagnetic_spherical);
    write_st_array("delta_shielding_paramagnetic.npy",
                   &MatchedAtomData::delta_shielding_paramagnetic_spherical);

    // 2026-04-26: WT and mutant shielding stored separately, so the
    // calibration pipeline can compute deltas at arbitrary irrep
    // recombinations rather than only the precomputed delta.
    write_st_array("wt_shielding_total.npy",
                   &MatchedAtomData::wt_shielding_total_spherical);
    write_st_array("wt_shielding_diamagnetic.npy",
                   &MatchedAtomData::wt_shielding_diamagnetic_spherical);
    write_st_array("wt_shielding_paramagnetic.npy",
                   &MatchedAtomData::wt_shielding_paramagnetic_spherical);
    write_st_array("mut_shielding_total.npy",
                   &MatchedAtomData::mut_shielding_total_spherical);
    write_st_array("mut_shielding_diamagnetic.npy",
                   &MatchedAtomData::mut_shielding_diamagnetic_spherical);
    write_st_array("mut_shielding_paramagnetic.npy",
                   &MatchedAtomData::mut_shielding_paramagnetic_spherical);

    // -------------------------------------------------------------------
    // mutation_match.npy (N,) int32 — per-WT-atom mutant atom_index.
    //
    // For each WT atom, the matched mutant atom_index, or -1 if unmatched
    // (mutation-site sidechain atoms removed at PHE→ALA, etc.). Atom
    // indices fit comfortably in int32 for our use case; -1 is the
    // sentinel matching the AtomLocator ostream convention.
    // -------------------------------------------------------------------
    {
        std::vector<int32_t> match(N, -1);
        for (size_t i = 0; i < N; ++i) {
            if (HasMatch(i)) {
                match[i] = static_cast<int32_t>(
                    matched_atoms_[wt_to_matched_[i]].mut_index);
            }
        }
        if (NpyWriter::WriteInt32(output_dir + "/mutation_match.npy",
                                   match.data(), N))
            written++;
    }

    // -------------------------------------------------------------------
    // mutation_removed_rings.npy (R, 3) int32 — rings removed at mutation
    // sites.
    //
    // Three columns per removed ring:
    //   0  residue_index    parent residue (in WT) of the removed ring
    //   1  ring_index       WT ring index (into protein.RingAt)
    //   2  ring_type        RingTypeIndex enum (PheBenzene/TyrPhenol/...)
    //
    // R = total rings removed across all mutation sites. Empty (no NPY
    // emitted) when there are no mutations at aromatic-ring residues.
    // -------------------------------------------------------------------
    {
        const Protein& wt_protein = conf.ProteinRef();
        size_t total_removed = 0;
        for (const auto& site : mutation_sites_)
            total_removed += site.wt_ring_indices.size();

        if (total_removed > 0) {
            std::vector<int32_t> data(total_removed * 3, 0);
            size_t row = 0;
            for (const auto& site : mutation_sites_) {
                for (size_t ri : site.wt_ring_indices) {
                    const Ring& ring = wt_protein.RingAt(ri);
                    int32_t* r = &data[row * 3];
                    r[0] = static_cast<int32_t>(site.residue_index);
                    r[1] = static_cast<int32_t>(ri);
                    r[2] = static_cast<int32_t>(ring.type_index);
                    row++;
                }
            }
            if (NpyWriter::WriteInt32(output_dir + "/mutation_removed_rings.npy",
                                       data.data(), total_removed * 3))
                written++;
        }
    }

    // delta_scalars: (N, 6) — [matched, delta_T0, nearest_ring_dist, delta_charge, delta_mopac_charge, match_dist]
    {
        std::vector<double> data(N * 6, 0.0);
        for (size_t i = 0; i < N; ++i) {
            if (HasMatch(i)) {
                const auto& m = matched_atoms_[wt_to_matched_[i]];
                data[i*6 + 0] = 1.0;
                data[i*6 + 1] = m.delta_shielding_total_spherical.T0;
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
